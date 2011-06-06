/*

mfcc~ - A mel Frequency Cepstral Analysis external.

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.


version 0.2.5, May 16, 2010

¥ 0.2.5 changed filterbank_multiply so that sum of power in each filter is divided by the number of points in the filter.
¥Ê0.2.4 added an ifndef M_PI for guaranteed windows compilation
¥ 0.2.3 adds a #define M_PI for windows compilation, and declares all functions except _setup static. also changed filter_power[] array in _filterbank_multiply to be resizable rather than fixed at declaration time.
¥ 0.2.2 is part of the update that ensures all function names are prepended by the external name (bfcc_ or bfcc_tilde_, etc).
¥ 0.2.1 adds flexible filterwidth memory allocation
¥ 0.2.0 implements mayer_realfft 
¥ 0.1.9 added normalization option

*/

#include "m_pd.h"
#include <math.h>
#include <limits.h>
#include <stdio.h>
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

static t_class *mfcc_tilde_class;

typedef struct filter
{
    float *filter;
} t_filter;

typedef struct indices
{
    int index[2];
} t_indices;

typedef struct _mfcc_tilde
{
    t_object x_obj;
    float sr;
    float n;
	float window;
	int window_function;
	int debug;
    int overlap;
    int power_spectrum;
    int normalize;
    double last_dsp_time;
    int num_filters;
	int *x_filterbank_sizes;
    float mel_spacing;
	t_filter *x_filterbank;
	t_indices *x_indices;
	t_sample *signal_R;
    float *mfcc;
    float *hann;
    float x_f;
    t_outlet *x_featureList;

} t_mfcc_tilde;




/* ---------------- utility functions for building filterbank ---------------------- */

static int mfcc_tilde_nearest_bin_index(float target, float *bin_freqs)
{
	int i, nearest;
	float after, before;
	
	i=0;
	before=0;
	after=0;
	nearest=0;
	
	while(bin_freqs[i] < target)
		i++;
	
	after = bin_freqs[i];
	
	if(i != 0)
		before = bin_freqs[i-1];
	
	if( fabs(before-target) < fabs(after-target) )
		nearest = i-1;
	else
		nearest = i;
		
		
	return(nearest);
}


static float mfcc_tilde_mel2freq(float mel)
{
	float freq;
	
	freq = 700 * (exp(mel/1127.01048) - 1);
	if(freq < 0)
		freq = 0;
		
	return(freq);
}


static void mfcc_tilde_make_filterbank(t_mfcc_tilde *x, float *mel_spaced_freqs)
{	
	int i, j, k, nearest, window_half;
	float *bin_freqs;
	
	window_half = x->window * 0.5;
	
	// create local memory
	bin_freqs = (float *)getbytes(0);
	bin_freqs = (float *)t_resizebytes(bin_freqs, 0, window_half * sizeof(float));

	nearest=0;	
 			
	// first, find the actual freq for each bin based on current window size
	for(i=0; i<window_half; i++)
		bin_freqs[i] = ((x->sr/x->overlap/x->window) * i);
	
	// finally, build the filterbank
	for(i=1; i<=x->num_filters; i++)
	{
		int start_index, peak_index, finish_index, filter_width;
		start_index = peak_index = finish_index = 0;
		
		start_index = mfcc_tilde_nearest_bin_index(mel_spaced_freqs[i-1], bin_freqs);
		peak_index = mfcc_tilde_nearest_bin_index(mel_spaced_freqs[i], bin_freqs);
		finish_index = mfcc_tilde_nearest_bin_index(mel_spaced_freqs[i+1], bin_freqs);

		// grab memory for this filter
		filter_width = finish_index - start_index + 1;
		x->x_filterbank[i-1].filter = (float *)t_resizebytes(x->x_filterbank[i-1].filter, 0, filter_width * sizeof(float));
		x->x_filterbank_sizes[i-1] = filter_width; // store the sizes for freeing memory later

		// initialize this filter
		for(j=0; j<filter_width; j++)
			x->x_filterbank[i-1].filter[j] = 0.0;
			
		// PROTECT AGAINST indices < 0 or > Nyquist
		if( start_index > x->window/2.0 )
			start_index = x->window/2.0;

		if( peak_index > x->window/2.0 )
			peak_index = x->window/2.0;

		if( finish_index > x->window/2.0 )
			finish_index = x->window/2.0;

		// SPECIAL CASE FOR DUPLICATE START/PEAK/FINISH
		if( (finish_index-start_index)==0 )
		{
			x->x_filterbank[i-1].filter[0] = mel_spaced_freqs[i]/bin_freqs[peak_index];
			if(	x->x_filterbank[i-1].filter[0]>1 )
				x->x_filterbank[i-1].filter[0]=1;
		}
		else if( (finish_index - start_index)==1 )
		{
			x->x_filterbank[i-1].filter[0] = mel_spaced_freqs[i]/bin_freqs[peak_index];
			if(	x->x_filterbank[i-1].filter[0]>1 )
				x->x_filterbank[i-1].filter[0]=1;

			x->x_filterbank[i-1].filter[1] = mel_spaced_freqs[i+1]/bin_freqs[finish_index];
			if(	x->x_filterbank[i-1].filter[1]>1 )
				x->x_filterbank[i-1].filter[1]=1;
		}
		else
		{	
			// UPWARD
			for(j=start_index, k=0; j<=peak_index; j++, k++)
			{
				// (bin_freqs(j)-start)/(peak-start);
				x->x_filterbank[i-1].filter[k] = (bin_freqs[j] - mel_spaced_freqs[i-1])/(mel_spaced_freqs[i] - mel_spaced_freqs[i-1]);

				// all references to x_filterbank or x_indices must use [i-1] since we start at 1 and end at num_filters
			
				if(   x->x_filterbank[i-1].filter[k] < 0   )
					x->x_filterbank[i-1].filter[k] = 0;

				if(   x->x_filterbank[i-1].filter[k] > 1   )
					x->x_filterbank[i-1].filter[k] = 1;
			};


			// DOWNWARD (k continues where it is...)
			for(j=(peak_index+1); j<finish_index; j++, k++)
			{
				x->x_filterbank[i-1].filter[k] = (mel_spaced_freqs[i+1] - bin_freqs[j])/(mel_spaced_freqs[i+1] - mel_spaced_freqs[i]);
		
				if(   x->x_filterbank[i-1].filter[k] < 0   )
					x->x_filterbank[i-1].filter[k] = 0;

				if(   x->x_filterbank[i-1].filter[k] > 1   )
					x->x_filterbank[i-1].filter[k] = 1;
			};
		};

		x->x_indices[i-1].index[0] = start_index;
		x->x_indices[i-1].index[1] = finish_index;
	};

	// free local memory
	t_freebytes(bin_freqs, (int)x->window * sizeof(float));

}

/* ---------------- END filterbank functions ---------------------- */




/* ---------------- dsp utility functions ---------------------- */

static void mfcc_tilde_hann(int n, t_sample *in1, float *hann)
{	
	while (n--)
    {
    	*in1 = *in1 * *hann;

    	in1++;
    	hann++;
    };
}

static void mfcc_tilde_realfft_unpack(int n, int n_half, t_sample *input, t_sample *imag)
{
	int i, j;
		
	imag[0]=input[0];  // DC
	
	for(i=(n-1), j=1; i>n_half; i--, j++)
		imag[j] = input[i];
}

static void mfcc_tilde_abs(int n, t_sample *in1, t_sample *in2)
{
	while (n--)
    {   
        *in1 = sqrt( (*in1 * *in1) + (*in2 * *in2) );
        in1++;
        in2++;
    };
}

static void mfcc_tilde_filterbank_multiply(t_sample *in1, int normalize, t_filter *x_filterbank, t_indices *x_indices, int num_filters)
{
	int i, j, k;
	float sum, sumsum, *filter_power;

	filter_power = (float *)getbytes(0);
	filter_power = (float *)t_resizebytes(filter_power, 0, num_filters * sizeof(float));

 	sumsum = 0;
	
	for(i=0; i<num_filters; i++)
    {
	   	sum = 0;
 
		for(j=x_indices[i].index[0], k=0; j< (int)(x_indices[i].index[1] + 1); j++, k++)	
	    	sum += in1[j] * x_filterbank[i].filter[k];
		
		sum /= k;
		filter_power[i] = sum;  // get the total power.  another weighting might be better.

 		sumsum += sum;  // normalize so power in all bands sums to 1
	};
	
	if(normalize)
	{
		// prevent divide by 0
		if(sumsum==0)
			sumsum=1;
		else
			sumsum = 1/sumsum; // take the reciprocal here to save a divide below
	}
	else
		sumsum=1;
		
	for(i=0; i<num_filters; i++)
		in1[i] = filter_power[i]*sumsum;

	t_freebytes(filter_power, num_filters * sizeof(float));
}

static void mfcc_tilde_compute_mfccs(float *mfcc, t_sample *in1, int num_filters)
{
	int i, k;
	float pi_over_nfilters;
	
	pi_over_nfilters = M_PI/num_filters; // save multiple divides below
	
	for(i=0; i<num_filters; i++)
    {
	   	mfcc[i] = 0;
 		
		for(k=0; k<num_filters; k++)
 	    	mfcc[i] += in1[k] * cos(i * (k+0.5) * pi_over_nfilters);  // DCT-II
	};
}

/* ---------------- END filterbank functions ---------------------- */




/* ------------------------ mfcc~ -------------------------------- */

static void mfcc_tilde_bang(t_mfcc_tilde *x)
{
	int i, j, window, window_half, bang_sample;
	t_atom *listOut;
	t_sample *signal_R, *signal_I;
	double current_time;

	window = x->window;
	window_half = window * 0.5;
	
	// create local memory
	listOut = (t_atom *)getbytes(0);
	signal_R = (t_sample *)getbytes(0);
	signal_I = (t_sample *)getbytes(0);
	listOut = (t_atom *)t_resizebytes(listOut, 0, x->num_filters * sizeof(t_atom));
	signal_R = (t_sample *)t_resizebytes(signal_R, 0, window*sizeof(t_sample));
	signal_I = (t_sample *)t_resizebytes(signal_I, 0, window_half*sizeof(t_sample));
	
	current_time = clock_gettimesince(x->last_dsp_time);
	bang_sample = (int)(((current_time/1000.0)*x->sr)+0.5); // round

	if (bang_sample < 0)
        bang_sample = 0;
    else if ( bang_sample >= x->n )
        bang_sample = x->n - 1;
            
	// construct analysis window using bang_sample as the end of the window
	for(i=0, j=bang_sample; i<window; i++, j++)
		signal_R[i] = x->signal_R[j];
	
	if(x->window_function==1)
		mfcc_tilde_hann(window, signal_R, x->hann);
	
	mayer_realfft(window, signal_R);
	mfcc_tilde_realfft_unpack(window, window_half, signal_R, signal_I);
	mfcc_tilde_abs(window_half, signal_R, signal_I);

	mfcc_tilde_filterbank_multiply(signal_R, x->normalize, x->x_filterbank, x->x_indices, x->num_filters);
	mfcc_tilde_compute_mfccs(x->mfcc, signal_R, x->num_filters);
		
		
	for(i=0; i<x->num_filters; i++)
		SETFLOAT(listOut+i, x->mfcc[i]);

	outlet_list(x->x_featureList, 0, x->num_filters, listOut);

	// free local memory
	t_freebytes(listOut, x->num_filters * sizeof(t_atom));
	t_freebytes(signal_R, window * sizeof(t_sample));
	t_freebytes(signal_I, window_half * sizeof(t_sample));
}


static void mfcc_tilde_create_filterbank(t_mfcc_tilde *x, t_floatarg ms)
{
	int i, j, old_num_filters, size_msf;
	float mel_sum, *mel_spaced_freqs;

	// check early for appropriate size of mel_spaced_freqs
	mel_sum = ms;
	for(i=1; mfcc_tilde_mel2freq(mel_sum)<((x->sr/x->overlap)*0.5); i++)
		mel_sum += ms;

	size_msf = i;

	// create local memory
	mel_spaced_freqs = (float *)getbytes(0);
	mel_spaced_freqs = (float *)t_resizebytes(mel_spaced_freqs, 0, size_msf*sizeof(float));
	
	old_num_filters = x->num_filters;
	
	if(ms >= 5.0 && ms <= 1000.0)
	{
	 	x->mel_spacing = ms;
	
		// find i+1 mel-spaced frequencies up to but not including Nyquist
		mel_sum = x->mel_spacing;
		mel_spaced_freqs[0] = 0;
		for(i=1; i<size_msf; i++)
		{
			mel_spaced_freqs[i] = mfcc_tilde_mel2freq(mel_sum);
			mel_sum += (int)x->mel_spacing;
		};
			
		// i-2 is the correct number of filters
		x->num_filters = i-2;
		
		x->mfcc = (float *)t_resizebytes(x->mfcc, old_num_filters * sizeof(float), x->num_filters * sizeof(float));

		// initialize mfcc
		for(i=0; i<x->num_filters; i++)
				x->mfcc[i] = 0.0;
				
	    // free memory for each filter
		for(i=0; i<old_num_filters; i++)
			t_freebytes(x->x_filterbank[i].filter, x->x_filterbank_sizes[i]*sizeof(float));

		x->x_filterbank = (t_filter *)t_resizebytes(x->x_filterbank, old_num_filters * sizeof(t_filter), x->num_filters * sizeof(t_filter));

		// getbytes for each new filter
		for(i=0; i<x->num_filters; i++)
			x->x_filterbank[i].filter = (float *)getbytes(0);

		x->x_indices = (t_indices *)t_resizebytes(x->x_indices, old_num_filters * sizeof(t_indices), x->num_filters * sizeof(t_indices));

		// initialize indices
		for(i=0; i<x->num_filters; i++)
			for(j=0; j<2; j++)
				x->x_indices[i].index[j] = 0;

		x->x_filterbank_sizes = (int *)t_resizebytes(x->x_filterbank_sizes, old_num_filters*sizeof(int), x->num_filters*sizeof(int));

		// initialize filterbank sizes
		for(i=0; i<x->num_filters; i++)
			x->x_filterbank_sizes[i] = 0;
				
		// do the dirty work in a separate function.
		mfcc_tilde_make_filterbank(x, mel_spaced_freqs);
			
		post("no. of filters: %i", x->num_filters);
		post("mel spacing: %f", x->mel_spacing);
    }
    else
		error("mel spacing must be between 5 and 1000 mels");

	// free local memory
	t_freebytes(mel_spaced_freqs, size_msf*sizeof(float));
}


static void mfcc_tilde_hat(t_mfcc_tilde *x, t_floatarg f, t_floatarg g)
{
	int i;
	
	for(i=0; i<(int)g; i++)
		post("val %i: %f", i, x->x_filterbank[(int)f].filter[i]);

}


static void mfcc_tilde_window(t_mfcc_tilde *x, t_floatarg f)
{
	int i, isPow2;
	
	isPow2 = (int)f && !( ((int)f-1) & (int)f );
	
	if( !isPow2 )
		error("requested window size is not a power of 2");
	else
	{
		x->signal_R = (t_sample *)t_resizebytes(x->signal_R, (x->window+x->n) * sizeof(t_sample), (f+x->n) * sizeof(t_sample));
		x->hann = (float *)t_resizebytes(x->hann, x->window * sizeof(float), f * sizeof(float));
	
		x->window = (int)f;
		
		for(i=0; i<x->window; i++)
			x->hann[i] = 0.5 * (1 - cos(2*M_PI*i/x->window));

		// init signal buffer
		for(i=0; i<(x->window+x->n); i++)
			x->signal_R[i] = 0.0;
			
		mfcc_tilde_create_filterbank(x, x->mel_spacing);
		
		post("window size: %i. overlap: %i. sampling rate: %i", (int)x->window, x->overlap, (int)(x->sr/x->overlap));
	}
}


static void mfcc_tilde_overlap(t_mfcc_tilde *x, t_floatarg f)
{
	if((int)f > 0)
		x->overlap = (int)f;
	else
		error("overlap must be at least 1");

    post("overlap: %i", x->overlap);
	mfcc_tilde_create_filterbank(x, x->mel_spacing);

}


// magnitude spectrum == 0, power spectrum == 1
static void mfcc_tilde_power_spectrum(t_mfcc_tilde *x, t_floatarg spec)
{
	x->power_spectrum = spec;
	
	if(x->power_spectrum)
		post("using power spectrum for MFCC computation.");
	else
		post("using magnitude spectrum for MFCC computation.");
}


static void mfcc_tilde_normalize(t_mfcc_tilde *x, t_floatarg norm)
{
	if(norm<0)
		x->normalize = 0;
	else if (norm>1)
		x->normalize = 1;
	else
		x->normalize = norm;
	
	if(x->normalize)
		post("spectrum normalization ON.");
	else
		post("spectrum normalization OFF.");
}


// this can be used to later to actually try different windows.
// Blackman and Hamming would be a start
static void mfcc_tilde_window_function(t_mfcc_tilde *x, t_floatarg f)
{
	x->window_function = (int)f;
	
	if((int)f==0)
		post("window function: rectangular");
	else
		post("window function: hann");
}


static void mfcc_tilde_debug(t_mfcc_tilde *x, t_floatarg f)
{
	if((int)f < 0)
		x->debug = 0;
	else if((int)f > 2)
		x->debug = 2;
	else
		x->debug = (int)f;

	post("debug mode: %i", (int)f);
}


static void *mfcc_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
    t_mfcc_tilde *x = (t_mfcc_tilde *)pd_new(mfcc_tilde_class);
	int i, isPow2;
	
	x->x_featureList = outlet_new(&x->x_obj, gensym("list"));
	s=s;
	
	x->mfcc = (float *)getbytes(0);
	x->signal_R = (t_sample *)getbytes(0);
	x->hann = (float *)getbytes(0);
	

	if(argc > 1)
	{
		x->window = atom_getfloat(argv);  // should perform a check for >64 && power of two
		isPow2 = (int)x->window && !( ((int)x->window-1) & (int)x->window );
		
		if(!isPow2)
		{
			error("requested window size is not a power of 2. default value of 1024 used instead");
			x->window = 1024;
		};

		x->mel_spacing = atom_getfloat(argv+1);  // should check for > 0

		if(x->mel_spacing <= 0)
		{
			error("mel spacing must be a positive real number. default value of 100 used instead");
			x->mel_spacing = 100;
		};
	}
	else if(argc > 0)
	{
		x->window = atom_getfloat(argv);
		isPow2 = (int)x->window && !( ((int)x->window-1) & (int)x->window );
		
		if(!isPow2)
		{
			error("requested window size is not a power of 2. default value of 1024 used instead");
			x->window = 1024;
		};
		
		x->mel_spacing = 100;
	}
	else
	{
		x->window = 1024;	
		x->mel_spacing = 100;
	}

	x->sr = 44100.0;
	x->n = 64.0;
	x->window_function = 1; // 1 is hann window, 0 is nothing
	x->power_spectrum = 0; // choose mag (0) or power (1) spec in the BFCC computation
	x->normalize = 1; // this is generally a good thing, but should be off for concatenative synth
	x->last_dsp_time = clock_getlogicaltime();
	x->overlap = 1;
	x->debug = 0;
    x->num_filters = 38;
    
	x->signal_R = (t_sample *)t_resizebytes(x->signal_R, 0, (x->window+x->n) * sizeof(t_sample));
	x->hann = (float *)t_resizebytes(x->hann, 0, x->window * sizeof(float));
	// filterbank and indices are resized in mfcc_tilde_create_filterbank

    // initialize signal buffer
    for(i=0; i<(x->window+x->n); i++)
			x->signal_R[i] = 0.0;

	// initialize hann window
 	for(i=0; i<x->window; i++)
 		x->hann[i] = 0.5 * (1 - cos(2*M_PI*i/x->window));

	x->x_filterbank_sizes = (int *)getbytes(0);
	x->x_filterbank_sizes = (int *)t_resizebytes(x->x_filterbank_sizes, 0, x->num_filters*sizeof(int));

    // initialize filterbank sizes
    for(i=0; i<x->num_filters; i++)
		x->x_filterbank_sizes[i] = 0;	

	// grab memory for filterbank
	x->x_filterbank = (t_filter *)getbytes(0);
	x->x_filterbank = (t_filter *)t_resizebytes(x->x_filterbank, 0, x->num_filters * sizeof(t_filter));

	// getbytes for each filter
	for(i=0; i<x->num_filters; i++)
		x->x_filterbank[i].filter = (float *)getbytes(0);

	// grab memory for indices
	x->x_indices = (t_indices *)getbytes(0);
	x->x_indices = (t_indices *)t_resizebytes(x->x_indices, 0, x->num_filters * sizeof(t_indices));
	
	mfcc_tilde_create_filterbank(x, x->mel_spacing);
	
    post("mfcc~: window size: %i", (int)x->window);
    
    return (x);
}


static t_int *mfcc_tilde_perform(t_int *w)
{
    int i, n;
	
    t_mfcc_tilde *x = (t_mfcc_tilde *)(w[1]);

    t_sample *in1 = (t_float *)(w[2]);
    n = w[3];
 			
 	// shift signal buffer contents back.
	for(i=0; i<x->window; i++)
		x->signal_R[i] = x->signal_R[i+n];
	
	// write new block to end of signal buffer.
	for(i=0; i<n; i++)
		x->signal_R[(int)x->window+i] = in1[i];
		
	x->last_dsp_time = clock_getlogicaltime();

    return (w+4);
}


static void mfcc_tilde_dsp(t_mfcc_tilde *x, t_signal **sp)
{
	int i;
	
	dsp_add(
		mfcc_tilde_perform,
		3,
		x,
		sp[0]->s_vec,
		sp[0]->s_n
	); 

// compare n to stored n and recalculate filterbank if different
	if( sp[0]->s_sr != x->sr || sp[0]->s_n != x->n )
	{
		x->signal_R = (t_sample *)t_resizebytes(x->signal_R, (x->window+x->n) * sizeof(t_sample), (x->window+sp[0]->s_n) * sizeof(t_sample));

		x->sr = sp[0]->s_sr;
		x->n = sp[0]->s_n;
		x->last_dsp_time = clock_getlogicaltime();

		// init signal buffer
		for(i=0; i<(x->window+x->n); i++)
			x->signal_R[i] = 0.0;
			
		mfcc_tilde_create_filterbank(x, x->mel_spacing);

    	post("mfcc~: window size: %i. overlap: %i. sampling rate: %i, block size: %i", (int)x->window, x->overlap, (int)(x->sr/x->overlap), (int)x->n);
	};

};


static void mfcc_tilde_free(t_mfcc_tilde *x)
{	
	int i;
	
	// free the input buffer memory
    t_freebytes(x->signal_R, (x->window+x->n)*sizeof(t_sample));

	// free the hann window memory
    t_freebytes(x->hann, x->num_filters*sizeof(float));

	// free the mfcc memory
    t_freebytes(x->mfcc, x->num_filters*sizeof(float));
       
    // free the filterbank memory
	for(i=0; i<x->num_filters; i++)
		t_freebytes(x->x_filterbank[i].filter, x->x_filterbank_sizes[i]*sizeof(float));
    
	t_freebytes(x->x_filterbank, x->num_filters*sizeof(t_filter));

	// free the indices memory
    t_freebytes(x->x_indices, x->num_filters*sizeof(t_indices));
}


void mfcc_tilde_setup(void)
{
    mfcc_tilde_class = 
    class_new(
    	gensym("mfcc~"),
    	(t_newmethod)mfcc_tilde_new,
    	(t_method)mfcc_tilde_free,
        sizeof(t_mfcc_tilde),
        CLASS_DEFAULT, 
        A_GIMME,
		0
    );

    CLASS_MAINSIGNALIN(mfcc_tilde_class, t_mfcc_tilde, x_f);

	class_addbang(mfcc_tilde_class, mfcc_tilde_bang);

	class_addmethod(
		mfcc_tilde_class, 
        (t_method)mfcc_tilde_create_filterbank,
		gensym("filterbank"),
        A_DEFFLOAT,
        A_DEFFLOAT,
		0
	);

	class_addmethod(
		mfcc_tilde_class, 
        (t_method)mfcc_tilde_hat,
		gensym("hat"),
        A_DEFFLOAT,
        A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		mfcc_tilde_class, 
        (t_method)mfcc_tilde_window,
		gensym("window"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		mfcc_tilde_class, 
        (t_method)mfcc_tilde_overlap,
		gensym("overlap"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		mfcc_tilde_class, 
        (t_method)mfcc_tilde_window_function,
		gensym("window_function"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		mfcc_tilde_class, 
        (t_method)mfcc_tilde_power_spectrum,
		gensym("power_spectrum"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		mfcc_tilde_class, 
        (t_method)mfcc_tilde_normalize, 
		gensym("normalize"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		mfcc_tilde_class, 
        (t_method)mfcc_tilde_debug,
		gensym("debug"),
		A_DEFFLOAT,
		0
	);
	
    class_addmethod(
    	mfcc_tilde_class,
    	(t_method)mfcc_tilde_dsp,
    	gensym("dsp"),
    	0
    );
}