/*

barkSpec~ - A Bark Frequency Spectrum Analysis external.

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.


version 0.2.5, May 16, 2010

¥ 0.2.5 changed filterbank_multiply so that sum of power in each filter is divided by the number of points in the filter.
¥Ê0.2.4 added an ifndef M_PI for guaranteed windows compilation
¥ 0.2.3 adds a #define M_PI for windows compilation, and declares all functions except _setup static
¥ 0.2.2 is part of the update that ensures all function names are prepended by the external name (bfcc_ or bfcc_tilde_, etc).
¥ 0.2.1 uses dynamic memory allocation for filter widths
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

static t_class *barkSpec_tilde_class;

typedef struct filter
{
    float *filter;
} t_filter;

typedef struct indices
{
    int index[2];
} t_indices;

typedef struct _barkSpec_tilde
{
    t_object x_obj;
    float sr;
    float n;
	float window;
	int power_spectrum;
	int normalize;
	int window_function;
	int debug;
    int overlap;
    double last_dsp_time;
    int num_filters;
	int *x_filterbank_sizes;
    float bark_spacing;
	t_filter *x_filterbank;
	t_indices *x_indices;
	t_sample *signal_R;
    float *blackman;
    float *cosine;
    float *hamming;
    float *hann;
    float x_f;
    t_outlet *x_featureList;

} t_barkSpec_tilde;




/* ---------------- utility functions for building filterbank ---------------------- */

static int barkSpec_tilde_nearest_bin_index(float target, float *bin_freqs)
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


static float barkSpec_tilde_bark2freq(float bark)
{
	float freq;
	
	freq = 1960/(26.81/(bark+0.53) - 1);
	if(freq < 0)
		freq = 0;
		
	return(freq);
}


static void barkSpec_tilde_make_filterbank(t_barkSpec_tilde *x, float *bark_spaced_freqs)
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
		
		start_index = barkSpec_tilde_nearest_bin_index(bark_spaced_freqs[i-1], bin_freqs);
		peak_index = barkSpec_tilde_nearest_bin_index(bark_spaced_freqs[i], bin_freqs);
		finish_index = barkSpec_tilde_nearest_bin_index(bark_spaced_freqs[i+1], bin_freqs);

		// grab memory for this filter
		filter_width = finish_index - start_index + 1;
		x->x_filterbank[i-1].filter = (float *)t_resizebytes(x->x_filterbank[i-1].filter, 0, filter_width * sizeof(float));
		x->x_filterbank_sizes[i-1] = filter_width; // store the sizes for freeing memory later

		// initialize this filter
		for(j=0; j<filter_width; j++)
			x->x_filterbank[i-1].filter[j] = 0.0;

		// PROTECT AGAINST indices < 0 or > Nyquist
		if( start_index > x->window * 0.5 )
			start_index = x->window * 0.5;

		if( peak_index > x->window * 0.5 )
			peak_index = x->window * 0.5;

		if( finish_index > x->window * 0.5 )
			finish_index = x->window * 0.5;

		// SPECIAL CASE FOR DUPLICATE START/PEAK/FINISH
		if( (finish_index-start_index)==0 )
		{
			x->x_filterbank[i-1].filter[0] = bark_spaced_freqs[i]/bin_freqs[peak_index];
			if(	x->x_filterbank[i-1].filter[0]>1 )
				x->x_filterbank[i-1].filter[0]=1;
		}
		else if( (finish_index - start_index)==1 )
		{
			x->x_filterbank[i-1].filter[0] = bark_spaced_freqs[i]/bin_freqs[peak_index];
			if(	x->x_filterbank[i-1].filter[0]>1 )
				x->x_filterbank[i-1].filter[0]=1;

			x->x_filterbank[i-1].filter[1] = bark_spaced_freqs[i+1]/bin_freqs[finish_index];
			if(	x->x_filterbank[i-1].filter[1]>1 )
				x->x_filterbank[i-1].filter[1]=1;
		}
		else
		{	
			// UPWARD
			for(j=start_index, k=0; j<=peak_index; j++, k++)
			{
				// (bin_freqs(j)-start)/(peak-start);
				x->x_filterbank[i-1].filter[k] = (bin_freqs[j] - bark_spaced_freqs[i-1])/(bark_spaced_freqs[i] - bark_spaced_freqs[i-1]);

				// all references to x_filterbank or x_indices must use [i-1] since we start at 1 and end at num_filters
			
				if(   x->x_filterbank[i-1].filter[k] < 0   )
					x->x_filterbank[i-1].filter[k] = 0;

				if(   x->x_filterbank[i-1].filter[k] > 1   )
					x->x_filterbank[i-1].filter[k] = 1;
			};


			// DOWNWARD (k continues where it is...)
			for(j=(peak_index+1); j<finish_index; j++, k++)
			{
				x->x_filterbank[i-1].filter[k] = (bark_spaced_freqs[i+1] - bin_freqs[j])/(bark_spaced_freqs[i+1] - bark_spaced_freqs[i]);
		
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

static void barkSpec_tilde_DSP_window(int n, t_sample *in, float *window_func)
{	
	while (n--)
    {
    	*in *= *window_func;

    	in++;
    	window_func++;
    };
}

static void barkSpec_tilde_realfft_unpack(int n, int n_half, t_sample *input, t_sample *imag)
{
	int i, j;
		
	imag[0]=0.0;  // DC
	
	for(i=(n-1), j=1; i>n_half; i--, j++)
		imag[j] = input[i];
}

static void barkSpec_tilde_abs(int n, t_sample *in1, t_sample *in2)
{
	while (n--)
    {   
        *in1 = sqrt( (*in1 * *in1) + (*in2 * *in2) );
        in1++;
        in2++;
    };
}

static void barkSpec_tilde_square(int n, t_sample *in)
{
	while (n--)
    {   
	    *in *= *in;	        
        in++;
    };
}

static void barkSpec_tilde_filterbank_multiply(t_sample *in, int normalize, t_filter *x_filterbank, t_indices *x_indices, int num_filters)
{
	int i, j, k;
	float sum, sumsum, *filter_power;

	// create local memory
	filter_power = (float *)getbytes(0);
	filter_power = (float *)t_resizebytes(filter_power, 0, num_filters * sizeof(float));

 	sumsum = 0;
	
	for(i=0; i<num_filters; i++)
    {
	   	sum = 0;
 
		for(j=x_indices[i].index[0], k=0; j< (int)(x_indices[i].index[1] + 1); j++, k++)	
	    	sum += in[j] * x_filterbank[i].filter[k];
		
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
		in[i] = filter_power[i]*sumsum;

	// free local memory
	t_freebytes(filter_power, num_filters * sizeof(float));
}

/* ---------------- END utility functions ---------------------- */




/* ------------------------ barkSpec~ -------------------------------- */

static void barkSpec_tilde_bang(t_barkSpec_tilde *x)
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
	
	switch(x->window_function)
	{
		case 0:
			barkSpec_tilde_DSP_window(window, signal_R, x->blackman);
			break;
		case 1:
			barkSpec_tilde_DSP_window(window, signal_R, x->cosine);
			break;
		case 2:
			barkSpec_tilde_DSP_window(window, signal_R, x->hamming);
			break;
		case 3:
			barkSpec_tilde_DSP_window(window, signal_R, x->hann);
			break;
		default:
			break;
	};
	
	mayer_realfft(window, signal_R);
	barkSpec_tilde_realfft_unpack(window, window_half, signal_R, signal_I);
	barkSpec_tilde_abs(window_half, signal_R, signal_I);

// power spectrum sometimes generates lower scores than magnitude. make it optional.
	if(x->power_spectrum)
		barkSpec_tilde_square(window_half, signal_R);

	barkSpec_tilde_filterbank_multiply(signal_R, x->normalize, x->x_filterbank, x->x_indices, x->num_filters);
		
	for(i=0; i<x->num_filters; i++)
		SETFLOAT(listOut+i, signal_R[i]);

	outlet_list(x->x_featureList, 0, x->num_filters, listOut);

	// free local memory
	t_freebytes(listOut, x->num_filters * sizeof(t_atom));
	t_freebytes(signal_R, window * sizeof(t_sample));
	t_freebytes(signal_I, window_half * sizeof(t_sample));
}


static void barkSpec_tilde_create_filterbank(t_barkSpec_tilde *x, t_floatarg bs, t_floatarg suppress_post)
{
	int i, j, old_num_filters, size_bsf;
	float bark_sum, *bark_spaced_freqs;
		
	// check early for appropriate size of mel_spaced_freqs
	bark_sum = bs;
	for(i=1; bark_sum < (24.0 + bs); i++)
		bark_sum += bs;

	size_bsf = i;

	// create local memory
	bark_spaced_freqs = (float *)getbytes(0);
	bark_spaced_freqs = (float *)t_resizebytes(bark_spaced_freqs, 0, size_bsf * sizeof(float));
	
	old_num_filters = x->num_filters;
	
	if(bs >= 0.02 && bs <= 10.0)
	{
	 	x->bark_spacing = bs;

		// find i+1 bark-spaced frequencies up to 24 barks
		bark_sum = x->bark_spacing;
		bark_spaced_freqs[0] = 0;
		for(i=1; i<size_bsf; i++)
		{
			bark_spaced_freqs[i] = barkSpec_tilde_bark2freq(bark_sum);
			bark_sum += x->bark_spacing;
		};
				
		// i-2 is the correct number of filters
		x->num_filters = i-2;

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
		barkSpec_tilde_make_filterbank(x, bark_spaced_freqs);

		if(!suppress_post)
		{
			post("bark spacing: %f", x->bark_spacing);
			post("no. of filters: %i", x->num_filters);
		}
    }
    else
		error("Bark spacing must be between 0.02 and 10.0 Barks.");

	// free local memory
	t_freebytes(bark_spaced_freqs, size_bsf * sizeof(float));
}


static void barkSpec_tilde_hat(t_barkSpec_tilde *x, t_floatarg f, t_floatarg g)
{
	int i;
	
	for(i=0; i<(int)g; i++)
		post("val %i: %f", i, x->x_filterbank[(int)f].filter[i]);

}


static void barkSpec_tilde_window(t_barkSpec_tilde *x, t_floatarg f)
{
	int i, isPow2;
	
	isPow2 = (int)f && !( ((int)f-1) & (int)f );
	
	if( !isPow2 )
		error("requested window size is not a power of 2.");
	else
	{
		x->signal_R = (t_sample *)t_resizebytes(x->signal_R, (x->window+x->n) * sizeof(t_sample), (f+x->n) * sizeof(t_sample));
		x->blackman = (float *)t_resizebytes(x->blackman, x->window * sizeof(float), f * sizeof(float));
		x->cosine = (float *)t_resizebytes(x->cosine, x->window * sizeof(float), f * sizeof(float));
		x->hamming = (float *)t_resizebytes(x->hamming, x->window * sizeof(float), f * sizeof(float));
		x->hann = (float *)t_resizebytes(x->hann, x->window * sizeof(float), f * sizeof(float));
	
		x->window = (int)f;

		// re-init window functions
		for(i=0; i<x->window; i++)
		{
			x->blackman[i] = 0.42 - (0.5 * cos(2*M_PI*i/x->window)) + (0.08 * cos(4*M_PI*i/x->window));
			x->cosine[i] = sin(M_PI*i/x->window);
			x->hamming[i] = 0.5 - (0.46 * cos(2*M_PI*i/x->window));
			x->hann[i] = 0.5 * (1 - cos(2*M_PI*i/x->window));
		};	

 		
		// init signal buffer
		for(i=0; i<(x->window+x->n); i++)
			x->signal_R[i] = 0.0;
			
		barkSpec_tilde_create_filterbank(x, x->bark_spacing, 0);
		
		post("window size: %i. overlap: %i. sampling rate: %i", (int)x->window, x->overlap, (int)(x->sr/x->overlap));
	}
}


static void barkSpec_tilde_overlap(t_barkSpec_tilde *x, t_floatarg f)
{
	if((int)f > 0)
		x->overlap = (int)f;
	else
		error("overlap must be at least 1.");

    post("overlap: %i", x->overlap);
	barkSpec_tilde_create_filterbank(x, x->bark_spacing, 0);

}


// this can be used to later to actually try different windows.
// Blackman and Hamming would be a start
static void barkSpec_tilde_window_function(t_barkSpec_tilde *x, t_floatarg f)
{
	x->window_function = (int)f;
	
	switch((int)f)
	{
		case 0:
			post("window function: blackman.");
			break;
		case 1:
			post("window function: cosine.");
			break;
		case 2:
			post("window function: hamming.");
			break;
		case 3:
			post("window function: hann.");
			break;
		default:
			break;
	};
}


// magnitude spectrum == 0, power spectrum == 1
static void barkSpec_tilde_power_spectrum(t_barkSpec_tilde *x, t_floatarg spec)
{
	x->power_spectrum = spec;
	
	if(x->power_spectrum)
		post("using power spectrum for BFCC computation.");
	else
		post("using magnitude spectrum for BFCC computation.");
}


static void barkSpec_tilde_normalize(t_barkSpec_tilde *x, t_floatarg norm)
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


static void barkSpec_tilde_debug(t_barkSpec_tilde *x, t_floatarg f)
{
	if((int)f < 0)
		x->debug = 0;
	else if((int)f > 2)
		x->debug = 2;
	else
		x->debug = (int)f;

	post("debug mode: %i", (int)f);
}


static void *barkSpec_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
    t_barkSpec_tilde *x = (t_barkSpec_tilde *)pd_new(barkSpec_tilde_class);
	int i, isPow2;
	s=s;
	
	x->x_featureList = outlet_new(&x->x_obj, gensym("list"));

	if(argc > 1)
	{
		x->window = atom_getfloat(argv);  // should perform a check for >64 && power of two
		isPow2 = (int)x->window && !( ((int)x->window-1) & (int)x->window );
		
		if(!isPow2)
		{
			error("requested window size is not a power of 2. default value of 1024 used instead.");
			x->window = 1024;
		};

		x->bark_spacing = atom_getfloat(argv+1);  // should check for > 0

		if(x->bark_spacing <= 0)
		{
			error("bark spacing must be a positive real number. default value of 0.5 used instead.");
			x->bark_spacing = 0.5;
		};
	}
	else if(argc > 0)
	{
		x->window = atom_getfloat(argv);
		isPow2 = (int)x->window && !( ((int)x->window-1) & (int)x->window );
		
		if(!isPow2)
		{
			error("requested window size is not a power of 2. default value of 1024 used instead.");
			x->window = 1024;
		};
		
		x->bark_spacing = 0.5;
	}
	else
	{
		x->window = 1024;	
		x->bark_spacing = 0.5;
	}

	x->sr = 44100.0;
	x->n = 64.0;
	x->window_function = 3; // 3 is hann window
	x->power_spectrum = 0; // choose mag (0) or power (1) spec in the BFCC computation
	x->normalize = 1; // this is generally a good thing, but should be off for concatenative synth
	x->last_dsp_time = clock_getlogicaltime();
	x->overlap = 1;
	x->debug = 0;
    x->num_filters = 47;
	
	x->signal_R = (t_sample *)getbytes(0);
	x->signal_R = (t_sample *)t_resizebytes(x->signal_R, 0, (x->window+x->n) * sizeof(t_sample));
 
     // initialize signal buffer
    for(i=0; i<(x->window+x->n); i++)
		x->signal_R[i] = 0.0;
		
  	x->blackman = (float *)getbytes(0);
	x->blackman = (float *)t_resizebytes(x->blackman, 0, x->window * sizeof(float));
 
 	x->cosine = (float *)getbytes(0);
	x->cosine = (float *)t_resizebytes(x->cosine, 0, x->window * sizeof(float));
 	
 	x->hamming = (float *)getbytes(0);
 	x->hamming = (float *)t_resizebytes(x->hamming, 0, x->window * sizeof(float));
	
	x->hann = (float *)getbytes(0);
	x->hann = (float *)t_resizebytes(x->hann, 0, x->window * sizeof(float));

 	// initialize blackman window
 	// alpha = 0.16. a0 = 0.42; a1 = 0.5; a2 = 0.08;
  	// initialize cosine window
 	// initialize hamming window
	// initialize hann window
 	for(i=0; i<x->window; i++)
 	{
 		x->blackman[i] = 0.42 - (0.5 * cos(2*M_PI*i/x->window)) + (0.08 * cos(4*M_PI*i/x->window));
  		x->cosine[i] = sin(M_PI*i/x->window);
  		x->hamming[i] = 0.5 - (0.46 * cos(2*M_PI*i/x->window));
 		x->hann[i] = 0.5 * (1 - cos(2*M_PI*i/x->window));
 	};		
		
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
		
	barkSpec_tilde_create_filterbank(x, x->bark_spacing, 0);
	
    post("barkSpec~: window size: %i.", (int)x->window);
    
    return (x);
}


static t_int *barkSpec_tilde_perform(t_int *w)
{
    int i, n;
	
    t_barkSpec_tilde *x = (t_barkSpec_tilde *)(w[1]);

    t_sample *in = (t_float *)(w[2]);
    n = w[3];
 			
 	// shift signal buffer contents back.
	for(i=0; i<(x->window-n); i++)
		x->signal_R[i] = x->signal_R[i+n];
	
	// write new block to end of signal buffer.
	for(i=0; i<n; i++)
		x->signal_R[(int)x->window-n+i] = in[i];
		
	x->last_dsp_time = clock_getlogicaltime();

    return (w+4);
}


static void barkSpec_tilde_dsp(t_barkSpec_tilde *x, t_signal **sp)
{
	int i;
	
	dsp_add(
		barkSpec_tilde_perform,
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

		// init signal buffer
		for(i=0; i<(x->window+x->n); i++)
			x->signal_R[i] = 0.0;
			
		barkSpec_tilde_create_filterbank(x, x->bark_spacing, 0);

    	post("barkSpec~: window size: %i. overlap: %i. sampling rate: %i, block size: %i", (int)x->window, x->overlap, (int)(x->sr/x->overlap), (int)x->n);
	};

};


static void barkSpec_tilde_free(t_barkSpec_tilde *x)
{	
	int i;
	
	// free the input buffer memory
    t_freebytes(x->signal_R, (x->window+x->n)*sizeof(t_sample));

	// free the blackman window memory
    t_freebytes(x->blackman, x->window*sizeof(float));

	// free the cosine window memory
    t_freebytes(x->cosine, x->window*sizeof(float));

	// free the hamming window memory
    t_freebytes(x->hamming, x->window*sizeof(float));

	// free the hann window memory
    t_freebytes(x->hann, x->window*sizeof(float));
       
    // free the filterbank memory
	for(i=0; i<x->num_filters; i++)
		t_freebytes(x->x_filterbank[i].filter, x->x_filterbank_sizes[i]*sizeof(float));
    
	t_freebytes(x->x_filterbank, x->num_filters*sizeof(t_filter));

	// free the indices memory
    t_freebytes(x->x_indices, x->num_filters*sizeof(t_indices));
}


void barkSpec_tilde_setup(void)
{
    barkSpec_tilde_class = 
    class_new(
    	gensym("barkSpec~"),
    	(t_newmethod)barkSpec_tilde_new,
    	(t_method)barkSpec_tilde_free,
        sizeof(t_barkSpec_tilde),
        CLASS_DEFAULT, 
        A_GIMME,
		0
    );

    CLASS_MAINSIGNALIN(barkSpec_tilde_class, t_barkSpec_tilde, x_f);

	class_addbang(barkSpec_tilde_class, barkSpec_tilde_bang);

	class_addmethod(
		barkSpec_tilde_class, 
        (t_method)barkSpec_tilde_create_filterbank,
		gensym("filterbank"),
        A_DEFFLOAT,
        A_DEFFLOAT,
		0
	);

	class_addmethod(
		barkSpec_tilde_class, 
        (t_method)barkSpec_tilde_hat,
		gensym("hat"),
        A_DEFFLOAT,
        A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		barkSpec_tilde_class, 
        (t_method)barkSpec_tilde_window,
		gensym("window"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		barkSpec_tilde_class, 
        (t_method)barkSpec_tilde_overlap,
		gensym("overlap"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		barkSpec_tilde_class, 
        (t_method)barkSpec_tilde_window_function,
		gensym("window_function"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		barkSpec_tilde_class, 
        (t_method)barkSpec_tilde_power_spectrum,
		gensym("power_spectrum"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		barkSpec_tilde_class, 
        (t_method)barkSpec_tilde_normalize, 
		gensym("normalize"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		barkSpec_tilde_class, 
        (t_method)barkSpec_tilde_debug,
		gensym("debug"),
		A_DEFFLOAT,
		0
	);
	
    class_addmethod(
    	barkSpec_tilde_class,
    	(t_method)barkSpec_tilde_dsp,
    	gensym("dsp"),
    	0
    );
}