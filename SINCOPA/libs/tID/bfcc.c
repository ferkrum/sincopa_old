/*

bfcc - A non-real-time Bark Frequency Cepstral Analysis external.

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.


version 0.0.6, May 16, 2010

¥ 0.0.6 changed filterbank_multiply so that sum of power in each filter is divided by the number of points in the filter.
¥Ê0.0.5 added an ifndef M_PI for guaranteed windows compilation
¥ 0.0.4 adds a #define M_PI for windows compilation, and declares all functions except _setup static
¥ 0.0.3 is part of the update that ensures all function names are prepended by the external name (bfcc_ or bfcc_tilde_, etc).

*/

#include "m_pd.h"
#include <math.h>
#include <limits.h>
#include <stdio.h>
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

static t_class *bfcc_class;

typedef struct filter
{
    float *filter;
} t_filter;

typedef struct indices
{
    int index[2];
} t_indices;

typedef struct _bfcc
{
    t_object x_obj;
    float sr;
	float window;
	int power_spectrum;
	int normalize;
	int window_function;
	int debug;
	int max_window_size;
	int *powers_of_two;
    int  pow_two_arr_size;
    int num_filters;
	int *x_filterbank_sizes;
    float bark_spacing;
	t_filter *x_filterbank;
	t_indices *x_indices;
	t_sample *signal_R;
    float *bfcc;
	t_word *x_vec;
	t_symbol *x_arrayname;
	int x_array_points;
    t_outlet *x_featureList;

} t_bfcc;




/* ---------------- utility functions for building filterbank ---------------------- */

static int bfcc_nearest_bin_index(float target, float *bin_freqs)
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


static float bfcc_bark2freq(float bark)
{
	float freq;
	
	freq = 1960/(26.81/(bark+0.53) - 1);
	if(freq < 0)
		freq = 0;
		
	return(freq);
}


static void bfcc_make_filterbank(t_bfcc *x, float *bark_spaced_freqs)
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
		bin_freqs[i] = ((x->sr/x->window) * i);
	

	// finally, build the filterbank
	for(i=1; i<=x->num_filters; i++)
	{
		int start_index, peak_index, finish_index, filter_width;
		start_index = peak_index = finish_index = 0;
		
		start_index = bfcc_nearest_bin_index(bark_spaced_freqs[i-1], bin_freqs);
		peak_index = bfcc_nearest_bin_index(bark_spaced_freqs[i], bin_freqs);
		finish_index = bfcc_nearest_bin_index(bark_spaced_freqs[i+1], bin_freqs);

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

static void bfcc_blackman_window(int n, t_sample *in)
{
	int i;
	
	for(i=0; i<n; i++, in++)
    	*in *= 0.42 - (0.5 * cos(2*M_PI*i/n)) + (0.08 * cos(4*M_PI*i/n));
}


static void bfcc_cosine_window(int n, t_sample *in)
{
	int i;

	for(i=0; i<n; i++, in++)
    	*in *= sin(M_PI*i/n);
}

static void bfcc_hamming_window(int n, t_sample *in)
{	
	int i;

	for(i=0; i<n; i++, in++)
    	*in *= 0.5 - (0.46 * cos(2*M_PI*i/n));
}

static void bfcc_hann_window(int n, t_sample *in)
{
	int i;

	for(i=0; i<n; i++, in++)
    	*in *= 0.5 * (1 - cos(2*M_PI*i/n));
}
 		
static void bfcc_realfft_unpack(int n, int n_half, t_sample *input, t_sample *imag)
{
	int i, j;
		
	imag[0]=0.0;  // DC
	
	for(i=(n-1), j=1; i>n_half; i--, j++)
		imag[j] = input[i];
}

static void bfcc_abs(int n, t_sample *in1, t_sample *in2)
{
	while (n--)
    {   
        *in1 = sqrt( (*in1 * *in1) + (*in2 * *in2) );
        in1++;
        in2++;
    };
}

static void bfcc_square(int n, t_sample *in)
{
	while (n--)
    {   
	    *in *= *in;	        
        in++;
    };
}

static void bfcc_filterbank_multiply(t_sample *in, int normalize, t_filter *x_filterbank, t_indices *x_indices, int num_filters)
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

static void bfcc_compute_bfccs(float *bfcc, t_sample *in, int num_filters)
{
	int i, k;
	float pi_over_nfilters;
	
	pi_over_nfilters = M_PI/num_filters; // save multiple divides below
	
	for(i=0; i<num_filters; i++)
    {
	   	bfcc[i] = 0;
 		
		for(k=0; k<num_filters; k++)
 	    	bfcc[i] += in[k] * cos(i * (k+0.5) * pi_over_nfilters);  // DCT-II
	};
}

/* ---------------- END filterbank functions ---------------------- */




/* ------------------------ bfcc -------------------------------- */

static void bfcc_create_filterbank(t_bfcc *x, t_floatarg bs, t_floatarg suppress_post)
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
			bark_spaced_freqs[i] = bfcc_bark2freq(bark_sum);
			bark_sum += x->bark_spacing;
		};
				
		// i-2 is the correct number of filters
		x->num_filters = i-2;

		x->bfcc = (float *)t_resizebytes(x->bfcc, old_num_filters * sizeof(float), x->num_filters * sizeof(float));

		// initialize bfcc
		for(i=0; i<x->num_filters; i++)
				x->bfcc[i] = 0.0;
				
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
		bfcc_make_filterbank(x, bark_spaced_freqs);
	
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


static void bfcc_analyze(t_bfcc *x, t_floatarg start, t_floatarg n)
{
	int i, j, old_window, window, window_half, start_samp, end_samp, length_samp;
	t_atom *listOut;
	t_sample *signal_I;
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayname, garray_class)))
        pd_error(x, "%s: no such array", x->x_arrayname->s_name);
    else if(!garray_getfloatwords(a, &x->x_array_points, &x->x_vec))
    	pd_error(x, "%s: bad template for bfcc", x->x_arrayname->s_name);
	else
	{

	start_samp = start;
	
	if(start_samp < 0)
		start_samp = 0;

	if(n)
		end_samp = start_samp + n-1;
	else
		end_samp = start_samp + x->window-1;
		
	if(end_samp > x->x_array_points)
		end_samp = x->x_array_points-1;

	length_samp = end_samp - start_samp + 1;

	if(end_samp <= start_samp)
	{
		error("bad range of samples.");
		return;
	}
		
	if(length_samp > x->powers_of_two[x->pow_two_arr_size-1])
	{
		post("WARNING: bfcc: window truncated because requested size is larger than the current max_window setting. Use the max_window method to allow larger windows.");
		length_samp = x->powers_of_two[x->pow_two_arr_size-1];
		window = length_samp;
		end_samp = start_samp + window - 1;
	}
	else
	{
		i=0;
		while(length_samp > x->powers_of_two[i])
			i++;

		window = x->powers_of_two[i];
	}

	window_half = window * 0.5;

	old_window = x->window;
	x->window = window;
	bfcc_create_filterbank(x, x->bark_spacing, 1);
	
	// create local memory
	listOut = (t_atom *)getbytes(0);
	signal_I = (t_sample *)getbytes(0);
	listOut = (t_atom *)t_resizebytes(listOut, 0, x->num_filters * sizeof(t_atom));
	x->signal_R = (t_sample *)t_resizebytes(x->signal_R, old_window*sizeof(t_sample), window*sizeof(t_sample));
	signal_I = (t_sample *)t_resizebytes(signal_I, 0, window_half*sizeof(t_sample));

	// construct analysis window
	for(i=0, j=start_samp; j<=end_samp; i++, j++)
		x->signal_R[i] = x->x_vec[j].w_float;

	// window first	
	switch(x->window_function)
	{
		case 0:
			bfcc_blackman_window(length_samp, x->signal_R);
			break;
		case 1:
			bfcc_cosine_window(length_samp, x->signal_R);
			break;
		case 2:
			bfcc_hamming_window(length_samp, x->signal_R);
			break;
		case 3:
			bfcc_hann_window(length_samp, x->signal_R);
			break;
		default:
			break;
	};

	// then zero pad the end
	for(; i<window; i++)
		x->signal_R[i] = 0.0;

	mayer_realfft(window, x->signal_R);
	bfcc_realfft_unpack(window, window_half, x->signal_R, signal_I);
	bfcc_abs(window_half, x->signal_R, signal_I);

// power spectrum sometimes generates lower scores than magnitude. make it optional.
	if(x->power_spectrum)
		bfcc_square(window_half, x->signal_R);

	bfcc_filterbank_multiply(x->signal_R, x->normalize, x->x_filterbank, x->x_indices, x->num_filters);
	bfcc_compute_bfccs(x->bfcc, x->signal_R, x->num_filters);
		
	for(i=0; i<x->num_filters; i++)
		SETFLOAT(listOut+i, x->bfcc[i]);

	outlet_list(x->x_featureList, 0, x->num_filters, listOut);

	// free local memory
	t_freebytes(listOut, x->num_filters * sizeof(t_atom));
	t_freebytes(signal_I, window_half * sizeof(t_sample));

	}
}


// analyze the whole damn array
static void bfcc_bang(t_bfcc *x)
{
	int i, window, start_samp, end_samp, length_samp;
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayname, garray_class)))
        pd_error(x, "%s: no such array", x->x_arrayname->s_name);
    else if(!garray_getfloatwords(a, &x->x_array_points, &x->x_vec))
    	pd_error(x, "%s: bad template for bfcc", x->x_arrayname->s_name);
	else
	{

	start_samp = 0;
	end_samp = x->x_array_points-1;

	length_samp = end_samp - start_samp + 1;

	if(length_samp > x->powers_of_two[x->pow_two_arr_size-1])
	{
		post("WARNING: bfcc: window truncated because requested size is larger than the current max_window setting. Use the max_window method to allow larger windows.");
		length_samp = x->powers_of_two[x->pow_two_arr_size-1];
		window = length_samp;
		end_samp = start_samp + window - 1;
	}
	else
	{
		i=0;
		while(length_samp > x->powers_of_two[i])
			i++;

		window = x->powers_of_two[i];
	}

	bfcc_analyze(x, start_samp, window);

	}
}


static void bfcc_set(t_bfcc *x, t_symbol *s)
{
	t_garray *a;
	
	if(!(a = (t_garray *)pd_findbyclass(s, garray_class)))
		pd_error(x, "%s: no such array", s->s_name);
	else if(!garray_getfloatwords(a, &x->x_array_points, &x->x_vec))
		pd_error(x, "%s: bad template for bfcc", s->s_name);
	else
	    x->x_arrayname = s;
}



static void bfcc_hat(t_bfcc *x, t_floatarg f, t_floatarg g)
{
	int i;
	
	for(i=0; i<(int)g; i++)
		post("val %i: %f", i, x->x_filterbank[(int)f].filter[i]);
}


static void bfcc_samplerate(t_bfcc *x, t_floatarg sr)
{
	if(sr<64)
		x->sr = 64;
	else
		x->sr = sr;

	post("samplerate: %i", (int)x->sr);
}


static void bfcc_max_window(t_bfcc *x, t_floatarg w)
{
	int i;
		
	if(w<64)
		x->max_window_size = 64;
	else
		x->max_window_size = w;

	x->powers_of_two = (int *)t_resizebytes(x->powers_of_two, x->pow_two_arr_size*sizeof(int), sizeof(int));

	x->powers_of_two[0] = 64; // must have at least this large of a window

	i=1;
	while(x->powers_of_two[i-1] < x->max_window_size)
	{
		x->powers_of_two = (int *)t_resizebytes(x->powers_of_two, (i)*sizeof(int), (i+1)*sizeof(int));
		x->powers_of_two[i] = pow(2, i+6); // +6 because we're starting at 2**6
		i++;
	}

	x->pow_two_arr_size = i;
	
	post("maximum window size: %i", x->max_window_size);
}


static void bfcc_window(t_bfcc *x, t_floatarg w)
{
	int i, isPow2;
	
	isPow2 = (int)w && !( ((int)w-1) & (int)w );
	
	if( !isPow2 )
		error("requested window size is not a power of 2.");
	else
	{
		x->signal_R = (t_sample *)t_resizebytes(x->signal_R, x->window * sizeof(t_sample), w * sizeof(t_sample));
	
		x->window = (int)w;
 		
		// init signal buffer
		for(i=0; i<x->window; i++)
			x->signal_R[i] = 0.0;
			
		bfcc_create_filterbank(x, x->bark_spacing, 0);
		
		post("window size: %i. sampling rate: %i", (int)x->window, (int)x->sr);
	}
}


static void bfcc_window_function(t_bfcc *x, t_floatarg f)
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
static void bfcc_power_spectrum(t_bfcc *x, t_floatarg spec)
{
	x->power_spectrum = spec;
	
	if(x->power_spectrum)
		post("using power spectrum for BFCC computation.");
	else
		post("using magnitude spectrum for BFCC computation.");
}


static void bfcc_normalize(t_bfcc *x, t_floatarg norm)
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


static void bfcc_debug(t_bfcc *x, t_floatarg f)
{
	if((int)f < 0)
		x->debug = 0;
	else if((int)f > 2)
		x->debug = 2;
	else
		x->debug = (int)f;

	post("debug mode: %i", (int)f);
}


static void *bfcc_new(t_symbol *s, t_floatarg bark_spacing)
{
    t_bfcc *x = (t_bfcc *)pd_new(bfcc_class);
	int i;
	t_garray *a;

	x->x_featureList = outlet_new(&x->x_obj, gensym("list"));

	if(bark_spacing)
	{
		x->x_arrayname = s;

	    if(!(a = (t_garray *)pd_findbyclass(x->x_arrayname, garray_class)))
	        pd_error(x, "%s: no such array", x->x_arrayname->s_name);
	    else if(!garray_getfloatwords(a, &x->x_array_points, &x->x_vec))
	    	pd_error(x, "%s: bad template for bfcc", x->x_arrayname->s_name);

		x->bark_spacing = bark_spacing;  // should check for > 0

		if(x->bark_spacing <= 0)
		{
			error("bark spacing must be a positive real number. default value of 0.5 used instead.");
			x->bark_spacing = 0.5;
		};
	}
	else if(s)
	{
		x->x_arrayname = s;

	    if(!(a = (t_garray *)pd_findbyclass(x->x_arrayname, garray_class)))
	        pd_error(x, "%s: no such array", x->x_arrayname->s_name);
	    else if(!garray_getfloatwords(a, &x->x_array_points, &x->x_vec))
	    	pd_error(x, "%s: bad template for bfcc", x->x_arrayname->s_name);

		x->bark_spacing = 0.5;
	}
	else
	{
		error("bfcc: no array specified.");
		x->bark_spacing = 0.5;
	}

	x->sr = 44100.0;
	x->window = 1024.0;
	x->window_function = 3; // 3 is hann window
	x->power_spectrum = 0; // choose mag (0) or power (1) spec in the BFCC computation
	x->normalize = 1; // this is generally a good thing, but should be off for concatenative synth
	x->debug = 0;
    x->num_filters = 47;


	x->max_window_size = x->sr * 10; // 10 second maximum;
	x->powers_of_two = (int *)getbytes(0);
	x->powers_of_two = (int *)t_resizebytes(x->powers_of_two, 0, sizeof(int));


	x->powers_of_two[0] = 64; // must have at least this large of a window

	i=1;
	while(x->powers_of_two[i-1] < x->max_window_size)
	{
		x->powers_of_two = (int *)t_resizebytes(x->powers_of_two, (i)*sizeof(int), (i+1)*sizeof(int));
		x->powers_of_two[i] = pow(2, i+6); // +6 because we're starting at 2**6
		i++;
	}

	x->pow_two_arr_size = i;

	x->signal_R = (t_sample *)getbytes(0);
	x->signal_R = (t_sample *)t_resizebytes(x->signal_R, 0, x->window * sizeof(t_sample));

	// initialize signal_R
    for(i=0; i<x->window; i++)
		x->signal_R[i] = 0.0;
		
	x->bfcc = (float *)getbytes(0);
	x->bfcc = (float *)t_resizebytes(x->bfcc, 0, x->num_filters*sizeof(float));
		
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
		
	bfcc_create_filterbank(x, x->bark_spacing, 0);
    
    return (x);
}


static void bfcc_free(t_bfcc *x)
{
	int i;

	// free the input buffer memory
    t_freebytes(x->signal_R, x->window*sizeof(t_sample));
    
	// free the bfcc memory
    t_freebytes(x->bfcc, x->num_filters*sizeof(float));
       
    // free the filterbank memory
	for(i=0; i<x->num_filters; i++)
		t_freebytes(x->x_filterbank[i].filter, x->x_filterbank_sizes[i]*sizeof(float));
    
	t_freebytes(x->x_filterbank, x->num_filters*sizeof(t_filter));

	// free the indices memory
    t_freebytes(x->x_indices, x->num_filters*sizeof(t_indices));

	// free the powers of two table
    t_freebytes(x->powers_of_two, x->pow_two_arr_size*sizeof(int));
}


void bfcc_setup(void)
{
    bfcc_class = 
    class_new(
    	gensym("bfcc"),
    	(t_newmethod)bfcc_new,
    	(t_method)bfcc_free,
        sizeof(t_bfcc),
        CLASS_DEFAULT, 
        A_DEFSYM,
		A_DEFFLOAT,
		0
    );

	class_addbang(bfcc_class, bfcc_bang);

	class_addmethod(
		bfcc_class, 
        (t_method)bfcc_analyze,
		gensym("analyze"),
        A_DEFFLOAT,
        A_DEFFLOAT,
		0
	);

	class_addmethod(
		bfcc_class,
		(t_method)bfcc_set,
		gensym("set"),
		A_SYMBOL,
		0
	);
        
	class_addmethod(
		bfcc_class, 
        (t_method)bfcc_create_filterbank,
		gensym("filterbank"),
        A_DEFFLOAT,
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		bfcc_class, 
        (t_method)bfcc_hat,
		gensym("hat"),
        A_DEFFLOAT,
        A_DEFFLOAT,
		0
	);

	class_addmethod(
		bfcc_class, 
        (t_method)bfcc_window,
		gensym("window"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		bfcc_class, 
        (t_method)bfcc_samplerate,
		gensym("samplerate"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		bfcc_class, 
        (t_method)bfcc_max_window,
		gensym("max_window"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		bfcc_class, 
        (t_method)bfcc_window_function,
		gensym("window_function"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		bfcc_class, 
        (t_method)bfcc_power_spectrum,
		gensym("power_spectrum"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		bfcc_class, 
        (t_method)bfcc_normalize, 
		gensym("normalize"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		bfcc_class, 
        (t_method)bfcc_debug,
		gensym("debug"),
		A_DEFFLOAT,
		0
	);
}