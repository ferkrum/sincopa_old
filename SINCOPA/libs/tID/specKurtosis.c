/*

specKurtosis - A non-real-time spectral kurtosis analysis external.

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

version 0.0.3, May 15, 2010

• 0.0.3 added an ifndef M_PI for guaranteed windows compilation
• 0.0.2 adds a #define M_PI for windows compilation, and declares all functions except _setup static

*/

#include "m_pd.h"
#include <math.h>
#include <limits.h>
#include <stdio.h>
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

static t_class *specKurtosis_class;

typedef struct _specKurtosis
{
    t_object x_obj;
    float sr;
	float window;
	int window_function;
	int max_window_size;
	int *powers_of_two;
    int  pow_two_arr_size;
	t_sample *signal_R;
	t_word *x_vec;
	t_symbol *x_arrayname;
	int x_array_points;
    t_outlet *x_kurtosis;

} t_specKurtosis;


/* ---------------- dsp utility functions ---------------------- */

static void specKurtosis_blackman_window(int n, t_sample *in)
{
	int i;
	
	for(i=0; i<n; i++, in++)
    	*in *= 0.42 - (0.5 * cos(2*M_PI*i/n)) + (0.08 * cos(4*M_PI*i/n));
}


static void specKurtosis_cosine_window(int n, t_sample *in)
{
	int i;

	for(i=0; i<n; i++, in++)
    	*in *= sin(M_PI*i/n);
}

static void specKurtosis_hamming_window(int n, t_sample *in)
{	
	int i;

	for(i=0; i<n; i++, in++)
    	*in *= 0.5 - (0.46 * cos(2*M_PI*i/n));
}

static void specKurtosis_hann_window(int n, t_sample *in)
{
	int i;

	for(i=0; i<n; i++, in++)
    	*in *= 0.5 * (1 - cos(2*M_PI*i/n));
}
 		
static void specKurtosis_realfft_unpack(int n, int n_half, t_sample *input, t_sample *imag)
{
	int i, j;
		
	imag[0]=0.0;  // DC
	
	for(i=(n-1), j=1; i>n_half; i--, j++)
		imag[j] = input[i];
}

static void specKurtosis_abs(int n, t_sample *in1, t_sample *in2)
{
	while (n--)
    {   
        *in1 = sqrt((*in1 * *in1) + (*in2 * *in2));
        in1++;
        in2++;
    };
}

/* ---------------- END utility functions ---------------------- */




/* ------------------------ specKurtosis -------------------------------- */

static void specKurtosis_analyze(t_specKurtosis *x, t_floatarg start, t_floatarg n)
{
	int i, j, old_window, window, window_half, start_samp, end_samp, length_samp;
    float dividend, divisor, centroid, std, kurtosis, *bin_freqs;
	t_sample *signal_I;
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayname, garray_class)))
        pd_error(x, "%s: no such array", x->x_arrayname->s_name);
    else if(!garray_getfloatwords(a, &x->x_array_points, &x->x_vec))
    	pd_error(x, "%s: bad template for bfcc", x->x_arrayname->s_name);
	else
	{

	bin_freqs = (float *)getbytes(0);

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
		post("WARNING: specKurtosis: window truncated because requested size is larger than the current max_window setting. Use the max_window method to allow larger windows.");
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

	bin_freqs = (float *)t_resizebytes(bin_freqs, 0, window_half * sizeof(float));

	// freqs for each bin based on current window size
	for(i=0; i<window_half; i++)
		bin_freqs[i] = (x->sr/x->window) * i;
		
	// create local memory
	signal_I = (t_sample *)getbytes(0);
	signal_I = (t_sample *)t_resizebytes(signal_I, 0, window_half*sizeof(t_sample));
	
	// resize signal_R
	x->signal_R = (t_sample *)t_resizebytes(x->signal_R, old_window*sizeof(t_sample), window*sizeof(t_sample));

	// construct analysis window
	for(i=0, j=start_samp; j<=end_samp; i++, j++)
		x->signal_R[i] = x->x_vec[j].w_float;

	// window first	
	switch(x->window_function)
	{
		case 0:
			specKurtosis_blackman_window(length_samp, x->signal_R);
			break;
		case 1:
			specKurtosis_cosine_window(length_samp, x->signal_R);
			break;
		case 2:
			specKurtosis_hamming_window(length_samp, x->signal_R);
			break;
		case 3:
			specKurtosis_hann_window(length_samp, x->signal_R);
			break;
		default:
			break;
	};

	// then zero pad the end
	for(; i<window; i++)
		x->signal_R[i] = 0.0;

	mayer_realfft(window, x->signal_R);
	specKurtosis_realfft_unpack(window, window_half, x->signal_R, signal_I);
	specKurtosis_abs(window_half, x->signal_R, signal_I);

	dividend=0;
	divisor=0;
	centroid=0;
	
	for(i=0; i<window_half; i++)
	{
		dividend += x->signal_R[i] * bin_freqs[i];  // weight by bin freq
		divisor += x->signal_R[i];
	}
	
	if(divisor==0) divisor = 1; // don't divide by zero
	
	centroid = dividend/divisor;

	dividend=0;
	std=0;

	for(i=0; i<window_half; i++)
		dividend += pow((bin_freqs[i] - centroid), 2) * x->signal_R[i];

	std = sqrt(dividend/divisor);
	std = pow(std, 4);
	if(std==0) std = 1; // don't divide by zero
	
	dividend=0;
	kurtosis=0;

	for(i=0; i<window_half; i++)
		dividend += pow((bin_freqs[i] - centroid), 4) * x->signal_R[i];

	kurtosis = (dividend/divisor)/std;

	outlet_float(x->x_kurtosis, kurtosis);

	// free local memory
	t_freebytes(signal_I, window_half * sizeof(t_sample));
	}
}


// analyze the whole damn array
static void specKurtosis_bang(t_specKurtosis *x)
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
		post("WARNING: specKurtosis: window truncated because requested size is larger than the current max_window setting. Use the max_window method to allow larger windows.");
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

	specKurtosis_analyze(x, start_samp, window);

	}
}


static void specKurtosis_set(t_specKurtosis *x, t_symbol *s)
{
	t_garray *a;
	
	if(!(a = (t_garray *)pd_findbyclass(s, garray_class)))
		pd_error(x, "%s: no such array", s->s_name);
	else if(!garray_getfloatwords(a, &x->x_array_points, &x->x_vec))
		pd_error(x, "%s: bad template for specKurtosis", s->s_name);
	else
	    x->x_arrayname = s;
}


static void specKurtosis_samplerate(t_specKurtosis *x, t_floatarg sr)
{
	if(sr<64)
		x->sr = 64;
	else
		x->sr = sr;

	post("samplerate: %i", (int)x->sr);
}


static void specKurtosis_max_window(t_specKurtosis *x, t_floatarg w)
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


static void specKurtosis_window(t_specKurtosis *x, t_floatarg w)
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
					
		post("window size: %i. sampling rate: %i", (int)x->window, (int)x->sr);
	}
}


static void specKurtosis_window_function(t_specKurtosis *x, t_floatarg f)
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


static void *specKurtosis_new(t_symbol *s)
{
    t_specKurtosis *x = (t_specKurtosis *)pd_new(specKurtosis_class);
	int i;
	t_garray *a;

	x->x_kurtosis = outlet_new(&x->x_obj, &s_float);

	if(s)
	{
		x->x_arrayname = s;

	    if(!(a = (t_garray *)pd_findbyclass(x->x_arrayname, garray_class)))
	        pd_error(x, "%s: no such array", x->x_arrayname->s_name);
	    else if(!garray_getfloatwords(a, &x->x_array_points, &x->x_vec))
	    	pd_error(x, "%s: bad template for bfcc", x->x_arrayname->s_name);
	}
	else
		error("specKurtosis: no array specified.");

	x->sr = 44100.0;
	x->window = 1024.0;
	x->window_function = 3; // 3 is hann window

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
	
    return (x);
}


static void specKurtosis_free(t_specKurtosis *x)
{
	// free the input buffer memory
    t_freebytes(x->signal_R, x->window*sizeof(t_sample));

	// free the powers of two table
    t_freebytes(x->powers_of_two, x->pow_two_arr_size*sizeof(int));
}


void specKurtosis_setup(void)
{
    specKurtosis_class = 
    class_new(
    	gensym("specKurtosis"),
    	(t_newmethod)specKurtosis_new,
    	(t_method)specKurtosis_free,
        sizeof(t_specKurtosis),
        CLASS_DEFAULT, 
        A_DEFSYM,
		0
    );

	class_addbang(specKurtosis_class, specKurtosis_bang);

	class_addmethod(
		specKurtosis_class, 
        (t_method)specKurtosis_analyze,
		gensym("analyze"),
        A_DEFFLOAT,
        A_DEFFLOAT,
		0
	);

	class_addmethod(
		specKurtosis_class,
		(t_method)specKurtosis_set,
		gensym("set"),
		A_SYMBOL,
		0
	);

	class_addmethod(
		specKurtosis_class, 
        (t_method)specKurtosis_window,
		gensym("window"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		specKurtosis_class, 
        (t_method)specKurtosis_samplerate,
		gensym("samplerate"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		specKurtosis_class, 
        (t_method)specKurtosis_max_window,
		gensym("max_window"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		specKurtosis_class, 
        (t_method)specKurtosis_window_function,
		gensym("window_function"),
		A_DEFFLOAT,
		0
	);
}