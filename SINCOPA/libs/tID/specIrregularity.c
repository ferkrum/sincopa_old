/*

specIrregularity - A non-real-time spectral irregularity analysis external.

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

version 0.0.2, May 15, 2010

¥Ê0.0.2 added an ifndef M_PI for guaranteed windows compilation

*/

#include "m_pd.h"
#include <math.h>
#include <limits.h>
#include <stdio.h>
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

static t_class *specIrregularity_class;

typedef struct _specIrregularity
{
    t_object x_obj;
    float sr;
	float window;
	int window_function;
	int max_window_size;
	int *powers_of_two;
    int  pow_two_arr_size;
    int algorithm;
	t_sample *signal_R;
	t_word *x_vec;
	t_symbol *x_arrayname;
	int x_array_points;
    t_outlet *x_irregularity;

} t_specIrregularity;


/* ---------------- dsp utility functions ---------------------- */

static void specIrregularity_blackman_window(int n, t_sample *in)
{
	int i;
	
	for(i=0; i<n; i++, in++)
    	*in *= 0.42 - (0.5 * cos(2*M_PI*i/n)) + (0.08 * cos(4*M_PI*i/n));
}


static void specIrregularity_cosine_window(int n, t_sample *in)
{
	int i;

	for(i=0; i<n; i++, in++)
    	*in *= sin(M_PI*i/n);
}

static void specIrregularity_hamming_window(int n, t_sample *in)
{	
	int i;

	for(i=0; i<n; i++, in++)
    	*in *= 0.5 - (0.46 * cos(2*M_PI*i/n));
}

static void specIrregularity_hann_window(int n, t_sample *in)
{
	int i;

	for(i=0; i<n; i++, in++)
    	*in *= 0.5 * (1 - cos(2*M_PI*i/n));
}
 		
static void specIrregularity_realfft_unpack(int n, int n_half, t_sample *input, t_sample *imag)
{
	int i, j;
		
	imag[0]=0.0;  // DC
	
	for(i=(n-1), j=1; i>n_half; i--, j++)
		imag[j] = input[i];
}

static void specIrregularity_abs(int n, t_sample *in1, t_sample *in2)
{
	while (n--)
    {   
        *in1 = sqrt( (*in1 * *in1) + (*in2 * *in2) );
        in1++;
        in2++;
    };
}

/* ---------------- END utility functions ---------------------- */




/* ------------------------ specIrregularity -------------------------------- */

static void specIrregularity_analyze(t_specIrregularity *x, t_floatarg start, t_floatarg n)
{
	int i, j, old_window, window, window_half, start_samp, end_samp, length_samp;
    float divisor, irregularity;
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
		post("WARNING: specIrregularity: window truncated because requested size is larger than the current max_window setting. Use the max_window method to allow larger windows.");
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
			specIrregularity_blackman_window(length_samp, x->signal_R);
			break;
		case 1:
			specIrregularity_cosine_window(length_samp, x->signal_R);
			break;
		case 2:
			specIrregularity_hamming_window(length_samp, x->signal_R);
			break;
		case 3:
			specIrregularity_hann_window(length_samp, x->signal_R);
			break;
		default:
			break;
	};

	// then zero pad the end
	for(; i<window; i++)
		x->signal_R[i] = 0.0;

	mayer_realfft(window, x->signal_R);
	specIrregularity_realfft_unpack(window, window_half, x->signal_R, signal_I);
	specIrregularity_abs(window_half, x->signal_R, signal_I);

	divisor=0.0;
	irregularity=0.0;
	
	if(x->algorithm)
	{
		// Krimphoff
		for(i=1; i<window_half-1; i++)
		{
			float local_avg;
			local_avg = 0;
			
			for(j=0; j<3; j++)
				local_avg += x->signal_R[i-1+j];
			
			local_avg *= 0.333333333333;
				
			irregularity += fabs(x->signal_R[i] - local_avg);
			//irregularity = log10(irregularity);
		}
	}
	else
	{
		// Jensen
		for(i=0; i<window_half; i++)
		{
			if(i==window_half-1)
				irregularity += x->signal_R[i] * x->signal_R[i];
			else
				irregularity += pow(x->signal_R[i] - x->signal_R[i+1], 2);
			
			divisor += x->signal_R[i] * x->signal_R[i];
		}
		
		irregularity /= divisor;
	}
		
	outlet_float(x->x_irregularity, irregularity);

	// free local memory
	t_freebytes(signal_I, window_half * sizeof(t_sample));

	}
}


// analyze the whole damn array
static void specIrregularity_bang(t_specIrregularity *x)
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
		post("WARNING: specIrregularity: window truncated because requested size is larger than the current max_window setting. Use the max_window method to allow larger windows.");
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

	specIrregularity_analyze(x, start_samp, window);

	}
}


static void specIrregularity_set(t_specIrregularity *x, t_symbol *s)
{
	t_garray *a;
	
	if(!(a = (t_garray *)pd_findbyclass(s, garray_class)))
		pd_error(x, "%s: no such array", s->s_name);
	else if(!garray_getfloatwords(a, &x->x_array_points, &x->x_vec))
		pd_error(x, "%s: bad template for specIrregularity", s->s_name);
	else
	    x->x_arrayname = s;
}


static void specIrregularity_samplerate(t_specIrregularity *x, t_floatarg sr)
{
	if(sr<64)
		x->sr = 64;
	else
		x->sr = sr;

	post("samplerate: %i", (int)x->sr);
}


static void specIrregularity_max_window(t_specIrregularity *x, t_floatarg w)
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


static void specIrregularity_window(t_specIrregularity *x, t_floatarg w)
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


static void specIrregularity_window_function(t_specIrregularity *x, t_floatarg f)
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

static void specIrregularity_algorithm(t_specIrregularity *x, t_floatarg a)
{
	if(a > 1.0)
		a = 1.0;
	else if(a < 0.0)
		a = 0.0;
	
	x->algorithm = (int)a;
	
	switch(x->algorithm)
	{
		case 0:
			post("Jensen irregularity.");
			break;
		case 1:
			post("Krimphoff irregularity.");
			break;
		default:
			break;
	};
}

static void *specIrregularity_new(t_symbol *s)
{
    t_specIrregularity *x = (t_specIrregularity *)pd_new(specIrregularity_class);
	int i;
	t_garray *a;

	x->x_irregularity = outlet_new(&x->x_obj, &s_float);

	if(s)
	{
		x->x_arrayname = s;

	    if(!(a = (t_garray *)pd_findbyclass(x->x_arrayname, garray_class)))
	        pd_error(x, "%s: no such array", x->x_arrayname->s_name);
	    else if(!garray_getfloatwords(a, &x->x_array_points, &x->x_vec))
	    	pd_error(x, "%s: bad template for bfcc", x->x_arrayname->s_name);
	}
	else
		error("specIrregularity: no array specified.");

	x->algorithm = 0; // Jensen by default. 1 for Krimphoff
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


static void specIrregularity_free(t_specIrregularity *x)
{
	// free the input buffer memory
    t_freebytes(x->signal_R, x->window*sizeof(t_sample));

	// free the powers of two table
    t_freebytes(x->powers_of_two, x->pow_two_arr_size*sizeof(int));
}


void specIrregularity_setup(void)
{
    specIrregularity_class = 
    class_new(
    	gensym("specIrregularity"),
    	(t_newmethod)specIrregularity_new,
    	(t_method)specIrregularity_free,
        sizeof(t_specIrregularity),
        CLASS_DEFAULT, 
        A_DEFSYM,
		0
    );

	class_addbang(specIrregularity_class, specIrregularity_bang);

	class_addmethod(
		specIrregularity_class, 
        (t_method)specIrregularity_analyze,
		gensym("analyze"),
        A_DEFFLOAT,
        A_DEFFLOAT,
		0
	);

	class_addmethod(
		specIrregularity_class,
		(t_method)specIrregularity_set,
		gensym("set"),
		A_SYMBOL,
		0
	);

	class_addmethod(
		specIrregularity_class, 
        (t_method)specIrregularity_window,
		gensym("window"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		specIrregularity_class, 
        (t_method)specIrregularity_samplerate,
		gensym("samplerate"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		specIrregularity_class, 
        (t_method)specIrregularity_max_window,
		gensym("max_window"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		specIrregularity_class, 
        (t_method)specIrregularity_window_function,
		gensym("window_function"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		specIrregularity_class, 
        (t_method)specIrregularity_algorithm,
		gensym("algorithm"),
		A_DEFFLOAT,
		0
	);
}