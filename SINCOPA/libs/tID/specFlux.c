/*

specFlux - A non-real-time spectral flux analysis external.

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

version 0.0.3, May 15, 2010

¥Ê0.0.3 added an ifndef M_PI for guaranteed windows compilation
¥ 0.0.2 adds a #define M_PI for windows compilation, and declares all functions except _setup static

*/

#include "m_pd.h"
#include <math.h>
#include <limits.h>
#include <stdio.h>
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

static t_class *specFlux_class;

typedef struct _specFlux
{
    t_object x_obj;
    float sr;
	float window;
	int window_function;
	int max_window_size;
	int *powers_of_two;
    int  pow_two_arr_size;
    int separation;
    int squaredDiff;
	int normalize;
	t_word *x_vec;
	t_symbol *x_arrayname;
	int x_array_points;
    t_outlet *x_flux;
    t_outlet *x_fluxList;

} t_specFlux;


/* ---------------- dsp utility functions ---------------------- */

static void specFlux_blackman_window(int n, t_sample *in)
{
	int i;
	
	for(i=0; i<n; i++, in++)
    	*in *= 0.42 - (0.5 * cos(2*M_PI*i/n)) + (0.08 * cos(4*M_PI*i/n));
}


static void specFlux_cosine_window(int n, t_sample *in)
{
	int i;

	for(i=0; i<n; i++, in++)
    	*in *= sin(M_PI*i/n);
}

static void specFlux_hamming_window(int n, t_sample *in)
{	
	int i;

	for(i=0; i<n; i++, in++)
    	*in *= 0.5 - (0.46 * cos(2*M_PI*i/n));
}

static void specFlux_hann_window(int n, t_sample *in)
{
	int i;

	for(i=0; i<n; i++, in++)
    	*in *= 0.5 * (1 - cos(2*M_PI*i/n));
}
 		
static void specFlux_realfft_unpack(int n, int n_half, t_sample *input, t_sample *imag)
{
	int i, j;
		
	imag[0]=0.0;  // DC
	
	for(i=(n-1), j=1; i>n_half; i--, j++)
		imag[j] = input[i];
}

static void specFlux_abs(int n, t_sample *in1, t_sample *in2)
{
	while (n--)
    {   
        *in1 = sqrt( (*in1 * *in1) + (*in2 * *in2) );
        in1++;
        in2++;
    };
}

static void specFlux_tilde_normal(int n, t_sample *in)
{
	int i;
	float sum, norm_scalar;
	
	sum=0;
	norm_scalar=0;
	
	for(i=0; i<n; i++)
		sum += in[i];
	
	if(sum==0) sum=1;
	
	norm_scalar = 1/sum;
	
	for(i=0; i<n; i++)
		in[i] *= norm_scalar;
}

/* ---------------- END utility functions ---------------------- */




/* ------------------------ specFlux -------------------------------- */

static void specFlux_analyze(t_specFlux *x, t_floatarg start, t_floatarg n)
{
	int i, j, window, window_half, rear_start_samp, rear_end_samp, start_samp, end_samp, length_samp;
    float diff, val, flux;
    t_atom *listOut;
    t_sample *signal2_R, *signal2_I, *signal1_R, *signal1_I;
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
		post("WARNING: specFlux: window truncated because requested size is larger than the current max_window setting. Use the max_window method to allow larger windows.");
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
	x->window = window;


	// create local memory
	listOut = (t_atom *)getbytes(0);
	signal2_R = (t_sample *)getbytes(0);
	signal2_I = (t_sample *)getbytes(0);
	signal1_R = (t_sample *)getbytes(0);
	signal1_I = (t_sample *)getbytes(0);
	listOut = (t_atom *)t_resizebytes(listOut, 0, window_half * sizeof(t_atom));
	signal2_R = (t_sample *)t_resizebytes(signal2_R, 0, window*sizeof(t_sample));
	signal2_I = (t_sample *)t_resizebytes(signal2_I, 0, window_half*sizeof(t_sample));
	signal1_R = (t_sample *)t_resizebytes(signal1_R, 0, window*sizeof(t_sample));
	signal1_I = (t_sample *)t_resizebytes(signal1_I, 0, window_half*sizeof(t_sample));


	rear_start_samp = start_samp - x->separation;
	if(rear_start_samp < 0)
		rear_start_samp = 0;

	rear_end_samp = end_samp - x->separation;
	if(rear_end_samp >  x->x_array_points-1)
		rear_end_samp =  x->x_array_points-1;

	// construct the rear analysis window
	for(i=0, j=rear_start_samp; j<=rear_end_samp; i++, j++)
		signal1_R[i] = x->x_vec[j].w_float;

	// window first	
	switch(x->window_function)
	{
		case 0:
			specFlux_blackman_window(length_samp, signal1_R);
			break;
		case 1:
			specFlux_cosine_window(length_samp, signal1_R);
			break;
		case 2:
			specFlux_hamming_window(length_samp, signal1_R);
			break;
		case 3:
			specFlux_hann_window(length_samp, signal1_R);
			break;
		default:
			break;
	};

	// then zero pad the end
	for(; i<window; i++)
		signal1_R[i] = 0.0;

	mayer_realfft(window, signal1_R);
	specFlux_realfft_unpack(window, window_half, signal1_R, signal1_I);
	specFlux_abs(window_half, signal1_R, signal1_I);
	
	if(x->normalize)
		specFlux_tilde_normal(window_half, signal1_R);

	// construct the forward analysis window
	for(i=0, j=start_samp; j<=end_samp; i++, j++)
		signal2_R[i] = x->x_vec[j].w_float;

	// window first	
	switch(x->window_function)
	{
		case 0:
			specFlux_blackman_window(length_samp, signal2_R);
			break;
		case 1:
			specFlux_cosine_window(length_samp, signal2_R);
			break;
		case 2:
			specFlux_hamming_window(length_samp, signal2_R);
			break;
		case 3:
			specFlux_hann_window(length_samp, signal2_R);
			break;
		default:
			break;
	};

	// then zero pad the end
	for(; i<window; i++)
		signal2_R[i] = 0.0;

	mayer_realfft(window, signal2_R);
	specFlux_realfft_unpack(window, window_half, signal2_R, signal2_I);
	specFlux_abs(window_half, signal2_R, signal2_I);
	
	if(x->normalize)
		specFlux_tilde_normal(window_half, signal2_R);
	
	flux=0;

	for(i=0; i<window_half; i++)
	{
		diff = signal2_R[i] - signal1_R[i];
		
		if(x->squaredDiff)
			val = diff*diff;
		else
			val = fabs(diff);
			
		SETFLOAT(listOut+i, diff);
		flux += val;
	}
		
 	outlet_list(x->x_fluxList, 0, window_half, listOut);
	outlet_float(x->x_flux, flux);

	// free local memory
	t_freebytes(listOut, window_half * sizeof(t_atom));
	t_freebytes(signal1_R, window * sizeof(t_sample));
	t_freebytes(signal1_I, window_half * sizeof(t_sample));
	t_freebytes(signal2_R, window * sizeof(t_sample));
	t_freebytes(signal2_I, window_half * sizeof(t_sample));

	}
}


// analyze the whole damn array
static void specFlux_bang(t_specFlux *x)
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
		post("WARNING: specFlux: window truncated because requested size is larger than the current max_window setting. Use the max_window method to allow larger windows.");
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

	specFlux_analyze(x, start_samp, window);

	}
}


static void specFlux_set(t_specFlux *x, t_symbol *s)
{
	t_garray *a;
	
	if(!(a = (t_garray *)pd_findbyclass(s, garray_class)))
		pd_error(x, "%s: no such array", s->s_name);
	else if(!garray_getfloatwords(a, &x->x_array_points, &x->x_vec))
		pd_error(x, "%s: bad template for specFlux", s->s_name);
	else
	    x->x_arrayname = s;
}


static void specFlux_samplerate(t_specFlux *x, t_floatarg sr)
{
	if(sr<64)
		x->sr = 64;
	else
		x->sr = sr;

	post("samplerate: %i", (int)x->sr);
}


static void specFlux_max_window(t_specFlux *x, t_floatarg w)
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


static void specFlux_window(t_specFlux *x, t_floatarg w)
{
	int isPow2;
	
	isPow2 = (int)w && !( ((int)w-1) & (int)w );
	
	if( !isPow2 )
		error("requested window size is not a power of 2.");
	else
	{
		x->window = (int)w;
		post("window size: %i. sampling rate: %i", (int)x->window, (int)x->sr);
	}
}


static void specFlux_window_function(t_specFlux *x, t_floatarg f)
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


static void specFlux_separation(t_specFlux *x, t_floatarg s)
{
	x->separation = (int)s;

	if(x->separation > x->window)
	{
		error("analysis windows cannot be more than %i samples apart", x->window);
		x->separation = 128;
	}
	else if(x->separation < 0)
	{
		error("frame separation must be > 0");
		x->separation = 128;
	};
		
    post("frame separation: %i", x->separation);

}


static void specFlux_squaredDiff(t_specFlux *x, t_floatarg sd)
{
	if(sd<0)
		x->squaredDiff = 0;
	else if (sd>1)
		x->squaredDiff = 1;
	else
		x->squaredDiff = sd;
	
	if(x->squaredDiff)
		post("difference**2");
	else
		post("|difference|");
}


static void specFlux_normalize(t_specFlux *x, t_floatarg norm)
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


static void *specFlux_new(t_symbol *s, t_floatarg separation)
{
    t_specFlux *x = (t_specFlux *)pd_new(specFlux_class);
	int i;
	t_garray *a;

	x->x_flux = outlet_new(&x->x_obj, &s_float);
	x->x_fluxList = outlet_new(&x->x_obj, gensym("list"));

	if(separation)
	{
		x->x_arrayname = s;

	    if(!(a = (t_garray *)pd_findbyclass(x->x_arrayname, garray_class)))
	        pd_error(x, "%s: no such array", x->x_arrayname->s_name);
	    else if(!garray_getfloatwords(a, &x->x_array_points, &x->x_vec))
	    	pd_error(x, "%s: bad template for specBrightness", x->x_arrayname->s_name);

		x->separation = separation;
		
		if(x->separation > 1024)
		{
			error("analysis windows cannot by more than the window size");
			x->separation = 128;
		}
		else if(x->separation < 0)
		{
			error("frame separation must be > 0");
			x->separation = 128;
		};

	}
	else if(s)
	{
		x->x_arrayname = s;

	    if(!(a = (t_garray *)pd_findbyclass(x->x_arrayname, garray_class)))
	        pd_error(x, "%s: no such array", x->x_arrayname->s_name);
	    else if(!garray_getfloatwords(a, &x->x_array_points, &x->x_vec))
	    	pd_error(x, "%s: bad template for specBrightness", x->x_arrayname->s_name);

		x->separation = 128;
	}
	else
	{
		error("specFlux: no array specified.");
		x->separation = 128;
	}

	x->sr = 44100.0;
	x->window = 1024.0;
	x->squaredDiff = 0; // absolute value by default
	x->normalize = 0;
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

    post("specFlux: window size: %i. separation: %i", (int)x->window, (int)x->separation);

    return (x);
}


static void specFlux_free(t_specFlux *x)
{
	// free the powers of two table
    t_freebytes(x->powers_of_two, x->pow_two_arr_size*sizeof(int));
}


void specFlux_setup(void)
{
    specFlux_class = 
    class_new(
    	gensym("specFlux"),
    	(t_newmethod)specFlux_new,
    	(t_method)specFlux_free,
        sizeof(t_specFlux),
        CLASS_DEFAULT, 
        A_DEFSYM,
		A_DEFFLOAT,
		0
    );

	class_addbang(specFlux_class, specFlux_bang);

	class_addmethod(
		specFlux_class, 
        (t_method)specFlux_analyze,
		gensym("analyze"),
        A_DEFFLOAT,
        A_DEFFLOAT,
		0
	);

	class_addmethod(
		specFlux_class,
		(t_method)specFlux_set,
		gensym("set"),
		A_SYMBOL,
		0
	);

	class_addmethod(
		specFlux_class, 
        (t_method)specFlux_window,
		gensym("window"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		specFlux_class, 
        (t_method)specFlux_samplerate,
		gensym("samplerate"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		specFlux_class, 
        (t_method)specFlux_max_window,
		gensym("max_window"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		specFlux_class, 
        (t_method)specFlux_window_function,
		gensym("window_function"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		specFlux_class, 
        (t_method)specFlux_squaredDiff, 
		gensym("squaredDiff"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		specFlux_class, 
        (t_method)specFlux_separation,
		gensym("separation"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		specFlux_class, 
        (t_method)specFlux_normalize,
		gensym("normalize"),
		A_DEFFLOAT,
		0
	);
}