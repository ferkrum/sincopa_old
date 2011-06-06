/*

zeroCrossing - A non-real-time zero crossing rate analysis external.

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.


version 0.0.3, March 8, 2010

¥ 0.0.3 declares all functions except _setup static
¥ 0.0.2 is part of the update that ensures all function names are prepended by the external name (bfcc_ or bfcc_tilde_, etc).

*/

#include "m_pd.h"
#include <math.h>
#include <limits.h>
#include <stdio.h>


static t_class *zeroCrossing_class;

typedef struct _zeroCrossing
{
    t_object x_obj;
	float window;
	int window_function;
    float freq_boundary;
	t_sample *signal_R;
	t_word *x_vec;
	t_symbol *x_arrayname;
	int x_array_points;
    t_outlet *x_crossings;

} t_zeroCrossing;


/* ------------------------ utility functions -------------------------------- */

static int zeroCrossing_signum(t_sample sample)
{
	int sign, crossings;
	
	sign=0;
	crossings = 0;
	
	if(sample>0)
		sign = 1;
	else if(sample<0)
		sign = -1;
	else
		sign = 0;
	
	return(sign);
}


static float zeroCrossing_zero_crossing_rate(int n, t_sample *in1)
{
	int i;
	float crossings;
	
	crossings = 0.0;
	
	for(i=1; i<n; i++)
		crossings += fabs(zeroCrossing_signum(in1[i]) - zeroCrossing_signum(in1[i-1]));
	
	crossings /= 2;
		
	return(crossings);
}
/* ------------------------ end utility functions -------------------------------- */


/* ------------------------ zeroCrossing -------------------------------- */

static void zeroCrossing_analyze(t_zeroCrossing *x, t_floatarg start, t_floatarg n)
{
	int i, j, old_window, window, start_samp, end_samp, crossings;
	
	
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

	window = end_samp - start_samp + 1;

	if(end_samp <= start_samp)
	{
		error("bad range of samples.");
		return;
	}
		
	old_window = x->window;
	x->window = window;

	// resize signal_R
	x->signal_R = (t_sample *)t_resizebytes(x->signal_R, old_window*sizeof(t_sample), window*sizeof(t_sample));

	// construct analysis window
	for(i=0, j=start_samp; j<=end_samp; i++, j++)
		x->signal_R[i] = x->x_vec[j].w_float;

	// then zero pad the end
	for(; i<window; i++)
		x->signal_R[i] = 0.0;

	crossings=0;
	crossings = zeroCrossing_zero_crossing_rate(window, x->signal_R);

	outlet_float(x->x_crossings, crossings);
	
	}
}


// analyze the whole damn array
static void zeroCrossing_bang(t_zeroCrossing *x)
{
	int window, start_samp, end_samp;
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayname, garray_class)))
        pd_error(x, "%s: no such array", x->x_arrayname->s_name);
    else if(!garray_getfloatwords(a, &x->x_array_points, &x->x_vec))
    	pd_error(x, "%s: bad template for bfcc", x->x_arrayname->s_name);
	else
	{

	start_samp = 0;
	end_samp = x->x_array_points-1;

	window = end_samp - start_samp + 1;

	zeroCrossing_analyze(x, start_samp, window);
	
	}
}


static void zeroCrossing_set(t_zeroCrossing *x, t_symbol *s)
{
	t_garray *a;
	
	if(!(a = (t_garray *)pd_findbyclass(s, garray_class)))
		pd_error(x, "%s: no such array", s->s_name);
	else if(!garray_getfloatwords(a, &x->x_array_points, &x->x_vec))
		pd_error(x, "%s: bad template for zeroCrossing", s->s_name);
	else
	    x->x_arrayname = s;
}


static void zeroCrossing_window(t_zeroCrossing *x, t_floatarg w)
{
	int i;
	
	if(w < 64)
	{
		error("minimum window size is 64 samples.");
		w = 64;
	};

	x->signal_R = (t_sample *)t_resizebytes(x->signal_R, x->window * sizeof(t_sample), w * sizeof(t_sample));
	x->window = (int)w;

	// init signal buffer
	for(i=0; i<x->window; i++)
		x->signal_R[i] = 0.0;
				
	post("window size: %i.", (int)x->window);
}


static void zeroCrossing_window_function(t_zeroCrossing *x, t_floatarg f)
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


static void *zeroCrossing_new(t_symbol *s)
{
    t_zeroCrossing *x = (t_zeroCrossing *)pd_new(zeroCrossing_class);
	int i;
	t_garray *a;

	x->x_crossings = outlet_new(&x->x_obj, &s_float);

	if(s)
	{
		x->x_arrayname = s;

	    if(!(a = (t_garray *)pd_findbyclass(x->x_arrayname, garray_class)))
	        pd_error(x, "%s: no such array", x->x_arrayname->s_name);
	    else if(!garray_getfloatwords(a, &x->x_array_points, &x->x_vec))
	    	pd_error(x, "%s: bad template for bfcc", x->x_arrayname->s_name);
	}
	else
		error("zeroCrossing: no array specified.");

	x->window = 1024.0;
	x->window_function = 3; // 3 is hann window

	x->signal_R = (t_sample *)getbytes(0);
	x->signal_R = (t_sample *)t_resizebytes(x->signal_R, 0, x->window * sizeof(t_sample));

	// initialize signal_R
    for(i=0; i<x->window; i++)
		x->signal_R[i] = 0.0;

    return (x);
}


static void zeroCrossing_free(t_zeroCrossing *x)
{
	// free the input buffer memory
    t_freebytes(x->signal_R, x->window*sizeof(t_sample));
}


void zeroCrossing_setup(void)
{
    zeroCrossing_class = 
    class_new(
    	gensym("zeroCrossing"),
    	(t_newmethod)zeroCrossing_new,
    	(t_method)zeroCrossing_free,
        sizeof(t_zeroCrossing),
        CLASS_DEFAULT, 
        A_DEFSYM,
		0
    );

	class_addbang(zeroCrossing_class, zeroCrossing_bang);

	class_addmethod(
		zeroCrossing_class, 
        (t_method)zeroCrossing_analyze,
		gensym("analyze"),
        A_DEFFLOAT,
        A_DEFFLOAT,
		0
	);

	class_addmethod(
		zeroCrossing_class,
		(t_method)zeroCrossing_set,
		gensym("set"),
		A_SYMBOL,
		0
	);

	class_addmethod(
		zeroCrossing_class, 
        (t_method)zeroCrossing_window,
		gensym("window"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		zeroCrossing_class, 
        (t_method)zeroCrossing_window_function,
		gensym("window_function"),
		A_DEFFLOAT,
		0
	);
}