/*

zeroCrossing~

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

version 0.2.1, March 8, 2010

¥ 0.2.1 declares all functions except _setup static
¥ 0.2.0 is part of the update that ensures all function names are prepended by the external name (bfcc_ or bfcc_tilde_, etc).
¥Ê0.1.9 removes power of 2 restriction for window size
¥ 0.1.8 goes with a normal help file name.

*/

#include "m_pd.h"
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>


static t_class *zeroCrossing_tilde_class;

typedef struct _zeroCrossing_tilde
{
    t_object x_obj;
    float sr;
    float n;
    int overlap;
    int window;
    double last_dsp_time;
    t_sample *signal_R;
    float x_f;
    t_outlet *x_crossings;

} t_zeroCrossing_tilde;



/* ------------------------ utility functions -------------------------------- */

static int zeroCrossing_tilde_signum(t_sample sample)
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


static float zeroCrossing_tilde_zero_crossing_rate(int n, t_sample *in1)
{
	int i;
	float crossings;
	
	crossings = 0.0;
	
	for(i=1; i<n; i++)
		crossings += fabs(zeroCrossing_tilde_signum(in1[i]) - zeroCrossing_tilde_signum(in1[i-1]));
	
	crossings /= 2;
		
	return(crossings);
}
/* ------------------------ end utility functions -------------------------------- */



/* ------------------------ zeroCrossing~ -------------------------------- */

static void zeroCrossing_tilde_bang(t_zeroCrossing_tilde *x)
{
    int i, j, window, bang_sample;
    float crossings;
    t_sample *signal_R;
	double current_time;

    window = x->window;

	signal_R = (t_sample *)getbytes(0);
	signal_R = (t_sample *)t_resizebytes(signal_R, 0, window * sizeof(t_sample));
    
	current_time = clock_gettimesince(x->last_dsp_time);
	bang_sample = (int)(((current_time/1000.0)*x->sr)+0.5); // round

	if (bang_sample < 0)
        bang_sample = 0;
    else if ( bang_sample >= x->n )
        bang_sample = x->n - 1;
            
	// construct analysis window using bang_sample as the end of the window
	for(i=0, j=bang_sample; i<window; i++, j++)
		signal_R[i] = x->signal_R[j];

	
	crossings=0;
	
	crossings = zeroCrossing_tilde_zero_crossing_rate(window, signal_R);
		
	outlet_float(x->x_crossings, crossings);

	// free local memory
    t_freebytes(signal_R, window*sizeof(t_sample));
}


static void zeroCrossing_tilde_window(t_zeroCrossing_tilde *x, t_floatarg w)
{
	int i;

	if(w < 64)
	{
		error("minimum window size is 64 samples.");
		w = 64;
	};
		

	x->signal_R = (t_sample *)t_resizebytes(x->signal_R, (x->window+x->n) * sizeof(t_sample), (w+x->n) * sizeof(t_sample));
	x->window = (int)w;

	// init signal buffer
	for(i=0; i<(x->window+x->n); i++)
		x->signal_R[i] = 0.0;
				
	post("window size: %i. overlap: %i. sampling rate: %i", (int)x->window, x->overlap, (int)(x->sr/x->overlap));
}


static void zeroCrossing_tilde_overlap(t_zeroCrossing_tilde *x, t_floatarg f)
{
	x->overlap = (int)f;

    post("overlap: %i", x->overlap);
}


static void *zeroCrossing_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
    t_zeroCrossing_tilde *x = (t_zeroCrossing_tilde *)pd_new(zeroCrossing_tilde_class);
	int i;
	
	x->x_crossings = outlet_new(&x->x_obj, &s_float);
	s=s;
	
	x->signal_R = (t_sample *)getbytes(0);

	if(argc > 0)
	{
		x->window = atom_getfloat(argv);
		
		if(x->window < 64)
		{
			error("minimum window size is 64 samples.");
			x->window = 64;
		};
	}
	else
		x->window = 1024;	
	
	
	x->sr = 44100.0;
	x->n = 64.0;
	x->overlap = 1;
	x->last_dsp_time = clock_getlogicaltime();

	x->signal_R = (t_sample *)t_resizebytes(x->signal_R, 0, (x->window+x->n) * sizeof(t_sample));

 	for(i=0; i<(x->window+x->n); i++)
		x->signal_R[i] = 0.0;
		
    post("zeroCrossing~: window size: %i", (int)x->window);
    
    return (x);
}


static t_int *zeroCrossing_tilde_perform(t_int *w)
{
    int i, n;

    t_zeroCrossing_tilde *x = (t_zeroCrossing_tilde *)(w[1]);

    t_sample *in = (t_float *)(w[2]);
    n = w[3];
 			
 	// shift signal buffer contents back.
	for(i=0; i<x->window; i++)
		x->signal_R[i] = x->signal_R[i+n];
	
	// write new block to end of signal buffer.
	for(i=0; i<n; i++)
		x->signal_R[(int)x->window+i] = in[i];
		
	x->last_dsp_time = clock_getlogicaltime();

    return (w+4);
}


static void zeroCrossing_tilde_dsp(t_zeroCrossing_tilde *x, t_signal **sp)
{
	int i;
	
	dsp_add(
		zeroCrossing_tilde_perform,
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
			
    	post("zeroCrossing~: window size: %i. overlap: %i. sampling rate: %i, block size: %i", (int)x->window, x->overlap, (int)(x->sr/x->overlap), (int)x->n);
	};
};

static void zeroCrossing_tilde_free(t_zeroCrossing_tilde *x)
{	
	// free the input buffer memory
    t_freebytes(x->signal_R, (x->window+x->n)*sizeof(t_sample));
}

void zeroCrossing_tilde_setup(void)
{
    zeroCrossing_tilde_class = 
    class_new(
    	gensym("zeroCrossing~"),
    	(t_newmethod)zeroCrossing_tilde_new,
    	(t_method)zeroCrossing_tilde_free,
        sizeof(t_zeroCrossing_tilde),
        CLASS_DEFAULT, 
        A_GIMME,
		0
    );

    CLASS_MAINSIGNALIN(zeroCrossing_tilde_class, t_zeroCrossing_tilde, x_f);

	class_addbang(zeroCrossing_tilde_class, zeroCrossing_tilde_bang);
	
	class_addmethod(
		zeroCrossing_tilde_class, 
        (t_method)zeroCrossing_tilde_window,
		gensym("window"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		zeroCrossing_tilde_class, 
        (t_method)zeroCrossing_tilde_overlap,
		gensym("overlap"),
		A_DEFFLOAT,
		0
	);
	
    class_addmethod(
    	zeroCrossing_tilde_class,
    	(t_method)zeroCrossing_tilde_dsp,
    	gensym("dsp"),
    	0
    );
}