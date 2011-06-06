/*

specFlux~

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.


version 0.2.2, May 15, 2010

¥Ê0.2.2 added an ifndef M_PI for guaranteed windows compilation
¥ 0.2.1 implements mayer_realfft
¥ 0.2.0 making squared difference vs. absolute value of difference an option. Also making a second outlet, so that it reports the full list of squared differences out the left, and the sum of the squared differences at the right.
¥ 0.1.9 fixed bug: normalization infrastructure was there, but it wasn't actually implemented in the _bang function! Also, trying fabs(diff) instead of squared difference to see if a larger number will be more useful.
¥ 0.1.8 added normalization option

*/

#include "m_pd.h"
#include <math.h>
#include <limits.h>
#include <stdio.h>
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

static t_class *specFlux_tilde_class;

typedef struct _specFlux_tilde
{
    t_object x_obj;
    float sr;
    float n;
    int overlap;
    int squaredDiff;
    int normalize;
    int window;
    double last_dsp_time;
    int separation;
    t_sample *signal_R;
    float *hann;
    float x_f;
    t_outlet *x_fluxList;
    t_outlet *x_flux;

} t_specFlux_tilde;



/* ------------------------ spectrum functions -------------------------------- */

void specFlux_tilde_hann(int n, t_sample *in, float *hann)
{	
	while (n--)
    {
    	*in = *in * *hann;

    	in++;
    	hann++;
    };
}

void specFlux_tilde_realfft_unpack(int n, int n_half, t_sample *input, t_sample *imag)
{
	int i, j;
		
	imag[0]=0.0;  // DC
	
	for(i=(n-1), j=1; i>n_half; i--, j++)
		imag[j] = input[i];
}

void specFlux_tilde_abs(int n, t_sample *in1, t_sample *in2)
{
	while (n--)
    {   
        *in1 = sqrt( (*in1 * *in1) + (*in2 * *in2) );
        in1++;
        in2++;
    };
}

void specFlux_tilde_normal(int n, t_sample *in)
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
/* ------------------------ end spectrum functions -------------------------------- */



/* ------------------------ specFlux~ -------------------------------- */

void specFlux_tilde_bang(t_specFlux_tilde *x)
{
    int i, j, window, window_half, bang_sample;
    float diff, val, flux;
    t_atom *listOut;
    t_sample *signal2_R, *signal2_I, *signal1_R, *signal1_I;
	double current_time;

    window = x->window;
    window_half = window*0.5;

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
		
	current_time = clock_gettimesince(x->last_dsp_time);
	bang_sample = (int)(((current_time/1000.0)*x->sr)+0.5); // round

	if (bang_sample < 0)
        bang_sample = 0;
    else if ( bang_sample >= x->n )
        bang_sample = x->n - 1;
            
	// construct newest analysis window
	for(i=0, j=bang_sample; i<window; i++, j++)
		signal2_R[i] = x->signal_R[window + j];

	// construct older analysis window x->separation frames earlier
	for(i=0, j=bang_sample; i<window; i++, j++)
		signal1_R[i] = x->signal_R[window - x->separation + j];
	
	specFlux_tilde_hann(window, signal1_R, x->hann);
	mayer_realfft(window, signal1_R);
	specFlux_tilde_realfft_unpack(window, window_half, signal1_R, signal1_I);
	specFlux_tilde_abs(window_half, signal1_R, signal1_I);
	
	if(x->normalize)
		specFlux_tilde_normal(window_half, signal1_R);

	specFlux_tilde_hann(window, signal2_R, x->hann);
	mayer_realfft(window, signal2_R);
	specFlux_tilde_realfft_unpack(window, window_half, signal2_R, signal2_I);
	specFlux_tilde_abs(window_half, signal2_R, signal2_I);
	
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


void specFlux_tilde_window(t_specFlux_tilde *x, t_floatarg f)
{
	int i, isPow2;
	
	isPow2 = (int)f && !( ((int)f-1) & (int)f );
	
	if( !isPow2 )
		error("requested window size is not a power of 2");
	else
	{
		x->signal_R = (t_sample *)t_resizebytes(x->signal_R, ((x->window*2)+x->n) * sizeof(t_sample), (f*2+x->n) * sizeof(t_sample));
		x->hann = (float *)t_resizebytes(x->hann, x->window * sizeof(float), f * sizeof(float));
		x->window = (int)f;
		
		for(i=0; i<x->window; i++)
			x->hann[i] = 0.5 * (1 - cos(2*M_PI*i/x->window));

		// init signal buffer
		for(i=0; i<((x->window*2)+x->n); i++)
			x->signal_R[i] = 0.0;
					
		post("window size: %i. overlap: %i. sampling rate: %i", (int)x->window, x->overlap, (int)(x->sr/x->overlap));
	}
}


void specFlux_tilde_separation(t_specFlux_tilde *x, t_floatarg f)
{
	x->separation = (int)f;

	if(x->separation > x->window)
	{
		error("analysis windows cannot be more than %i frames apart", x->window);
		x->separation = 128;
	}
	else if(x->separation < 0)
	{
		error("frame separation must be > 0");
		x->separation = 128;
	};
		
    post("frame separation: %i", x->separation);

}


void specFlux_tilde_squaredDiff(t_specFlux_tilde *x, t_floatarg sd)
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


void specFlux_tilde_overlap(t_specFlux_tilde *x, t_floatarg f)
{
	x->overlap = (int)f;

    post("overlap: %i", x->overlap);

}


void specFlux_tilde_normalize(t_specFlux_tilde *x, t_floatarg norm)
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


void *specFlux_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
    t_specFlux_tilde *x = (t_specFlux_tilde *)pd_new(specFlux_tilde_class);
	int i, isPow2;
	s=s;
	
	x->x_flux = outlet_new(&x->x_obj, &s_float);
	x->x_fluxList = outlet_new(&x->x_obj, gensym("list"));

	x->signal_R = (t_sample *)getbytes(0);
	x->hann = (float *)getbytes(0);
	
	if(argc > 1)
	{
		x->window = atom_getfloat(argv);
		isPow2 = (int)x->window && !( ((int)x->window-1) & (int)x->window );
		
		if(!isPow2)
		{
			error("requested window size is not a power of 2. default value of 1024 used instead");
			x->window = 1024;
		};
		
		x->separation = atom_getfloat(argv+1);
		
		if(x->separation > x->window)
		{
			error("analysis windows cannot be more than %i frames apart", x->window);
			x->separation = x->window;
		}
		else if(x->separation < 0)
		{
			error("frame separation must be > 0");
			x->separation = x->window*0.125;
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
		
		x->separation = x->window*0.125;
	}
	else
	{
		x->window = 1024;
		x->separation = x->window*0.125;
	}
	
	
	x->sr = 44100.0;
	x->n = 64.0;
	x->overlap = 1;
	x->squaredDiff = 0; // absolute value by default
	x->normalize = 1;
	x->last_dsp_time = clock_getlogicaltime();

	x->signal_R = (t_sample *)t_resizebytes(x->signal_R, 0, ((x->window*2)+x->n) * sizeof(t_sample));
	x->hann = (float *)t_resizebytes(x->hann, 0, x->window * sizeof(float));
	
	// initialize hann window
 	for(i=0; i<x->window; i++)
 		x->hann[i] = 0.5 * (1 - cos(2*M_PI*i/x->window));
 
	for(i=0; i<((x->window*2)+x->n); i++)
		x->signal_R[i] = 0.0;
			
    post("specFlux~: window size: %i. separation: %i", (int)x->window, (int)x->separation);
    
    return (x);
}


t_int *specFlux_tilde_perform(t_int *w)
{
    int i, n;

    t_specFlux_tilde *x = (t_specFlux_tilde *)(w[1]);

    t_sample *in = (t_float *)(w[2]);
    n = w[3];
 			
 	// shift signal buffer contents back.
	for(i=0; i<(x->window*2); i++)
		x->signal_R[i] = x->signal_R[i+n];
	
	// write new block to end of signal buffer.
	for(i=0; i<n; i++)
		x->signal_R[(int)(x->window*2)+i] = in[i];
		
	x->last_dsp_time = clock_getlogicaltime();

    return (w+4);
}


void specFlux_tilde_dsp(t_specFlux_tilde *x, t_signal **sp)
{
	int i;
	
	dsp_add(
		specFlux_tilde_perform,
		3,
		x,
		sp[0]->s_vec,
		sp[0]->s_n
	); 

// compare n to stored n and recalculate filterbank if different
	if( sp[0]->s_sr != x->sr || sp[0]->s_n != x->n )
	{
		x->signal_R = (t_sample *)t_resizebytes(x->signal_R, ((x->window*2)+x->n) * sizeof(t_sample), ((x->window*2)+sp[0]->s_n) * sizeof(t_sample));

		x->sr = sp[0]->s_sr;
		x->n = sp[0]->s_n;
		x->last_dsp_time = clock_getlogicaltime();

		// init signal buffer
		for(i=0; i<((x->window*2)+x->n); i++)
			x->signal_R[i] = 0.0;
			
    	post("specFlux~: window size: %i. overlap: %i. sampling rate: %i, block size: %i", (int)x->window, x->overlap, (int)(x->sr/x->overlap), (int)x->n);
	};
};

static void specFlux_tilde_free(t_specFlux_tilde *x)
{
	// free the input buffer memory
    t_freebytes(x->signal_R, ((x->window*2)+x->n)*sizeof(t_sample));

	// free the hann window memory
    t_freebytes(x->hann, x->window*sizeof(float));
}

void specFlux_tilde_setup(void)
{
    specFlux_tilde_class = 
    class_new(
    	gensym("specFlux~"),
    	(t_newmethod)specFlux_tilde_new,
    	(t_method)specFlux_tilde_free,
        sizeof(t_specFlux_tilde),
        CLASS_DEFAULT, 
        A_GIMME,
		0
    );

    CLASS_MAINSIGNALIN(specFlux_tilde_class, t_specFlux_tilde, x_f);

	class_addbang(specFlux_tilde_class, specFlux_tilde_bang);
	
	class_addmethod(
		specFlux_tilde_class, 
        (t_method)specFlux_tilde_window,
		gensym("window"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		specFlux_tilde_class, 
        (t_method)specFlux_tilde_separation,
		gensym("separation"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		specFlux_tilde_class, 
        (t_method)specFlux_tilde_overlap,
		gensym("overlap"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		specFlux_tilde_class, 
        (t_method)specFlux_tilde_squaredDiff, 
		gensym("squaredDiff"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		specFlux_tilde_class, 
        (t_method)specFlux_tilde_normalize, 
		gensym("normalize"),
		A_DEFFLOAT,
		0
	);

    class_addmethod(
    	specFlux_tilde_class,
    	(t_method)specFlux_tilde_dsp,
    	gensym("dsp"),
    	0
    );
}