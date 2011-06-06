/*

specIrregularity~

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

static t_class *specIrregularity_tilde_class;

typedef struct _specIrregularity_tilde
{
    t_object x_obj;
    float sr;
    float n;
    int overlap;
    int normalize;
    int window;
    int algorithm;
    double last_dsp_time;
    t_sample *signal_R;
    float *hann;
    float x_f;
    t_outlet *x_irregularity;

} t_specIrregularity_tilde;



/* ------------------------ spectrum functions -------------------------------- */

static void specIrregularity_tilde_hann(int n, t_sample *in, float *hann)
{	
	while (n--)
    {
    	*in = *in * *hann;

    	in++;
    	hann++;
    };
}

static void specIrregularity_tilde_realfft_unpack(int n, int n_half, t_sample *input, t_sample *imag)
{
	int i, j;
		
	imag[0]=0.0;  // DC
	
	for(i=(n-1), j=1; i>n_half; i--, j++)
		imag[j] = input[i];
}

static void specIrregularity_tilde_abs(int n, t_sample *in1, t_sample *in2)
{
	while (n--)
    {   
        *in1 = sqrt( (*in1 * *in1) + (*in2 * *in2) );
        in1++;
        in2++;
    };
}

static void specIrregularity_tilde_normal(int n, t_sample *in)
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



/* ------------------------ specIrregularity~ -------------------------------- */

static void specIrregularity_tilde_bang(t_specIrregularity_tilde *x)
{
    int i, j, window, window_half, bang_sample;
    float divisor, irregularity;
    t_sample *signal_R, *signal_I;
	double current_time;

    window = x->window;
    window_half = window*0.5;

	// create local memory
	signal_R = (t_sample *)getbytes(0);
	signal_I = (t_sample *)getbytes(0);
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
	
	specIrregularity_tilde_hann(window, signal_R, x->hann);
	mayer_realfft(window, signal_R);
	specIrregularity_tilde_realfft_unpack(window, window_half, signal_R, signal_I);
	specIrregularity_tilde_abs(window_half, signal_R, signal_I);
	
// 	if(x->normalize)
// 		specIrregularity_tilde_normal(window_half, signal_R);

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
				local_avg += signal_R[i-1+j];
			
			local_avg *= 0.333333333333;
				
			irregularity += fabs(signal_R[i] - local_avg);
			//irregularity = log10(irregularity);
		}
	}
	else
	{
		// Jensen
		for(i=0; i<window_half; i++)
		{
			if(i==window_half-1)
				irregularity += signal_R[i] * signal_R[i];
			else
				irregularity += pow(signal_R[i] - signal_R[i+1], 2);
			
			divisor += signal_R[i] * signal_R[i];
		}
		
		if(divisor==0)
			divisor = 1.0;
			
		irregularity /= divisor;
	}
		
	outlet_float(x->x_irregularity, irregularity);

	// free local memory
	t_freebytes(signal_R, window * sizeof(t_sample));
	t_freebytes(signal_I, window_half * sizeof(t_sample));
}


static void specIrregularity_tilde_window(t_specIrregularity_tilde *x, t_floatarg f)
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
					
		post("window size: %i. overlap: %i. sampling rate: %i", (int)x->window, x->overlap, (int)(x->sr/x->overlap));
	}
}


static void specIrregularity_tilde_overlap(t_specIrregularity_tilde *x, t_floatarg f)
{
	x->overlap = (int)f;

    post("overlap: %i", x->overlap);

}


static void specIrregularity_tilde_normalize(t_specIrregularity_tilde *x, t_floatarg norm)
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

static void specIrregularity_tilde_algorithm(t_specIrregularity_tilde *x, t_floatarg a)
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

static void *specIrregularity_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
    t_specIrregularity_tilde *x = (t_specIrregularity_tilde *)pd_new(specIrregularity_tilde_class);
	int i, isPow2;
	s=s;
	
	x->x_irregularity = outlet_new(&x->x_obj, &s_float);

	x->signal_R = (t_sample *)getbytes(0);
	x->hann = (float *)getbytes(0);
	
	if(argc > 0)
	{
		x->window = atom_getfloat(argv);
		isPow2 = (int)x->window && !( ((int)x->window-1) & (int)x->window );
		
		if(!isPow2)
		{
			error("requested window size is not a power of 2. default value of 1024 used instead");
			x->window = 1024;
		};
	}
	else
		x->window = 1024;	
	
	
	x->sr = 44100.0;
	x->n = 64.0;
	x->overlap = 1;
	x->normalize = 1;
	x->last_dsp_time = clock_getlogicaltime();

	x->signal_R = (t_sample *)t_resizebytes(x->signal_R, 0, (x->window+x->n) * sizeof(t_sample));
	x->hann = (float *)t_resizebytes(x->hann, 0, x->window * sizeof(float));
	
	// initialize hann window
 	for(i=0; i<x->window; i++)
 		x->hann[i] = 0.5 * (1 - cos(2*M_PI*i/x->window));

 	for(i=0; i<(x->window+x->n); i++)
		x->signal_R[i] = 0.0;
		
    post("specIrregularity~: window size: %i", (int)x->window);
    
    return (x);
}


static t_int *specIrregularity_tilde_perform(t_int *w)
{
    int i, n;

    t_specIrregularity_tilde *x = (t_specIrregularity_tilde *)(w[1]);

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


static void specIrregularity_tilde_dsp(t_specIrregularity_tilde *x, t_signal **sp)
{
	int i;
	
	dsp_add(
		specIrregularity_tilde_perform,
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
			
    	post("specIrregularity~: window size: %i. overlap: %i. sampling rate: %i, block size: %i", (int)x->window, x->overlap, (int)(x->sr/x->overlap), (int)x->n);
	};
};

static void specIrregularity_tilde_free(t_specIrregularity_tilde *x)
{
	// free the input buffer memory
    t_freebytes(x->signal_R, (x->window+x->n)*sizeof(t_sample));

	// free the hann window memory
    t_freebytes(x->hann, x->window*sizeof(float));
}

void specIrregularity_tilde_setup(void)
{
    specIrregularity_tilde_class = 
    class_new(
    	gensym("specIrregularity~"),
    	(t_newmethod)specIrregularity_tilde_new,
    	(t_method)specIrregularity_tilde_free,
        sizeof(t_specIrregularity_tilde),
        CLASS_DEFAULT, 
        A_GIMME,
		0
    );

    CLASS_MAINSIGNALIN(specIrregularity_tilde_class, t_specIrregularity_tilde, x_f);

	class_addbang(specIrregularity_tilde_class, specIrregularity_tilde_bang);
	
	class_addmethod(
		specIrregularity_tilde_class, 
        (t_method)specIrregularity_tilde_window,
		gensym("window"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		specIrregularity_tilde_class, 
        (t_method)specIrregularity_tilde_overlap,
		gensym("overlap"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		specIrregularity_tilde_class, 
        (t_method)specIrregularity_tilde_normalize, 
		gensym("normalize"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		specIrregularity_tilde_class, 
        (t_method)specIrregularity_tilde_algorithm,
		gensym("algorithm"),
		A_DEFFLOAT,
		0
	);
	
    class_addmethod(
    	specIrregularity_tilde_class,
    	(t_method)specIrregularity_tilde_dsp,
    	gensym("dsp"),
    	0
    );
}