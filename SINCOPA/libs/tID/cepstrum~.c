/*

cepstrum~ - A Cepstral Analysis external.

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.


version 0.2.4, May 15, 2010

��0.2.4 added an ifndef M_PI for guaranteed windows compilation
��0.2.3 goes with the standard cepstrum algorithm, and no longer zeroes out the negative frequencies.
� 0.2.2 fixes a bug in the cepstrum computation. the imaginary buffer was not zeroed out before the final IFT. this caused a wild ripple in the upper coefficients at times.
� 0.2.1 adds a #define M_PI for windows compilation, and declares all functions except _setup static
� 0.2.0 implements mayer_realfft
� 0.1.8 added normalization option

*/

#include "m_pd.h"
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

static t_class *cepstrum_tilde_class;

typedef struct _cepstrum_tilde
{
    t_object x_obj;
    float sr;
    float n;
    int overlap;
	int normalize;
    float window;
    double last_dsp_time;
    t_sample *signal_R;
    float *hann;
    float x_f;
    t_outlet *x_featureList;

} t_cepstrum_tilde;



/* ------------------------ spectrum functions -------------------------------- */

static void cepstrum_tilde_hann(int n, t_sample *in, float *hann)
{	
	while (n--)
    {
    	*in = *in * *hann;

    	in++;
    	hann++;
    };
}

static void cepstrum_tilde_realfft_unpack(int n, int n_half, t_sample *input, t_sample *imag)
{
	int i, j;
		
	imag[0]=0.0;  // DC
	
	for(i=(n-1), j=1; i>n_half; i--, j++)
		imag[j] = input[i];
}

static void cepstrum_tilde_abs(int n, t_sample *in1, t_sample *in2)
{
	while (n--)
    {   
        *in1 = sqrt( (*in1 * *in1) + (*in2 * *in2) );
        in1++;
        in2++;
    };
}

static void cepstrum_tilde_normal(int n, t_sample *in)
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

static void cepstrum_tilde_log(int n, t_sample *spectrum)
{
	while (n--)
    {   
    // if to protect against log(0);
    	if(*spectrum==0)
    		*spectrum = 0;
    	else
	        *spectrum = log(*spectrum);
	        
        spectrum++;
    };
}

static void cepstrum_tilde_realifft(int n, int n_half, t_sample *real, t_sample *imag)
{
    int i, j;
    float n_recip;
    
	for(i=(n_half+1), j=(n_half-1); i<n; i++, j--)
		real[i] = imag[j];
			
    mayer_realifft(n, real);
    
    n_recip = 1.0/n;
    
// normalize by 1/N, since mayer apparently doesn't
	for(i=0; i<n; i++)
		real[i] *= n_recip;
}

/* ------------------------ end spectrum functions -------------------------------- */



/* ------------------------ cepstrum~ -------------------------------- */

static void cepstrum_tilde_bang(t_cepstrum_tilde *x)
{
    int i, j, window, window_half, bang_sample;
   	t_atom *listOut;
    t_sample *signal_R, *signal_I;
	double current_time;

    window = x->window;
    window_half = window*0.5;

   	signal_R = (t_sample *)getbytes(0);
	signal_I = (t_sample *)getbytes(0);
	listOut = (t_atom *)getbytes(0);
	signal_R = (t_sample *)t_resizebytes(signal_R, 0, window*sizeof(t_sample));
	signal_I = (t_sample *)t_resizebytes(signal_I, 0, window_half*sizeof(t_sample));
	listOut = (t_atom *)t_resizebytes(listOut, 0, ((int)x->window*0.5)*sizeof(t_atom));
	
	current_time = clock_gettimesince(x->last_dsp_time);
	bang_sample = (int)(((current_time/1000.0)*x->sr)+0.5); // round

	if (bang_sample < 0)
        bang_sample = 0;
    else if ( bang_sample >= x->n )
        bang_sample = x->n - 1;
            
	// construct analysis window using bang_sample as the end of the window
	for(i=0, j=bang_sample; i<window; i++, j++)
		signal_R[i] = x->signal_R[j];
	
	cepstrum_tilde_hann(window, signal_R, x->hann);
	mayer_realfft(window, signal_R);
	cepstrum_tilde_realfft_unpack(window, window_half, signal_R, signal_I);
	cepstrum_tilde_abs(window_half, signal_R, signal_I);
	
	if(x->normalize)
		cepstrum_tilde_normal(window_half, signal_R);
	
	cepstrum_tilde_log(window_half, signal_R);

	// fill out the negative frequencies
	for(i=0; i<window_half; i++)
		signal_I[i] = signal_R[i];
//		signal_I[i] = 0.0; // or we could zero them out...
		
	cepstrum_tilde_realifft(window, window_half, signal_R, signal_I);

				
	for(i=0; i<window_half; i++)
		SETFLOAT(listOut+i, signal_R[i]);

	outlet_list(x->x_featureList, 0, window_half, listOut);

	t_freebytes(signal_R, window * sizeof(t_sample));
	t_freebytes(signal_I, window_half * sizeof(t_sample));
	t_freebytes(listOut, window_half * sizeof(t_atom));
}

static void cepstrum_tilde_window(t_cepstrum_tilde *x, t_floatarg f)
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
		
		// re-calc hann window
		for(i=0; i<x->window; i++)
			x->hann[i] = 0.5 * (1 - cos(2*M_PI*i/x->window));

		// init signal buffer
		for(i=0; i<(x->window+x->n); i++)
			x->signal_R[i] = 0.0;
		
		post("window size: %i. overlap: %i. sampling rate: %i", (int)x->window, x->overlap, (int)(x->sr/x->overlap));
	}
}


static void cepstrum_tilde_overlap(t_cepstrum_tilde *x, t_floatarg f)
{	
	x->overlap = (int)f;
		
    post("overlap: %i", x->overlap);

}


static void cepstrum_tilde_normalize(t_cepstrum_tilde *x, t_floatarg norm)
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


static void *cepstrum_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
    t_cepstrum_tilde *x = (t_cepstrum_tilde *)pd_new(cepstrum_tilde_class);
	int i, isPow2;
	s=s;
	
	x->x_featureList = outlet_new(&x->x_obj, gensym("list"));

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

	x->sr = 44100;
	x->n = 64.0;
	x->overlap = 1;
	x->normalize = 0; // normalization doesn't seem to matter at all
	x->last_dsp_time = clock_getlogicaltime();

	x->signal_R = (t_sample *)t_resizebytes(x->signal_R, 0, (x->window+x->n) * sizeof(t_sample));
	x->hann = (float *)t_resizebytes(x->hann, 0, x->window * sizeof(float));
	
	// initialize hann window
 	for(i=0; i<x->window; i++)
 		x->hann[i] = 0.5 * (1 - cos(2*M_PI*i/x->window));

 	for(i=0; i<(x->window+x->n); i++)
		x->signal_R[i] = 0.0;
		
    post("cepstrum~: window size: %i", (int)x->window);
    
    return (x);
}


static t_int *cepstrum_tilde_perform(t_int *w)
{
    int i, n;

    t_cepstrum_tilde *x = (t_cepstrum_tilde *)(w[1]);

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


static void cepstrum_tilde_dsp(t_cepstrum_tilde *x, t_signal **sp)
{
	int i;
	
	dsp_add(
		cepstrum_tilde_perform,
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
			
    	post("cepstrum~: window size: %i. overlap: %i. sampling rate: %i, block size: %i", (int)x->window, x->overlap, (int)(x->sr/x->overlap), (int)x->n);
	};
};

static void cepstrum_tilde_free(t_cepstrum_tilde *x)
{	
	// free the input buffer memory
    t_freebytes(x->signal_R, (x->window+x->n)*sizeof(t_sample));

	// free the hann window memory
    t_freebytes(x->hann, x->window*sizeof(float));
};

void cepstrum_tilde_setup(void)
{
    cepstrum_tilde_class = 
    class_new(
    	gensym("cepstrum~"),
    	(t_newmethod)cepstrum_tilde_new,
    	(t_method)cepstrum_tilde_free,
        sizeof(t_cepstrum_tilde),
        CLASS_DEFAULT, 
        A_GIMME,
		0
    );

    CLASS_MAINSIGNALIN(cepstrum_tilde_class, t_cepstrum_tilde, x_f);

	class_addbang(cepstrum_tilde_class, cepstrum_tilde_bang);
	
	class_addmethod(
		cepstrum_tilde_class, 
        (t_method)cepstrum_tilde_window,
		gensym("window"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		cepstrum_tilde_class, 
        (t_method)cepstrum_tilde_overlap,
		gensym("overlap"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		cepstrum_tilde_class, 
        (t_method)cepstrum_tilde_normalize, 
		gensym("normalize"),
		A_DEFFLOAT,
		0
	);
	
    class_addmethod(
    	cepstrum_tilde_class,
    	(t_method)cepstrum_tilde_dsp,
    	gensym("dsp"),
    	0
    );
}