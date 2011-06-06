/*

std - Calculate standard deviation of a list of numbers.

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.


version 0.0.2, March 8, 2010

• 0.0.2 declares all functions except _setup static

*/

#include "m_pd.h"
#include <math.h>
#include <stdio.h>

static t_class *std_class;

typedef struct _std
{
    t_object x_obj;
    t_outlet *std;
    
} t_std;


/* ------------------------ std -------------------------------- */

static void std_calculate(t_std *x, t_symbol *s, int argc, t_atom *argv)
{
	float n, sum, mean, std, *input;
    int i;

	n = argc;
	
	if(n <= 1)
		error("std: too few elements in list.");
	else
	{
		// create local memory
		input = (float *)getbytes(0);
		input = (float *)t_resizebytes(input, 0, n*sizeof(float));

		for(i=0; i<n; i++)
			input[i] = atom_getfloat(argv+i);
	
		sum = 0.0;
	
		for(i=0; i<n; i++)
			sum += input[i];
	
		mean = sum/n;
		
		sum = 0.0;
		
		// center & square the data
		for(i=0; i<n; i++)
		{
			input[i] -= mean;
			input[i] *= input[i];
			sum += input[i];
		}
		
		std = sum/(n-1.0);
		std = sqrt(std);
			
		outlet_float(x->std, std);

		// free local memory
		t_freebytes(input, n*sizeof(float));
	}	
}

static void *std_new(void)
{	
    t_std *x = (t_std *)pd_new(std_class);
    x->std = outlet_new(&x->x_obj, &s_float);
    
    return (x);
}

void std_setup(void)
{
    std_class = 
    class_new(
    	gensym("std"),
    	(t_newmethod)std_new,
    	0,
        sizeof(t_std),
        CLASS_DEFAULT, 
		0
    );

	class_addlist(
		std_class,
		(t_method)std_calculate
	);
}