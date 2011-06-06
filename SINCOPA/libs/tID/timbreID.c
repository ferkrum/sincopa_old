/*

timbreID - A generic classification external.

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.


version 0.3.5, May 16, 2010

¥ part of the update acknowledging the change to one tID folder full of sources and binaries. This update also included filter power averaging in barkSpec, bfcc, and mfcc.
¥Êremoved cosine similarity distance metric
¥ 0.3.3 reflects the general package update, which includes specIrregularity(~), and fixes a bug in cepstrum(~). also, changed timbreID's default KNN setting to 1.  It should really only be anything different if clustering is done.
¥ 0.3.3 is based on the windows source, where the _read functions were changed to avoid using fdopen.  also corrected the use of a resizable array for semicolon in the _read_cluster_text function. also added static to all functions except _setup.
¥ 0.3.2 no changes to [timbreID]. Added RT and NRT spectral spread externs, plus mean and standard deviation externs for summarizing spectro-temporal features of various lengths
¥ 0.3.1G this gets around to adding the cluster read/write functions. F is cleaned up version of D. 0.3.1G does work with read/write clusters.
¥ 0.3.1F skips the different approach in E (didn't take that very far), and fixed D to work.  The problem in _compute_cluster was after clustering and in memory cleanup. Shifting elements of x->cluster_members and x->cluster_member_lengths was done in the wrong order and memory was resized incorrectly. F cleans up D, and adds the cluster_write functions.
¥ 0.3.1D all that remains is cluster write/read.
¥ 0.3.1C extends 0.3.1B to shoot for making x->cluster_members memory dynamically allocated as well.  C works up to and including unclustering - next will be the read/write functions.
¥ 0.3.1 finally makes feature database and cluster member memory dynamically allocated. Set up the _read and _write functions so that the header includes the entire instance_feature_lengths array. Might as well keep all that information so that mixed instance databases won't be a problem in the future. BUT - write/read_text do not store/retrieve this information. More importantly, x->feature_length is simply set to the first element in x->instance_feature_lengths in the _read function. So, the repercussions of mixed length features in the database has not been thought through. A common feature length is still assumed.
¥ 0.3.0 adds a method for finding the worst match, and two others for outputting the maxes and mins of all features in the database (as long as it has been normalized). Also - found a bug in the concat_id chain. timbreID_clear never set x->neighborhood back to 0, which caused crashes after training/clearing/training a few times.
¥ 0.2.9 removes some warnings about unitialized variables. the entire timbreID set of externs was updated on this date as well to include non real-time versions of all feature externs. 0.3.0 will finally implement dynamic memory allocation for feature database instances, rather than using MAXFEATURELENGTH for all (very wasteful).
¥ 0.2.8 feature externs now use mayer_realfft, switched to timbreID-help.pd rather than help-timbreID.pd for the help patch name.
¥Ê0.2.7C Added cluster_membership function
¥ 0.2.7 Added ability to output a instance's feature list from the 4th outlet.
¥ 0.2.6 Added ability to write and read cluster instances in binary and text formats.
¥ 0.2.5 Added correlation as a distance metric. Fixed a bug in compute_cluster: it now throws an error if you ask for 0 clusters.  Also slightly tweaked the clustering algorithm itself.
¥ 0.2.4 Fixed major bug where timbreID_train wasn't updating the size of x->feature_input if the feature length changes. This caused crashes in timbreID_id after training with a feature length different than the default of 47 points. Fixed bug where requesting compute_cluster with a different number of clusters when already clustered causes crash. Fixed error in _write: j was used uninitialized. Distance metric choice is now used in compute_cluster. Previously, it was squared_euclid only. Added method for outputting cluster member lists, or the entire cluster list. 
¥ 0.2.3 adds binary file output (.timid now default), and an option for text output.
¥ 0.2.2 adds MATLAB .mat file output.
¥ 0.2.0 adds cosine distance function (doesn't account for attribute normalization yet though!).
¥ 0.1.9 adds ARFF function to export data for WEKA.

¥ fixed bug where the distances of the final elements in a compute_order result were sometimes 
larger than INT_MAX, so that they were all index 0 (init value). Using FLT_MAX instead now, which 
should be big enough to be safe.

¥ Accounted for normalization flag in the compute_variance function (requiring a new timbreID_mean 
function as well).

*/

#include "m_pd.h"
#include <math.h>
#include <limits.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXCLUSTERMEMS 8192

static t_class *timbreID_class;
    	
typedef struct instance
{
    float *instance;
} t_instance;

typedef struct member
{
    int *member;
} t_member;

typedef struct knn_info
{
    float dist;
	float safe_dist;
    int idx;
    int cluster;
} t_knn_info;

typedef struct norm_data
{
    float max;
    float min;
	float denominator;
} t_norm_data;


typedef struct _timbreID
{
    t_object x_obj;
    t_instance *instances;
    t_member *cluster_members;
    int *cluster_member_lengths;
    t_knn_info *knn_dists_idxs;
    t_norm_data *norm_data;
	t_float *feature_input;
    int *instance_cluster_membership;
    int *instance_feature_lengths;
    int feature_length;
    int num_instances;
    int num_instr;
    int dist_metric;
    int k;
    int normalize;
    int relative_ordering;
    // must deal with resizing attribute_order and weights next...
    int *attribute_order;
    float *weights;
    
    int reorient_flag;
    int neighborhood;
    int search_center;
    int prev_match;
    int max_matches;
    float jump_prob;

    int attributelo;
    int attributehi;
    t_canvas *x_canvas;
    t_outlet *id;
    t_outlet *nearest_dist;
    t_outlet *confidence;
    t_outlet *x_orderList;
    t_outlet *x_featureList;
} t_timbreID;


/* ---------------- utility functions ---------------------- */

static void timbreID_sort_knn_info(int k, int num_instances, int prev_match, t_knn_info *list)
{
	int i, j, top_i, *top_matches;
	
	top_matches = (int *)getbytes(0);
	top_matches = (int *)t_resizebytes(top_matches, 0, k * sizeof(int));

	for(i=0; i<k; i++)
	{
		float max_best;
			
		max_best = FLT_MAX;
		top_i = 0;
			
		for(j=0; j<num_instances; j++)
			if(list[j].dist < max_best)
				if(list[j].idx != prev_match) // doesn't include previous match - this is good
				{
					max_best = list[j].dist;
					top_i = j;
				};
		
		list[top_i].dist = FLT_MAX;
	
		top_matches[i] = list[top_i].idx;
	}

	for(i=0; i<k; i++)
	{
		t_knn_info tmp;

		tmp = list[i];
		list[i] = list[top_matches[i]];
		list[top_matches[i]] = tmp;
	}

	// free memory
	t_freebytes(top_matches, k*sizeof(int));

	// now, the list passed to the function will have the first k elements in order, 
	// and these elements have the lowest distances in the whole list.
}

static void timbreID_sort_float(int n, float *list)
{
	int i, j, flag;
	
	for(i=0; i<n; i++)
	{
		flag = 0;
		
		for(j=0; j<(n-1); j++)
		{
			if(list[j] > list[j+1])
			{
				float tmp;
				
				flag = 1;

				tmp = list[j+1];
				list[j+1] = list[j];
				list[j] = tmp;
			}
		}

		if(flag==0)
			break;
	}
}

static float timbreID_mean(int num_rows, int column, t_instance *instances, int normal_flag, t_norm_data *norm_data)
{
	int i;
	float avg, min, denominator;
	
	avg=0;
	
	for(i=0; i<num_rows; i++)
	{
		if(instances[i].instance[column] == FLT_MAX)
			continue;
		else
		{
			if(normal_flag)
			{
				min = norm_data[column].min;
				denominator = norm_data[column].denominator;
			
				avg += (instances[i].instance[column] - min) * denominator;
			}
			else
				avg += instances[i].instance[column];
		}
	}
	
	avg /= num_rows;
	
	return(avg);
}

static float timbreID_squared_euclid(t_timbreID *x, float *v1, float *v2)
{
	int i;
	float min, max, norm_denominator, sum, dist;
	
	sum=dist=max=0;
	
	for(i=x->attributelo; i<= x->attributehi; i++)
	{
		if(x->normalize)
		{
			if(v1[x->attribute_order[i]] < x->norm_data[x->attribute_order[i]].min)
				min = v1[x->attribute_order[i]];
			else
				min = x->norm_data[x->attribute_order[i]].min;

			if(v1[x->attribute_order[i]] > x->norm_data[x->attribute_order[i]].max)
			{
				max = v1[x->attribute_order[i]];

				if(max <= min) // don't divide by zero
				{
					max=2.0;
					min=1.0;
				};
				
				norm_denominator = 1.0/(max-min);
			}
			else
			{					
				max = x->norm_data[x->attribute_order[i]].max;
				norm_denominator = x->norm_data[x->attribute_order[i]].denominator;
			}	
			
			dist = (  (v1[x->attribute_order[i]] - min) * norm_denominator  ) - (  (v2[x->attribute_order[i]] - min) * norm_denominator  );
			sum += dist*dist*x->weights[x->attribute_order[i]];
		}
		else
		{
			dist = v1[x->attribute_order[i]] - v2[x->attribute_order[i]];
			sum += dist*dist*x->weights[x->attribute_order[i]];
		}
	}
	
	return(sum);
}

static float timbreID_manhattan(t_timbreID *x, float *v1, float *v2)
{
	int i;
	float min, max, norm_denominator, sum, dist;
	
	sum=dist=0;

	for(i=x->attributelo; i<= x->attributehi; i++)
	{
		if(x->normalize)
		{
			if(v1[x->attribute_order[i]] < x->norm_data[x->attribute_order[i]].min)
				min = v1[x->attribute_order[i]];
			else
				min = x->norm_data[x->attribute_order[i]].min;

			if(v1[x->attribute_order[i]] > x->norm_data[x->attribute_order[i]].max)
			{
				max = v1[x->attribute_order[i]];

				if(max <= min) // don't divide by zero
				{
					max=2.0;
					min=1.0;
				};
				
				norm_denominator = 1.0/(max-min);
			}
			else
			{					
				max = x->norm_data[x->attribute_order[i]].max;
				norm_denominator = x->norm_data[x->attribute_order[i]].denominator;
			}	
			
			dist = (  (v1[x->attribute_order[i]] - min) * norm_denominator  ) - (  (v2[x->attribute_order[i]] - min) * norm_denominator  );
		}
		else
			dist = v1[x->attribute_order[i]] - v2[x->attribute_order[i]];

		sum += fabs(dist) * x->weights[x->attribute_order[i]];
	}
	
	return(sum);
}

static float timbreID_correlation(t_timbreID *x, float *v1, float *v2)
{
	int i, j, vecLen, vecLenM1;
	float min, max, norm_denominator, vecLenRecip, vecLenM1Recip;
	float mean1, mean2, std1, std2, covariance, correlation, *vec1Centered, *vec2Centered;

	mean1=mean2=std1=std2=covariance=correlation=0;
	min=max=norm_denominator=1;

	vecLen = x->attributehi - x->attributelo + 1;
	vecLenM1 = vecLen-1;
	if(vecLen <= 0) vecLen = 1;
	if(vecLenM1 <=0) vecLenM1 = 1;
	vecLenRecip = 1.0/(float)vecLen;
	vecLenM1Recip = 1.0/(float)vecLenM1;

	vec1Centered = (float *)getbytes(0);
	vec2Centered = (float *)getbytes(0);
	vec1Centered = (float *)t_resizebytes(vec1Centered, 0, vecLen * sizeof(float));
	vec2Centered = (float *)t_resizebytes(vec2Centered, 0, vecLen * sizeof(float));
	
	for(i=x->attributelo; i<= x->attributehi; i++)
	{
		if(x->normalize)
		{
			if(v1[x->attribute_order[i]] < x->norm_data[x->attribute_order[i]].min)
				min = v1[x->attribute_order[i]];
			else
				min = x->norm_data[x->attribute_order[i]].min;

			if(v1[x->attribute_order[i]] > x->norm_data[x->attribute_order[i]].max)
			{
				if(max <= min) // don't divide by zero
				{
					max=2.0;
					min=1.0;
				};
				
				max = v1[x->attribute_order[i]];
				norm_denominator = 1.0/(max-min);
			}
			else
			{					
				max = x->norm_data[x->attribute_order[i]].max;
				norm_denominator = x->norm_data[x->attribute_order[i]].denominator;
			}	
			
			mean1 += (v1[x->attribute_order[i]] - min) * norm_denominator;
			mean2 += (v2[x->attribute_order[i]] - min) * norm_denominator;
		}
		else
		{
			mean1 += v1[x->attribute_order[i]];
			mean2 += v2[x->attribute_order[i]];
		}
	};

	mean1 *= vecLenRecip;
	mean2 *= vecLenRecip;


	// min, max, and norm_denominator have already been established
	for(i=x->attributelo, j=0; i<= x->attributehi; i++, j++)
	{	
		if(x->normalize)
		{
			vec1Centered[j] = ((v1[x->attribute_order[i]] -  min) * norm_denominator) - mean1;
			vec2Centered[j] = ((v2[x->attribute_order[i]] -  min) * norm_denominator) - mean2;
		}
		else
		{
			vec1Centered[j] = v1[x->attribute_order[i]] - mean1;
			vec2Centered[j] = v2[x->attribute_order[i]] - mean2;
		}
	};

	for(i=0; i<vecLen; i++)
	{
		std1 += vec1Centered[i]*vec1Centered[i];
		std2 += vec2Centered[i]*vec2Centered[i];
	};
	
	std1 *= vecLenM1Recip;
	std2 *= vecLenM1Recip;

	// take sqrt to convert variance to standard deviation
	std1 = sqrt(std1);
	std2 = sqrt(std2);
	
	if(std1==0 || std2==0) std1=std2=0; // don't divide by zero below

	for(i=0; i<vecLen; i++)
		covariance += vec1Centered[i]*vec2Centered[i];

	// covariance usually averaged via N-1, not N
	covariance *= vecLenM1Recip;
	
	// Pearson correlation coefficient
	correlation = covariance/(std1*std2);

	// bash to the 0-2 range, then flip sign so that lower is better. this keeps things consistent with other distance metrics.
	correlation += 1;
	correlation *= -1;
	
	// normally in -1 to 1 range, where higher is more simlar
	// bash it to the 0 to 2 range, where lower is more similar
	//correlation = 2-(correlation+1);

	// free memory
	t_freebytes(vec1Centered, vecLen*sizeof(float));
	t_freebytes(vec2Centered, vecLen*sizeof(float));

	return(correlation);
}

/* ---------------- END utility functions ---------------------- */




/* ------------------------ timbreID -------------------------------- */

static void timbreID_train(t_timbreID *x, t_symbol *s, int argc, t_atom *argv)
{
	int i, instance_ind, list_length;

	instance_ind = x->num_instances;
	list_length = argc;
	s=s; // to get rid of 'unused variable' warning

	x->instances = (t_instance *)t_resizebytes(x->instances, x->num_instances * sizeof(t_instance), (x->num_instances+1) * sizeof(t_instance));
	x->instance_feature_lengths = (int *)t_resizebytes(x->instance_feature_lengths, x->num_instances * sizeof(int), (x->num_instances+1) * sizeof(int));
	x->knn_dists_idxs = (t_knn_info *)t_resizebytes(x->knn_dists_idxs, x->num_instances * sizeof(t_knn_info), (x->num_instances+1) * sizeof(t_knn_info));
	x->instance_cluster_membership = (int *)t_resizebytes(x->instance_cluster_membership, x->num_instances * sizeof(int), (x->num_instances+1) * sizeof(int));
	x->cluster_members = (t_member *)t_resizebytes(x->cluster_members, x->num_instances * sizeof(t_member), (x->num_instances+1) * sizeof(t_member));
	x->cluster_member_lengths = (int *)t_resizebytes(x->cluster_member_lengths, x->num_instances * sizeof(int), (x->num_instances+1) * sizeof(int));
	
	x->instance_cluster_membership[instance_ind] = instance_ind;
	x->instance_feature_lengths[instance_ind] = list_length;
	x->cluster_member_lengths[instance_ind] = 2; // 2 because we're unclustered to start, and each instance has a cluster with itself as a member, plus -1 as the 2nd element to terminate the list

	x->instances[instance_ind].instance = (float *)getbytes(0);
	x->instances[instance_ind].instance = (float *)t_resizebytes(x->instances[instance_ind].instance, 0, list_length * sizeof(float));

	x->cluster_members[instance_ind].member = (int *)getbytes(0);
	x->cluster_members[instance_ind].member = (int *)t_resizebytes(x->cluster_members[instance_ind].member, 0, 2 * sizeof(int));
	
	// init new cluster_members
	x->cluster_members[instance_ind].member[0] = instance_ind; // first member of the cluster is the instance index
	x->cluster_members[instance_ind].member[1] = -1;
	
	x->num_instances++;
	x->num_instr++;
	x->neighborhood++;
	
	if(x->feature_length != list_length)
	{
		x->feature_input = (t_float *)t_resizebytes(x->feature_input, x->feature_length * sizeof(t_float), list_length * sizeof(t_float));
		x->attribute_order = (int *)t_resizebytes(x->attribute_order, x->feature_length * sizeof(int), list_length * sizeof(int));
		x->weights = (float *)t_resizebytes(x->weights, x->feature_length * sizeof(float), list_length * sizeof(float));
		x->feature_length = list_length;

		// initialize attribute_order
		for(i=0; i<x->feature_length; i++)
			x->attribute_order[i] = i;
	
		// initialize weights
		for(i=0; i<x->feature_length; i++)
			x->weights[i] = 1.0;
		
		x->attributelo = 0;
		x->attributehi = x->feature_length-1;
		post("feature length: %i.", x->feature_length);
		post("attribute range: %i through %i.", x->attributelo, x->attributehi);
	};
			
	for(i=0; i<list_length; i++)
		x->instances[instance_ind].instance[i] = atom_getfloat(argv+i);
	
	outlet_float(x->id, instance_ind); // output received feedback here rather than a post (in case of hi-speed training)
}


static void timbreID_id(t_timbreID *x, t_symbol *s, int argc, t_atom *argv)
{
	float sum, best, second_best;
	float distance_output, confidence;
    int i, list_length, id, *votes, top_vote;
	s=s; // to get rid of 'unused variable' warning

	if(x->num_instances)
	{
		votes = (int *)getbytes(0);
		votes = (int *)t_resizebytes(votes, 0, x->num_instr * sizeof(int));
		
		// init votes to 0
		for(i=0; i<x->num_instr; i++)
			votes[i] = 0;
			
		// init cluster info to instance idx
		for(i=0; i<x->num_instances; i++)
			x->knn_dists_idxs[i].cluster = i;
			
		list_length = argc;
		distance_output = 0;
		confidence = 0;
		
		if(x->feature_length != list_length)
		{
			x->feature_input = (t_float *)t_resizebytes(x->feature_input, x->feature_length * sizeof(t_float), list_length * sizeof(t_float));
			x->feature_length = list_length;
			x->attributelo = 0;
			x->attributehi = x->feature_length-1;
			post("feature length: %i.", x->feature_length);
			post("attribute range: %i through %i.", x->attributelo, x->attributehi);
		};
			
		for(i=0; i<x->feature_length; i++)
			x->feature_input[i] = atom_getfloat(argv+i);		
			
		id = 0;
		best = FLT_MAX;
		
		for(i=0; i<x->num_instances; i++)
		{
			sum = 0;		

			switch(x->dist_metric)
			{
				case 0:
					sum = timbreID_squared_euclid(x, x->feature_input, x->instances[i].instance);
					break;
				case 1:
					sum = timbreID_manhattan(x, x->feature_input, x->instances[i].instance);
					break;
				case 2:
					sum = timbreID_correlation(x, x->feature_input, x->instances[i].instance);
					break;
				default:
					break;
			};
		
			x->knn_dists_idxs[i].dist = x->knn_dists_idxs[i].safe_dist = sum; // store the distance
			x->knn_dists_idxs[i].idx = i; // store the idx
		};
	
		// a reduced sort, so that the first k elements in knn_dists_idxs will be the lowest distances in the list, and in order to boot.	
		timbreID_sort_knn_info(x->k, x->num_instances, -1, x->knn_dists_idxs); // pass a prev_match value of -1, since it's unused here
			
		// store instance's cluster id
		for(i=0; i<x->k; i++)
			x->knn_dists_idxs[i].cluster = x->instance_cluster_membership[x->knn_dists_idxs[i].idx];
		
		// vote
		for(i=0; i<x->k; i++)
			votes[x->knn_dists_idxs[i].cluster]++;
	
		top_vote = -1;
		for(i=0; i<x->num_instr; i++)
			if(votes[i] > top_vote)
			{
				top_vote = votes[i];
				id = i; // store cluster id of winner
			};
		
		// in case of a tie, pick the shortest distance
		if(top_vote <= (x->k*0.5))
			id = x->knn_dists_idxs[0].cluster;
	
		for(i=0; i<x->k; i++)
			if(x->knn_dists_idxs[i].cluster==id)
			{
				best = x->knn_dists_idxs[i].safe_dist;
				break;
			};
		
		second_best = FLT_MAX;
		
		for(i=0; i<x->k; i++)
			if(x->knn_dists_idxs[i].cluster!=id)
			{
				second_best = x->knn_dists_idxs[i].safe_dist;
				break;
			};

		// if no second best assignment is made (because all K items belong to same cluster), make 2nd best the 2nd in list
		if(second_best==FLT_MAX)
			second_best = x->knn_dists_idxs[1].safe_dist;
		
		if(second_best<=0 || second_best==FLT_MAX)
			confidence = 0;
		else
			confidence = 1-(best/second_best);
		
		distance_output = best;
	
		// free memory
		t_freebytes(votes, x->num_instr*sizeof(int));
		
		outlet_float(x->confidence, confidence);	
		outlet_float(x->nearest_dist, distance_output);
		outlet_float(x->id, id);
    }
    else
    	error("timbreID: no training instances have been loaded. cannot perform ID.");
	
}


static void timbreID_worst_match(t_timbreID *x, t_symbol *s, int argc, t_atom *argv)
{
	float sum, worst;
	float distance_output;
    int i, list_length, id;
	s=s; // to get rid of 'unused variable' warning

	if(x->num_instances)
	{
		list_length = argc;
		distance_output = 0;
			
		if(x->feature_length != list_length)
		{
			x->feature_input = (t_float *)t_resizebytes(x->feature_input, x->feature_length * sizeof(t_float), list_length * sizeof(t_float));
			x->feature_length = list_length;
			x->attributelo = 0;
			x->attributehi = x->feature_length-1;
			post("feature length: %i.", x->feature_length);
			post("attribute range: %i through %i.", x->attributelo, x->attributehi);
		};
			
		for(i=0; i<argc; i++)
			x->feature_input[i] = atom_getfloat(argv+i);		
			
		id = 0;
		worst = 0;
		
		for(i=0; i<x->num_instances; i++)
		{
			sum = 0;		

			switch(x->dist_metric)
			{
				case 0:
					sum = timbreID_squared_euclid(x, x->feature_input, x->instances[i].instance);
					break;
				case 1:
					sum = timbreID_manhattan(x, x->feature_input, x->instances[i].instance);
					break;
				case 2:
					sum = timbreID_correlation(x, x->feature_input, x->instances[i].instance);
					break;
				default:
					break;
			};
		
			if(sum > worst)
			{
				worst = sum;
				id = i;
				post("updated worst: %i, %f", id, worst);
			}
		};
	
		id = x->instance_cluster_membership[id];
		
		distance_output = worst;
		
		outlet_float(x->confidence, 0);	
		outlet_float(x->nearest_dist, distance_output);
		outlet_float(x->id, id);
    }
    else
    	error("timbreID: no training instances have been loaded. cannot perform worst match.");
	
}


static void timbreID_forget(t_timbreID *x)
{	
	if(x->num_instances > 0)
	{
		// free the instance
		t_freebytes(x->instances[x->num_instances-1].instance, x->instance_feature_lengths[x->num_instances-1]*sizeof(float));		
		x->instances = (t_instance *)t_resizebytes(x->instances, x->num_instances * sizeof(t_instance), (x->num_instances-1) * sizeof(t_instance));
		
		// shrink instance_feature_lengths
		x->instance_feature_lengths = (int *)t_resizebytes(x->instance_feature_lengths, x->num_instances * sizeof(int), (x->num_instances-1) * sizeof(int));
		
		x->knn_dists_idxs = (t_knn_info *)t_resizebytes(x->knn_dists_idxs, x->num_instances * sizeof(t_knn_info), (x->num_instances-1) * sizeof(t_knn_info));
		x->instance_cluster_membership = (int *)t_resizebytes(x->instance_cluster_membership, x->num_instances * sizeof(int), (x->num_instances-1) * sizeof(int));

		// should probably do a double-check here that the database isn't clustered. if it is clustered, we need to find
		// which cluster the instance belongs to and remove that instance from the cluster's member list.
		
		// shrink cluster_members
		t_freebytes(x->cluster_members[x->num_instances-1].member, x->cluster_member_lengths[x->num_instances-1]*sizeof(int));		
		x->cluster_members = (t_member *)t_resizebytes(x->cluster_members, x->num_instances * sizeof(t_member), (x->num_instances-1) * sizeof(t_member));
		
		// shrink instance_feature_lengths
		x->cluster_member_lengths = (int *)t_resizebytes(x->cluster_member_lengths, x->num_instances * sizeof(int), (x->num_instances-1) * sizeof(int));
			
		x->num_instances--;
		x->num_instr--;
	
		post("forgot last instance. instances 0 through %i remain.", x->num_instances-1);
	}
	else
		error("timbreID: nothing to forget.");
}



//************************************* concatenative synthesis functions

static void timbreID_concat_id(t_timbreID *x, t_symbol *s, int argc, t_atom *argv)
{
	float sum, best;
	float distance_output;
    int i, list_length, half_neighborhood, mod_index, id;
	s=s; // to get rid of 'unused variable' warning

	if(x->num_instances)
	{
		// init cluster info to instance idx
		for(i=0; i<x->num_instances; i++)
			x->knn_dists_idxs[i].idx = i;
	
		list_length = argc;
		distance_output = 0;
		half_neighborhood = x->neighborhood/2;
		
		if(x->feature_length != list_length)
		{
			x->feature_input = (t_float *)t_resizebytes(x->feature_input, x->feature_length * sizeof(t_float), list_length * sizeof(t_float));
			x->feature_length = list_length;
			x->attributelo = 0;
			x->attributehi = x->feature_length-1;
			post("feature length: %i.", x->feature_length);
			post("attribute range: %i through %i.", x->attributelo, x->attributehi);
		};
			
		for(i=0; i<argc; i++)
			x->feature_input[i] = atom_getfloat(argv+i);		
			
		id = 0;
		best = FLT_MAX;
		
		for(i=(x->search_center - half_neighborhood); i<(x->search_center + half_neighborhood + 0); i++)
		{
			if( i<0 )
				mod_index = x->num_instances + i;  // wraps in reverse to end of table.  + i because i is neg.
			else
				mod_index = i%x->num_instances;

			sum = 0;
	
			switch(x->dist_metric)
			{
				case 0:
					sum = timbreID_squared_euclid(x, x->feature_input, x->instances[mod_index].instance);
					break;
				case 1:
					sum = timbreID_manhattan(x, x->feature_input, x->instances[mod_index].instance);
					break;
				case 2:
					sum = timbreID_correlation(x, x->feature_input, x->instances[mod_index].instance);
					break;
				default:
					break;
			};

				
			x->knn_dists_idxs[mod_index].dist = x->knn_dists_idxs[mod_index].safe_dist = sum; // store the distance
			x->knn_dists_idxs[mod_index].idx = mod_index; // store the idx
		};
	
		
		// a reduced sort, so that the first max_matches elements in knn_dists_idxs will be the lowest distances in the list, and in order to boot.
		// pass x->prev_match to make sure we don't output the same match two times in a row (to prevent one grain being played back several
		// times in sequence.
		timbreID_sort_knn_info(x->max_matches, x->num_instances, x->prev_match, &x->knn_dists_idxs[0]);

		if(x->prev_match == -1)
			id = x->knn_dists_idxs[0].idx;
		else
		{
			for(i=0, best=FLT_MAX; i<x->max_matches; i++)
			{	
				float dist;
				dist=sum=0.0;

				switch(x->dist_metric)
				{
					case 0:
						sum = timbreID_squared_euclid(x, x->instances[x->prev_match].instance, x->instances[x->knn_dists_idxs[i].idx].instance);
						break;
					case 1:
						sum = timbreID_manhattan(x, x->instances[x->prev_match].instance, x->instances[x->knn_dists_idxs[i].idx].instance);
						break;
					case 2:
						sum = timbreID_correlation(x, x->instances[x->prev_match].instance, x->instances[x->knn_dists_idxs[i].idx].instance);
						break;
					default:
						break;
				};

				if( sum < best)
				{
					best = sum;
					id = x->knn_dists_idxs[i].idx;
				};

				sum = 0;
			};
		};
		
		if( rand() < (RAND_MAX*(1-x->jump_prob)) )
		{
			if( x->reorient_flag == 1 )
				x->search_center = id;
		}
		else
		{
			id += rand();
			id = (int)id%x->num_instances;
			x->search_center = id;
		};
		

		x->prev_match = id;
		distance_output = best;
	
		outlet_float(x->nearest_dist, distance_output);
		outlet_float(x->id, id);
    }
    else
    	error("timbreID: no training instances have been loaded. cannot perform ID.");
	
}


static void timbreID_concat_neighborhood(t_timbreID *x, t_floatarg n)
{
	if( (int)n > x->num_instances)
		x->neighborhood = x->num_instances;
	else if( (int)n < 1 )
		x->neighborhood = 1;
	else
		x->neighborhood = (int)n;
}


static void timbreID_concat_jump_prob(t_timbreID *x, t_floatarg jp)
{
	if( jp<0 )
		x->jump_prob = 0;
	else if( jp>1 )
		x->jump_prob = 1;
	else
		x->jump_prob = jp;
}


static void timbreID_concat_reorient(t_timbreID *x, t_floatarg r)
{
	if( r<0 )
		x->reorient_flag = 0;
	else if( r>1 )
		x->reorient_flag = 1;
	else
		x->reorient_flag = (int)r;
}


static void timbreID_concat_search_center(t_timbreID *x, t_floatarg sc)
{
	if( (int)sc < 0 )
		x->search_center = 0;
	else if( (int)sc > x->num_instances )
		x->search_center = x->num_instances;
	else
		x->search_center = (int)sc;
}


static void timbreID_concat_max_matches(t_timbreID *x, t_floatarg mm)
{	
	if( (int)mm < 1 )
		x->max_matches = 1;
	else if( (int)mm > 50 )
		x->max_matches = 50;
	else
		x->max_matches = (int)mm;
}

//************************************* END concatenative synthesis functions




static void timbreID_knn(t_timbreID *x, t_floatarg k)
{		
	if(k<1.0)
		post("k must be greater than zero.");
	else if(k>x->num_instances)
		post("k must be less than the total number of instances.");
	else
	{
		x->k = (int)k;
		post("searching %i neighbors for KNN.", x->k);
	}
}


static void timbreID_normalize(t_timbreID *x, t_floatarg n)
{	
	int i, j;
	float *attribute_column;
	
	// create local memory
	attribute_column = (float *)getbytes(0);
	attribute_column = (float *)t_resizebytes(attribute_column, 0, x->num_instances * sizeof(float));
	
	if(n<=0)
	{
		// free memory	
		t_freebytes(x->norm_data, x->feature_length*sizeof(t_norm_data));
		
		x->normalize=0;
		post("feature attribute normalization OFF.");
	}
	else
	{
		if(x->num_instances)
		{
			// create memory;
			x->norm_data = (t_norm_data *)getbytes(0);
			x->norm_data = (t_norm_data *)t_resizebytes(x->norm_data, 0, x->feature_length * sizeof(t_norm_data));
	
			// j for columns (attributes), i for rows (instances)
			for(j=0; j<x->feature_length; j++)
			{
				for(i=0; i<x->num_instances; i++)
					attribute_column[i] = x->instances[i].instance[j];
		
				timbreID_sort_float(x->num_instances, &attribute_column[0]);
				
				x->norm_data[j].min = attribute_column[0];
				x->norm_data[j].max = attribute_column[x->num_instances-1];
				
				// don't divide by zero
				if(x->norm_data[j].max <= x->norm_data[j].min)
				{
					x->norm_data[j].max = 2.0;
					x->norm_data[j].min = 1.0;
				};
				
				x->norm_data[j].denominator = 1.0/(x->norm_data[j].max - x->norm_data[j].min);
			};
			
			x->normalize=1;
			post("feature attribute normalization ON.");
		}
		else
			error("timbreID: no training instances have been loaded. cannot calculate normalization terms.");
	}

	// free local memory	
	t_freebytes(attribute_column, x->num_instances*sizeof(float));

}


static void timbreID_manual_cluster(t_timbreID *x, t_floatarg num_instr, t_floatarg cluster_idx, t_floatarg low, t_floatarg hi)
{
	int i, j, cluster_idx_i, low_i, hi_i, num_members;

	cluster_idx_i = cluster_idx;
	low_i = low;
	hi_i = hi;
	num_members = hi_i - low_i + 1;

	if(x->num_instances < (int)num_instr)
		error("timbreID: not enough instances to cluster.");
	else
	{
		// only change memory size if x->num_instr hasn't been updated to be equal to num_instr
		if(x->num_instr != num_instr)
		{
			x->cluster_members = (t_member *)t_resizebytes(x->cluster_members, x->num_instr * sizeof(t_member), num_instr * sizeof(t_member));
			x->cluster_member_lengths = (int *)t_resizebytes(x->cluster_member_lengths, x->num_instr * sizeof(int), num_instr * sizeof(int));
		
			x->num_instr = num_instr;
		};

		// free the old memory for this cluster member list
		t_freebytes(x->cluster_members[cluster_idx_i].member, x->cluster_member_lengths[cluster_idx_i]*sizeof(int));

		// update the size of the list
		x->cluster_member_lengths[cluster_idx_i] = num_members+1; // +1 for the terminating -1

		x->cluster_members[cluster_idx_i].member = (int *)getbytes(0);
		x->cluster_members[cluster_idx_i].member = (int *)t_resizebytes(x->cluster_members[cluster_idx_i].member, 0, x->cluster_member_lengths[cluster_idx_i] * sizeof(int));
		

		for(i=low_i; i<=hi_i; i++)
			x->instance_cluster_membership[i] = cluster_idx_i;

		for(i=low_i, j=0; i<=hi_i; i++, j++)
			x->cluster_members[cluster_idx_i].member[j] = i;
	
		// terminate with -1
		x->cluster_members[cluster_idx_i].member[j] = -1;
			
		post("cluster %i contains instances %i through %i.", cluster_idx_i, low_i, hi_i);
	};
}


static void timbreID_compute_cluster(t_timbreID *x, t_floatarg num_instr)
{
	int i, j, k, numInstances, numInstancesM1, num_pairs, num_cluster_members1, num_cluster_members2, num_cluster_members_sum, num_clusters, *min_dist_idx;
	t_instance *cluster_data;
	float *pair_dists;
	float min_dist, num_cluster_members1_recip;
	t_atom *listOut;

	if(x->num_instances < num_instr)
		error("timbreID: not enough instances to cluster.");
	else if(x->num_instr != x->num_instances)
		error("timbreID: instances already clustered. uncluster first.");
	else if(num_instr == 0)
		error("timbreID: cannot create 0 clusters.");
	else
	{
	x->num_instr = num_instr;
	numInstances = x->num_instances;
	numInstancesM1 = numInstances-1;
	num_pairs = (numInstances*numInstancesM1) * 0.5;
	num_clusters = numInstances;
	num_cluster_members1 = 0;
	num_cluster_members2 = 0;
	num_cluster_members1_recip = 1;
	i=j=k=0;

	// create local memory
	min_dist_idx = (int *)getbytes(0);
	cluster_data = (t_instance *)getbytes(0);
	pair_dists = (float *)getbytes(0);
	listOut = (t_atom *)getbytes(0);

	min_dist_idx = (int *)t_resizebytes(min_dist_idx, 0, 2 * sizeof(int));
	cluster_data = (t_instance *)t_resizebytes(cluster_data, 0, numInstances * sizeof(t_instance));
	pair_dists = (float *)t_resizebytes(pair_dists, 0, num_pairs * sizeof(float));
	listOut = (t_atom *)t_resizebytes(listOut, 0, numInstances * sizeof(t_atom));
	
	for(i=0; i<numInstances; i++)
	{	
		x->cluster_members[i].member[0] = i; // first member of the cluster is the instance index
		x->cluster_members[i].member[1] = -1;
	}

	// copy x->instances into a safe local copy: cluster_data
	for(i=0; i<numInstances; i++)
	{
		cluster_data[i].instance = (float *)getbytes(0);
		cluster_data[i].instance = (float *)t_resizebytes(cluster_data[i].instance, 0, x->feature_length * sizeof(float));

		for(j=0; j<x->feature_length; j++)
			cluster_data[i].instance[j] = x->instances[i].instance[j];
	}

	
	while(num_clusters > x->num_instr)
	{		
		min_dist = FLT_MAX;
		
		// init min_dist_idx
		for(i=0; i<2; i++)
			min_dist_idx[i] = -1;
		
		// init pair distances 
		for(i=0; i<num_pairs; i++)
			pair_dists[i] = FLT_MAX;
			
		
		// get distances between all possible pairs in cluster_data
		for(i=0, k=0; i<numInstancesM1; i++)
		{
			if( cluster_data[i].instance[0] != -9999 ) // if this is true, the data hasn't been clustered yet.
			{
				for(j=1; j<numInstances; j++)
				{	
					if( (i+j) < numInstances )
					{
						if( cluster_data[i+j].instance[0] != -9999 )
						{
							switch(x->dist_metric)
							{
								case 0:
									pair_dists[k] = timbreID_squared_euclid(x, cluster_data[i].instance, cluster_data[i+j].instance);
									break;
								case 1:
									pair_dists[k] = timbreID_manhattan(x, cluster_data[i].instance, cluster_data[i+j].instance);
									break;
								case 2:
									pair_dists[k] = timbreID_correlation(x, cluster_data[i].instance, cluster_data[i+j].instance);
									break;
								default:
									break;
							};
	
							num_cluster_members1 = x->cluster_member_lengths[i]-1; // -1 because the list is terminated with -1
							num_cluster_members2 = x->cluster_member_lengths[i+j]-1;
								
							// definition of Ward's linkage from MATLAB linkage doc
							// pair_dists[k] is already squared euclidean distance
							
							num_cluster_members_sum = num_cluster_members1 + num_cluster_members2;
							
							if(num_cluster_members_sum > 0)
								pair_dists[k] = num_cluster_members1*num_cluster_members2 * (pair_dists[k]/(num_cluster_members1+num_cluster_members2));
							else
								pair_dists[k] = FLT_MAX;
							
							if(pair_dists[k]<min_dist)
							{
								min_dist=pair_dists[k];
								min_dist_idx[0]=i;
								min_dist_idx[1]=i+j;
							};
										
							k++; // increment pair_dists index if something was actually written to it.
						};
					}
					else
						break;
				}
			}
		};
		
		// we've found the smallest distance between cluster_data elements and stored it 
		// in min_dist. we've store the cluster_data indices of the two elements in 
		// min_dist_idx[0] and min_dist_idx[1].
		
		// set i to the index for storing the new member(s) of the cluster.
		i = x->cluster_member_lengths[min_dist_idx[0]]-1;
		
		// actually store the new member(s).
		j=0;
		while(x->cluster_members[min_dist_idx[1]].member[j] != -1)
		{
			// make some more memory for the new member(s)
			x->cluster_members[min_dist_idx[0]].member = (int *)t_resizebytes(x->cluster_members[min_dist_idx[0]].member, x->cluster_member_lengths[min_dist_idx[0]] * sizeof(int), (x->cluster_member_lengths[min_dist_idx[0]]+1) * sizeof(int));
			x->cluster_member_lengths[min_dist_idx[0]]++; // remember to update this member list's length

			x->cluster_members[min_dist_idx[0]].member[i++] = x->cluster_members[min_dist_idx[1]].member[j++];
		}

		i = x->cluster_member_lengths[min_dist_idx[0]]-1;
		x->cluster_members[min_dist_idx[0]].member[i] = -1; // terminate

		num_cluster_members1 = x->cluster_member_lengths[min_dist_idx[0]]-1;

		if(num_cluster_members1 > 0)
			num_cluster_members1_recip = 1.0/(float)num_cluster_members1;
		else
			num_cluster_members1_recip = 1.0;

		// resize the usurped cluster's cluster list memory, and update its size to 1
		x->cluster_members[min_dist_idx[1]].member = (int *)t_resizebytes(x->cluster_members[min_dist_idx[1]].member, x->cluster_member_lengths[min_dist_idx[1]] * sizeof(int), sizeof(int));
		x->cluster_members[min_dist_idx[1]].member[0] = -1;
		x->cluster_member_lengths[min_dist_idx[1]] = 1;

		// grab the first original instance for this cluster index
		for(i=0; i<x->feature_length; i++)
			cluster_data[min_dist_idx[0]].instance[i] = x->instances[min_dist_idx[0]].instance[i];
				
		// sum the original instances of the cluster members to compute centroid below
		for(i=1; i<num_cluster_members1; i++)
			for(j=0; j<x->feature_length; j++)
				cluster_data[min_dist_idx[0]].instance[j] += x->instances[  x->cluster_members[min_dist_idx[0]].member[i]  ].instance[j];

		// compute centroid
		for(i=0; i<x->feature_length; i++)
			cluster_data[min_dist_idx[0]].instance[i] *= num_cluster_members1_recip;

		// write -9999 to the first element in the nearest neighbor's instance to indicate it's now vacant.
		// this is all that's needed since all previous members were averaged and stored here.
		cluster_data[min_dist_idx[1]].instance[0] = -9999;
		
		num_clusters--;
	};

	// since the indices of the clusters have gaps from the process,
	// shift the cluster_members arrays that actually have content (!= -1)
	// to the head of cluster_members.  this will produce indices from 0 through num_instr-1.
	for(i=0, k=0; i<numInstances; i++)
		if( x->cluster_members[i].member[0] != -1)
		{
			// resize this member list
 			x->cluster_members[k].member = (int *)t_resizebytes(x->cluster_members[k].member, x->cluster_member_lengths[k] * sizeof(int), x->cluster_member_lengths[i] * sizeof(int));

			for(j=0; j<x->cluster_member_lengths[i]; j++)
				x->cluster_members[k].member[j] = x->cluster_members[i].member[j];
			
			// shift the list length info back
 			x->cluster_member_lengths[k] = x->cluster_member_lengths[i];
 						
			k++;
		};

	// free the excess cluster_members memory
	for(i=x->num_instr; i<numInstances; i++)
		t_freebytes(x->cluster_members[i].member, x->cluster_member_lengths[i]*sizeof(int));
		
	// resize cluster_members so it is only x->num_instr big
	x->cluster_members = (t_member *)t_resizebytes(x->cluster_members, numInstances * sizeof(t_member), x->num_instr * sizeof(t_member));

	// resize cluster_member_lengths so it is only x->num_instr big
	x->cluster_member_lengths = (int *)t_resizebytes(x->cluster_member_lengths, numInstances * sizeof(int), x->num_instr * sizeof(int));
		
	for(i=0, k=0; i<x->num_instr; i++)
		for(j=0; j<(x->cluster_member_lengths[i]-1); j++)
		{
			x->instance_cluster_membership[x->cluster_members[i].member[j]] = i;
			SETFLOAT(listOut+k, x->cluster_members[i].member[j]);
			k++;
		};
		
	outlet_list(x->x_orderList, 0, x->num_instances, listOut);
	
	// free memory
	t_freebytes(min_dist_idx, 2*sizeof(int));

	// free the database memory
	for(i=0; i<numInstances; i++)
		t_freebytes(cluster_data[i].instance, x->feature_length*sizeof(float));
		
	t_freebytes(cluster_data, numInstances*sizeof(t_instance));

	t_freebytes(pair_dists, num_pairs*sizeof(float));
	t_freebytes(listOut, numInstances*sizeof(t_atom));
		
	post("instances clustered.");

	} // end of main if/else	
}


static void timbreID_uncluster(t_timbreID *x)
{
	int i;

	// free each x->cluster_members list's memory
	for(i=0; i<x->num_instr; i++)
		t_freebytes(x->cluster_members[i].member, x->cluster_member_lengths[i]*sizeof(int));

	x->cluster_members = (t_member *)t_resizebytes(x->cluster_members, x->num_instr * sizeof(t_member), x->num_instances * sizeof(t_member));

	for(i=0; i<x->num_instances; i++)
	{
		x->cluster_members[i].member = (int *)getbytes(0);
		x->cluster_members[i].member = (int *)t_resizebytes(x->cluster_members[i].member, 0, 2 * sizeof(int));
	}

	// expand size of cluster_member_lengths again
	x->cluster_member_lengths = (int *)t_resizebytes(x->cluster_member_lengths, x->num_instr * sizeof(int), x->num_instances * sizeof(int));

	x->num_instr = x->num_instances;
	
	for(i=0; i<x->num_instances; i++)
	{
		x->instance_cluster_membership[i]=i; // init membership to index
		x->cluster_members[i].member[0]=i; // first member of the cluster is the instance index
		x->cluster_members[i].member[1] = -1;

		x->cluster_member_lengths[i] = 2;
	}
	
    post("instances unclustered.");
}


static void timbreID_compute_variance(t_timbreID *x)
{
	int i, j;
	float max, *attribute_var;
	t_instance *meanCentered;
	
	if(x->num_instances > 0)
	{
		// create local memory
		attribute_var = (float *)getbytes(0);
		meanCentered = (t_instance *)getbytes(0);
	
		attribute_var = (float *)t_resizebytes(attribute_var, 0, x->feature_length * sizeof(float));
		meanCentered = (t_instance *)t_resizebytes(meanCentered, 0, x->num_instances * sizeof(t_instance));
	
		for(i=0; i<x->num_instances; i++)
		{
			meanCentered[i].instance = (float *)getbytes(0);
			meanCentered[i].instance = (float *)t_resizebytes(meanCentered[i].instance, 0, x->feature_length * sizeof(float));
		}
		
		// init mean centered
		for(i=0; i<x->num_instances; i++)
			for(j=0; j<x->feature_length; j++)
				meanCentered[i].instance[j] = 0.0;
		
		// get the mean of each attribute
		// mean() checks for FLT_MAX and doesn't include those rows
		for(i=0; i<x->feature_length; i++)
			attribute_var[i] = timbreID_mean(x->num_instances, i, x->instances, x->normalize, x->norm_data);
		
		// center the data and write the matrix B
		for(i=0; i<x->num_instances; i++)
			for(j=0; j<x->feature_length; j++)
			{
				if(x->normalize)
					meanCentered[i].instance[j] = ((x->instances[i].instance[j] - x->norm_data[j].min) * x->norm_data[j].denominator) - attribute_var[j];
				else
					meanCentered[i].instance[j] = x->instances[i].instance[j] - attribute_var[j];				
			}
			
	// 	// variance is calculated as: sum(B(:,1).^2)/(M-1) for the first attribute
	// 	// run process by matrix columns rather than rows, hence the j, i order
		for(j=0; j<x->feature_length; j++)
		{		
			attribute_var[j] = 0;
			
			for(i=0; i<x->num_instances; i++)
				if(x->instances[i].instance[0] == FLT_MAX)
					continue;
				else
					attribute_var[j] += meanCentered[i].instance[j] * meanCentered[i].instance[j];
			
			if((x->num_instances-1) > 0)
				attribute_var[j] /= x->num_instances-1;
		}
	
	// 	for(i=0; i<x->feature_length; i++)
	// 		post("attribute variance %i: %f.", i, attribute_var[i]);
		
		// sort attribute_order by largest variances: find max in attribute_var,
		// replace it with -99999, find next max.
		for(i=0; i<x->feature_length; i++)
		{
			max=0.0;
			for(j=0; j<x->feature_length; j++)
			{   
				if(attribute_var[j] > max)
				{
					max = attribute_var[j];
					x->attribute_order[i] = j;
				}
			};
			
			attribute_var[x->attribute_order[i]] = -99999.0;
		};
		
		// free local memory
		t_freebytes(attribute_var, x->feature_length*sizeof(float));

		// free the meanCentered memory
		for(i=0; i<x->num_instances; i++)
			t_freebytes(meanCentered[i].instance, x->feature_length*sizeof(float));
		
		t_freebytes(meanCentered, x->num_instances*sizeof(t_instance));
	
		post("attributes ordered by variance.");
	}
	else
		error("timreID: no instances for variance computation.");

}


static void timbreID_clusters_list(t_timbreID *x)
{
	int i, j, k;
	t_atom *listOut;

	// create local memory
	listOut = (t_atom *)getbytes(0);
	listOut = (t_atom *)t_resizebytes(listOut, 0, x->num_instances * sizeof(t_atom));

	for(i=0, k=0; i<x->num_instr; i++)
		for(j=0; j<(x->cluster_member_lengths[i]-1); j++, k++) // -1 because it's terminated by -1
			SETFLOAT(listOut+k, x->cluster_members[i].member[j]);

	outlet_list(x->x_orderList, 0, x->num_instances, listOut);
	
	// free local memory
	t_freebytes(listOut, x->num_instances*sizeof(t_atom));
}


static void timbreID_cluster_list(t_timbreID *x, t_floatarg idx)
{
	int i, idx_i;
	t_atom *listOut;

	idx_i = (int)idx;

	if(idx_i >= x->num_instr || idx_i < 0)
		error("timbreID: cluster %i does not exist.", idx_i);
	else
	{
		// create local memory
		listOut = (t_atom *)getbytes(0);

		for(i=0; i<(x->cluster_member_lengths[idx_i]-1); i++)
		{
			listOut = (t_atom *)t_resizebytes(listOut, i * sizeof(t_atom), (i+1) * sizeof(t_atom));
			SETFLOAT(listOut+i, x->cluster_members[idx_i].member[i]);
		};

		outlet_list(x->x_orderList, 0, i, listOut);
		
		// free local memory
		t_freebytes(listOut, i*sizeof(t_atom));
	}
}


static void timbreID_cluster_membership(t_timbreID *x, t_floatarg idx)
{
	int idx_i;
	t_atom *listOut;

	idx_i = (int)idx;

	if(idx_i >= x->num_instances || idx_i < 0)
		error("timbreID: instance %i does not exist.", idx_i);
	else
	{
		// create local memory for a single element
		listOut = (t_atom *)getbytes(0);
		listOut = (t_atom *)t_resizebytes(listOut, 0, sizeof(t_atom));

		SETFLOAT(listOut, x->instance_cluster_membership[idx_i]);
		outlet_list(x->x_orderList, 0, 1, listOut);
		
		// free local memory
		t_freebytes(listOut, sizeof(t_atom));
	}
}


static void timbreID_compute_order(t_timbreID *x, t_floatarg reference)
{
	int i, j, smallIdx, ref;
	float smallest, sum;
	t_instance *instances;
	t_atom *listOut;

	// create local memory
	instances = (t_instance *)getbytes(0);
	listOut = (t_atom *)getbytes(0);

	instances = (t_instance *)t_resizebytes(instances, 0, x->num_instances * sizeof(t_instance));
	listOut = (t_atom *)t_resizebytes(listOut, 0, x->num_instances * sizeof(t_atom));

	for(i=0; i<x->num_instances; i++)
	{
		instances[i].instance = (float *)getbytes(0);
		instances[i].instance = (float *)t_resizebytes(instances[i].instance, 0, x->feature_length * sizeof(float));
	}
	
	if(reference >= x->num_instances)
		ref = x->num_instances-1;
	else if(reference < 0)
		ref = 0;
	else
	    ref = reference;

	// make a local copy of instances so they can be abused
    for(i=0; i<x->num_instances; i++)
    	for(j=0; j<x->feature_length; j++)
    		instances[i].instance[j] = x->instances[i].instance[j];
    		
    
    for(i=0; i<x->num_instances; i++)
    {
		smallest = FLT_MAX;
		smallIdx = 0;		
		
		for(j=0; j<x->num_instances; j++)
		{
			sum = 0;	

			// break out of this for iteration early if this instance slot has already been used.
			if(instances[j].instance[0] == FLT_MAX)
				continue;

			switch(x->dist_metric)
			{
				case 0:
					sum = timbreID_squared_euclid(x, x->instances[ref].instance, instances[j].instance);
					break;
				case 1:
					sum = timbreID_manhattan(x, x->instances[ref].instance, instances[j].instance);
					break;
				case 2:
					sum = timbreID_correlation(x, x->instances[ref].instance, instances[j].instance);
					break;
				default:
					break;
			};
	
			if(sum<smallest)
			{
				smallest = sum;
				smallIdx = j;
			};
			
		};

		SETFLOAT(listOut+i, smallIdx); // store the best from this round;

		if(x->relative_ordering)
			ref = smallIdx; // reorient search to nearest match;
		
		// set this instance to something huge so it will never be chosen as a good match
		for(j=0; j<x->feature_length; j++)
			instances[smallIdx].instance[j] = FLT_MAX;

	};

	outlet_list(x->x_orderList, 0, x->num_instances, listOut);

	// free local memory
	for(i=0; i<x->num_instances; i++)
		t_freebytes(instances[i].instance, x->feature_length*sizeof(float));
		
	t_freebytes(instances, x->num_instances*sizeof(t_instance));
	t_freebytes(listOut, x->num_instances*sizeof(t_atom));
}


static void timbreID_relative_ordering(t_timbreID *x, t_floatarg rel)
{
	if(rel<0)
		x->relative_ordering = 0;
	else if (rel>1)
		x->relative_ordering = 1;
	else
		x->relative_ordering = rel;
	
	if(x->relative_ordering)
		post("relative ordering ON.");
	else
		post("relative ordering OFF.");
}


static void timbreID_dist_metric(t_timbreID *x, t_floatarg f)
{		
	if((int)f < 0)
		x->dist_metric = 0;
	if((int)f > 3)
		x->dist_metric = 3;
	else
		x->dist_metric = (int)f;

	switch(x->dist_metric)
	{
		case 0:
			post("distance metric: EUCLIDEAN.");
			break;
		case 1:
			post("distance metric: MANHATTAN (taxicab distance).");
			break;
		case 2:
			post("distance metric: PEARSON CORRELATION COEFF.");
			break;
		default:
			break;
	};
}



static void timbreID_weights(t_timbreID *x, t_symbol *s, int argc, t_atom *argv)
{
	int i;
	s=s; // to get rid of 'unused variable' warning
	
	if(argc > x->feature_length)
	{
		post("WARNING: weights list longer than current feature length");
		argc = x->feature_length;
	}
	else
		post("%i weights received.", argc);

	
	for(i=0; i<argc; i++)
		x->weights[i] = atom_getfloat(argv+i);
		
	// if only the first few of a long feature vector are specified, fill in the rest with 1.0
	for(i=argc; i<x->feature_length; i++)
		x->weights[i] = 1.0;
}


static void timbreID_attributes(t_timbreID *x, t_symbol *s, int argc, t_atom *argv)
{
	int i;
	s=s; // to get rid of 'unused variable' warning
	
	if(argc > x->feature_length)
	{
		post("WARNING: attribute list longer than timbreID's current feature length");
		argc = x->feature_length;
	}
	else
		post("attribute list received.");

	
	for(i=0; i<argc; i++)
		x->attribute_order[i] = atom_getfloat(argv+i);
		
	// fill any remainder with attribute 0
	for(i=argc; i<x->feature_length; i++)
		x->attribute_order[i] = 0;
}


static void timbreID_attribute_range(t_timbreID *x, t_floatarg lo, t_floatarg hi)
{

	if( (int)lo < x->feature_length)
		x->attributelo = (int)lo;
	else
		x->attributelo = x->feature_length-1;
	
	
	if( (int)hi < x->feature_length)
		x->attributehi = (int)hi;
	else
		x->attributehi = x->feature_length-1;
	
	
	if( (int)lo > (int)hi)
		post("WARNING: low attribute > high attribute.  Correct this before further use.");
		
    post("attribute range: %i through %i.", x->attributelo, x->attributehi);
    	
}


static void timbreID_order_attributes(t_timbreID *x)
{
	int i;

	// initialize attribute_order
	for(i=0; i<x->feature_length; i++)
		x->attribute_order[i] = i;
		
	post("attribute order initialized.");
}


static void timbreID_print_instance(t_timbreID *x, t_floatarg e, t_floatarg f, t_floatarg g)
{
	int i;
   

	if( (int)e >= x->num_instances || (int)e < 0)
	{
		error("instance %i does not exist", (int)e);
	}
	else
	{
		post("T%i = [", (int)e);
	
		for(i=f; i<(g+1); i++)
		{
			
			if(i != g)
				post("%f, ", x->instances[(int)e].instance[i]);
			else
				post("%f", x->instances[(int)e].instance[i]);
		};
				
	
		post("]");
		post("");
    };
    	
}


static void timbreID_feature_list(t_timbreID *x, t_floatarg idx)
{
	int i, idx_i;
	t_atom *listOut;

	idx_i = (int)idx;

	if(idx_i >= x->num_instances || idx_i < 0)
		error("timbreID: instance %i does not exist.", idx_i);
	else
	{
		// create local memory
		listOut = (t_atom *)getbytes(0);
		listOut = (t_atom *)t_resizebytes(listOut, 0, x->feature_length * sizeof(t_atom));
		
		for(i=0; i<x->feature_length; i++)
		{
			if(x->normalize)
			{
				if( x->norm_data[i].max <= x->norm_data[i].min )
				{
					x->norm_data[i].max = 2.0;
					x->norm_data[i].min = 1.0;
				}
				
				SETFLOAT(listOut+i, (x->instances[idx_i].instance[i] - x->norm_data[i].min)/(x->norm_data[i].max - x->norm_data[i].min));	
			}
			else
				SETFLOAT(listOut+i, x->instances[idx_i].instance[i]);
		}
		
		outlet_list(x->x_featureList, 0, x->feature_length, listOut);
		
		// free local memory
		t_freebytes(listOut, x->feature_length*sizeof(t_atom));
	}
}


static void timbreID_max_values(t_timbreID *x)
{
	int i;
	t_atom *listOut;

	if(x->normalize)
	{
		// create local memory
		listOut = (t_atom *)getbytes(0);
		listOut = (t_atom *)t_resizebytes(listOut, 0, x->feature_length * sizeof(t_atom));
		
		for(i=0; i<x->feature_length; i++)
			SETFLOAT(listOut+i, x->norm_data[i].max);
		
		outlet_list(x->x_featureList, 0, x->feature_length, listOut);
		
		// free local memory
		t_freebytes(listOut, x->feature_length*sizeof(t_atom));
	}
	else
		error("timbreID: feature database not normalized yet");
}


static void timbreID_min_values(t_timbreID *x)
{
	int i;
	t_atom *listOut;

	if(x->normalize)
	{
		// create local memory
		listOut = (t_atom *)getbytes(0);
		listOut = (t_atom *)t_resizebytes(listOut, 0, x->feature_length * sizeof(t_atom));
		
		for(i=0; i<x->feature_length; i++)
			SETFLOAT(listOut+i, x->norm_data[i].min);
		
		outlet_list(x->x_featureList, 0, x->feature_length, listOut);
		
		// free local memory
		t_freebytes(listOut, x->feature_length*sizeof(t_atom));
	}
	else
		error("timbreID: feature database not normalized yet");
}


static void timbreID_clear(t_timbreID *x)
{
	int i;
	
	// free the database memory
	for(i=0; i<x->num_instances; i++)
		t_freebytes(x->instances[i].instance, x->instance_feature_lengths[i]*sizeof(float));
	
	x->instances = (t_instance *)t_resizebytes(x->instances, x->num_instances * sizeof(t_instance), 0);

	x->instance_feature_lengths = (int *)t_resizebytes(x->instance_feature_lengths, x->num_instances * sizeof(int), 0);
	
	x->knn_dists_idxs = (t_knn_info *)t_resizebytes(x->knn_dists_idxs, x->num_instances * sizeof(t_knn_info), 0);
	x->instance_cluster_membership = (int *)t_resizebytes(x->instance_cluster_membership, x->num_instances * sizeof(int), 0);
	x->cluster_members = (t_member *)t_resizebytes(x->cluster_members, x->num_instr * sizeof(t_member), 0);
	x->num_instances = 0;
	x->neighborhood = 0;
	x->num_instr = 0;
	
    post("all instances cleared.");
}


static void timbreID_write(t_timbreID *x, t_symbol *s)
{
	FILE *fd;
	int i, *header;
    float *fp;
    char *filename = s->s_name;
    char *buf;

	// create local memory
	buf = (char *)getbytes(0);
	buf = (char *)t_resizebytes(buf, 0, MAXPDSTRING * sizeof(char));

	header = (int *)getbytes(0);
	header = (int *)t_resizebytes(header, 0, (x->num_instances+1) * sizeof(int)); // record the size of each instance's feature, plus 1 for num_instances


    canvas_makefilename(x->x_canvas, filename, buf, MAXPDSTRING);

	fd = fopen(buf, "wb");
	
    if(!fd)
    {
        error("%s: couldn't create", buf);
        return;
    }		

	header[0] = x->num_instances;
	for(i=0; i<x->num_instances; i++)
		header[i+1] = x->instance_feature_lengths[i];
		
	fwrite(header, sizeof(int), x->num_instances+1, fd);

    for(i=0; i<x->num_instances; i++)
    {    
		fp = x->instances[i].instance;
		fwrite(fp, sizeof(float), x->instance_feature_lengths[i], fd);
   	};
   	
   	
    post("wrote %i instances to file: %s.", x->num_instances, buf);
    
    fclose(fd);

    // free memory
    t_freebytes(buf, MAXPDSTRING*sizeof(char));
    t_freebytes(header, (x->num_instances+1)*sizeof(int));
}


static void timbreID_read(t_timbreID *x, t_symbol *s)
{
	FILE *fd;
    float *fp;
	int i;
    char *filename = s->s_name;
    char *buf;

	// create local memory
	buf = (char *)getbytes(0);
	buf = (char *)t_resizebytes(buf, 0, MAXPDSTRING * sizeof(char));

    canvas_makefilename(x->x_canvas, filename, buf, MAXPDSTRING);

    // erase old instances & clusters and resize to 0.
    
	// free the database memory
	for(i=0; i<x->num_instances; i++)
		t_freebytes(x->instances[i].instance, x->instance_feature_lengths[i]*sizeof(float));
		
    x->instances = (t_instance *)t_resizebytes(x->instances, x->num_instances * sizeof(t_instance), 0);
	x->instance_feature_lengths = (int *)t_resizebytes(x->instance_feature_lengths, x->num_instances * sizeof(int), 0);
	x->knn_dists_idxs = (t_knn_info *)t_resizebytes(x->knn_dists_idxs, x->num_instances * sizeof(t_knn_info), 0);
	x->instance_cluster_membership = (int *)t_resizebytes(x->instance_cluster_membership, x->num_instances * sizeof(int), 0);

	for(i=0; i<x->num_instr; i++)
		t_freebytes(x->cluster_members[i].member, x->cluster_member_lengths[i]*sizeof(int));

	x->cluster_members = (t_member *)t_resizebytes(x->cluster_members, x->num_instr * sizeof(t_member), 0);
	x->cluster_member_lengths = (int *)t_resizebytes(x->cluster_member_lengths, x->num_instr * sizeof(int), 0);

	x->feature_input = (t_float *)t_resizebytes(x->feature_input, x->feature_length * sizeof(t_float), 0);   
	x->attribute_order = (int *)t_resizebytes(x->attribute_order, x->feature_length * sizeof(int), 0);
	x->weights = (float *)t_resizebytes(x->weights, x->feature_length * sizeof(float), 0);

    fd = fopen(buf, "rb");

    if (!fd)
    {
        post("%s: open failed", buf);
        return;
    }

	fread(&x->num_instances, sizeof(int), 1, fd);
	
	x->instance_feature_lengths = (int *)t_resizebytes(x->instance_feature_lengths, 0, x->num_instances * sizeof(int));

	fread(x->instance_feature_lengths, sizeof(int), x->num_instances, fd);
	
	// should search for the min and max instance sizes and store them both.
	// could resize relevant memory according to max, but read point limits according to min.
	// for now, just assume they're all the same.
	x->feature_length = x->instance_feature_lengths[0];	
	x->neighborhood = x->num_instances;
	x->num_instr = x->num_instances;
	x->attributelo = 0;
	x->attributehi = x->feature_length-1;
		
    // resize instances & cluster_members to num_instances
    x->instances = (t_instance *)t_resizebytes(x->instances, 0, x->num_instances*sizeof(t_instance));

	for(i=0; i<x->num_instances; i++)
	{
		x->instances[i].instance = (float *)getbytes(0);
		x->instances[i].instance = (float *)t_resizebytes(x->instances[i].instance, 0, x->instance_feature_lengths[i] * sizeof(float));
	}
	
	x->knn_dists_idxs = (t_knn_info *)t_resizebytes(x->knn_dists_idxs, 0, x->num_instances * sizeof(t_knn_info));
	x->instance_cluster_membership = (int *)t_resizebytes(x->instance_cluster_membership, 0, x->num_instances * sizeof(int));
	x->cluster_members = (t_member *)t_resizebytes(x->cluster_members, 0, x->num_instances * sizeof(t_member));
	x->cluster_member_lengths = (int *)t_resizebytes(x->cluster_member_lengths, 0, x->num_instances * sizeof(int));
	x->feature_input = (t_float *)t_resizebytes(x->feature_input, 0, x->feature_length * sizeof(t_float));
	x->attribute_order = (int *)t_resizebytes(x->attribute_order, 0, x->feature_length * sizeof(int));
	x->weights = (float *)t_resizebytes(x->weights, 0, x->feature_length * sizeof(float));

	// initialize attribute_order
	for(i=0; i<x->feature_length; i++)
		x->attribute_order[i] = i;

	// initialize weights
	for(i=0; i<x->feature_length; i++)
		x->weights[i] = 1.0;
			
	// initialize feature input buffer
	for(i=0; i<x->feature_length; i++)
		x->feature_input[i] = 0.0;
	
	for(i=0; i<x->num_instances; i++)
	{
		x->cluster_members[i].member = (int *)getbytes(0);
		x->cluster_members[i].member = (int *)t_resizebytes(x->cluster_members[i].member, 0, 2 * sizeof(int));

		x->instance_cluster_membership[i]=i; // init membership to index
		x->cluster_members[i].member[0]=i; // first member of the cluster is the instance index
		x->cluster_members[i].member[1] = -1;

		x->cluster_member_lengths[i] = 2;
	}


	// finally, read in the instance data
    for(i=0; i<x->num_instances; i++)
    {
		fp = x->instances[i].instance;
		fread(fp, sizeof(float), x->instance_feature_lengths[i], fd);
    };
    
    
    post("read %i instances from file: %s.\n", x->num_instances, buf);    
    post("feature length: %i.", x->feature_length);
	post("attribute range: %i through %i.", x->attributelo, x->attributehi);
    
    fclose(fd);

    // free memory
    t_freebytes(buf, MAXPDSTRING*sizeof(char));
}


static void timbreID_write_text(t_timbreID *x, t_symbol *s)
{
	FILE *fd;
    int i, j, *header;
    float *fp;
    char *filename = s->s_name;
    char *buf;

	// create local memory
	buf = (char *)getbytes(0);
	buf = (char *)t_resizebytes(buf, 0, MAXPDSTRING * sizeof(char));

	header = (int *)getbytes(0);
	header = (int *)t_resizebytes(header, 0, 2 * sizeof(int));
	
	j=0; // to keep track of no. of instances written.
    canvas_makefilename(x->x_canvas, filename, buf, MAXPDSTRING);

	fd = fopen(buf, "w");
	
    if(!fd)
    {
        error("%s: couldn't create", buf);
        return;
    }		

	// unlike the binary write/read, here we just assume a common feature length, and never look at instance_feature_lengths...for now.
	header[0] = x->num_instances;
	header[1] = x->feature_length;

	for(i=0; i<2; i++)
		fprintf(fd, "%i ", header[i]);

	fprintf(fd, "\n\n");

    for(i=0; i<x->num_instances; i++)
    {    
		fp = x->instances[i].instance;

		j=0;

		// only write actual values, not the FLT_MAX placeholders
		while(j<x->instance_feature_lengths[i])
		{
			fprintf(fd, "%6.20f ", *fp++);
			j++;
		};
		
		fprintf(fd, "\n\n");
   	};
    
    post("wrote %i instances to file: %s.", x->num_instances, buf);

    // free memory
    t_freebytes(buf, MAXPDSTRING*sizeof(char));
    t_freebytes(header, 2*sizeof(int));
    
    fclose(fd);
}


static void timbreID_read_text(t_timbreID *x, t_symbol *s)
{

    FILE *fd;

    int i, j, *header;

    float *fp;
    char *filename = s->s_name;
    char *buf;

	// create local memory
	buf = (char *)getbytes(0);
	buf = (char *)t_resizebytes(buf, 0, MAXPDSTRING * sizeof(char));

	header = (int *)getbytes(0);
	header = (int *)t_resizebytes(header, 0, 2 * sizeof(int));
	
    canvas_makefilename(x->x_canvas, filename, buf, MAXPDSTRING);

	// free the database memory
	for(i=0; i<x->num_instances; i++)
		t_freebytes(x->instances[i].instance, x->instance_feature_lengths[i]*sizeof(float));
		
    x->instances = (t_instance *)t_resizebytes(x->instances, x->num_instances * sizeof(t_instance), 0);
	x->instance_feature_lengths = (int *)t_resizebytes(x->instance_feature_lengths, x->num_instances * sizeof(int), 0);
	x->knn_dists_idxs = (t_knn_info *)t_resizebytes(x->knn_dists_idxs, x->num_instances * sizeof(t_knn_info), 0);
	x->instance_cluster_membership = (int *)t_resizebytes(x->instance_cluster_membership, x->num_instances * sizeof(int), 0);

	for(i=0; i<x->num_instr; i++)
		t_freebytes(x->cluster_members[i].member, x->cluster_member_lengths[i]*sizeof(int));

	x->cluster_members = (t_member *)t_resizebytes(x->cluster_members, x->num_instr * sizeof(t_member), 0);
	x->cluster_member_lengths = (int *)t_resizebytes(x->cluster_member_lengths, x->num_instr * sizeof(int), 0);

	x->feature_input = (t_float *)t_resizebytes(x->feature_input, x->feature_length * sizeof(t_float), 0);   
	x->attribute_order = (int *)t_resizebytes(x->attribute_order, x->feature_length * sizeof(int), 0);
	x->weights = (float *)t_resizebytes(x->weights, x->feature_length * sizeof(float), 0);

    fd = fopen(buf, "r");

    if (!fd)
    {
        post("%s: open failed", buf);
        return;
    }

	//// unlike the binary write/read, here we just assume a common feature length, and never look at instance_feature_lengths...for now.
	for(i=0; i<2; i++, header++)
		fscanf(fd, "%i", header);
	
	// in reverse order due to ptr arithmetic
	x->feature_length = *(--header);
	x->num_instances = *(--header);

	x->neighborhood = x->num_instances;
	x->num_instr = x->num_instances;
	x->attributelo = 0;
	x->attributehi = x->feature_length-1;
		
    // resize instances & cluster_members to num_instances
    x->instances = (t_instance *)t_resizebytes(x->instances, 0, x->num_instances*sizeof(t_instance));

	for(i=0; i<x->num_instances; i++)
	{
		x->instances[i].instance = (float *)getbytes(0);
		x->instances[i].instance = (float *)t_resizebytes(x->instances[i].instance, 0, x->feature_length * sizeof(float));
	}

	x->instance_feature_lengths = (int *)t_resizebytes(x->instance_feature_lengths, 0, x->num_instances * sizeof(int));

	x->knn_dists_idxs = (t_knn_info *)t_resizebytes(x->knn_dists_idxs, 0, x->num_instances * sizeof(t_knn_info));
	x->instance_cluster_membership = (int *)t_resizebytes(x->instance_cluster_membership, 0, x->num_instances * sizeof(int));
	x->cluster_members = (t_member *)t_resizebytes(x->cluster_members, 0, x->num_instances * sizeof(t_member));
	x->cluster_member_lengths = (int *)t_resizebytes(x->cluster_member_lengths, 0, x->num_instances * sizeof(int));
	x->feature_input = (t_float *)t_resizebytes(x->feature_input, 0, x->feature_length * sizeof(t_float));
	x->attribute_order = (int *)t_resizebytes(x->attribute_order, 0, x->feature_length * sizeof(int));
	x->weights = (float *)t_resizebytes(x->weights, 0, x->feature_length * sizeof(float));

	// initialize instance sizes
	for(i=0; i<x->num_instances; i++)
		x->instance_feature_lengths[i] = x->feature_length;
		
	// initialize attribute_order
	for(i=0; i<x->feature_length; i++)
		x->attribute_order[i] = i;

	// initialize weights
	for(i=0; i<x->feature_length; i++)
		x->weights[i] = 1.0;
		
	// initialize feature input buffer
	for(i=0; i<x->feature_length; i++)
		x->feature_input[i] = 0.0;

	for(i=0; i<x->num_instances; i++)
	{
		x->cluster_members[i].member = (int *)getbytes(0);
		x->cluster_members[i].member = (int *)t_resizebytes(x->cluster_members[i].member, 0, 2 * sizeof(int));

		x->instance_cluster_membership[i]=i; // init membership to index
		x->cluster_members[i].member[0]=i; // first member of the cluster is the instance index
		x->cluster_members[i].member[1] = -1;

		x->cluster_member_lengths[i] = 2;
	}


    for(i=0; i<x->num_instances; i++)
    {
		fp = x->instances[i].instance;

		for(j=0; j<x->feature_length; j++, fp++)
			fscanf(fd, "%f", fp);
    };
    	
    post("read %i instances from file: %s.\n", x->num_instances, buf);    
    post("feature length: %i.", x->feature_length);
	post("attribute range: %i through %i.", x->attributelo, x->attributehi);
    
    fclose(fd);
    
    // free memory
    t_freebytes(buf, MAXPDSTRING*sizeof(char));
    t_freebytes(header, 2*sizeof(int));
}


static void timbreID_ARFF(t_timbreID *x, t_symbol *s, int argc, t_atom *argv)
{
	FILE *fd;
    int i, j, features_written, att_range_low, att_range_hi;
    float *fp;
	t_symbol *filename_symbol, *relation_symbol, *att_symbol;
    char *buf, *filename, *relation, *att_name;
	
	s=s;
	
	att_range_low = 0;
	att_range_hi = -1;
	att_symbol = 0;
	att_name = 0;

	filename_symbol = atom_getsymbol(argv);
	filename = filename_symbol->s_name;

	// create local memory
	buf = (char *)getbytes(0);
	buf = (char *)t_resizebytes(buf, 0, MAXPDSTRING * sizeof(char));

    canvas_makefilename(x->x_canvas, filename, buf, MAXPDSTRING);

	fd = fopen(buf, "w");
	
    if(!fd)
    {
        error("%s: couldn't create", buf);
        return;
    }		

	if(argc>1)
	{
		relation_symbol = atom_getsymbol(argv+1);
		relation = relation_symbol->s_name;
	}
	else
		relation = "relation";

	fprintf(fd, "@RELATION %s\n\n\n", relation);

	if(argc>2)
	{
	    for(i=2; i<argc; i++)
		{

			switch((i-2)%3)
			{
				case 0:
					att_range_low = atom_getfloat(argv+i);
					break;
				case 1:
					att_range_hi = atom_getfloat(argv+i);
					break;
				case 2:
					att_symbol = atom_getsymbol(argv+i);
					att_name = att_symbol->s_name;
					for(j=0; j<=att_range_hi-att_range_low; j++)
						fprintf(fd, "@ATTRIBUTE %s-%i NUMERIC\n", att_name, j);
					break;
				default:
					break;
			}

		}

		// in case the argument list was incomplete
		for(j=0; j<(x->feature_length-1-att_range_hi); j++)
			fprintf(fd, "@ATTRIBUTE undefined-attribute-%i NUMERIC\n", j);

	}
	else
	{
		for(i=0; i<x->feature_length; i++)
			fprintf(fd, "@ATTRIBUTE undefined-attribute-%i NUMERIC\n", i);
	}


	fprintf(fd, "\n\n");
	fprintf(fd, "@DATA\n\n");

    for(i=0; i<x->num_instances; i++)
    {    
		fp = x->instances[i].instance;

		features_written=0; // to keep track of each instances no. of features written.

		while(1)
		{
			if(features_written++ == (x->instance_feature_lengths[i]-1))
			{
				fprintf(fd, "%6.20f", *fp++);
				break;
			}
			else
				fprintf(fd, "%6.20f, ", *fp++);
		};
		
		fprintf(fd, "\n");
   	};
    
    post("wrote %i instances to file: %s.", x->num_instances, buf);
    
    fclose(fd);

    // free memory
    t_freebytes(buf, MAXPDSTRING*sizeof(char));
}


static void timbreID_MATLAB(t_timbreID *x, t_symbol *file_symbol, t_symbol *var_symbol)
{
	FILE *fd;
    int i, features_written;
    float *fp;
    char *buf, *filename, *varname;

	filename = file_symbol->s_name;
	varname = var_symbol->s_name;

	// create local memory
	buf = (char *)getbytes(0);
	buf = (char *)t_resizebytes(buf, 0, MAXPDSTRING * sizeof(char));

    canvas_makefilename(x->x_canvas, filename, buf, MAXPDSTRING);

	fd = fopen(buf, "w");
	
    if(!fd)
    {
        error("%s: couldn't create", buf);
        return;
    }		

	fprintf(fd, "%% name: %s\n", varname);
	fprintf(fd, "%% type: matrix\n");
	fprintf(fd, "%% rows: %i\n", x->num_instances);
	fprintf(fd, "%% columns: %i\n\n", x->feature_length);
		
    for(i=0; i<x->num_instances; i++)
    {    
		fp = x->instances[i].instance;

		features_written=0; // to keep track of each instances no. of features written.

		while(1)
		{
			if(features_written++ == (x->instance_feature_lengths[i]-1))
			{
				fprintf(fd, "%6.20f", *fp++);
				break;
			}
			else
				fprintf(fd, "%6.20f, ", *fp++);
		};
		
		fprintf(fd, "\n");
   	};
    
    post("wrote %i instances to file: %s.", x->num_instances, buf);
    
    fclose(fd);

    // free memory
    t_freebytes(buf, MAXPDSTRING*sizeof(char));
}


static void timbreID_write_clusters(t_timbreID *x, t_symbol *s)
{
	FILE *fd;
	int i, *fp;
    char *filename = s->s_name;
    char *buf;

	// create local memory
	buf = (char *)getbytes(0);
	buf = (char *)t_resizebytes(buf, 0, MAXPDSTRING * sizeof(char));

    canvas_makefilename(x->x_canvas, filename, buf, MAXPDSTRING);

	fd = fopen(buf, "wb");
	
    if(!fd)
    {
        error("%s: couldn't create", buf);
        return;
    }		

	// write a header indicating the number of clusters (x->num_instr)
	fwrite(&x->num_instr, sizeof(int), 1, fd);

	fp = x->instance_cluster_membership;
    fwrite(fp, sizeof(int), x->num_instances, fd);

	fp = x->cluster_member_lengths;
    fwrite(fp, sizeof(int), x->num_instr, fd);

    for(i=0; i<x->num_instr; i++)
    {
		fp = x->cluster_members[i].member;
		fwrite(fp, sizeof(int), x->cluster_member_lengths[i], fd);
   	};
   	
   	
    post("wrote %i clusters to file: %s.", x->num_instr, buf);
    
    fclose(fd);

    // free memory
    t_freebytes(buf, MAXPDSTRING*sizeof(char));
}


static void timbreID_read_clusters(t_timbreID *x, t_symbol *s)
{
	FILE *fd;
	int i, *fp, header;
    char *filename = s->s_name;
    char *buf;

	// create local memory
	buf = (char *)getbytes(0);
	buf = (char *)t_resizebytes(buf, 0, MAXPDSTRING * sizeof(char));

    canvas_makefilename(x->x_canvas, filename, buf, MAXPDSTRING);	

	// free current x->cluster_members memory for each list
	for(i=0; i<x->num_instr; i++)
		t_freebytes(x->cluster_members[i].member, x->cluster_member_lengths[i]*sizeof(int));


	fd = fopen(buf, "rb");
	
    if(!fd)
    {
        error("%s: couldn't create", buf);
        return;
    }	

	// read header indicating number of instruments
	fread(&header, sizeof(int), 1, fd);

	x->cluster_member_lengths = (int *)t_resizebytes(x->cluster_member_lengths, x->num_instr * sizeof(int), header * sizeof(int));
	x->cluster_members = (t_member *)t_resizebytes(x->cluster_members, x->num_instr * sizeof(t_member), header * sizeof(t_member));

	x->num_instr = header;

	for(i=0; i<x->num_instr; i++)
		x->cluster_members[i].member = (int *)getbytes(0);

	fp = x->instance_cluster_membership;
    fread(fp, sizeof(int), x->num_instances, fd);

	fp = x->cluster_member_lengths;
    fread(fp, sizeof(int), x->num_instr, fd);

    for(i=0; i<x->num_instr; i++)
    {
		x->cluster_members[i].member = (int *)t_resizebytes(x->cluster_members[i].member, 0, x->cluster_member_lengths[i] * sizeof(int));

		fp = x->cluster_members[i].member;
		fread(fp, sizeof(int), x->cluster_member_lengths[i], fd);
    };
    
    post("read %i clusters from file: %s.\n", x->num_instr, buf);    

	// send the cluster list out the 4th outlet
	timbreID_clusters_list(x);


    fclose(fd);
    
    // free memory
    t_freebytes(buf, MAXPDSTRING*sizeof(char));
}


static void timbreID_write_clusters_text(t_timbreID *x, t_symbol *s)
{
	FILE *fd;
	int i, j, *fp;
    char *filename = s->s_name;
    char *buf;

	// create local memory
	buf = (char *)getbytes(0);
	buf = (char *)t_resizebytes(buf, 0, MAXPDSTRING * sizeof(char));

    canvas_makefilename(x->x_canvas, filename, buf, MAXPDSTRING);

	fd = fopen(buf, "w");
	
    if(!fd)
    {
        error("%s: couldn't create", buf);
        return;
    }		

	// write a header indicating number of clusters
	fprintf(fd, "%i\n\n", x->num_instr);
		
    for(i=0; i<x->num_instr; i++)
    {
		fp = x->cluster_members[i].member;
		
		j=0;
		
		while(x->cluster_members[i].member[j] != -1)
		{
			fprintf(fd, "%i ", *fp++);
			j++;
		};

		fprintf(fd, ";\n");
   	};
   	
    post("wrote %i clusters to file: %s.", x->num_instr, buf);
    
    fclose(fd);

    // free memory
    t_freebytes(buf, MAXPDSTRING*sizeof(char));
}


static void timbreID_read_clusters_text(t_timbreID *x, t_symbol *s)
{
	FILE *fd;
    int i, j, *fp, header;
    char *filename = s->s_name;
    char *buf, semicolon;

	// create local memory
	buf = (char *)getbytes(0);
	buf = (char *)t_resizebytes(buf, 0, MAXPDSTRING * sizeof(char));

    canvas_makefilename(x->x_canvas, filename, buf, MAXPDSTRING);

	// free current x->cluster_members memory for each list
	for(i=0; i<x->num_instr; i++)
		t_freebytes(x->cluster_members[i].member, x->cluster_member_lengths[i]*sizeof(int));


	fd = fopen(buf, "r");
	
    if(!fd)
    {
        error("%s: couldn't create", buf);
        return;
    }
		
	
	// read header indicating number of instruments
	fscanf(fd, "%i", &header);

	x->cluster_member_lengths = (int *)t_resizebytes(x->cluster_member_lengths, x->num_instr * sizeof(int), header * sizeof(int));
	x->cluster_members = (t_member *)t_resizebytes(x->cluster_members, x->num_instr * sizeof(t_member), header * sizeof(t_member));
	
	x->num_instr = header;

    for(i=0; i<x->num_instr; i++)
	{
		// can't be bigger than this, and memory can't be added incrementally during the read. have to get a big block here, then shrink later.
    	x->cluster_member_lengths[i] = x->num_instances;
		x->cluster_members[i].member = (int *)getbytes(0);
		x->cluster_members[i].member = (int *)t_resizebytes(x->cluster_members[i].member, 0, x->cluster_member_lengths[i]*sizeof(int));
	}
	
    for(i=0; i<x->num_instr; i++)
    { 	
		fp = x->cluster_members[i].member;
		j=0;
		
 		while(fscanf(fd, "%i", fp))
 		{
 			j++;
 			fp++;
 		}
 		
 		fscanf(fd, "%c", &semicolon);
 		
 		// terminate the list with -1
 		*fp = -1;

		// shrink off the excess
		x->cluster_members[i].member = (int *)t_resizebytes(x->cluster_members[i].member, x->cluster_member_lengths[i] * sizeof(int), (j+1) * sizeof(int));
		x->cluster_member_lengths[i] = (j+1);
    };
    
    for(i=0; i<x->num_instr; i++)
    {
    	int idx;
    	
    	j=0;
    	
    	while(x->cluster_members[i].member[j] != -1)
    	{
    		idx = x->cluster_members[i].member[j];
    		x->instance_cluster_membership[idx] = i;
    		j++;
    	};
    };
    	
    post("read %i clusters from file: %s.\n", x->num_instr, buf);    
    
    fclose(fd);
    
    // free memory
    t_freebytes(buf, MAXPDSTRING*sizeof(char));
}


static void *timbreID_new(void)
{	
    t_timbreID *x = (t_timbreID *)pd_new(timbreID_class);
    x->id = outlet_new(&x->x_obj, &s_float);
    x->nearest_dist = outlet_new(&x->x_obj, &s_float);
    x->confidence = outlet_new(&x->x_obj, &s_float);
	x->x_orderList = outlet_new(&x->x_obj, gensym("list"));
	x->x_featureList = outlet_new(&x->x_obj, gensym("list"));
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("list"), gensym("id"));
	inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("list"), gensym("concat_id"));

    x->instances = (t_instance *)getbytes(0);
    x->instance_feature_lengths = (int *)getbytes(0);
    x->cluster_members = (t_member *)getbytes(0);
    x->cluster_member_lengths = (int *)getbytes(0);
    x->instance_cluster_membership = (int *)getbytes(0);
    x->weights = (float *)getbytes(0);
    x->attribute_order = (int *)getbytes(0);
    x->knn_dists_idxs = (t_knn_info *)getbytes(0);
    x->feature_input = (t_float *)getbytes(0);
	
	x->feature_length = 0;
	x->num_instr=0;
    x->num_instances = 0;
    x->dist_metric = 0;  // euclidean distance by default
    x->k = 1;
    x->normalize = 0;
    x->relative_ordering = 1;

	x->prev_match = -1;
	x->max_matches = 3;
    x->reorient_flag = 0;
    x->neighborhood = 0;
    x->search_center = 0;
    x->jump_prob = 0.00;
		
    x->x_canvas = canvas_getcurrent();
    
	post("timbreID version 0.3.5");
    return (x);
}


static void timbreID_free(t_timbreID *x)
{
	int i;

	// free the database memory
	for(i=0; i<x->num_instances; i++)
		t_freebytes(x->instances[i].instance, x->instance_feature_lengths[i]*sizeof(float));
	
	t_freebytes(x->instances, x->num_instances*sizeof(t_instance));

	// free the cluster_members memory
	for(i=0; i<x->num_instr; i++)
		t_freebytes(x->cluster_members[i].member, x->cluster_member_lengths[i]*sizeof(int));
		
	t_freebytes(x->cluster_members, x->num_instr*sizeof(t_member));

	// free cluster member lengths
	t_freebytes(x->cluster_member_lengths, x->num_instr*sizeof(int));

	t_freebytes(x->knn_dists_idxs, x->num_instances*sizeof(t_knn_info));
	t_freebytes(x->instance_cluster_membership, x->num_instances*sizeof(int));
	t_freebytes(x->feature_input, x->feature_length*sizeof(t_float));
}


void timbreID_setup(void)
{
    timbreID_class = 
    class_new(
    	gensym("timbreID"),
    	(t_newmethod)timbreID_new,
    	(t_method)timbreID_free,
        sizeof(t_timbreID),
        CLASS_DEFAULT, 
		0
    );

	class_addlist(
		timbreID_class,
		(t_method)timbreID_train
	);

	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_id,
		gensym("id"),
        A_GIMME,
		0
	);

	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_worst_match,
		gensym("worst_match"),
        A_GIMME,
		0
	);
	
	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_forget, 
		gensym("forget"),
		0
	);

	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_concat_id,
		gensym("concat_id"),
        A_GIMME,
		0
	);

	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_concat_neighborhood,
		gensym("neighborhood"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_concat_jump_prob,
		gensym("jump_prob"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_concat_reorient,
		gensym("reorient"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_concat_search_center,
		gensym("search_center"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_concat_max_matches,
		gensym("max_matches"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_knn, 
		gensym("knn"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_normalize, 
		gensym("normalize"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_manual_cluster, 
		gensym("manual_cluster"),
		A_DEFFLOAT,
		A_DEFFLOAT,
		A_DEFFLOAT,
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_compute_cluster, 
		gensym("cluster"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_uncluster, 
		gensym("uncluster"),
		0
	);

	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_clusters_list, 
		gensym("clusters_list"),
		0
	);

	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_cluster_list, 
		gensym("cluster_list"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_cluster_membership, 
		gensym("cluster_membership"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_compute_order, 
		gensym("order"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_relative_ordering, 
		gensym("relative_ordering"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_compute_variance, 
		gensym("variance"),
		0
	);
	
	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_dist_metric,
		gensym("dist_metric"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_weights,
		gensym("weights"),
		A_GIMME,
		0
	);

	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_attributes,
		gensym("attributes"),
        A_GIMME,
		0
	);
	
	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_attribute_range,
		gensym("attribute_range"),
        A_DEFFLOAT,
        A_DEFFLOAT,
		0
	);

	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_order_attributes, 
		gensym("order_attributes"),
		0
	);

	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_print_instance,
		gensym("print_instance"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        A_DEFFLOAT,
		0
	);

	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_feature_list, 
		gensym("feature_list"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_max_values, 
		gensym("max_values"),
		0
	);
	
	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_min_values, 
		gensym("min_values"),
		0
	);
	
	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_clear,
		gensym("clear"),
		0
	);

	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_write,
		gensym("write"),
		A_SYMBOL,
		0
	);

	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_read,
		gensym("read"),
		A_SYMBOL,
		0
	);
	
	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_write_text,
		gensym("write_text"),
		A_SYMBOL,
		0
	);
	
	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_read_text,
		gensym("read_text"),
		A_SYMBOL,
		0
	);

	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_ARFF,
		gensym("ARFF"),
		A_GIMME,
		0
	);

	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_MATLAB,
		gensym("MATLAB"),
		A_SYMBOL,
		A_SYMBOL,
		0
	);

	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_write_clusters,
		gensym("write_clusters"),
		A_SYMBOL,
		0
	);
	
	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_read_clusters,
		gensym("read_clusters"),
		A_SYMBOL,
		0
	);
	
	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_write_clusters_text,
		gensym("write_clusters_text"),
		A_SYMBOL,
		0
	);

	class_addmethod(
		timbreID_class, 
        (t_method)timbreID_read_clusters_text,
		gensym("read_clusters_text"),
		A_SYMBOL,
		0
	);
}
