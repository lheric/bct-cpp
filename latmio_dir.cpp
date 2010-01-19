#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

/*
 * Returns a 'latticized' graph R, with equivalent degree sequence to the original 
 * weighted directed graph G, and with preserved connectedness 
 * (hence the input graph must be connected).
 *
 * Each edge is rewired (on average) ITER times. The strength distributions 
 * are not preserved for weighted graphs.
 *
 * Rewiring algorithm: Maslov and Sneppen (2002) Science 296:910
 * Latticizing algorithm: Sporns and Zwi (2004); Neuroinformatics 2:145
 */

gsl_matrix* bct::latmio_dir(const gsl_matrix* m, int iters) {
	//It seems D could also be loaded from a file. How to consider that?
	gsl_matrix* R = copy(m);
	int n = m->size1;
	gsl_matrix* D = zeros(n);
	gsl_vector* v1 = sequence(1, 1, n-1); 
	gsl_vector* v2 = sequence(n-1, -1, 1); 
	gsl_matrix* u1 = concatenate_columns(v1, v2); 
	gsl_vector* v3 = min(u1); 
	gsl_vector* u = concatenate(0.0, v3); 
	gsl_vector_free(v1);
	gsl_vector_free(v2);
	gsl_vector_free(v3);
	gsl_matrix_free(u1);
	int upper_lim = (int)ceil((n-1)/2);
	for(int v = 0;v <= upper_lim;v++) {
		gsl_vector* index1 = sequence((v+1), 1, (n-1)); 
		gsl_vector* index2 = sequence(0, 1, v); 
		gsl_vector* u_seg1 = ordinal_index(u, index1); 
		gsl_vector* u_seg2 = ordinal_index(u, index2); 
		gsl_vector* seg_splice = concatenate(u_seg1, u_seg2); 
		gsl_vector_view bottom_row = gsl_matrix_row(D, (n-1-v));
		gsl_vector_memcpy(&bottom_row.vector, seg_splice);
		gsl_vector_view top_row = gsl_matrix_row(D, v);
		gsl_vector* rev_seg_splice = reverse(seg_splice); 
		gsl_vector_memcpy(&top_row.vector, rev_seg_splice);
		gsl_vector_free(index1);
		gsl_vector_free(index2);
		gsl_vector_free(u_seg1);
		gsl_vector_free(u_seg2);
		gsl_vector_free(seg_splice);
		gsl_vector_free(rev_seg_splice);
		
	}
	gsl_vector_free(u);
	
	//[i j]=find(R);
	gsl_matrix* R_ij = find_ij(R); 
	gsl_vector_view i = gsl_matrix_column(R_ij, 0); 
	gsl_vector_view j = gsl_matrix_column(R_ij, 1); 
	int K = R_ij->size1;
	int tot_iters = iters;
	tot_iters *= K;
	int rewire = 0;
	srand(time(0));
	for(int iter = 1;iter <= tot_iters;iter++) {
		while(1) {
			rewire = 1;
			int a,b,c,d;
			int e1, e2;
			while(1) {
				e1 = ceil((K-1) * (((double)rand())/((double)RAND_MAX)));
				e2 = ceil((K-1) * (((double)rand())/((double)RAND_MAX)));
				while(e2 == e1) {
					e2 = ceil((K-1) * (((double)rand())/((double)RAND_MAX)));
				}
				a = gsl_vector_get(&i.vector, e1);
				b = gsl_vector_get(&j.vector, e1);
				c = gsl_vector_get(&i.vector, e2);
				d = gsl_vector_get(&j.vector, e2);
				//if all(a~=[c d]) && all(b~=[c d]);
				if(a != c && a != d && b!=c && b!=d) {  //all four vertices must be different
					break;
				}
			}
	        
	        //if ~(R(a,d) || R(c,b))
	        //rewiring condition
	        if(!(gsl_matrix_get(R, a, d) || gsl_matrix_get(R, c, b))) {
	        	//if (D(a,b)+D(c,d))>=(D(a,d)+D(c,b))
	        	double val1 = gsl_matrix_get(D, a, b);
	        	double val2 = gsl_matrix_get(D, c, d);
	        	double val3 = gsl_matrix_get(D, a, d);
	        	double val4 = gsl_matrix_get(D, c, b);
	        	if((val1 + val2) >= (val3 + val4)) { //lattice condition
	        		//R(a,d)=R(a,b); R(a,b)=0;
	    			gsl_matrix_set(R, a, d, gsl_matrix_get(R, a, b));
	    			gsl_matrix_set(R, a, b, 0.0);
                	//R(c,b)=R(c,d); R(c,d)=0;
                	gsl_matrix_set(R, c, b, gsl_matrix_get(R, c, d));
	    			gsl_matrix_set(R, c, d, 0.0);
	    			gsl_vector_set(&j.vector, e1, d);
	    			gsl_vector_set(&j.vector, e2, b);
	    			break;
		    	}
			}
		}
	}
	gsl_matrix_free(R_ij);
	gsl_matrix_free(D);
	return R;
}      	
	        					
	        				
	        	
	            
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
		
	
				










