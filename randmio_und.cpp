#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

/*
 * Returns a 'randomized' graph R, with equivalent degree sequence to the original 
 * weighted undirected graph G, and with preserved connectedness 
 * (hence the input graph must be connected).
 *
 * Each edge is rewired (on average) ITER times. The strength distributions 
 * are not preserved for weighted graphs.
 *
 * Rewiring algorithm: Maslov and Sneppen (2002) Science 296:910
 */

gsl_matrix* bct::randmio_und(const gsl_matrix* m, const int iters) {
	//[i j]=find(tril(R));
	gsl_matrix* R = copy(m);
	int n = m->size1;
	gsl_matrix* tril_R = tril(R);
	gsl_matrix* R_ij = find_ij(tril_R); 
	gsl_matrix_free(tril_R);
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
			double val = ((double)rand())/((double)RAND_MAX);
			//flip edge c-d with 50% probability to explore all potential rewirings
			if(val > 0.5) {
	            gsl_vector_set(&i.vector, e2, d);
	            gsl_vector_set(&j.vector, e2, c);
	            c = gsl_vector_get(&i.vector, e2); 
	            d = gsl_vector_get(&j.vector, e2);
	        }
	        
	        //if ~(R(a,d) || R(c,b))
	        //rewiring condition
	        if(!(gsl_matrix_get(R, a, d) || gsl_matrix_get(R, c, b))) {
        		//R(a,d)=R(a,b); R(a,b)=0;
    			gsl_matrix_set(R, a, d, gsl_matrix_get(R, a, b));
    			gsl_matrix_set(R, a, b, 0.0);
	    		//R(d,a)=R(b,a); R(b,a)=0;
	    		gsl_matrix_set(R, d, a, gsl_matrix_get(R, b, a));
    			gsl_matrix_set(R, b, a, 0.0);
            	//R(c,b)=R(c,d); R(c,d)=0;
            	gsl_matrix_set(R, c, b, gsl_matrix_get(R, c, d));
    			gsl_matrix_set(R, c, d, 0.0);
            	//R(b,c)=R(d,c); R(d,c)=0;
            	gsl_matrix_set(R, b, c, gsl_matrix_get(R, d, c));
    			gsl_matrix_set(R, d, c, 0.0);
    			gsl_vector_set(&j.vector, e1, d);
    			gsl_vector_set(&j.vector, e2, b);
    			break;
			}
		}
	}
	gsl_matrix_free(R_ij);
	return R;
}      	
	        					
	        				
	        	
	            
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
		
	
				










