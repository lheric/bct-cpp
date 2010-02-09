#include "bct.h"
#include <ctime>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>

/*
 * Returns a randomized graph with equivalent degree sequence to the original
 * weighted directed graph.  On average, each edge is rewired ITER times.  Out-
 * strength is preserved for weighted graphs, while in-strength is not.
 */
gsl_matrix* bct::randmio_dir(const gsl_matrix* R, int ITER) {
	if (R->size1 != R->size2) throw size_exception();
	
	gsl_rng_default_seed = std::time(NULL);
	gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);
	
	// [i j]=find(R);
	gsl_matrix* R_ij = find_ij(R);
	gsl_vector* i = gsl_matrix_alloc(R->size1);
	gsl_vector* j = gsl_matrix_alloc(R->size1);
	gsl_matrix_get_col(i, R_ij, 0);
	gsl_matrix_get_col(i, R_ij, 1);
	gsl_matrix_free(R_ij);
	
	// K=length(i);
	int K = length(i);
	
	// ITER=K*ITER;
	ITER = K * ITER;
	
	gsl_matrix* _R = copy(R);
	
	// for iter=1:ITER
	for (int iter = 1; iter <= ITER; iter++) {
		
		// while 1
		while (true) {
			
			// while 1
			while (true) {
				
				// e1=ceil(K*rand);
				int e1 = gsl_rng_uniform_int(rng, K);
				
				// e2=ceil(K*rand);
				int e2 = gsl_rng_uniform_int(rng, K);
				
				// while (e2==e1),
				while (e2 == e1) {
					
					// e2=ceil(K*rand);
					e2 = gsl_rng_uniform_int(rng, K);
				}
				
				// a=i(e1); b=j(e1);
				int a = (int)gsl_vector_get(i, e1);
				int b = (int)gsl_vector_get(j, e1);
				
				// c=i(e2); d=j(e2);
				int c = (int)gsl_vector_get(i, e2);
				int d = (int)gsl_vector_get(j, e2);
				
				// if all(a~=[c d]) && all(b~=[c d]);
				if (a != c && a != d && b != c && b != d) {
					
					// break
					break;
				}
			}
			
			// if ~(R(a,d) || R(c,b))
			if (fp_zero(gsl_matrix_get(_R, a, d)) && fp_zero(gsl_matrix_get(_R, c, b))) {
				
				// R(a,d)=R(a,b); R(a,b)=0;
				gsl_matrix_set(_R, a, d, gsl_matrix_get(_R, a, b));
				gsl_matrix_set(_R, a, b, 0.0);
				
				// R(c,b)=R(c,d); R(c,d)=0;
				gsl_matrix_set(_R, c, b, gsl_matrix_get(_R, c, d));
				gsl_matrix_set(_R, c, d, 0.0);
				
				// j(e1) = d;
				gsl_vector_set(j, e1, (double)d);
				
				// j(e2) = b;
				gsl_vector_set(j, e2, (double)b);
				
				// break;
				break;
			}
		}
	}
	
	gsl_vector_free(i);
	gsl_vector_free(j);
	return _R;
}
