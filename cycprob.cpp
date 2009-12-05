#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*
 * cycprob_fcyc computes the proportion of paths that are cycels.
 * cycprob_pcyc computes the probability that non-cyclic path of length 'q-1' 
 * can be extended to form a cycle of length 'q',  for each path length 'q'.
 */
 
gsl_vector* bct::cycprob_fcyc(gsl_matrix** pq, int len_pq) {
	gsl_matrix* fcyc_m = zeros(1, len_pq);
	gsl_vector* fcyc = to_vector(fcyc_m);
	gsl_matrix_free(fcyc_m);
	for(int q = 0;q < len_pq;q++) {
		gsl_vector* first_sum = sum(pq[q]);
		double total_paths = sum(first_sum);
		if(total_paths > 0.0) {
			gsl_vector_view diagonal = gsl_matrix_diagonal(pq[q]);
			double total_cycles = sum(&diagonal.vector);
			double frac_cycles = total_cycles/total_paths;
			gsl_vector_set(fcyc, q, frac_cycles);
		}
		else {
			gsl_vector_set(fcyc, q, 0.0);
		}
		gsl_vector_free(first_sum);
	}
	return fcyc;
}

gsl_vector* bct::cycprob_pcyc(gsl_matrix** pq, int len_pq) {
	gsl_matrix* pcyc_m = zeros(1, len_pq);
	gsl_vector* pcyc = to_vector(pcyc_m);
	gsl_matrix_free(pcyc_m);
	for(int q = 1;q < len_pq;q++) {
		gsl_vector* first_sum1 = sum(pq[q-1]);
		double total_paths = sum(first_sum1);
		gsl_vector_view diagonal_1 = gsl_matrix_diagonal(pq[q-1]);
		double total_cycles1 = sum(&diagonal_1.vector);
		if((total_paths - total_cycles1) > 0) {
			gsl_vector_view diagonal_0 = gsl_matrix_diagonal(pq[q]);
			double total_cycles0 = sum(&diagonal_0.vector);
			double prob_cycles = total_cycles0/(total_paths - total_cycles1);
			gsl_vector_set(pcyc, q, prob_cycles);
		}
		else {
			gsl_vector_set(pcyc, q, 0.0);
		}
		gsl_vector_free(first_sum1);
	}
	return pcyc;
}
	
	
	
	
	
	
	
	
	
