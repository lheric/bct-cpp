#include "bct.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>

/*
 * Computes the normalized average shortest path length for a weighted graph.
 */
double bct::normalized_path_length_wei(const gsl_matrix* m, double wmax) {
	if (safe_mode) check_status(m, WEIGHTED, "norm_avr_shortest_path_length_wei");
	if (m->size1 != m->size2) throw size_exception();
	int N = m->size1;
	gsl_matrix* inv_m = invert_elements(m);
	gsl_matrix* D = distance_wei(inv_m);
	gsl_matrix_free(inv_m);
	double dmax = (double)N / wmax;
	gsl_matrix* D_gt_dmax = compare_elements(D, fp_greater, dmax);
	logical_index_assign(D, D_gt_dmax, dmax);
	gsl_matrix_free(D_gt_dmax);
	double sum = 0.0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			sum += gsl_matrix_get(D, i, j);
		}
	}
	gsl_matrix_free(D);
	double dmin = 1.0 / wmax;
	return ((sum / (double)(N * N - N)) - dmin) / (dmax - dmin);
}
