#include "bct.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>

/*
 * Computes the normalized average shortest path length for a binary graph.
 */
double bct::norm_avr_shortest_path_length_bin(const gsl_matrix* m) {
	if (safe_mode) check_status(m, BINARY, "norm_avr_shortest_path_length_bin");
	if (m->size1 != m->size2) throw size_exception();
	int N = m->size1;
	gsl_matrix* D = distance_bin(m);
	gsl_matrix* D_eq_inf = compare_elements(D, fp_equal, GSL_POSINF);
	logical_index_assign(D, D_eq_inf, (double)N);
	gsl_matrix_free(D_eq_inf);
	double sum = 0.0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			sum += gsl_matrix_get(D, i, j);
		}
	}
	gsl_matrix_free(D);
	return ((sum / (double)(N * N - N)) - 1.0) / (double)(N - 1);
}
