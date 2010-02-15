#include "bct.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>

/*
 * Computes the normalized average shortest path length for a weighted graph.
 */
double bct::norm_avr_shortest_path_length_wei(const gsl_matrix* m, double max_weight) {
	if (safe_mode) check_status(m, WEIGHTED, "norm_avr_shortest_path_length_wei");
	if (m->size1 != m->size2) throw size_exception();
	int N = m->size1;
	gsl_matrix* _m = gsl_matrix_calloc(N, N);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double value = gsl_matrix_get(m, i, j);
			if (fp_nonzero(value)) {
				gsl_matrix_set(_m, i, j, 1.0 / value);
			}
		}
	}
	gsl_matrix* D = distance_wei(_m);
	gsl_matrix_free(_m);
	double max_distance = (double)N / max_weight;
	gsl_matrix* D_gt_max_distance = compare_elements(D, fp_greater, max_distance);
	logical_index_assign(D, D_gt_max_distance, max_distance);
	gsl_matrix_free(D_gt_max_distance);
	double sum = 0.0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			sum += gsl_matrix_get(D, i, j);
		}
	}
	gsl_matrix_free(D);
	double min_distance = 1.0 / max_weight;
	return ((sum / (double)(N * N - N)) - min_distance) / (max_distance - min_distance);
}
