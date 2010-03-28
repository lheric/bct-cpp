#include "bct.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <iostream>

/*
 * Computes the normalized average shortest path length for a weighted graph.
 * 
 * Assumes it has been given distances rather than weights.
 */
double bct::normalized_path_length_wei(const gsl_matrix* m, double wmax) {
	if (safe_mode) check_status(m, SQUARE | WEIGHTED, "normalized_path_length_wei");
	int N = m->size1;
	gsl_matrix* D = distance_wei(m);
	double dmin = 1.0 / wmax;
	double dmax = (double) N / wmax;
	double sum = 0.0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double d = gsl_matrix_get(D, i, j);
			if (d < dmax)
				sum += d;
			else
				sum += dmax;
		}
	}
	gsl_matrix_free(D);
	double npl = ((sum / (double)(N*(N - 1))) - dmin) / (dmax - dmin);
	if (npl < -1.0e-5)
		std::cerr << "npl too negative (" << npl << ")" << std::endl;
	return fabs(npl);
}

// double bct::normalized_path_length_wei(const gsl_matrix* m, double wmax=1.0) {
// 	if (safe_mode) check_status(m, SQUARE | WEIGHTED, "normalized_path_length_wei");
// 	int N = m->size1;
// 	gsl_matrix* inv_m = invert_elements(m);
// 	gsl_matrix* D = distance_wei(inv_m);
// 	gsl_matrix_free(inv_m);
// 	double dmax = (double)N / wmax;
// 	gsl_matrix* D_gt_dmax = compare_elements(D, fp_greater, dmax);
// 	logical_index_assign(D, D_gt_dmax, dmax);
// 	gsl_matrix_free(D_gt_dmax);
// 	double sum = 0.0;
// 	for (int i = 0; i < N; i++) {
// 		for (int j = 0; j < N; j++) {
// 			sum += gsl_matrix_get(D, i, j);
// 		}
// 	}
// 	gsl_matrix_free(D);
// 	double dmin = 1.0 / wmax;
// 	return ((sum / (double)(N*(N - 1))) - dmin) / (dmax - dmin);
// }
