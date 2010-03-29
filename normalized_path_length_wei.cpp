#include "bct.h"
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <iostream>

/*
 * Computes the normalized average shortest path length for a weighted graph.
 * Assumes edges are distances rather than weights.
 */
double bct::normalized_path_length_wei(const gsl_matrix* m, double wmax) {
	if (safe_mode) check_status(m, SQUARE | WEIGHTED, "normalized_path_length_wei");
	int N = m->size1;
	gsl_matrix* D = distance_wei(m);
	double dmin = 1.0 / wmax;
	double dmax = (double)N / wmax;
	double sum = 0.0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double d = gsl_matrix_get(D, i, j);
			sum += (d < dmax) ? d : dmax;
		}
	}
	gsl_matrix_free(D);
	double npl = ((sum / (double)(N * (N - 1))) - dmin) / (dmax - dmin);
	if (fp_negative(npl)) {
		std::cerr << "Negative NPL (" << npl << ")." << std::endl;
	}
	return std::abs(npl);
}
