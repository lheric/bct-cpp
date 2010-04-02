#include "bct.h"
#include <cmath>
#include <gsl/gsl_matrix.h>
#include <iostream>

/*
 * Given a distance matrix, computes the normalized shortest path length.
 */
double bct::normalized_path_length(const gsl_matrix* D, double wmax) {
	int N = D->size1;
	double dmin = 1.0 / wmax;
	double dmax = (double)N / wmax;
	double sum = 0.0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double d = gsl_matrix_get(D, i, j);
			sum += (d < dmax) ? d : dmax;
		}
	}
	double npl = ((sum / (double)(N * (N - 1))) - dmin) / (dmax - dmin);
	if (npl < 0.0) {
		std::cerr << "Negative normalized path length (" << npl << ")." << std::endl;
	}
	return std::abs(npl);
}
