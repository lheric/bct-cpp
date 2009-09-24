#include "bct.h"
#include <gsl/gsl_matrix.h>

/*
 * Returns a binary copy of the given matrix.
 */
gsl_matrix* bct::binary_matrix(const gsl_matrix* m) {
	gsl_matrix* bm = gsl_matrix_calloc(m->size1, m->size2);
	for (int i = 0; i < m->size1; i++) {
		for (int j = 0; j < m->size2; j++) {
			if (abs(gsl_matrix_get(m, i, j)) > EPSILON) {
				gsl_matrix_set(bm, i, j, 1.0);
			}
		}
	}
	return bm;
}
