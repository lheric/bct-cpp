#include "bct.h"
#include <gsl/gsl_vector.h>

/*
 * Returns a binary copy of the given vector.
 */
gsl_vector* bct::binary_vector(const gsl_vector* v) {
	gsl_vector* bv = gsl_vector_calloc(v->size);
	for (int i = 0; i < v->size; i++) {
		if (abs(gsl_vector_get(v, i)) > EPSILON) {
			gsl_vector_set(bv, i, 1.0);
		}
	}
	return bv;
}
