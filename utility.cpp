#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*
 * Catches GSL errors and throws BCT exceptions.
 */
void bct::gsl_error_handler(const char* reason, const char* file, int line, int gsl_errno) {
	if (gsl_errno == GSL_ENOMEM) {
		throw out_of_memory_exception();
	} else {
		throw gsl_exception();
	}
}

/*
 * Overloaded convenience function for freeing GSL vectors and matrices.
 */
void bct::gsl_free(gsl_vector* v) { gsl_vector_free(v); }
void bct::gsl_free(gsl_matrix* m) { gsl_matrix_free(m); }

/*
 * Initializes the BCT library for external use.
 */
void bct::init() {
	gsl_set_error_handler(gsl_error_handler);
}

/* M
 * Splice two vectors into one
 */
gsl_vector* bct::splice(const gsl_vector* v1, const gsl_vector* v2) {
	int spliceVindex = 0;
	gsl_vector* spliceV = gsl_vector_alloc(v1->size + v2->size);
	for(int i = 0;i < v1->size;i++)
		gsl_vector_set(spliceV, spliceVindex++, gsl_vector_get(v1, i));
	for(int i = 0;i < v2->size;i++)
		gsl_vector_set(spliceV, spliceVindex++, gsl_vector_get(v2, i));
	return spliceV;
}
