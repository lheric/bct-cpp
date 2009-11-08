#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/* M
 * Returns a vector of indices of the elements in m that satisfy a condition
 * given by a comparison with cmprVal. Presently, the comparsion operators are
 * coded in the cmprFlag parameter, as follows:
 * 0 -> 'greater than operator, >'
 */
gsl_matrix* bct::find(const gsl_matrix* m, int cmprFlag, double cmprVal) {
	gsl_matrix* indices = gsl_matrix_alloc(2, m->size1 * m->size2);
	int size = 0;
	if(cmprFlag==0) { //means, perform a '>'
		for(int i = 0;i < m->size1;i++) {
			for(int j = 0;j < m->size2;j++) {
				double matrixVal = gsl_matrix_get(m, i, j);
				if(matrixVal > cmprVal) {
					gsl_matrix_set(indices, 0, size, i);
					gsl_matrix_set(indices, 1, size++, j);
				}
			}
		}
		gsl_matrix* trimIndices = gsl_matrix_alloc(2, size);
		gsl_matrix_view trim = gsl_matrix_submatrix(indices, 0, 0, 2, size); //the 4th and 5th params DO NOT represent indices
		gsl_matrix_memcpy(trimIndices, &trim.matrix);
		gsl_matrix_free(indices);
		return trimIndices;
	}
}

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
