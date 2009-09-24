#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*
 * Computes the in-degree for each node in a directed binary matrix.  Weights
 * are discarded.  Results are returned in a vector where each element is the
 * in-degree of the corresponding node.
 */
gsl_vector* bct::degrees_dir_in(const gsl_matrix* m) {

	// Calculated the same way as degrees_und
	return degrees_und(m);
}

/*
 * Computes the out-degree for each node in a directed binary matrix.  Weights
 * are discarded.  Results are returned in a vector where each element is the
 * out-degree of the corresponding node.
 */
gsl_vector* bct::degrees_dir_out(const gsl_matrix* m) {
	gsl_matrix* bm = binary_matrix(m);
	gsl_vector* out_degrees = gsl_vector_calloc(bm->size1);
	for (int i = 0; i < bm->size2; i++) {
		gsl_vector_const_view column = gsl_matrix_const_column(bm, i);
		gsl_vector_add(out_degrees, &column.vector);
	}
	gsl_matrix_free(bm);
	return out_degrees;
}

/*
 * Computes the degree for each node in a directed binary matrix.  Weights are
 * discarded.  Results are returned in a vector where each element is the degree
 * of the corresponding node.
 */
gsl_vector* bct::degrees_dir(const gsl_matrix* m) {
	gsl_vector* in_degrees = degrees_dir_in(m);
	gsl_vector* out_degrees = degrees_dir_out(m);
	gsl_vector_add(in_degrees, out_degrees);
	gsl_vector_free(out_degrees);
	return in_degrees;
}
