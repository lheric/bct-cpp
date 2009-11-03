#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*
 * Computes the degree for each node in a directed binary matrix.  Weights are
 * discarded.  Results are returned in a vector where each element is the degree
 * of the corresponding node.  In- and out-degree may be obtained by passing
 * pointers to appropriately sized vectors for the last two arguments.
 */
gsl_vector* bct::degrees_dir(const gsl_matrix* m, gsl_vector* in_degrees, gsl_vector* out_degrees) {
	if (safe_mode) check_status(m, DIRECTED, "degrees_dir");
	if (m->size1 != m->size2) {
		throw size_exception();
	}
	
	// Allocate storage for in- and out-degree vectors if necessary
	bool free_in_degrees = false;
	bool free_out_degrees = false;
	if (in_degrees == NULL) {
		free_in_degrees = true;
		in_degrees = gsl_vector_alloc(m->size2);
	}
	if (out_degrees == NULL) {
		free_out_degrees = true;
		out_degrees = gsl_vector_alloc(m->size1);
	}
	
	// Calculate in-degree
	for (int i = 0; i < m->size2; i++) {
		gsl_vector_const_view column = gsl_matrix_const_column(m, i);
		gsl_vector_set(in_degrees, i, nnz(&column.vector));
	}
	
	// Calculate out-degree
	for (int i = 0; i < m->size1; i++) {
		gsl_vector_const_view row = gsl_matrix_const_row(m, i);
		gsl_vector_set(out_degrees, i, nnz(&row.vector));
	}
	
	// Sum in- and out-degree to get total degree
	gsl_vector* degrees = gsl_vector_calloc(m->size1);
	gsl_vector_add(degrees, in_degrees);
	gsl_vector_add(degrees, out_degrees);
	
	// Free in- and out-degree vectors if necessary
	if (free_in_degrees) {
		gsl_vector_free(in_degrees);
	}
	if (free_out_degrees) {
		gsl_vector_free(out_degrees);
	}
	
	return degrees;
}
