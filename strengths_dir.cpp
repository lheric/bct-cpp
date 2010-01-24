#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*
 * Computes the strength for each node in a directed matrix.
 */
gsl_vector* bct::strengths_dir(const gsl_matrix* m, gsl_vector** is, gsl_vector** os) {
	if (safe_mode) check_status(m, DIRECTED, "strengths_dir");
	if (m->size1 != m->size2) throw size_exception();
	
	// is = sum(CIJ,1);
	gsl_vector* _is = sum(m, 1);
	
	// os = sum(CIJ,2);
	gsl_vector* _os = sum(m, 2);
	
	// str = is+os;
	gsl_vector* str = copy(_is);
	gsl_vector_add(str, _os);
	
	if (is != NULL) *is = _is; else gsl_vector_free(_is);
	if (os != NULL) *os = _os; else gsl_vector_free(_os);
	return str;
}
