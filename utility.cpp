#include "bct.h"
#include <cmath>
#include <cstring>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <sstream>
#include <vector>

/*
 * Catches GSL errors and throws BCT exceptions.
 */
void bct::gsl_error_handler(const char* reason, const char* file, int line, int gsl_errno) {
	std::stringstream what;
	what << reason << " in " << file << ", line " << line << ".";
	throw bct_exception(what.str());
}

/*
 * Overloaded convenience function for freeing GSL vectors and matrices.
 */
void bct::gsl_free(gsl_vector* v) { gsl_vector_free(v); }
void bct::gsl_free(gsl_matrix* m) { gsl_matrix_free(m); }
void bct::gsl_free(std::vector<gsl_matrix*>& m) {
	for (int i = 0; i < (int)m.size(); i++) {
		if (m[i] != NULL) {
			gsl_matrix_free(m[i]);
			m[i] = NULL;
		}
	}
}

/*
 * Initializes the BCT library for external use.
 */
void bct::init() {
	gsl_set_error_handler(gsl_error_handler);
}

/*
 * Returns the number of edges in a directed graph.
 */
int bct::number_of_edges_dir(const gsl_matrix* m) {
	return nnz(m);
}

/*
 * Returns the number of edges in an undirected graph.
 */
int bct::number_of_edges_und(const gsl_matrix* m) {
	gsl_matrix* triu_m = triu(m);
	int ret = nnz(triu_m);
	gsl_matrix_free(triu_m);
	return ret;
}

/*
 * Returns the number of nodes in a graph.
 */
int bct::number_of_nodes(const gsl_matrix* m) {
	return (int)m->size1;
}

// TODO: These belong in matlab/functions.cpp
// Move them and remove unnecessary #includes when SWIG-accessible

double bct::mean(const gsl_vector* v, const char* opt) {
	if (std::strcmp(opt, "a") == 0) {
		double sum = 0.0;
		for (int i = 0; i < (int)v->size; i++) {
			sum += gsl_vector_get(v, i);
		}
		return sum / (double)v->size;
	} else if (std::strcmp(opt, "g") == 0) {
		double product = 1.0;
		for (int i = 0; i < (int)v->size; i++) {
			product *= gsl_vector_get(v, i);
		}
		return std::pow(product, 1.0 / (double)v->size);
	} else if (std::strcmp(opt, "h") == 0) {
		double sum = 0.0;
		for (int i = 0; i < (int)v->size; i++) {
			sum += 1.0 / gsl_vector_get(v, i);
		}
		return (double)v->size / sum;
	} else {
		// TODO: Error message?
		return 0.0;
	}
}

gsl_vector* bct::mean(const gsl_matrix* m, int dim, const char* opt) {
	if (dim == 1) {
		gsl_vector* mean_v = gsl_vector_alloc(m->size2);
		for (int i = 0; i < (int)m->size2; i++) {
			gsl_vector_const_view m_col_i = gsl_matrix_const_column(m, i);
			double value = mean(&m_col_i.vector, opt);
			gsl_vector_set(mean_v, i, value);
		}
		return mean_v;
	} else if (dim == 2) {
		gsl_vector* mean_v = gsl_vector_alloc(m->size1);
		for (int i = 0; i < (int)m->size1; i++) {
			gsl_vector_const_view m_row_i = gsl_matrix_const_row(m, i);
			double value = mean(&m_row_i.vector, opt);
			gsl_vector_set(mean_v, i, value);
		}
		return mean_v;
	} else {
		return NULL;
	}
}
