#include "bct.h"
#include <cmath>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*
 * Computes the clustering coefficient for a weighted undirected graph.  Results
 * are returned in a vector where each element is the clustering coefficient of
 * the corresponding node.
 */
gsl_vector* bct::clustering_coef_wu(const gsl_matrix* m) {
	if (safe_mode) check_status(m, WEIGHTED | UNDIRECTED, "clustering_coef_wu");
	if (m->size1 != m->size2) {
		throw size_exception();
	}
	
	// K=sum(W~=0,2)
	gsl_matrix* bm = binary(m);
	gsl_vector* K = sum(bm, 2);
	
	// cyc3=diag((W.^(1/3))^3)
	gsl_matrix* cyc3m = gsl_matrix_alloc(m->size1, m->size2);
	for (int i = 0; i < m->size1; i++) {
		for (int j = 0; j < m->size2; j++) {
			gsl_matrix_set(cyc3m, i, j, std::pow(gsl_matrix_get(m, i, j), 1.0 / 3.0));
		}
	}
	gsl_matrix* cyc3m2 = gsl_matrix_alloc(m->size1, m->size2);
	gsl_matrix* cyc3m3 = gsl_matrix_alloc(m->size1, m->size2);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, cyc3m, cyc3m, 1.0, cyc3m2);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, cyc3m2, cyc3m, 1.0, cyc3m3);
	gsl_vector* cyc3 = gsl_vector_alloc(m->size1);
	gsl_vector_view cyc3v = gsl_matrix_diagonal(cyc3m3);
	gsl_vector_memcpy(cyc3, &cyc3v.vector);
	
	// C=cyc3./(K.*(K-1))
	gsl_vector* K2 = gsl_vector_alloc(m->size1);
	for (int i = 0; i < K2->size; i++) {
		double Kval = gsl_vector_get(K, i);
		gsl_vector_set(K2, i, Kval * (Kval - 1));
	}
	gsl_vector_div(cyc3, K2);
	
	// Set clustering coefficient to 0 for nodes with no 3-cycles
	for (int i = 0; i < m->size1; i++) {
		if (is_zero(gsl_vector_get(&cyc3v.vector, i))) {
			gsl_vector_set(cyc3, i, 0.0);
		}
	}
	
	gsl_vector_free(K2);
	gsl_matrix_free(cyc3m3);
	gsl_matrix_free(cyc3m2);
	gsl_matrix_free(cyc3m);
	gsl_vector_free(K);
	gsl_matrix_free(bm);
	
	return cyc3;
}
