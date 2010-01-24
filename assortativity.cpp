#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

double assortativity(const gsl_vector*, const gsl_matrix*);

/*
 * Computes assortativity for a directed matrix.  Connection weights are
 * ignored.
 */
double bct::assortativity_dir(const gsl_matrix* m) {
	if (safe_mode) check_status(m, DIRECTED, "assortativity_dir");
	
	// [id,od,deg] = degrees_dir(CIJ);
	gsl_vector* deg = degrees_dir(m);
	
	// [i,j] = find(CIJ>0);
	gsl_matrix* m_gt_0 = compare_elements(m, fp_greater, 0.0);
	gsl_matrix* ij = find_ij(m_gt_0);
	gsl_matrix_free(m_gt_0);
	
	double ret = assortativity(deg, ij);
	gsl_vector_free(deg);
	gsl_matrix_free(ij);
	return ret;
}

/*
 * Computes assortativity for an undirected matrix.  Connection weights are
 * ignored.
 */
double bct::assortativity_und(const gsl_matrix* m) {
	if (safe_mode) check_status(m, UNDIRECTED, "assortativity_und");
	
	// [deg] = degrees_und(m);
	gsl_vector* deg = degrees_und(m);
	
	// [i,j] = find(triu(CIJ,1)>0);
	gsl_matrix* triu_m = triu(m, 1);
	gsl_matrix* triu_m_gt_0 = compare_elements(triu_m, fp_greater, 0.0);
	gsl_matrix_free(triu_m);
	gsl_matrix* ij = find_ij(triu_m_gt_0);
	gsl_matrix_free(triu_m_gt_0);
	
	double ret = assortativity(deg, ij);
	gsl_vector_free(deg);
	gsl_matrix_free(ij);
	return ret;
}

double assortativity(const gsl_vector* deg, const gsl_matrix* ij) {
	using namespace bct;
	
	gsl_vector_const_view i = gsl_matrix_const_column(ij, 0);
	gsl_vector_const_view j = gsl_matrix_const_column(ij, 1);
	gsl_vector* degi = gsl_vector_alloc(ij->size1);
	gsl_vector* degj = gsl_vector_alloc(ij->size1);		
	
	// K = length(i);
	int K = ij->size1;
	
	// for k=1:K
	for (int k = 0; k < K; k++) {
		
		// degi(k) = deg(i(k));
		int i_k = (int)gsl_vector_get(&i.vector, k);
		gsl_vector_set(degi, k, gsl_vector_get(deg, i_k));
		
		// degj(k) = deg(j(k));
		int j_k = (int)gsl_vector_get(&j.vector, k);
		gsl_vector_set(degj, k, gsl_vector_get(deg, j_k));
	}
	
	// r = (sum(degi.*degj)/K - (sum(0.5*(degi+degj))/K)^2)/(sum(0.5*(degi.^2+degj.^2))/K - (sum(0.5*(degi+degj))/K)^2);
	
	gsl_vector* temp1 = copy(degi);
	gsl_vector_memcpy(temp1, degi);
	gsl_vector_mul(temp1, degj);
	double r1 = sum(temp1) / (double)K;
	
	gsl_vector_memcpy(temp1, degi);
	gsl_vector_add(temp1, degj);
	gsl_vector_scale(temp1, 0.5);
	double r2 = sum(temp1) / (double)K;
	r2 *= r2;
	
	gsl_vector_free(temp1);
	temp1 = pow_elements(degi, 2);
	gsl_vector* temp2 = pow_elements(degj, 2);
	gsl_vector_add(temp1, temp2);
	gsl_vector_scale(temp1, 0.5);
	double r3 = sum(temp1) / (double)K;
	
	gsl_vector_memcpy(temp1, degi);
	gsl_vector_add(temp1, degj);
	gsl_vector_scale(temp1, 0.5);
	double r4 = sum(temp1) / (double)K;
	r4 *= r4;
	
	gsl_vector_free(temp1);
	gsl_vector_free(temp2);
	return (r1 - r2) / (r3 - r4);
}
