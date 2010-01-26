#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

double assortativity(const gsl_vector*, const gsl_matrix*);

/*
 * Computes assortativity for a directed graph.  Connection weights are ignored.
 */
double bct::assortativity_dir(const gsl_matrix* CIJ) {
	if (safe_mode) check_status(CIJ, DIRECTED, "assortativity_dir");
	
	// [id,od,deg] = degrees_dir(CIJ);
	gsl_vector* deg = degrees_dir(CIJ);
	
	// [i,j] = find(CIJ>0);
	gsl_matrix* m_gt_0 = compare_elements(CIJ, fp_greater, 0.0);
	gsl_matrix* ij = find_ij(m_gt_0);
	gsl_matrix_free(m_gt_0);
	
	double ret = assortativity(deg, ij);
	gsl_vector_free(deg);
	gsl_matrix_free(ij);
	return ret;
}

/*
 * Computes assortativity for an undirected graph.  Connection weights are
 * ignored.
 */
double bct::assortativity_und(const gsl_matrix* CIJ) {
	if (safe_mode) check_status(CIJ, UNDIRECTED, "assortativity_und");
	
	// [deg] = degrees_und(m);
	gsl_vector* deg = degrees_und(CIJ);
	
	// [i,j] = find(triu(CIJ,1)>0);
	gsl_matrix* triu_m = triu(CIJ, 1);
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
	
	gsl_vector* product = copy(degi);
	gsl_vector_mul(product, degj);
	double r1 = sum(product) / (double)K;
	gsl_vector_free(product);
	
	gsl_vector* sum_scaled = copy(degi);
	gsl_vector_add(sum_scaled, degj);
	gsl_vector_scale(sum_scaled, 0.5);
	double r2 = sum(sum_scaled) / (double)K;
	gsl_vector_free(sum_scaled);
	r2 *= r2;
	
	gsl_vector* sumsq_scaled = pow_elements(degi, 2);
	gsl_vector* degjsq = pow_elements(degj, 2);
	gsl_vector_add(sumsq_scaled, degjsq);
	gsl_vector_free(degjsq);
	gsl_vector_scale(sumsq_scaled, 0.5);
	double r3 = sum(sumsq_scaled) / (double)K;
	gsl_vector_free(sumsq_scaled);
	
	gsl_vector_free(degi);
	gsl_vector_free(degj);
	
	return (r1 - r2) / (r3 - r2);
}
