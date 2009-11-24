#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <limits>

/*
 * Computes the betweenness centrality for a binary graph.  Results are returned
 * in a vector where each element is the betweenness centrality of the
 * corresponding node.
 */
gsl_vector* bct::betweenness_bin(const gsl_matrix* m) {
	if (safe_mode) check_status(m, BINARY, "betweenness_bin");
	if (m->size1 != m->size2) {
		throw size_exception();
	}
	
	// d=1;
	int d = 1;
	
	// NPd=G;
	gsl_matrix* npd = copy(m);
	
	// NSPd=NPd;
	gsl_matrix* nspd = copy(m);
	
	// NSP=NSPd; NSP(I)=1;
	gsl_matrix* nsp = copy(m);
	gsl_vector_view nsp_diagonal = gsl_matrix_diagonal(nsp);
	gsl_vector_set_all(&nsp_diagonal.vector, 1.0);
	
	// L=NSPd; L(I)=1;
	gsl_matrix* l = copy(nsp);
	
	// while find(NSPd,1);
	while (nnz(nspd) > 0) {
		gsl_matrix* temp;
		
		// d=d+1;
		d++;
		
		// NPd=NPd*G;
		temp = mul(npd, m);
		gsl_matrix_free(npd);
		npd = temp;
		
		// NSPd=NPd.*(L==0);
		gsl_matrix* l_equal_0 = compare_elements(l, cmp_equal, 0.0);
		gsl_matrix_free(nspd);
		nspd = copy(npd);
		gsl_matrix_mul_elements(nspd, l_equal_0);
		
		// NSP=NSP+NSPd;
		gsl_matrix_add(nsp, nspd);
		
		// L=L+d.*(NSPd~=0);
		gsl_matrix* nspd_not_equal_0 = compare_elements(nspd, cmp_not_equal, 0.0);
		gsl_matrix_scale(nspd_not_equal_0, d);
		gsl_matrix_add(l, nspd_not_equal_0);
		
		gsl_matrix_free(nspd_not_equal_0);
		gsl_matrix_free(l_equal_0);
	}
	
	// L(~L)=inf; L(I)=0;
	gsl_matrix* not_l = compare_elements(l, cmp_equal, 0.0);
	logical_index_assign(l, not_l, std::numeric_limits<double>::max());
	gsl_vector_view l_diagonal = gsl_matrix_diagonal(l);
	gsl_vector_set_zero(&l_diagonal.vector);
	
	// NSP(~NSP) = 1;
	gsl_matrix* not_nsp = compare_elements(nsp, cmp_equal, 0.0);
	logical_index_assign(nsp, not_nsp, 1.0);
	
	// Gt=G.';
	gsl_matrix* transpose_m = gsl_matrix_alloc(m->size1, m->size2);
	gsl_matrix_transpose_memcpy(transpose_m, m);
	
	// DP=zeros(n);
	gsl_matrix* dp = zeros(m->size1);
	
	// diam=d-1;
	int diam = d - 1;
	
	// for d=diam:-1:2
	for (int d = diam; d >= 2; d--) {
		
		// DPd1=(((L==d).*(1+DP)./NSP)*Gt).*((L==(d-1)).*NSP);
		gsl_matrix* l_equal_d = compare_elements(l, cmp_equal, d);
		gsl_matrix* dp_plus_1 = copy(dp);
		gsl_matrix_add_constant(dp_plus_1, 1.0);
		gsl_matrix_mul_elements(l_equal_d, dp_plus_1);
		gsl_matrix_div_elements(l_equal_d, nsp);
		gsl_matrix* dpd1 = mul(l_equal_d, transpose_m);
		gsl_matrix* l_equal_d_less_1 = compare_elements(l, cmp_equal, d - 1.0);
		gsl_matrix_mul_elements(l_equal_d_less_1, nsp);
		gsl_matrix_mul_elements(dpd1, l_equal_d_less_1);
		gsl_matrix_add(dp, dpd1);
		
		gsl_matrix_free(l_equal_d_less_1);
		gsl_matrix_free(dpd1);
		gsl_matrix_free(dp_plus_1);
		gsl_matrix_free(l_equal_d);
	}
	
	// BC=sum(DP,1);
	gsl_vector* betweenness = sum(dp);
	
	gsl_matrix_free(dp);
	gsl_matrix_free(transpose_m);
	gsl_matrix_free(not_nsp);
	gsl_matrix_free(not_l);
	gsl_matrix_free(nsp);
	gsl_matrix_free(nspd);
	gsl_matrix_free(npd);
	
	return betweenness;
}
