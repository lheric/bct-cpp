#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h> 
#include <cstdio>

/*
 * Computes the Global/local efficiency for binary undirected graph G.
 * Global efficiency: outputs the inverse distance matrix: the mean of this
 * matrix (excluding main diagonal) is equivalent to the global efficiency.
 * Local efficiency:  outputs individual nodal local efficiency.
 * For directed networks, local efficiency works with the out-degree.
 */
 
gsl_matrix* bct::efficiency(const gsl_matrix* CIJ, int* local_ind) {
	if (safe_mode) check_status(CIJ, BINARY | UNDIRECTED, "efficiency");
	
	gsl_matrix* m =copy(CIJ);
	if(local_ind != NULL) { //compute local efficiency
		int N = m->size1;
		gsl_matrix* E = zeros(N, 1);
		for(int node=0;node < N;node++) {
			gsl_vector_view neighbors_row = gsl_matrix_row(m, node);
			gsl_vector* neighbors = find(&neighbors_row.vector);
			int degree = neighbors->size;
			if(degree >= 2) {
				gsl_matrix* m_neighbors = index(m, neighbors, neighbors);
				gsl_matrix* eff = bct::distance_inv(m_neighbors);
				gsl_vector* eff_v = to_vector(eff);
				double factor = (double)1/((degree*degree)-degree);
				gsl_vector_scale(eff_v, factor);
				double efficiency = sum(eff_v);
				gsl_matrix_set(E, node, 0, efficiency);
			}
		}
		return E;
	}
	else {
		gsl_matrix* E = distance_inv(m);
		return E;
	}
}

gsl_matrix* bct::distance_inv(gsl_matrix* g) {
	gsl_matrix* D = eye(g->size1);
	int n = 1;
	gsl_matrix* npath = copy(g);
	gsl_matrix* L = binary(npath);
	
	gsl_vector* L_non_zero;
	while((L_non_zero = find(L, 1)) != NULL) {
		gsl_vector_free(L_non_zero);	
		// D=D+n.*L;
		gsl_matrix_scale(L, n);
		gsl_matrix_add(D, L);
		n++;
		gsl_matrix* npath_temp = mul(npath, g);
		gsl_matrix_free(npath);
		npath = npath_temp;
		gsl_matrix* L_temp = binary(npath);
		gsl_matrix* D_zero_ind = compare_elements(D, cmp_equal, 0.0);
		gsl_matrix_mul_elements(L_temp, D_zero_ind);
		gsl_matrix_free(L);
		L = L_temp;
		gsl_matrix_free(D_zero_ind);
	}
	
	gsl_matrix* notD_ind = compare_elements(D, cmp_equal, 0.0);
	
	logical_index_assign(D, notD_ind, GSL_POSINF);

	//invert D. This is not the same as inverse of a matrix
	gsl_matrix* div_D = ones(D->size1);
	gsl_matrix_div_elements(div_D, D);
	gsl_matrix_memcpy(D, div_D);
	gsl_matrix* g_eye = eye(g->size1);
    gsl_matrix_sub(D, g_eye);

    return D;
}
    
    
				
			
	
	
	
	
	
	
	
	
	
	
	
	
	
	
			
