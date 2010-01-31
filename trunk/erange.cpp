#include "bct.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*
 * Computes the range for each edge (i.e., the shortest path length from i to j
 * after edge i-j has been removed from the graph.  Optional returns:
 * - eta:    Average range for the entire graph.
 * - eshort: Logical matrix indicated edges that are shortcuts.
 * - fs:     Fraction of shortcuts in the graph.
 */
gsl_matrix* bct::erange(const gsl_matrix* m, double* eta, gsl_matrix* eshort, double* fs) {
	if (m->size1 != m->size2) {
		throw size_exception();
	}
	
	// K = length(nonzeros(CIJ));
	int k = nnz(m);
	
	// Erange = zeros(N,N);
	gsl_matrix* erange = gsl_matrix_calloc(m->size1, m->size2);
	
	// [i,j] = find(CIJ==1);
	gsl_matrix* m_equal_1 = compare_elements(m, fp_equal, 1.0);
	gsl_matrix* m_equal_1_indices = find_ij(m_equal_1);
	gsl_matrix_free(m_equal_1);	
	
	// for c=1:length(i)
	for (int c = 0; c < (int)m_equal_1_indices->size1; c++) {
		
		// CIJcut = CIJ;
		gsl_matrix* m_cut = copy(m);
		
		// CIJcut(i(c),j(c)) = 0;
		int i = (int)gsl_matrix_get(m_equal_1_indices, c, 0);
		int j = (int)gsl_matrix_get(m_equal_1_indices, c, 1);
		gsl_matrix_set(m_cut, i, j, 0.0);
		
		// [R,D] = reachdist(CIJcut);
		gsl_matrix* d = reachdist(m_cut);
		
		// Erange(i(c),j(c)) = D(i(c),j(c))
		gsl_matrix_set(erange, i, j, gsl_matrix_get(d, i, j));
		gsl_matrix_free(d);
		gsl_matrix_free(m_cut);
	}
	
	gsl_matrix_free(m_equal_1_indices);
	
	// eta = sum(Erange((Erange>0)&(Erange<Inf)))/length(Erange((Erange>0)&(Erange<Inf)));
	if (eta != NULL) {
		gsl_matrix* erange_positive = compare_elements(erange, fp_greater, 0.0);
		gsl_matrix* erange_finite = compare_elements(erange, fp_less, GSL_POSINF);
		gsl_matrix* erange_positive_finite = logical_and(erange_positive, erange_finite);
		gsl_vector* erange_indexed = logical_index(erange, erange_positive_finite);
		*eta = sum(erange_indexed) / (double)erange_indexed->size;
		gsl_vector_free(erange_indexed);
		gsl_matrix_free(erange_positive_finite);
		gsl_matrix_free(erange_finite);
		gsl_matrix_free(erange_positive);
	}
	
	// Eshort = Erange>2;
	if (eshort != NULL || fs != NULL) {
		if (eshort != NULL) {
			gsl_matrix_free(eshort);
		}
		eshort = compare_elements(erange, fp_greater, 2.0);
	}
	
	// fs = length(nonzeros(Eshort))/K;
	if (fs != NULL) {
		*fs = (double)nnz(eshort) / (double)k;
	}
	
	return erange;
}
