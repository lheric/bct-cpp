#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>

gsl_matrix* bct::breadthdist(const gsl_matrix* m, gsl_matrix* ret_reachability) {
	if (m->size1 != m->size2) {
		throw size_exception();
	}
	
	int N = m->size1;
	
	gsl_matrix* D = zeros(N);
	for(int i = 0;i < N;i++) {
		gsl_vector* D_ret = breadth(m, i+1);
		gsl_vector_view D_row = gsl_matrix_row(D, i);
		gsl_vector_memcpy(&D_row.vector, D_ret);
		gsl_vector_free(D_ret);
	}
	
	gsl_matrix* D_ind = compare_elements(D, cmp_equal, 0.0);
	logical_index_assign(D, D_ind, GSL_POSINF);
	
	gsl_matrix* R = compare_elements(D, cmp_not_equal, GSL_POSINF);
	
	if(ret_reachability != NULL) {
		ret_reachability = gsl_matrix_alloc(R->size1, R->size2);
		gsl_matrix_memcpy(ret_reachability, R);
		gsl_matrix_free(R);
	}
	
	return D;
}
	
	
	
	
