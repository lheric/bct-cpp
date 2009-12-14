#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>

 /* Original comments:
 * This function yields the reachability matrix and the distance matrix
 * based on the power of the adjacency matrix - this will execute a lot
 * faster for matrices with low average distance between vertices.
 * Another way to get the reachability matrix and the distance matrix uses 
 * breadth-first search (see 'breadthdist.m').  'reachdist' seems a 
 * little faster most of the time.  However, 'breadthdist' uses less memory 
 * in many cases.
 *
 * Olaf Sporns, Indiana University, 2002/2007/2008
 */

gsl_matrix* bct::reachdist(gsl_matrix *m, gsl_matrix* ret_R) {
	gsl_matrix* R = gsl_matrix_alloc(m->size1, m->size2);
	gsl_matrix* D = gsl_matrix_alloc(m->size1, m->size2);
    gsl_matrix_memcpy(R, m);
    gsl_matrix_memcpy(D, m);
	
	int powr = 2;
	int N = m->size1;
	
	gsl_matrix* CIJpwr = gsl_matrix_alloc(m->size1, m->size2);
	gsl_matrix_memcpy(CIJpwr, m);
	
	gsl_vector* id = sum(m, 1);
	gsl_vector* od = sum(m, 2);
	gsl_vector* id_0_ind = compare_elements(id, cmp_equal, 0.0);
	gsl_vector* id_0 = find(id_0_ind);
	gsl_vector* od_0_ind = compare_elements(od, cmp_equal, 0.0);
	gsl_vector* od_0 = find(od_0_ind);
	
	gsl_vector* row = find(id);
	gsl_vector* col = find(od);
	
	int stop = 0;
	while(!stop) {
		gsl_matrix* CIJpwr_temp = mul(CIJpwr, m);
		gsl_matrix_free(CIJpwr);
		CIJpwr = CIJpwr_temp;
		
		//R = double(R | ((CIJpwr)~=0));
		gsl_matrix* CIJpwr_bin = to_binary(CIJpwr);
		gsl_matrix* R_temp = logical_or(R, CIJpwr_bin);
		gsl_matrix_free(R);
		R = R_temp;
		gsl_matrix_free(CIJpwr_bin);
		gsl_matrix_add(D, R);

		//if ((powr<=N)&&(~isempty(nonzeros(R(row,col)==0)))) 
		gsl_matrix* R_indxd = ordinal_index(R, row, col);
		gsl_matrix* zeros_in_Rindxd = compare_elements(R_indxd, cmp_equal, 0.0);
		int nnz_value = nnz(zeros_in_Rindxd);
		gsl_matrix_free(R_indxd);
		gsl_matrix_free(zeros_in_Rindxd);
		if(powr <= N && nnz_value > 0) {
			powr++;
		}
		else {
			stop = 1;
		}
	}
	
	gsl_matrix_scale(D, -1.0);
	gsl_matrix_add_constant(D, (powr+1));
	gsl_matrix* D_ind = compare_elements(D, cmp_equal, (double)(N+2));
	logical_index_assign(D, D_ind, GSL_POSINF); //assigns the value 'inf', meaning positive infinity
	if(id_0 != NULL) {
		gsl_vector* all_rows = sequence(0, N-1);
		ordinal_index_assign(D, all_rows, id_0, GSL_POSINF);
		gsl_vector_free(all_rows);
	}
	if(od_0 != NULL) {
		gsl_vector* all_columns = sequence(0, N-1);
		ordinal_index_assign(D, od_0, all_columns, GSL_POSINF);
		gsl_vector_free(all_columns);
	}
	
	//return R
	if(ret_R != NULL) {
		ret_R = gsl_matrix_alloc(R->size1, R->size2);
		gsl_matrix_memcpy(ret_R, R);
		gsl_matrix_free(R);
	}
	
	return D;
}

