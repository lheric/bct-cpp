#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/* 
/* input:
 *           CIJ  = connection/adjacency matrix
 * output:
 *           Min  = matching index for incoming connections
 *           Mout = matching index for outgoing connections
 *           Mall = matching index for all connections
 *
 * Does not use self- or cross connections for comparison.
 * Does not use connections that are not present in BOTH i and j.
 * All output matrices are calculated for upper triangular only (symmetrical).
 */

//Compare all (incoming+outgoing) connections
gsl_matrix* bct::matching_ind(const gsl_matrix* m) {
	int N, ncon;
	double matchindex;
	gsl_matrix* Mall;
	gsl_vector* use;
	
	N = m->size1;
	Mall = gsl_matrix_calloc(N, N);
	for(int i = 0;i < (N-1);i++) {
		for(int j = (i+1);j < N;j++) {
			//c1 and r1 can be moved outside the inner 'for' loop
			gsl_vector_const_view c1 = gsl_matrix_const_column(m, i);
			gsl_vector_const_view c2 = gsl_matrix_const_column(m, j);
			gsl_vector_const_view r1 = gsl_matrix_const_row(m, i);
			gsl_vector_const_view r2 = gsl_matrix_const_row(m, j);
			gsl_vector* t1 = concatenate(&c1.vector, &r1.vector);
			gsl_vector* t2 = concatenate(&c2.vector, &r2.vector);
			use = logical_or(t1, t2);
			gsl_vector_set(use, i, 0);
			gsl_vector_set(use, j, 0);
			gsl_vector_set(use, i+N, 0);
			gsl_vector_set(use, j+N, 0);			
			ncon = sum(logical_index(t1, use)) + sum(logical_index(t2, use));
			if(ncon == 0)
				gsl_matrix_set(Mall, i, j, 0.0);
			else
			{
				//Improve this by storing the logical_index results, also calculated above?
				matchindex = ((double)sum(logical_and(logical_index(t1, use), logical_index(t2, use))))/ncon;
				matchindex *= 2.0;
				gsl_matrix_set(Mall, i, j, matchindex);
			}
			gsl_vector_free(t1);
			gsl_vector_free(t2);
		}
	}
	gsl_vector_free(use);
	return Mall;
}

//Compare the incoming connections only
gsl_matrix* bct::matching_ind_in(const gsl_matrix* m) {
	int N, ncon;
	double matchindex;
	gsl_matrix* Min;
	gsl_vector* use;
	
	N = m->size1;
	Min = gsl_matrix_calloc(N, N);
	for(int i = 0;i < (N-1);i++) {
		for(int j = (i+1);j < N;j++) {
			gsl_vector_const_view c1 = gsl_matrix_const_column(m, i);
			gsl_vector_const_view c2 = gsl_matrix_const_column(m, j);
			use = logical_not(logical_and(logical_not(&c1.vector), logical_not(&c2.vector)));
			gsl_vector_set(use, i, 0);
			gsl_vector_set(use, j, 0);
			ncon = sum(logical_index(&c1.vector, use)) + sum(logical_index(&c2.vector, use));
			if(ncon == 0)
				gsl_matrix_set(Min, i, j, 0.0);
			else
			{
				//Improve this by storing the logical_index results, also calculated above?
				matchindex = ((double)sum(logical_and(logical_index(&c1.vector, use), logical_index(&c2.vector, use))))/ncon;
				matchindex *= 2.0;
				gsl_matrix_set(Min, i, j, matchindex);
			}
		}
	}
	gsl_vector_free(use);
	return Min;	
}

//Compare the outgoing connections only
gsl_matrix* bct::matching_ind_out(const gsl_matrix* m) {
	int N, ncon;
	double matchindex;
	gsl_matrix* Mout;
	gsl_vector* use;
	
	N = m->size1;
	Mout = gsl_matrix_calloc(N, N);
	for(int i = 0;i < (N-1);i++) {
		for(int j = (i+1);j < N;j++) {
			gsl_vector_const_view r1 = gsl_matrix_const_row(m, i);
			gsl_vector_const_view r2 = gsl_matrix_const_row(m, j);
			use = logical_not(logical_and(logical_not(&r1.vector), logical_not(&r2.vector)));
			gsl_vector_set(use, i, 0);
			gsl_vector_set(use, j, 0);
			ncon = sum(logical_index(&r1.vector, use)) + sum(logical_index(&r2.vector, use));
			if(ncon == 0)
				gsl_matrix_set(Mout, i, j, 0.0);
			else
			{
				//Improve this by storing the logical_index results, also calculated above?
				matchindex = ((double)sum(logical_and(logical_index(&r1.vector, use), logical_index(&r2.vector, use))))/ncon;
				matchindex *= 2.0;
				gsl_matrix_set(Mout, i, j, matchindex);
			}
		}
	}
	gsl_vector_free(use);
	return Mout;
}
		
	

		
		
	
