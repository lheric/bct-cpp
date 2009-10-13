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

bct::matching_ind_outputs matching_ind(const gsl_matrix* m) {
	int N, ncon;
	double matchindex;
	gsl_matrix* Min;
	gsl_matrix* Mout;
	gsl_matrix* Mall;
	gsl_vector* use;

	N = m->size1;
	
	//Compare the incoming connections only
	Min = gsl_matrix_calloc(N, N);
	for(int i = 0;i < (N-1);i++)
		for(int j = (i+1);j < N;j++) {
			gsl_vector_const_view c1 = gsl_matrix_const_column(m, i);
			gsl_vector_const_view c2 = gsl_matrix_const_column(m, j);
			use = bct::vectorNot(bct::vectorAnd(bct::vectorNot(&c1.vector), bct::vectorNot(&c2.vector)));
			gsl_vector_set(use, i, 0);
			gsl_vector_set(use, j, 0);
			ncon = bct::sum(bct::pickCells(&c1.vector, use)) + bct::sum(bct::pickCells(&c2.vector, use));
			if(ncon == 0)
				gsl_matrix_set(Min, i, j, 0.0);
			else
			{
				//Improve this by storing the pickCells results, also calculated above?
				matchindex = ((double)bct::sum(bct::vectorAnd(bct::pickCells(&c1.vector, use), bct::pickCells(&c2.vector, use))))/ncon;
				matchindex *= 2.0;
				gsl_matrix_set(Min, i, j, matchindex);
			}
		}
	
	//Compare the outgoing connections only
	Mout = gsl_matrix_calloc(N, N);
	for(int i = 0;i < (N-1);i++)
		for(int j = (i+1);j < N;j++) {
			gsl_vector_const_view r1 = gsl_matrix_const_row(m, i);
			gsl_vector_const_view r2 = gsl_matrix_const_row(m, j);
			use = bct::vectorNot(bct::vectorAnd(bct::vectorNot(&r1.vector), bct::vectorNot(&r2.vector)));
			gsl_vector_set(use, i, 0);
			gsl_vector_set(use, j, 0);
			ncon = bct::sum(bct::pickCells(&r1.vector, use)) + bct::sum(bct::pickCells(&r2.vector, use));
			if(ncon == 0)
				gsl_matrix_set(Mout, i, j, 0.0);
			else
			{
				//Improve this by storing the pickCells results, also calculated above?
				matchindex = ((double)bct::sum(bct::vectorAnd(bct::pickCells(&r1.vector, use), bct::pickCells(&r2.vector, use))))/ncon;
				matchindex *= 2.0;
				gsl_matrix_set(Mout, i, j, matchindex);
			}
		}
	
	//Compare all (incoming+outgoing) connections
	Mall = gsl_matrix_calloc(N, N);
	for(int i = 0;i < (N-1);i++)
		for(int j = (i+1);j < N;j++) {
			gsl_vector_const_view c1 = gsl_matrix_const_column(m, i);
			gsl_vector_const_view c2 = gsl_matrix_const_column(m, j);
			gsl_vector_const_view r1 = gsl_matrix_const_row(m, i);
			gsl_vector_const_view r2 = gsl_matrix_const_row(m, j);
			gsl_vector* t1 = bct::splice(&c1.vector, &r1.vector);
			gsl_vector* t2 = bct::splice(&c2.vector, &r2.vector);
			use = bct::vectorNot(bct::vectorAnd(bct::vectorNot(t1), bct::vectorNot(t2)));
			gsl_vector_set(use, i, 0);
			gsl_vector_set(use, j, 0);
			gsl_vector_set(use, i+N, 0);
			gsl_vector_set(use, j+N, 0);			
			ncon = bct::sum(bct::pickCells(t1, use)) + bct::sum(bct::pickCells(t2, use));
			if(ncon == 0)
				gsl_matrix_set(Mall, i, j, 0.0);
			else
			{
				//Improve this by storing the pickCells results, also calculated above?
				matchindex = ((double)bct::sum(bct::vectorAnd(bct::pickCells(t1, use), bct::pickCells(t2, use))))/ncon;
				matchindex *= 2.0;
				gsl_matrix_set(Mall, i, j, matchindex);
			}
		}
	bct::matching_ind_outputs matching_inds;
	matching_inds.Min = Min;
	matching_inds.Mout = Mout;
	matching_inds.Mall = Mall;
	return matching_inds;
}
		
	

		
		
	
