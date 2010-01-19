#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <math.h>
#include <cstdio>

/* Original comments:
 * Distance matrix for weighted networks. The input matrix must be a mapping 
 * from weight to distance. For instance, in a weighted correlation network, higher 
 * correlations are more naturally interpreted as shorter distances. Consequently, 
 * in this case, the input matrix should be some inverse of the connectivity matrix.
 */
 
gsl_matrix* bct::distance_wei(const gsl_matrix* m) {
	if (safe_mode) check_status(m, WEIGHTED, "distance_wei");
	int N = m->size1;
	gsl_matrix* D = zeros(N);
	gsl_matrix* d_eye = eye(N);
	gsl_matrix* not_d_eye = logical_not(d_eye);
	logical_index_assign(D, not_d_eye, GSL_POSINF);
	gsl_matrix_free(d_eye);
	gsl_matrix_free(not_d_eye);

	for(int u = 0;u < N;u++) {
		gsl_matrix* S_m = ones(1, N);
		gsl_vector* S = to_vector(S_m);
		gsl_matrix_free(S_m);
		gsl_matrix* G1 = copy(m);	
		gsl_vector* V = gsl_vector_alloc(1);
		gsl_vector_set(V, 0, u);
		while(true) {
			ordinal_index_assign(S, V, 0.0);
			//G1(:,V)=0;
			gsl_vector* all_rows = sequence(0, N-1);
			ordinal_index_assign(G1, all_rows, V, 0.0);
			for(int i = 0;i < V->size;i++) {
				int v = gsl_vector_get(V, i);
				//W=find(G1(v,:))				
				gsl_vector_view G1_row_v = gsl_matrix_row(G1, v);
				gsl_vector* W = find(&G1_row_v.vector);
				//D(u,W)=min([D(u,W);D(u,v)+G1(v,W)]);
				if(W != NULL) {
					double dist_uv = gsl_matrix_get(D, u, v);
					gsl_vector* G1_v_indxd = ordinal_index(&G1_row_v.vector, W);
					gsl_vector_add_constant(G1_v_indxd, dist_uv);
					gsl_vector_view D_row_u = gsl_matrix_row(D, u);
					gsl_vector* D_u_indxd = ordinal_index(&D_row_u.vector, W);
					gsl_matrix* concat_DG = concatenate_columns(D_u_indxd, G1_v_indxd);
					gsl_vector* min_DG = min(concat_DG);
					gsl_matrix* min_DG_m = to_row_matrix(min_DG);
					gsl_vector* row_index = gsl_vector_alloc(1);
					gsl_vector_set(row_index, 0, u);
					ordinal_index_assign(D, row_index, W, min_DG_m);
					//free resources
					gsl_vector_free(G1_v_indxd);
					gsl_vector_free(D_u_indxd);
					gsl_vector_free(row_index);
					gsl_vector_free(min_DG);
					gsl_matrix_free(min_DG_m);
					gsl_matrix_free(concat_DG);
					gsl_vector_free(W);				
				}
			}
			//minD=min(D(u,S));
			gsl_vector* row_index = gsl_vector_alloc(1);
			gsl_vector_set(row_index, 0, u);
			gsl_matrix* D_S_indxd = ord_log_index(D, row_index, S);
			double minD = -1;
			if(D_S_indxd != NULL) {
				minD = gsl_matrix_min(D_S_indxd);
			}
			 //if isempty(minD)||isinf(minD)
			 if((minD == -1) || gsl_isinf(minD)) {
			 	break;
			 }
			 // V=find(D(u,:)==minD);
			 gsl_vector_view D_row_u = gsl_matrix_row(D, u);
			 gsl_vector* D_minD_ind = compare_elements(&D_row_u.vector, fp_equal, minD);
			 gsl_vector_free(V);
			 V = find(D_minD_ind);
			 gsl_vector_free(row_index);
			 gsl_vector_free(D_minD_ind);
			 gsl_matrix_free(D_S_indxd);
		}
        gsl_vector_free(S);		
        gsl_matrix_free(G1);
	}
	return D;
}
			
				
				
				
				
				
				
				
				
				
				
				
				
				
				
		
		
		
