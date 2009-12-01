#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>

#define WHITE 0
#define GRAY 1
#define BLACK 2

/* Original comments:
 * note: breadth-first search tree does not contain all paths 
 * (or all shortest paths), but allows the determination of at least one 
 * path with minimum distance.
 * the entire graph is explored, starting from source vertex 'source'
 *
 * Olaf Sporns, Indiana University, 2002/2007/2008
 */

gsl_vector* bct::breadth(const gsl_matrix* CIJ, int source, gsl_vector* ret_distance) {
	if (CIJ->size1 != CIJ->size2) {
		throw size_exception();
	}
	
	gsl_matrix* m =copy(CIJ);
	
	int N = m->size1;
	gsl_matrix* color_m = zeros(1, N);
	gsl_matrix* distance_m = yens(1, N, GSL_POSINF);
	gsl_matrix* branch_m = zeros(1, N);
	
	//Convert to vectors for comfy computation
	gsl_vector* color = to_vector(color_m);
	gsl_vector* distance = to_vector(distance_m);
	gsl_vector* branch = to_vector(branch_m);
	gsl_matrix_free(color_m);
	gsl_matrix_free(distance_m);
	gsl_matrix_free(branch_m);
	
	source = source - 1; //convert node range to C++ range, that is from (1.N) to (0,N-1)
	
	gsl_vector_set(color, source, GRAY);
	gsl_vector_set(distance, source, 0.0);
	gsl_vector_set(branch, source, -1.0);
	
	gsl_vector* Q = gsl_vector_alloc(1);
	gsl_vector_set(Q, 0, source);
	
	int start_index = 0;
	int end_index = 0;
	
	while(start_index <= end_index) {
		int u = gsl_vector_get(Q, start_index);
		int dist_u = gsl_vector_get(distance, u);
		gsl_vector_view u_row = gsl_matrix_row(m, u);
		gsl_vector* ns = find(&u_row.vector);
		for(int i = 0;i< ns->size;i++) {
			int v = gsl_vector_get(ns, i);
			int dist_v = gsl_vector_get(distance, v);
			if(dist_v == 0) {
				dist_v = dist_u + 1;
				gsl_vector_set(distance, v, dist_v);
			}
			int color_v = gsl_vector_get(color, v);
			if(color_v == WHITE) {
				color_v = GRAY;
				gsl_vector_set(color, v, color_v);
				int dist_v = gsl_vector_get(distance, v);
				dist_v = dist_u + 1;
				gsl_vector_set(distance, v, dist_v);
				gsl_vector_set(branch, v, u+1.0);  //Convert 'u' back to node range (1,N)
				gsl_vector* Q_cat = concatenate(Q, v);
				gsl_vector_free(Q);
				Q = Q_cat;
				end_index++;
			}
		}
		gsl_vector_set(color, u, BLACK);
		start_index++;
	}
	
	if(ret_distance == NULL) {
		ret_distance = gsl_vector_alloc(distance->size);
		gsl_vector_memcpy(ret_distance, distance);
		gsl_vector_free(distance);
	}
	return branch;
}
				
				
	
	
	
	
	
	
	
	
	
	
	
	
