#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/* Computes the joint degree distribution of the matrix m. The outputs are
 * in the structure jdegree_outputs, as follows:
 * J = matrix in which the value of each element (u,v) corresponds to the 
 * number of nodes that have u outgoing connections and v incoming connections.
 * J_od = number of vertices with od>id.
 * J_id = number of vertices with id>od.
 * J_bl = number of vertices with id=od.
 */
 
gsl_matrix* bct::jdegree(const gsl_matrix* m) {
	gsl_matrix* bin_m;
	gsl_matrix* J;
	gsl_vector* id;
	gsl_vector* od;
	int N, max_degree;
	int J_od, J_id, J_bl;
	double max_id, max_od;

	bin_m = binary(m);
	N = m->size1;
	id = sum(bin_m,1);
	od = sum(bin_m,2);
	max_id = gsl_vector_max(id);
	max_od = gsl_vector_max(od);
	if(max_id > max_od)
		max_degree = (int)max_id;
	else
		max_degree = (int)max_od;
	//calloc initializes the matrix to zeros
	J = gsl_matrix_calloc(max_degree+1, max_degree+1);
	for(int i = 0;i < N;i++) {
		int x = (int)gsl_vector_get(id, i);
		int y = (int)gsl_vector_get(od, i);
		//In jdegree.m, x and y are incremented by 1 because matlab matrices'
		//indices start from 1 and not 0		
		gsl_matrix_set(J, x, y, (int)gsl_matrix_get(J, x, y)+1); 
	}
	return J;
}

gsl_vector* bct::jdegree_id_od_bl(gsl_matrix* J) {
	gsl_vector* jdegree_others = gsl_vector_alloc(3);
	int J_id;
	J_id = (int)sum(sum(tril(J,-1)));
	gsl_vector_set(jdegree_others, 0, J_id);
	
	int J_od;
	J_od = (int)sum(sum(triu(J,1)));
	gsl_vector_set(jdegree_others, 1, J_od);
	
	int J_bl;
	gsl_vector_view diagonal = gsl_matrix_diagonal(J);
	J_bl = (int)sum(&diagonal.vector);
	gsl_vector_set(jdegree_others, 2, J_bl);

	return jdegree_others;
}


	
	
	
