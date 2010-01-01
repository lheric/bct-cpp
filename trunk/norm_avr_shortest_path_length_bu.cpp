#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*
 *Computes the normalized average shortest path length for undirected & binary network
 *
 */

double bct::norm_avr_shortest_path_length_bu(const gsl_matrix* m) {
	
	if (safe_mode) check_status(m, BINARY | UNDIRECTED, "norm_avr_shortest_path_length_bu");
	if (m->size1 != m->size2) {
		throw size_exception();
	}
	//gsl_matrix* D = gsl_matrix_alloc(m->size1,m->size2);
	gsl_matrix* D = bct::distance_bin(m);
	//gsl_matrix* D = gsl_matrix_alloc(m->size1,m->size2);
	//gsl_matrix_memcpy(D,m);
	
	double avr_dis_per_node = 0;
	double sum;
	double norm_dis;
	
	for ( int i = 0; i < D->size1; i++ ) {
		sum = 0;
		for ( int j = 0; j < D->size2; j++ ) {
			if (gsl_matrix_get(D,i,j) == GSL_POSINF)
				gsl_matrix_set(D,i,j,D->size1);
			
			sum += gsl_matrix_get(D,i,j);
		}
		avr_dis_per_node += sum/(D->size2-1);
	}
	
	norm_dis = (avr_dis_per_node/D->size1 - 1) / (D->size1 - 1);
	gsl_matrix_free(D);
	return norm_dis;
}
			
			
		
		

