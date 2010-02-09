#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*
 *Computes the normalized average shortest path length for directed & weighted network
 *
 */

double bct::norm_avr_shortest_path_length_wd(gsl_matrix* m) {
	
	if (safe_mode) check_status(m, WEIGHTED | DIRECTED, "norm_avr_shortest_path_length_wd");
	if (m->size1 != m->size2) {
		throw size_exception();
	}
	
	//Mappting from weights to distances
	double wei_edge;
	double wei_max = 8;
	for ( int i = 0; i < int(m->size1); i++ ) {
		for ( int j = 0; j < int(m->size2); j++ ) {
			wei_edge = gsl_matrix_get(m,i,j);
			if (wei_edge) {
				gsl_matrix_set(m,i,j,1/wei_edge);
			}
		}
	}
	
	
	
	
	//gsl_matrix* D = gsl_matrix_alloc(m->size1,m->size2);
	gsl_matrix* D = bct::distance_wei(m);
	//gsl_matrix* D = gsl_matrix_alloc(m->size1,m->size2);
	//gsl_matrix_memcpy(D,m);
	
	double avr_dis_per_node = 0;
	double sum;
	double norm_dis;
	double max_dis = D->size1*wei_max; //defined the largest value for element in D matrix
	
	for ( int i = 0; i < (int)D->size1; i++ ) {
		sum = 0;
		for ( int j = 0; j < (int)D->size2; j++ ) {
			if (gsl_matrix_get(D,i,j) == GSL_POSINF or gsl_matrix_get(D,i,j) > max_dis)//deal with infinite or invalid distances
				gsl_matrix_set(D,i,j,max_dis);
						
			sum += gsl_matrix_get(D,i,j);
		}
		avr_dis_per_node += sum/(D->size2-1);
	}
	
	norm_dis = (avr_dis_per_node/D->size1 - 1) / (max_dis - 1);
	gsl_matrix_free(D);
	return norm_dis;
}
			
			
		
		

