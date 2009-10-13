#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/* Assortativity coefficient. Essentially, the assortativity a correlation 
 * coefficient for the degrees of linked nodes. A positive assortativity coefficient 
 * indicates that nodes tend to link to other nodes with the same or similar degree.  
 * The function accepts weighted networks, but the connection weights are ignored.
 */
 
double bct::assortativity(const gsl_matrix* m, int flag) {
	gsl_vector* degi;
	gsl_vector* degj;
	int K;
	
	if(flag == 0) {
		gsl_vector* deg = bct::degrees_und(m);
		gsl_matrix* ij = bct::find(bct::triu(m, 1), 0, 0); //the second parameter translates to '>'
		gsl_vector_view i = gsl_matrix_row(ij, 0);
		gsl_vector_view j = gsl_matrix_row(ij, 1);
		K = ij->size2;
		for(int k = 0;k < K;k++) {
			gsl_vector_set(degi, gsl_vector_get(deg, ((int)gsl_vector_get(&i.vector, k))), k);
			gsl_vector_set(degj, gsl_vector_get(deg, ((int)gsl_vector_get(&j.vector, k))), k);
		}
	}
	
	if(flag == 1) {
		gsl_vector* deg = bct::degrees_dir(m);
		gsl_matrix* ij = bct::find(m, 0, 0); //the second parameter translates to '>'
		gsl_vector_view i = gsl_matrix_row(ij, 0);
		gsl_vector_view j = gsl_matrix_row(ij, 1);
		K = ij->size2;
		for(int k = 0;k < K;k++) {
			gsl_vector_set(degi, gsl_vector_get(deg, ((int)gsl_vector_get(&i.vector, k))), k);
			gsl_vector_set(degj, gsl_vector_get(deg, ((int)gsl_vector_get(&j.vector, k))), k);
		}
	}
	
	//r = (sum(degi.*degj)/K - (sum(0.5*(degi+degj))/K)^2)/(sum(0.5*(degi.^2+degj.^2))/K - (sum(0.5*(degi+degj))/K)^2);
	gsl_vector* degtemp = gsl_vector_calloc(degi->size);
	gsl_vector_memcpy(degtemp, degi);
	gsl_vector_mul(degtemp, degj);
	double term1 = bct::sum(degtemp)/(double)K;
	
	gsl_vector_memcpy(degtemp, degi);
	gsl_vector_add(degtemp, degj);
	gsl_vector_scale(degtemp, (double)0.5);
	double term2 = bct::sum(degtemp)/(double)K;
	term2 *= term2;
	
	double termA = term1 - term2;
	
	gsl_vector_memcpy(degtemp, degi);
	gsl_vector* degtemp2 = gsl_vector_calloc(degj->size);
	gsl_vector_memcpy(degtemp2, degj);
	gsl_vector_add(degtemp, degtemp2);
	gsl_vector_scale(degtemp, (double)0.5);
	double term3 = bct::sum(degtemp)/(double)K;
	
	gsl_vector_memcpy(degtemp, degi);
	gsl_vector_add(degtemp, degj);
	gsl_vector_scale(degtemp, (double)0.5);
	double term4 = bct::sum(degtemp)/(double)K;
	term4 *= term4;
	
	double termB = term3 - term4;
	
	double r = (double)termA/(double)termB;
	
	return r;
}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
