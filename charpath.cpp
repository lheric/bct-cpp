#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <math.h>

/*
 * Characteristic path length is calculated as the global mean of the
 * distance matrix D, not taking into account any 'Infs' but including the
 * distances on the main diagonal.
 */

double bct::charpath_lambda(const gsl_matrix* D) {
	//lambda = sum(sum(D(D~=Inf)))/length(nonzeros(D~=Inf));
	gsl_matrix* D_not_inf = compare_elements(D, cmp_not_equal, GSL_POSINF);
	gsl_vector* D_indxd = logical_index(D, D_not_inf);
	double total_dist = sum(D_indxd);
	double num_connected = nnz(D_not_inf);
	double lambda = (double)(total_dist/num_connected);
	gsl_matrix_free(D_not_inf);
	gsl_vector_free(D_indxd);
	return lambda;
}

gsl_vector* bct::charpath_ecc(const gsl_matrix* D, double *radius, double *diameter) {
	//ecc = max(D.*(D~=Inf),[],2);	
	gsl_matrix* D_not_inf = compare_elements(D, cmp_not_equal, GSL_POSINF);
	gsl_vector* D_indxd = logical_index(D, D_not_inf);
	gsl_matrix* D_inf_nullified = zeros(D->size1, D->size2);
	logical_index_assign(D_inf_nullified, D_not_inf, D_indxd);
	gsl_vector* ecc = max(D_inf_nullified, 2);
	gsl_matrix_free(D_not_inf);
	gsl_vector_free(D_indxd);
	//radius = min(ecc);
	if(radius == NULL) {
		radius = new double;
	}
	*radius = gsl_vector_min(ecc);
	//diameter = max(ecc);
	if(diameter == NULL) {
		diameter = new double;
	}
	*diameter = gsl_vector_max(ecc);
	return ecc;
}
