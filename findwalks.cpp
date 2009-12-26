#include "bct.h"
#include "matlab.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*
 * Uses the powers of the adjacency matrix to produce numbers of walks
 * Note that Wq grows very quickly for larger N,K,q.
 * Note: Weights are discarded.
 */

gsl_vector* bct::findwalks(gsl_matrix* m, gsl_matrix** ret_Wq, double* ret_twalk) {
	gsl_matrix* CIJ = to_binary(m);
	int N = CIJ->size1;
	//Wq = zeros(N,N,N);
	gsl_matrix* Wq[N];
	for(int i = 0;i < N;i++) {
		Wq[i] = gsl_matrix_alloc(N, N);  //'calloc' in original matlab source but is actually not necessary
	}
	gsl_matrix* CIJpwr = copy(CIJ);
	//Wq(:,:,1) = CIJ;
	Wq[0] = copy(CIJ);
	for(int q = 1;q < N;q++) {
		gsl_matrix* CIJpwr_temp = mul(CIJpwr, CIJ);
		gsl_matrix_free(CIJpwr);
		CIJpwr = CIJpwr_temp;
		Wq[q] = copy(CIJpwr);
	}
	//twalk = sum(sum(sum(Wq)));
	gsl_matrix* first_sum = NULL;
	for(int i=0;i < N;i++) {
		gsl_vector* temp_sum = sum(Wq[i]);
		gsl_matrix* concat_sum = concatenate_columns(first_sum, temp_sum);
		if(first_sum != NULL) {
			gsl_matrix_free(first_sum);
		}
		gsl_vector_free(temp_sum);
		first_sum = concat_sum;
	}
	gsl_vector* second_sum = sum(first_sum, 2);
	double twalk = sum(second_sum);
	gsl_matrix_free(first_sum);
	//wlq = reshape(sum(sum(Wq)),1,N);
	gsl_vector* wlq = second_sum;
	//assign ret_Wq
	if(ret_Wq == NULL) {
			ret_Wq = new gsl_matrix* [N];
	}
	for(int i=0;i < N;i++) {
		ret_Wq[i] = gsl_matrix_alloc(N, N);
	}
	for(int i=0;i < N;i++) {	
		gsl_matrix_memcpy(ret_Wq[i], Wq[i]);
	}
	for(int i=0;i < N;i++) {	
		gsl_matrix_free(Wq[i]);
	}
	//assign ret_twalk
	if(ret_twalk == NULL) {
		ret_twalk = new double;
	}
	*ret_twalk = twalk;
	return wlq;
}
	
	
	
