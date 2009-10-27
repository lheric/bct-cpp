#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <cassert>

bool bct::safe_mode = false;

/* M
 * Returns a vector of indices of the elements in m that satisfy a condition
 * given by a comparison with cmprVal. Presently, the comparsion operators are
 * coded in the cmprFlag parameter, as follows:
 * 0 -> 'greater than operator, >'
 */
gsl_matrix* bct::find(const gsl_matrix* m, int cmprFlag, double cmprVal) {
	gsl_matrix* indices = gsl_matrix_alloc(2, m->size1 * m->size2);
	int size = 0;
	if(cmprFlag==0) { //means, perform a '>'
		for(int i = 0;i < m->size1;i++) {
			for(int j = 0;j < m->size2;j++) {
				double matrixVal = gsl_matrix_get(m, i, j);
				if(matrixVal > cmprVal) {
					gsl_matrix_set(indices, 0, size, i);
					gsl_matrix_set(indices, 1, size++, j);
				}
			}
		}
		gsl_matrix* trimIndices = gsl_matrix_alloc(2, size);
		gsl_matrix_view trim = gsl_matrix_submatrix(indices, 0, 0, 2, size); //the 3rd and 4th params DO NOT represent indices
		gsl_matrix_memcpy(trimIndices, &trim.matrix);
		gsl_matrix_free(indices);
		return trimIndices;
	}
}

/* M
 * Strip a vector by picking only those cells where there is a corresponding
 * number 1 in the 'pick vector'.
 */
gsl_vector* bct::pick_cells(const gsl_vector* srcV, const gsl_vector* pickV) {
		int stripVindex = 0;
		int nnzV = nnz(pickV);
		if(!(nnzV > 0)) { //Matlab returns an empty matrix in this case, here it needs to be handled differently
			gsl_vector* stripV = gsl_vector_alloc(1);
			gsl_vector_set(stripV, 0, 0);
			return stripV;
		}
		else {
			gsl_vector* stripV = gsl_vector_alloc(nnzV);
			for(int i = 0;i < srcV->size;i++) {
				//1 or 1.0, does it make a difference? 
				//Nevertheless, pickV is created by logical_not method and it sets integer values
				if(gsl_vector_get(pickV, i) == 1) 
					gsl_vector_set(stripV, stripVindex++, gsl_vector_get(srcV, i));
			}
			return stripV;
		}
}

/*
 * Turns safe mode on or off.
 */
void bct::set_safe_mode(bool safe_mode) {
	bct::safe_mode = safe_mode;
}

/* M
 * Splice two vectors into one
 */
gsl_vector* bct::splice(const gsl_vector* v1, const gsl_vector* v2) {
	int spliceVindex = 0;
	gsl_vector* spliceV = gsl_vector_alloc(v1->size + v2->size);
	for(int i = 0;i < v1->size;i++)
		gsl_vector_set(spliceV, spliceVindex++, gsl_vector_get(v1, i));
	for(int i = 0;i < v2->size;i++)
		gsl_vector_set(spliceV, spliceVindex++, gsl_vector_get(v2, i));
	return spliceV;
}
