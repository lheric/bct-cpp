#include "bct.h"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*
 * Returns a binary copy of the given vector.
 */
gsl_vector* bct::binary(const gsl_vector* v) {
	gsl_vector* bv = gsl_vector_calloc(v->size);
	for (int i = 0; i < v->size; i++) {
		if (std::abs(gsl_vector_get(v, i)) > EPSILON) {
			gsl_vector_set(bv, i, 1.0);
		}
	}
	return bv;
}

/*
 * Returns a binary copy of the given matrix.
 */
gsl_matrix* bct::binary(const gsl_matrix* m) {
	gsl_matrix* bm = gsl_matrix_calloc(m->size1, m->size2);
	for (int i = 0; i < m->size1; i++) {
		for (int j = 0; j < m->size2; j++) {
			if (std::abs(gsl_matrix_get(m, i, j)) > EPSILON) {
				gsl_matrix_set(bm, i, j, 1.0);
			}
		}
	}
	return bm;
}

/*
 * Returns a copy of the given matrix with all elements along the diagonal set
 * to zero.
 */
gsl_matrix* bct::zero_diagonal(const gsl_matrix* m) {
	gsl_matrix* zdm = gsl_matrix_alloc(m->size1, m->size2);
	gsl_matrix_memcpy(zdm, m);
	gsl_vector_view diagonal = gsl_matrix_diagonal(zdm);
	gsl_vector_set_zero(&diagonal.vector);
	return zdm;
}

/*
 * Returns the number of nonzero elements in the given vector.
 */
int bct::nnz(const gsl_vector* v) {
	int count = 0;
	for (int i = 0; i < v->size; i++) {
		if (std::abs(gsl_vector_get(v, i)) > EPSILON) {
			count++;
		}
	}
	return count;
}

/*
 * Returns the number of nonzero elements in the given matrix.
 */
int bct::nnz(const gsl_matrix* m) {
	int count = 0;
	for (int i = 0; i < m->size1; i++) {
		for (int j = 0; j < m->size2; j++) {
			if (std::abs(gsl_matrix_get(m, i, j)) > EPSILON) {
				count++;
			}
		}
	}
	return count;
}

/*
 * Returns a vector of indices of nonzero elements in the given vector.
 */
gsl_vector* bct::find(const gsl_vector* v) {
	int size = nnz(v);
	assert(size > 0);
	gsl_vector* nz = gsl_vector_alloc(size);
	int pos = 0;
	for (int i = 0; i < v->size; i++) {
		if (std::abs(gsl_vector_get(v, i)) > EPSILON) {
			gsl_vector_set(nz, pos, i);
			pos++;
		}
	}
	return nz;
}

/*
 * Returns a matrix of size (rows->size, columns->size).  Elements are copied
 * from the specified rows and columns of m.
 */
gsl_matrix* bct::submatrix(const gsl_matrix* m, const gsl_vector* rows, const gsl_vector* columns) {
	gsl_matrix* s = gsl_matrix_alloc(rows->size, columns->size);
	for (int i = 0; i < rows->size; i++) {
		for (int j = 0; j < columns->size; j++) {
			gsl_matrix_set(s, i, j, gsl_matrix_get(m, gsl_vector_get(rows, i), gsl_vector_get(columns, j)));
		}
	}
	return s;
}

/* M
 * Returns the sum of the elements of a matrix, row-wise or column-wise
 * it is the equivalent of the 'sum' function in matlab. If the second
 * parameter 'dim' is 1, the sum is a cumulative sum of row slices. If
 * dim is 2, the result is a cumulative sum of column slices. In any case,
 *  it is assumed that m is a square matrix and therefore, the size of
 * sum is just set to m->size1 because m->size1 = m->size2 in a square matrix
 */
gsl_vector* bct::sum(const gsl_matrix* m, int dim) {
	gsl_vector* sum = gsl_vector_calloc(m->size1);
	if (dim == 1)
	{
		for (int i = 0; i < m->size1; i++) {
			gsl_vector_const_view row = gsl_matrix_const_row(m, i);
			gsl_vector_add(sum,&row.vector);
		}
	}
	else if (dim == 2)
	{
		for (int i = 0; i < m->size1; i++) {
			gsl_vector_const_view column = gsl_matrix_const_column(m, i);
			gsl_vector_add(sum,&column.vector);
		}
	}
	return sum;
}

/* M 
 * Returns the sum of the elements of a vector. Type cast the result
 * at the receiving end of this function
 */
double bct::sum(const gsl_vector* v) {
	double sum = 0.0;
	for (int i = 0; i < v->size; i++) {
		sum += gsl_vector_get(v, i);
	}
	return sum;
}

/* M
 * equivalent of 'triu' in octave. Returns a copy of the matrix m
 * retaining only the values of the upper triangle of m (including
 * the diagonal) while setting all the other values to zero. A 
 * memcpy followed by setting lower triangle to zeros is possibly a
 * more efficient approach than setting all to zeros and then copying
 * cell by cell, all the upper triangle elements
 * NOTE: K can't be > m->size1. This exception needs to be handled.
 */
gsl_matrix* bct::triu(const gsl_matrix* m, int K) {
	gsl_matrix* triu_m = gsl_matrix_calloc(m->size1, m->size2);
	gsl_matrix_memcpy(triu_m, m);
	for(int i = 0;i < m->size1; i++)
		for(int j=0;j < (i+K); j++)
			gsl_matrix_set(triu_m, i, j, 0.0);
	return triu_m;
}

/* M
 * Equivalent of 'tril' in octave. Returns a copy of the matrix m
 * retaining only the values of the lower triangle of m (including
 * the diagonal) while setting all the other values to zero. A 
 * memcpy followed by setting upper triangle to zeros is possibly a
 * more efficient approach than setting all to zeros and then copying
 * cell by cell, all the lower triangle elements.
 * NOTE: K can't be > m->size1. This exception needs to be handled.
 */
gsl_matrix* bct::tril(const gsl_matrix* m, int K) {
	gsl_matrix* tril_m = gsl_matrix_calloc(m->size1, m->size2);
	gsl_matrix_memcpy(tril_m, m);
	for(int i = 0;i < m->size1; i++)
		for(int j=i+1+K;j < m->size2; j++)
			gsl_matrix_set(tril_m, i, j, 0.0);
	return tril_m;
}

/* M
 * Perform a logical NOT on the elements of a vector
 */
gsl_vector* vectorNot(const gsl_vector* v) {
	gsl_vector* notV = gsl_vector_calloc(v->size);
	for(int i = 0;i < v->size;i++)
		gsl_vector_set(notV, i, !(int)gsl_vector_get(v, i));
	return notV;
}

/* M
 * Perform a logical AND between the elements of 2 vectors
 */
gsl_vector* vectorAnd(const gsl_vector* v1, const gsl_vector* v2) {
	gsl_vector* andV = gsl_vector_calloc(v1->size);
	for(int i = 0;i < v1->size;i++) {
		int val1 = (int)gsl_vector_get(v1, i);
		int val2 = (int)gsl_vector_get(v2, i);
		gsl_vector_set(andV, i, val1 && val2);
	}
	return andV;
}

/* M
 * Strip a vector by picking only those cells where there is a corresponding
 * number 1 in the 'pick vector'.
 */
gsl_vector* pickCells(const gsl_vector* srcV, const gsl_vector* pickV) {
		int stripVindex = 0;
		int nnzV = bct::nnz(pickV);
		gsl_vector* stripV = gsl_vector_calloc(nnzV);
		for(int i = 0;i < srcV->size;i++) {
			//1 or 1.0, does it make a difference? 
			//Nevertheless, pickV is created vectorNot method and it sets the values to integer
			if(gsl_vector_get(pickV, i) == 1) 
				gsl_vector_set(stripV, stripVindex++, gsl_vector_get(srcV, i));
		}
		return stripV;
}

/* M
 * Splice two vectors into one
 */
gsl_vector* splice(const gsl_vector* v1, const gsl_vector* v2) {
	int spliceVindex = 0;
	gsl_vector* spliceV = gsl_vector_calloc(v1->size + v2->size);
	for(int i = 0;i < v1->size;i++)
		gsl_vector_set(spliceV, spliceVindex++, gsl_vector_get(v1, i));
	for(int i = 0;i < v2->size;i++)
		gsl_vector_set(spliceV, spliceVindex++, gsl_vector_get(v2, i));
	return spliceV;
}
		
/*
 * Prints a vector using the given format for each element.  This is only
 * provided for debugging purposes.  In other cases, use gsl_vector_fprintf.
 */
void bct::printf(const gsl_vector* v, const char* format) {
	for (int i = 0; i < v->size; i++) {
		std::printf(format, gsl_vector_get(v, i));
		std::printf(" ");
	}
	std::printf("\n");
}

/*
 * Prints a matrix using the given format for each element.  This is only
 * provided for debugging purposes.  In other cases, use gsl_matrix_fprintf.
 */
void bct::printf(const gsl_matrix* m, const char* format) {
	for (int i = 0; i < m->size1; i++) {
		for (int j = 0; j < m->size2; j++) {
			std::printf(format, gsl_matrix_get(m, i, j));
			std::printf(" ");
		}
		std::printf("\n");
	}
}



