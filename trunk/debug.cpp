#include "bct.h"
#include <cstdio>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

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
