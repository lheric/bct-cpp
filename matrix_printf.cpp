#include "bct.h"
#include <cstdio>
#include <gsl/gsl_matrix.h>

/*
 * Prints a matrix using the given format for each element.  This is only
 * provided for debugging purposes.  In other cases, use gsl_matrix_fprintf.
 */
void bct::matrix_printf(const gsl_matrix* m, const char* format) {
	for (int i = 0; i < m->size1; i++) {
		for (int j = 0; j < m->size2; j++) {
			std::printf(format, gsl_matrix_get(m, i, j));
			std::printf(" ");
		}
		std::printf("\n");
	}
}
