#include "bct.h"
#include <cstdio>
#include <gsl/gsl_vector.h>

/*
 * Prints a vector using the given format for each element.  This is only
 * provided for debugging purposes.  In other cases, use gsl_vector_fprintf.
 */
void bct::vector_printf(const gsl_vector* v, const char* format) {
	for (int i = 0; i < v->size; i++) {
		std::printf(format, gsl_vector_get(v, i));
		std::printf(" ");
	}
	std::printf("\n");
}
