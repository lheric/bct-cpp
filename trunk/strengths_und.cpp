#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

gsl_vector* bct::strengths_und(const gsl_matrix* m) {
	return sum(m);
}
