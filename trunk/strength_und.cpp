#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

gsl_vector* bct::strength_und(const gsl_matrix* m) {
	gsl_vector* sum = bct::sum(m,1);
	return sum;
}
