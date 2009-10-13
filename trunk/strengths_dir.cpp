#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

gsl_vector* bct::strengths_dir(const gsl_matrix* m) {
	gsl_vector* strength = gsl_vector_calloc(m->size1);
	gsl_vector* instrength = bct::sum(m, 1);
	gsl_vector* outstrength = bct::sum(m, 2);
	gsl_vector_add(strength, instrength);
	gsl_vector_add(strength, outstrength);
	return strength;
}
