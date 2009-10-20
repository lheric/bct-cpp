#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

gsl_vector* bct::strengths_dir(const gsl_matrix* m, gsl_vector* in_strength, gsl_vector* out_strength) {
	bool free_in_strength = false;
	bool free_out_strength = false;
	
	if (in_strength == NULL) {
		free_in_strength = true;
		in_strength = gsl_vector_alloc(m->size1);
	}
	if (out_strength == NULL) {
		free_out_strength = true;
		out_strength = gsl_vector_alloc(m->size1);
	}
	
	gsl_vector* strength = gsl_vector_calloc(m->size1);
	gsl_vector* instrength = sum(m, 1);
	gsl_vector* outstrength = sum(m, 2);
	gsl_vector_add(strength, instrength);
	gsl_vector_add(strength, outstrength);
	
	if (free_in_strength)
		gsl_vector_free (in_strength);
	
	if (free_out_strength)
		gsl_vector_free (out_strength);
		
	return strength;
}
