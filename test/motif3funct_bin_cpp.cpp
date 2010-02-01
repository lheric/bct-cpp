#include <bct/bct.h>
#include "bct_test.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <octave/oct.h>

DEFUN_DLD(motif3funct_bin_cpp, args, , "Wrapper for C++ function.") {
	bct::set_motif_convention(bct::SPORNS);
	if (args.length() != 1) {
		return octave_value_list();
	}
	Matrix W = args(0).matrix_value();
	if (!error_state) {
		gsl_matrix* W_gsl = bct_test::to_gslm(W);
		gsl_vector* f;
		gsl_matrix* F = bct::motif3funct_bin(W_gsl, &f);
		octave_value_list ret;
		ret(0) = octave_value(bct_test::from_gsl(f));
		ret(1) = octave_value(bct_test::from_gsl(F));
		gsl_matrix_free(W_gsl);
		gsl_vector_free(f);
		gsl_matrix_free(F);
		return ret;
	} else {
		return octave_value_list();
	}
}
