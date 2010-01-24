#include <bct/bct.h>
#include "bct_test.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <octave/oct.h>

DEFUN_DLD(motif3struct_bin_cpp, args, , "Wrapper for C++ function.") {
	if (args.length() == 0) {
		return octave_value_list();
	}
	Matrix m = args(0).matrix_value();
	if (!error_state) {
		gsl_matrix* gslm = bct_test::to_gsl(m);
		gsl_matrix* F;
		gsl_vector* f = bct::motif3struct_bin(gslm, &F);
		octave_value_list ret;
		ret(0) = octave_value(bct_test::from_gsl(f));
		ret(1) = octave_value(bct_test::from_gsl(F));
		gsl_vector_free(f);
		gsl_matrix_free(F);
		return ret;
	} else {
		return octave_value_list();
	}
}
