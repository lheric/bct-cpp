#include <bct/bct.h>
#include "bct_test.h"
#include <gsl/gsl_matrix.h>
#include <octave/oct.h>

DEFUN_DLD(motif3struct_wei_cpp, args, , "Wrapper for C++ function.") {
	bct::set_motif_convention(bct::SPORNS);
	if (args.length() != 1) {
		return octave_value_list();
	}
	Matrix m = args(0).matrix_value();
	if (!error_state) {
		gsl_matrix* gslm = bct_test::to_gsl(m);
		gsl_matrix* I;
		gsl_matrix* Q;
		gsl_matrix* F = bct::motif3struct_wei(gslm, &I, &Q);
		octave_value_list ret;
		ret(0) = octave_value(bct_test::from_gsl(I));
		ret(1) = octave_value(bct_test::from_gsl(Q));
		ret(2) = octave_value(bct_test::from_gsl(F));
		gsl_matrix_free(I);
		gsl_matrix_free(Q);
		gsl_matrix_free(F);
		return ret;
	} else {
		return octave_value_list();
	}
}
