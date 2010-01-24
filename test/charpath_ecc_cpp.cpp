#include "bct_test.h"

DEFUN_DLD(charpath_ecc_cpp, args, , "Wrapper for C++ function.") {
	if (args.length() == 0) {
		return octave_value_list();
	}
	Matrix m = args(0).matrix_value();
	double *radius = new double;
	double *diameter = new double;
	octave_value_list retvals;
	if (!error_state) {
		gsl_matrix* gslm = bct_test::to_gsl(m);
		gsl_vector* ret_v = bct::charpath_ecc(gslm, radius, diameter);
		retvals(0) = octave_value(bct_test::from_gsl(ret_v));
		retvals(1) = octave_value(*radius);
		retvals(2) = octave_value(*diameter);		
		return retvals;
	} else {
		return octave_value_list();
	}
}

