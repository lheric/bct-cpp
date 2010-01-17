#include "bct_test.h"
#include "bct_test.cpp"
#include <gsl/gsl_matrix.h>
#include <octave/oct.h>

DEFUN_DLD(randmio_dir_cpp, args, , "Wrapper for C++ function.") {
		if (args.length() == 0) {
			return octave_value_list();
		}
		Matrix m = args(0).matrix_value();
		int iters = args(1).int_value();
		//return values
		if (!error_state) {
			gsl_matrix* gslm = bct_test::to_gsl(m);
			gsl_matrix* R = bct::randmio_dir(gslm, iters);
			return octave_value(bct_test::from_gsl(R));
		} else {
			return octave_value_list();
		}
}
