#include "bct_test.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <octave/oct.h>

DEFUN_DLD(breadth_cpp, args, , "Wrapper for C++ function.") {
		if (args.length() == 0) {
			return octave_value_list();
		}
		Matrix m = args(0).matrix_value();
		int source = args(1).int_value();
		if (!error_state) {
			gsl_matrix* gslm = bct_test::to_gsl(m);
			return octave_value(bct_test::from_gsl(bct::breadth(gslm, source)));
		} else {
			return octave_value_list();
		}
}
