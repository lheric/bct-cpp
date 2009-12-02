#include "bct_test.cpp"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <octave/oct.h>

DEFUN_DLD(findpaths_cpp, args, , "Wrapper for C++ function.") {
		if (args.length() == 0) {
			return octave_value_list();
		}
		Matrix m = args(0).matrix_value();
		Array<double> sources = args(1).array_value();
		int path_len_max = args(2).int_value();
		int savepaths = args(3).int_value();
		if (!error_state) {
			gsl_matrix* gslm = bct_test::to_gsl(m);
			gsl_vector* gslv = bct_test::to_gsl(sources);
			return octave_value(bct_test::from_gsl(bct::findpaths(gslm, gslv, path_len_max, savepaths)));
		} else {
			return octave_value_list();
		}
}
