#include <bct/bct.h>
#include "bct_test.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <octave/oct.h>

DEFUN_DLD(breadth_cpp, args, , "Wrapper for C++ function.") {
	if (args.length() != 2) {
		return octave_value_list();
	}
	Matrix m = args(0).matrix_value();
	int source = args(1).int_value() - 1;
	if (!error_state) {
		gsl_matrix* gslm = bct_test::to_gsl(m);
		gsl_vector* branch;
		gsl_vector* distance = bct::breadth(gslm, source, &branch);
		for (int i = 0; i < branch->size; i++) {
			int value = (int)gsl_vector_get(branch, i);
			if (gsl_isinf(gsl_vector_get(distance, i)) == 0 && value != -1) {
				gsl_vector_set(branch, i, (double)value + 1.0);
			}
		}
		octave_value_list ret;
		ret(0) = octave_value(bct_test::from_gsl(distance));
		ret(1) = octave_value(bct_test::from_gsl(branch));
		gsl_vector_free(branch);
		gsl_vector_free(distance);
		return ret;
	} else {
		return octave_value_list();
	}
}
