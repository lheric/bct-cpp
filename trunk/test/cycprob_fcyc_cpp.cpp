#include <bct/bct.h>
#include "bct_test.h"

DEFUN_DLD(cycprob_fcyc_cpp, args, , "Wrapper for C++ function.") {
	if (args.length() == 0) {
		return octave_value_list();
	}
	NDArray arr = args(0).array_value();
	int gslm_len = args(1).int_value();
	if (!error_state) {
		gsl_matrix** gslm = bct_test::to_gsl(arr);
		gsl_vector* ret_v = bct::cycprob_fcyc(gslm, gslm_len);
		for(int i =0;i < gslm_len;i++) {
			gsl_matrix_free(gslm[i]);
		}
		return octave_value(bct_test::from_gsl(ret_v));
	} else {
		return octave_value_list();
	}
}

