#include <bct/bct.h>
#include "bct_test.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <octave/oct.h>

DEFUN_DLD(findwalks_cpp, args, , "Wrapper for C++ function.") {
		if (args.length() == 0) {
			return octave_value_list();
		}
		Matrix m = args(0).matrix_value();
		int N = args(1).int_value();
		//return values
		octave_value_list retvals;
		gsl_matrix** ret_Wq = new gsl_matrix* [N];
        double* ret_twalk = new double;
		if (!error_state) {
			gsl_matrix* gslm = bct_test::to_gsl(m);
			gsl_vector* wlq = bct::findwalks(gslm, ret_Wq, ret_twalk);
			retvals(0) = octave_value(bct_test::from_gsl(wlq));
			retvals(1) = octave_value(*ret_twalk);			
			retvals(2) = octave_value(bct_test::from_gsl(ret_Wq, N));
			for(int i = 0;i < N;i++) {
				gsl_matrix_free(ret_Wq[i]);
			}
			return retvals;
		} else {
			return octave_value_list();
		}
}
