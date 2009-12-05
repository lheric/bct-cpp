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
		//return values
		octave_value_list retvals;
		gsl_matrix** ret_pq = new gsl_matrix* [path_len_max];
        long int* ret_tpath = NULL;
        gsl_vector* ret_plq = NULL;
        int* ret_qstop = new int;        
        gsl_matrix* ret_util = NULL;
		if (!error_state) {
			gsl_matrix* gslm = bct_test::to_gsl(m);
			gsl_vector* gslv = bct_test::to_gsl(sources);
			gsl_matrix* all_paths = bct::findpaths(gslm, gslv, path_len_max, savepaths, ret_pq, ret_tpath, ret_plq, ret_qstop, ret_util);
			retvals(0) = octave_value(bct_test::from_gsl(all_paths));
			retvals(1) = octave_value(bct_test::from_gsl(ret_pq, *ret_qstop));
			retvals(2) = octave_value(*ret_qstop);
			for(int i = 0;i < path_len_max;i++) {
				gsl_matrix_free(ret_pq[i]);
			}
			return retvals;
		} else {
			return octave_value_list();
		}
}
