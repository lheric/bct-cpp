#include <bct/bct.h>
#include "bct_test.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <vector>

DEFUN_DLD(cycprob_pcyc_cpp, args, , "Wrapper for C++ function.") {
	if (args.length() != 1) {
		return octave_value_list();
	}
	NDArray Pq = args(0).array_value();
	if (!error_state) {
		std::vector<gsl_matrix*> Pq_gsl = bct_test::to_gsl(Pq, 1);
		gsl_vector* pcyc = bct::cycprob_pcyc(Pq_gsl);
		octave_value ret = octave_value(bct_test::from_gsl(pcyc, 1));
		for (int i = 1; i < (int)Pq_gsl.size(); i++) {
			gsl_matrix_free(Pq_gsl[i]);
		}
		gsl_vector_free(pcyc);
		return ret;
	} else {
		return octave_value_list();
	}
}
