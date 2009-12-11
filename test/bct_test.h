#ifndef BCT_TEST_H
#define BCT_TEST_H

#include <bct/bct.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <octave/oct.h>

namespace bct_test {
	gsl_matrix* to_gsl(const Matrix);
	gsl_vector* to_gsl(const Array<int>);
	Matrix from_gsl(const gsl_matrix*);
	Matrix from_gsl(const gsl_vector*);
	NDArray from_gsl(gsl_matrix**, int);
	gsl_matrix** to_gsl(const NDArray);
};

#define MATRIX_TO_SCALAR_FUNCTION(function_name) \
	DEFUN_DLD(function_name##_cpp, args, , "Wrapper for C++ function.") { \
		if (args.length() == 0) { \
			return octave_value_list(); \
		} \
		Matrix m = args(0).matrix_value(); \
		if (!error_state) { \
			gsl_matrix* gslm = bct_test::to_gsl(m); \
			return octave_value(bct::function_name(gslm)); \
		} else { \
			return octave_value_list(); \
		} \
	}

#define MATRIX_TO_MATRIX_FUNCTION(function_name) \
	DEFUN_DLD(function_name##_cpp, args, , "Wrapper for C++ function.") { \
		if (args.length() == 0) { \
			return octave_value_list(); \
		} \
		Matrix m = args(0).matrix_value(); \
		if (!error_state) { \
			gsl_matrix* gslm = bct_test::to_gsl(m); \
			return octave_value(bct_test::from_gsl(bct::function_name(gslm))); \
		} else { \
			return octave_value_list(); \
		} \
	}

#endif
