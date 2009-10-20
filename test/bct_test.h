#ifndef BCT_TEST_H
#define BCT_TEST_H

#include <bct.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <octave/oct.h>

namespace bct_test {
	gsl_matrix* to_gsl(const Matrix);
	Matrix from_gsl(const gsl_matrix*);
	Matrix from_gsl(const gsl_vector*);
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

#define MATRIX_TO_SCALAR_MULT_PARAMS_FUNCTION(function_name) \
	DEFUN_DLD(function_name##_cpp, args, , "Wrapper for C++ function.") { \
		if (args.length() == 0) { \
			return octave_value_list(); \
		} \
		int int_param; \
		Matrix m = args(0).matrix_value(); \
		int_param = args(1).int_value(); \
		if (!error_state) { \
			gsl_matrix* gslm = bct_test::to_gsl(m); \
			return octave_value(bct::function_name(gslm, int_param)); \
		} else { \
			return octave_value_list(); \
		} \
	}	

#endif
