#ifndef BCT_TEST_H
#define BCT_TEST_H

#include <bct.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <octave/oct.h>
#include <cstdio>

namespace bct_test {
	gsl_matrix* to_gsl(const Matrix);
	gsl_vector* to_gsl(const Array<int>);
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

#define FINDPATHS_FUNCTION(function_name) \
	DEFUN_DLD(function_name##_cpp, args, , "Wrapper for C++ function.") { \
		if (args.length() == 0) { \
			return octave_value_list(); \
		} \
		Matrix m = args(0).matrix_value(); \
		Array<double> sources = args(1).array_value(); \
		int path_len_max = args(2).int_value(); \
		int savepaths = args(3).int_value(); \
		if (!error_state) { \
			gsl_matrix* gslm = bct_test::to_gsl(m); \
			gsl_vector* gslv = bct_test::to_gsl(sources); \
			return octave_value(bct_test::from_gsl(bct::function_name(gslm, gslv, path_len_max, savepaths))); \
		} else { \
			return octave_value_list(); \
		} \
	}	

#define BREADTH_FUNCTION(function_name) \
	DEFUN_DLD(function_name##_cpp, args, , "Wrapper for C++ function.") { \
		if (args.length() == 0) { \
			return octave_value_list(); \
		} \
		Matrix m = args(0).matrix_value(); \
		int source = args(1).int_value(); \
		if (!error_state) { \
			gsl_matrix* gslm = bct_test::to_gsl(m); \
			return octave_value(bct_test::from_gsl(bct::function_name(gslm, source))); \
		} else { \
			return octave_value_list(); \
		} \
	}		

#endif
