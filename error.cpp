#include "bct.h"
#include <gsl/gsl_errno.h>

void bct::gsl_error_handler(const char* reason, const char* file, int line, int gsl_errno) {
	if (gsl_errno == GSL_ENOMEM) {
		throw out_of_memory_exception();
	} else {
		throw gsl_exception();
	}
}
