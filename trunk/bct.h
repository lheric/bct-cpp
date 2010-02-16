#ifndef BCT_H
#define BCT_H

#ifdef SWIG
#define BCT_SWIG
%module bct
#endif

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector.h>
#include "matlab/matlab.h"
#include <vector>

#ifdef BCT_SWIG
#include <Python.h>
#endif

namespace bct {
#ifndef BCT_SWIG
	using namespace matlab;
#endif
	
	class bct_exception { };
	class gsl_exception : public bct_exception { };
	class out_of_memory_exception : public bct_exception { };
	class size_exception : public bct_exception { };

	// Density, degree, and assortativity
	double assortativity_dir(const gsl_matrix*);
	double assortativity_und(const gsl_matrix*);
	gsl_vector* degrees_dir(const gsl_matrix*, gsl_vector** = NULL, gsl_vector** = NULL);
	gsl_vector* degrees_und(const gsl_matrix*);
	double density_dir(const gsl_matrix*);
	double density_und(const gsl_matrix*);
	gsl_matrix* jdegree(const gsl_matrix*);
	int jdegree_bl(const gsl_matrix*);
	int jdegree_id(const gsl_matrix*);
	int jdegree_od(const gsl_matrix*);
	gsl_matrix* matching_ind(const gsl_matrix*);
	gsl_matrix* matching_ind_in(const gsl_matrix*);
	gsl_matrix* matching_ind_out(const gsl_matrix*);
	gsl_vector* strengths_dir(const gsl_matrix*, gsl_vector** = NULL, gsl_vector** = NULL);
	gsl_vector* strengths_und(const gsl_matrix*);

	// Clustering
	gsl_vector* clustering_coef_bd(const gsl_matrix*);
	gsl_vector* clustering_coef_bu(const gsl_matrix*);
	gsl_vector* clustering_coef_wd(const gsl_matrix*);
	gsl_vector* clustering_coef_wu(const gsl_matrix*);
	gsl_vector* efficiency_local(const gsl_matrix*);

	// Paths, distances, and cycles
	gsl_vector* breadth(const gsl_matrix*, int, gsl_vector** = NULL);
	gsl_matrix* breadthdist(const gsl_matrix*, gsl_matrix** = NULL);
	gsl_vector* charpath_ecc(const gsl_matrix*, double* = NULL, double* = NULL);
	double charpath_lambda(const gsl_matrix*);
	gsl_vector* cycprob_fcyc(const std::vector<gsl_matrix*>&);
	gsl_vector* cycprob_pcyc(const std::vector<gsl_matrix*>&);
	gsl_matrix* distance_bin(const gsl_matrix*);
	gsl_matrix* distance_wei(const gsl_matrix*); 
	gsl_matrix* efficiency_global(const gsl_matrix*);
	std::vector<gsl_matrix*> findpaths(const gsl_matrix*, const gsl_vector*, int, gsl_vector** = NULL, int* = NULL, gsl_matrix** = NULL, gsl_matrix** = NULL);
	std::vector<gsl_matrix*> findwalks(const gsl_matrix*, gsl_vector** = NULL);
	double norm_avr_shortest_path_length_bin(const gsl_matrix*);
	double norm_avr_shortest_path_length_wei(const gsl_matrix*, double = 1.0);
	gsl_matrix* reachdist(const gsl_matrix*, gsl_matrix** = NULL);

	// Centrality
	gsl_vector* betweenness_bin(const gsl_matrix*);
	gsl_vector* betweenness_wei(const gsl_matrix*);
	gsl_matrix* edge_betweenness_bin(const gsl_matrix*, gsl_vector** = NULL);
	gsl_matrix* edge_betweenness_wei(const gsl_matrix*, gsl_vector** = NULL);
	gsl_matrix* erange(const gsl_matrix*, double* = NULL, gsl_matrix** = NULL, double* = NULL);

	// Motifs
	enum motif_mode_enum { MILO, SPORNS };
	extern motif_mode_enum motif_mode;
	motif_mode_enum get_motif_mode();
	void set_motif_mode(motif_mode_enum);
	gsl_matrix* motif3funct_bin(const gsl_matrix*, gsl_vector** = NULL);
	gsl_matrix* motif3funct_wei(const gsl_matrix*, gsl_matrix** = NULL, gsl_matrix** = NULL);
	gsl_matrix* motif3generate(gsl_vector** = NULL, gsl_vector** = NULL);
	gsl_matrix* motif3struct_bin(const gsl_matrix*, gsl_vector** = NULL);
	gsl_matrix* motif3struct_wei(const gsl_matrix*, gsl_matrix** = NULL, gsl_matrix** = NULL);
	gsl_matrix* motif4generate(gsl_vector** = NULL, gsl_vector** = NULL);

	// Modularity and community structure
	
	// Synthetic connection networks
	
	// Graph randomization
	gsl_matrix* latmio_dir(const gsl_matrix*, int);
	gsl_matrix* latmio_dir_connected(const gsl_matrix*, int);	
	gsl_matrix* latmio_und(const gsl_matrix*, int);
	gsl_matrix* latmio_und_connected(const gsl_matrix*, int);
	gsl_matrix* randmio_dir(const gsl_matrix*, int);
	gsl_matrix* randmio_dir_connected(const gsl_matrix*, int);
	gsl_matrix* randmio_und(const gsl_matrix*, int);
	gsl_matrix* randmio_und_connected(const gsl_matrix*, int);

	// Data sets
	extern double cat_all[95 * 95];
	extern double cat_ctx[52 * 52];
	extern double fve30[30 * 30];
	extern double fve32[32 * 32];
	extern double macaque47[47 * 47];
	extern double macaque71[71 * 71];
	
	// Matrix status checking
	enum status {
		UNDIRECTED = 1, DIRECTED = 2,
		BINARY = 4, WEIGHTED = 8,
		POSITIVE = 16, SIGNED = 32,
		NO_LOOPS = 64, LOOPS = 128
	};
	extern bool safe_mode;
	bool get_safe_mode();
	void set_safe_mode(bool);
	bool check_status(const gsl_matrix*, int, const char* = NULL);
	bool is_undirected(const gsl_matrix*);
	bool is_directed(const gsl_matrix*);
	bool is_binary(const gsl_matrix*);
	bool is_weighted(const gsl_matrix*);
	bool is_positive(const gsl_matrix*);
	bool is_signed(const gsl_matrix*);
	bool has_loops(const gsl_matrix*);
	bool has_no_loops(const gsl_matrix*);
	
	// Matrix conversion
	gsl_matrix* remove_loops(const gsl_matrix*);
	gsl_matrix* to_binary(const gsl_matrix*);
	gsl_matrix* to_positive(const gsl_matrix*);
	gsl_matrix* to_undirected(const gsl_matrix*);
	
	// Utility
	void gsl_error_handler(const char*, const char*, int, int);
	void gsl_free(gsl_vector*);
	void gsl_free(gsl_matrix*);
	void init();
	
	// Debugging
	void printf(const gsl_vector*, const char*);
	void printf(const gsl_matrix*, const char*);
	void printf(const gsl_permutation*, const char*);
	
	// SWIG
#ifdef BCT_SWIG
	PyObject* from_gsl(const gsl_vector*);
	PyObject* from_gsl(const gsl_matrix*);
	gsl_vector* to_gslv(PyObject*);
	gsl_matrix* to_gslm(PyObject*);
#endif
};

#endif
