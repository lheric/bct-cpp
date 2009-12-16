#ifndef BCT_H
#define BCT_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector.h>
#include "matlab/matlab.h"

namespace bct {
	using namespace matlab;
	
	class bct_exception { };
	class gsl_exception : public bct_exception { };
	class out_of_memory_exception : public bct_exception { };
	class size_exception : public bct_exception { };

	// Density, degree, and assortativity
	double assortativity_dir(const gsl_matrix*);
	double assortativity_und(const gsl_matrix*);
	double assortativity(const gsl_vector*, const gsl_vector*, int);
	gsl_vector* degrees_dir(const gsl_matrix*, gsl_vector* = NULL, gsl_vector* = NULL);
	gsl_vector* degrees_und(const gsl_matrix*);
	double density_dir(const gsl_matrix*);
	double density_und(const gsl_matrix*);
	gsl_matrix* jdegree(const gsl_matrix*); //returns the Joint distribution matrix
	double jdegree_id(gsl_matrix*); //returns other metrices from the Joint dstrbn matrix
	double jdegree_od(gsl_matrix*);
	double jdegree_bl(gsl_matrix*);
	gsl_matrix* matching_ind(const gsl_matrix*); //returns all indices
	gsl_matrix* matching_ind_in(const gsl_matrix*);
	gsl_matrix* matching_ind_out(const gsl_matrix*);
	gsl_vector* strengths_dir(const gsl_matrix*, gsl_vector* = NULL, gsl_vector* = NULL);
	gsl_vector* strengths_und(const gsl_matrix*);

	// Clustering
	gsl_vector* clustering_coef_bd(const gsl_matrix*);
	gsl_vector* clustering_coef_bu(const gsl_matrix*);
	gsl_vector* clustering_coef_wd(const gsl_matrix*);
	gsl_vector* clustering_coef_wu(const gsl_matrix*);

	// Paths, distances, and cycles
	gsl_vector* breadth(const gsl_matrix*, int, gsl_vector* = NULL);
	gsl_matrix* breadthdist(const gsl_matrix*, gsl_matrix* = NULL);
	gsl_vector* cycprob_fcyc(gsl_matrix**, int);
	gsl_vector* cycprob_pcyc(gsl_matrix**, int);
	gsl_matrix* distance_inv(gsl_matrix*); //called by efficiency_*		
	gsl_matrix* distance_wei(const gsl_matrix*); 
	gsl_matrix* efficiency_local(const gsl_matrix*);
	gsl_matrix* efficiency_global(const gsl_matrix*);
	gsl_matrix* findpaths(const gsl_matrix*, const gsl_vector*, int, int, \
						  gsl_matrix** = NULL, long int* = NULL , gsl_vector* = NULL, int* = NULL, gsl_matrix* = NULL);
	gsl_matrix* reachdist(gsl_matrix*, gsl_matrix* = NULL);

	// Centrality
	gsl_vector* betweenness_bin(const gsl_matrix*);
	gsl_vector* betweenness_wei(const gsl_matrix*);
	gsl_matrix* edge_betweenness_bin(const gsl_matrix*);
	gsl_matrix* edge_betweenness_wei(const gsl_matrix*);
	gsl_matrix* erange(const gsl_matrix*, double* = NULL, gsl_matrix* = NULL, double* = NULL);
	void node_and_edge_betweenness_bin(const gsl_matrix*, gsl_vector* = NULL, gsl_matrix* = NULL);
	void node_and_edge_betweenness_wei(const gsl_matrix*, gsl_vector* = NULL, gsl_matrix* = NULL);

	// Motifs
	gsl_matrix* motif3generate(gsl_vector* = NULL, gsl_vector* = NULL, gsl_vector* = NULL);

	// Modularity and community structure

	// Data sets
	extern double cat_all[];
	extern double cat_ctx[];
	extern double fve30[];
	extern double fve32[];
	extern double macaque47[];
	extern double macaque71[];
	
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
	gsl_matrix* permute_columns(const gsl_permutation*, const gsl_matrix*);
	gsl_matrix* permute_rows(const gsl_permutation*, const gsl_matrix*);
	gsl_matrix* yens(int, double);
	gsl_matrix* yens(int, int, double);
	
	// Debugging
	void printf(const gsl_vector*, const char*);
	void printf(const gsl_matrix*, const char*);
};

#endif
