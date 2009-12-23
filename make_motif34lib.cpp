#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_vector.h>

/*
 * Returns the three-node motif library.
 */
gsl_matrix* bct::motif3generate(gsl_vector* Mn, gsl_vector* ID, gsl_vector* N) {
	
	// n=0;
	int n = -1;
	
	// M=false(54,6);
	gsl_matrix* M = gsl_matrix_calloc(54, 6);
	
	// CL=zeros(54,6,'uint8');
	gsl_matrix* CL = zeros(54, 6);
	
	// cl=zeros(1,6,'uint8');
	gsl_vector* cl = gsl_vector_calloc(6);
	
	// for i=0:2^6-1
	for (int i = 0; i < 64; i++) {
		
		// m=dec2bin(i);
		// m=[num2str(zeros(1,6-length(m)), '%d') m];
		char m[7];
		dec2bin(i, 6, m);
		
		// G=str2num ([ ...
		// '0'   ' '  m(3)  ' '  m(5) ;
		// m(1)  ' '  '0'   ' '  m(6) ;
		// m(2)  ' '  m(4)  ' '  '0'   ]);
		gsl_matrix* G = gsl_matrix_calloc(3, 3);
		for (int j = 0, k = 0; j < 3; j++) {
			for (int i = 0; i < 3; i++) {
				if (i != j && m[k++] == '1') {
					gsl_matrix_set(G, i, j, 1.0);
				}
			}
		}
		
		// Ko=sum(G,2);
		gsl_vector* Ko = sum(G, 2);
		
		// Ki=sum(G,1).';
		gsl_vector* Ki = sum(G, 1);
		
		// if Ko+Ki,
		gsl_vector* Kboth = gsl_vector_alloc(Ko->size);
		gsl_vector_memcpy(Kboth, Ko);
		gsl_vector_add(Kboth, Ki);
		if (to_bool(Kboth)) {
			
			// n=n+1;
			n++;
			
			// cl(:)=sortrows([Ko Ki]).';
			gsl_matrix* Ko_Ki = concatenate_rows(Ko, Ki);
			gsl_matrix* Ko_Ki_sort = sortrows(Ko_Ki);
			gsl_matrix* Ko_Ki_transpose = gsl_matrix_alloc(2, Ko->size);
			gsl_matrix_transpose_memcpy(Ko_Ki_transpose, Ko_Ki_sort);
			gsl_vector_free(cl);
			cl = to_vector(Ko_Ki_transpose);
			gsl_matrix_free(Ko_Ki_transpose);
			gsl_matrix_free(Ko_Ki_sort);
			gsl_matrix_free(Ko_Ki);
			
			// CL(n,:)=cl;
			gsl_matrix_set_row(CL, n, cl);
			
			// M(n,:)=G([2:4 6:8]);
			double indices[] = { 1, 2, 3, 5, 6, 7 };
			gsl_vector_view indices_vv = gsl_vector_view_array(indices, 6);
			gsl_vector* G_portion = ordinal_index(G, &indices_vv.vector);
			gsl_matrix_set_row(M, n, G_portion);
			gsl_vector_free(G_portion);
		}
		gsl_vector_free(Kboth);
		gsl_vector_free(Ki);
		gsl_vector_free(Ko);
		gsl_matrix_free(G);
	}
	
	// [u1 u2 ID]=unique(CL,'rows');
	gsl_vector* temp_ID = gsl_vector_alloc(54);
	gsl_matrix* u1 = unique_rows(CL, "last", NULL, temp_ID);
	gsl_vector_add_constant(temp_ID, 1.0);
	gsl_matrix_free(u1);
	gsl_vector_free(cl);
	gsl_matrix_free(CL);

	// id_mika=  [1  3  4  6  7  8  11];
	int id_mika[] = { 1, 3, 4, 6, 7, 8, 11 };
	
	// id_olaf= -[3  6  1 11  4  7   8];
	int id_olaf[] = { 3, 6, 1, 11, 4, 7, 8};
	
	// %convert IDs into Sporns & Kotter classification
	for (int i_ID = 0; i_ID < temp_ID->size; i_ID++) {
		int ID_value = (int)gsl_vector_get(temp_ID, i_ID);
		for (int i_id = 0; i_id < 7; i_id++) {
			if (ID_value == id_mika[i_id]) {
				gsl_vector_set(temp_ID, i_ID, id_olaf[i_id]);
			}
		}
	}
	
	// [X ind]=sortrows(ID);
	gsl_permutation* ind = gsl_permutation_alloc(54);
	gsl_sort_vector_index(ind, temp_ID);  // TODO: Not stable!
	
	// ID=ID(ind,:);
	if (ID == NULL) {
		gsl_vector_free(temp_ID);
	} else {
		gsl_permute_vector(ind, temp_ID);
		gsl_vector_memcpy(ID, temp_ID);
		gsl_vector_free(temp_ID);
	}
	
	// M=M(ind,:);
	gsl_matrix* permute_M = permute_rows(ind, M);
	gsl_matrix_free(M);
	M = permute_M;
	gsl_permutation_free(ind);
	
	// N=sum(M,2);
	if (N != NULL) {
		gsl_vector* temp_N = sum(M, 2);
		gsl_vector_memcpy(N, temp_N);
		gsl_vector_free(temp_N);
	}
	
	// Mn=uint32(sum(repmat(10.^(5:-1:0),size(M,1),1).*M,2));
	if (Mn != NULL) {
		for (int i = 0; i < M->size1; i++) {
			int Mn_value = 0;
			for (int j = 0; j < M->size2; j++) {
				Mn_value *= 10;
				Mn_value += (int)gsl_matrix_get(M, i, j);
			}
			gsl_vector_set(Mn, i, Mn_value);
		}
	}
	
	return M;
}
