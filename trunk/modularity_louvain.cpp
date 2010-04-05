#define GatherLoopCountData 0
#define InnerLoopCountMax 30

#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <vector>

/*
 * Detects communities in an undirected graph via Louvain modularity.  While the
 * MATLAB version returns intermediate values for community numbering and
 * modularity, the C++ version returns only the final values.
 */
double bct::modularity_und_louvain(const gsl_matrix* W, gsl_vector** Ci) {
	if (safe_mode) check_status(W, SQUARE | UNDIRECTED, "modularity_und_louvain");

	double modularity;
	
	while (true) {
		bool success = bct::_modularity_und_louvain(W, Ci, &modularity);
		if (success)
			break;
	}
	
	return modularity;
}

bool bct::_modularity_und_louvain(const gsl_matrix* W, gsl_vector** Ci, double* modularity) {
#if GatherLoopCountData
	long outer_loop_count = 0;
	static FILE* file = NULL;
	if (!file) {
		file = fopen("louvain_loop_data.txt", "w");
		if (!file) {
			fprintf(stderr, "Unable to open louvain loop data file");
			exit(1);
		}
	}
#endif

	long inner_loop_count_total = 0;

	// n=length(W);
	int n = length(W);
	
	// s=sum(W(:));
	gsl_vector* sum_W = sum(W);
	double s = sum(sum_W);
	gsl_vector_free(sum_W);
	
	// h=1;
	int h = 0;
	
	// Ci{h}=1:n;
	std::vector<gsl_vector*> _Ci;
	_Ci.push_back(sequence_double(1, n));
	
	// Q{h}=-1;
	std::vector<double> Q;
	Q.push_back(-1.0);
	
	// n0=n;
	int n0 = n;
	
	gsl_matrix* _W = copy(W);
	
	// while true
	while (true) {
		
	#if GatherLoopCountData
		outer_loop_count++;
	#endif
	
		// K=sum(W);
		gsl_vector* K = sum(_W);
		
		// Km=K;
		gsl_vector* Km = copy(K);
		
		// Knm=W;
		gsl_matrix* Knm = copy(_W);
		
		// M=1:n;
		gsl_vector* M = sequence_double(0, n - 1);
		
		// Nm=ones(1,n);
		gsl_vector* Nm = ones_vector_double(n);
		
		// flag=true;
		bool flag = true;
		
		// while flag
		while (flag) {
			
			inner_loop_count_total++;
			if (inner_loop_count_total > InnerLoopCountMax) {
				return false;	// no success
			}

			// flag=false;
			flag = false;
			
			// for i=randperm(n)
			gsl_permutation* randperm_n = randperm(n);
			for (int i_randperm_n = 0; i_randperm_n < n; i_randperm_n++) {
				int i = gsl_permutation_get(randperm_n, i_randperm_n);
				
				// dQ=(Knm(i,:)-Knm(i,M(i))+W(i,i)) - K(i).*(Km-Km(M(i))+K(i))/s;
				gsl_vector* dQ1 = gsl_vector_alloc(Knm->size2);
				gsl_matrix_get_row(dQ1, Knm, i);
				int M_i = (int)gsl_vector_get(M, i);
				gsl_vector_add_constant(dQ1, -gsl_matrix_get(Knm, i, M_i));
				gsl_vector_add_constant(dQ1, gsl_matrix_get(_W, i, i));
				gsl_vector* dQ2 = copy(Km);
				gsl_vector_add_constant(dQ2, -gsl_vector_get(Km, M_i));
				gsl_vector_add_constant(dQ2, gsl_vector_get(K, i));
				gsl_vector_scale(dQ2, gsl_vector_get(K, i) / s);
				gsl_vector* dQ = dQ1;
				gsl_vector_sub(dQ, dQ2);
				gsl_vector_free(dQ2);
				
				// dQ(M(i))=0;
				gsl_vector_set(dQ, M_i, 0.0);
				
				// max_dQ=max(dQ);
				double max_dQ = max(dQ);
				
				// if max_dQ>0;
				if (max_dQ > 0.0) {
					
					// j=find(dQ==max_dQ,1);
					gsl_vector* dQ_eq_max_dQ = compare_elements(dQ, fp_equal, max_dQ);
					gsl_vector* j_v = find(dQ_eq_max_dQ, 1);
					gsl_vector_free(dQ_eq_max_dQ);
					int j = (int)gsl_vector_get(j_v, 0);
					gsl_vector_free(j_v);
					
					// Knm(:,j)=Knm(:,j)+W(:,i);
					gsl_vector_view Knm_col_j = gsl_matrix_column(Knm, j);
					gsl_vector_view _W_col_i = gsl_matrix_column(_W, i);
					gsl_vector_add(&Knm_col_j.vector, &_W_col_i.vector);
					
					// Knm(:,M(i))=Knm(:,M(i))-W(:,i);
					gsl_vector_view Knm_col_M_i = gsl_matrix_column(Knm, M_i);
					gsl_vector_sub(&Knm_col_M_i.vector, &_W_col_i.vector);
					
					// Km(j)=Km(j)+K(i);
					gsl_vector_set(Km, j, gsl_vector_get(Km, j) + gsl_vector_get(K, i));
					
					// Km(M(i))=Km(M(i))-K(i);
					gsl_vector_set(Km, M_i, gsl_vector_get(Km, M_i) - gsl_vector_get(K, i));
					
					// Nm(j)=Nm(j)+1;
					gsl_vector_set(Nm, j, gsl_vector_get(Nm, j) + 1.0);
					
					// Nm(M(i))=Nm(M(i))-1;
					gsl_vector_set(Nm, M_i, gsl_vector_get(Nm, M_i) - 1.0);
					
					// M(i)=j;
					gsl_vector_set(M, i, (double)j);
					
					// flag=true;
					flag = true;
				}
				
				gsl_vector_free(dQ);
			}
			
			gsl_permutation_free(randperm_n);
		}
		
		gsl_vector_free(K);
		gsl_vector_free(Km);
		gsl_matrix_free(Knm);
		gsl_vector_free(Nm);
		
		// [x x M1]=unique(M);
		gsl_vector* x1;
		gsl_vector* M1;
		gsl_vector* x2 = unique(M, "last", &x1, &M1);
		gsl_vector_free(M);
		gsl_vector_free(x1);
		gsl_vector_free(x2);
		
		// h=h+1;
		h++;
		
		// Ci{h}=zeros(1,n0);
		_Ci.push_back(zeros_vector_double(n0));
		
		// for i=1:n
		for (int i = 0; i < n; i++) {
			
			// Ci{h}(Ci{h-1}==i)=M1(i);
			gsl_vector* _Ci_h_sub_1_eq_i_add_1 = compare_elements(_Ci[h - 1], fp_equal, (double)(i + 1));
			logical_index_assign(_Ci[h], _Ci_h_sub_1_eq_i_add_1, gsl_vector_get(M1, i) + 1.0);
			gsl_vector_free(_Ci_h_sub_1_eq_i_add_1);
		}
		
		// n=max(M1);
		n = (int)max(M1) + 1;
		
		// W1=zeros(n);
		gsl_matrix* _W1 = gsl_matrix_alloc(n, n);
		
		// for i=1:n
		for (int i = 0; i < n; i++) {
			
			// for j=1:n
			for (int j = 0; j < n; j++) {
				
				// w=sum(sum(W(M1==i,M1==j)));
				gsl_vector* M1_eq_i = compare_elements(M1, fp_equal, (double)i);
				gsl_vector* M1_eq_j = compare_elements(M1, fp_equal, (double)j);
				gsl_matrix* _W_idx = logical_index(_W, M1_eq_i, M1_eq_j);
				gsl_vector_free(M1_eq_i);
				gsl_vector_free(M1_eq_j);
				gsl_vector* sum__W_idx = sum(_W_idx);
				gsl_matrix_free(_W_idx);
				double w = sum(sum__W_idx);
				gsl_vector_free(sum__W_idx);
				
				// W1(i,j)=w;
				gsl_matrix_set(_W1, i, j, w);
				
				// W1(j,i)=w;
				gsl_matrix_set(_W1, j, i, w);
			}
		}
		
		gsl_vector_free(M1);
		
		// W=W1;
		gsl_matrix_free(_W);
		_W = _W1;
		
		// Q{h}=sum(diag(W))/s-sum(sum((W/s)^2));
		gsl_vector* diag__W = diag(_W);
		double sum_diag__W = sum(diag__W);
		gsl_vector_free(diag__W);
		gsl_matrix* _W_div_s = copy(_W);
		gsl_matrix_scale(_W_div_s, 1.0 / s);
		gsl_matrix* _W_div_s_pow_2 = pow(_W_div_s, 2);
		gsl_matrix_free(_W_div_s);
		gsl_vector* sum__W_div_s_pow_2 = sum(_W_div_s_pow_2);
		gsl_matrix_free(_W_div_s_pow_2);
		double sum_sum__W_div_s_pow_2 = sum(sum__W_div_s_pow_2);
		gsl_vector_free(sum__W_div_s_pow_2);
		Q.push_back(sum_diag__W / s - sum_sum__W_div_s_pow_2);
		
		// if Q{h}-Q{h-1}<=eps
		if (fp_less_or_equal(Q[h] - Q[h - 1], epsilon_double)) {
			
			// break
			break;
		}
	}
				
#if GatherLoopCountData
	fprintf(file, "%ld %ld\n", outer_loop_count, inner_loop_count_total);
	fflush(file);
#endif

	gsl_matrix_free(_W);
	
	// Ci([1 end])=[];
	for (int i = 0; i < (int)_Ci.size(); i++) {
		if (i != (int)_Ci.size() - 2) {
			gsl_vector_free(_Ci[i]);
		}
	}
	if (Ci == NULL) {
		gsl_vector_free(_Ci[_Ci.size() - 2]);
	} else {
		*Ci = _Ci[_Ci.size() - 2];
	} 
	
	// Q([1 end])=[];
	*modularity = Q[Q.size() - 2];
	
	return true;	// success
}