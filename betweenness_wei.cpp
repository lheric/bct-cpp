#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <limits>

/*
 * Computes the betweenness centrality for a weighted graph.  Results are
 * returned in a vector where each element is the betweenness centrality of the
 * corresponding node.
 */
gsl_vector* bct::betweenness_wei(const gsl_matrix* m) {
	if (safe_mode) check_status(m, WEIGHTED, "betweenness_wei");
	if (m->size1 != m->size2) {
		throw size_exception();
	}
	
	// BC=zeros(n,1);
	gsl_vector* betweenness = gsl_vector_calloc(m->size1);
	
	// for u=1:n
	for (int u = 0; u < m->size1; u++) {
		
		// D=inf(1,n); D(u) = 0;
		gsl_vector* d = gsl_vector_alloc(m->size1);
		gsl_vector_set_all(d, std::numeric_limits<double>::max());
		gsl_vector_set(d, u, 0.0);
		
		// NP=zeros(1,n); NP(u)=1;
		gsl_vector* np = gsl_vector_calloc(m->size1);
		gsl_vector_set(np, u, 1.0);
		
		// S=true(1,n);
		gsl_vector* s = gsl_vector_alloc(m->size1);
		gsl_vector_set_all(s, 1.0);
		
		// P=false(n);
		gsl_matrix* p = gsl_matrix_calloc(m->size1, m->size2);
		
		// Q=zeros(1,n); q=n;
		gsl_vector* Q = gsl_vector_calloc(m->size1);
		int q = m->size1 - 1;
		
		// G1=G;
		gsl_matrix* copy_m = copy(m);
		
		// V=u;
		gsl_vector* V = gsl_vector_alloc(1);
		gsl_vector_set(V, 0, u);
		
		// while 1
		while (true) {
			
			// S(V)=0;
			index_assign(s, V, 0.0);
			
			// G1(:,V)=0;
			for (int V_index = 0; V_index < V->size; V_index++) {
				int v = (int)gsl_vector_get(V, V_index);
				gsl_vector_view copy_m_column = gsl_matrix_column(copy_m, v);
				gsl_vector_set_zero(&copy_m_column.vector);
			}
			
			// for v=V
			for (int V_index = 0; V_index < V->size; V_index++) {
				int v = (int)gsl_vector_get(V, V_index);
				
				// Q(q)=v; q=q-1;
				gsl_vector_set(Q, q--, v);
				
				// W=find(G1(v,:));
				gsl_vector_view copy_m_row = gsl_matrix_row(copy_m, v);
				gsl_vector* W = find(&copy_m_row.vector);
				
				// for w=W
				for (int W_index = 0; W != NULL && W_index < W->size; W_index++) {
					int w = (int)gsl_vector_get(W, W_index);
					
					// Duw=D(v)+G1(v,w);
					int duw = (int)std::numeric_limits<double>::max();
					if (gsl_vector_get(d, v) < std::numeric_limits<double>::max()) {
						duw = (int)(gsl_vector_get(d, v) + gsl_matrix_get(copy_m, v, w));
					}
					
					// if Duw<D(w)
					if (duw < gsl_vector_get(d, w)) {
						
						// D(w)=Duw;
						gsl_vector_set(d, w, duw);
						
						// NP(w)=NP(v);
						gsl_vector_set(np, w, gsl_vector_get(np, v));
						
						// P(w,:)=0;
						gsl_vector_view p_row = gsl_matrix_row(p, w);
						gsl_vector_set_zero(&p_row.vector);
						
						// P(w,v)=1;
						gsl_matrix_set(p, w, v, 1.0);
						
					// elseif Duw==D(w)
					} else if (duw == (int)gsl_vector_get(d, w)) {
						
						// NP(w)=NP(w)+NP(v);
						double npw = gsl_vector_get(np, w);
						double npv = gsl_vector_get(np, v);
						gsl_vector_set(np, w, npw + npv);
						
						// P(w,v)=1;
						gsl_matrix_set(p, w, v, 1.0);
					}
				}
				if (W != NULL) {
					gsl_vector_free(W);
				}
			}
			
			// if isempty(minD), break
			if (nnz(s) == 0) {
				break;
			} else {
			
				// minD=min(D(S))
				gsl_vector* d_s = logical_index(d, s);
				double min_d = gsl_vector_min(d_s);
				gsl_vector_free(d_s);
				
				// elseif isinf(minD),
				if (fp_equal(min_d, std::numeric_limits<double>::max())) {
					
					// Q(1:q)=find(isinf(D)); break
					gsl_vector* isinf_d = compare_elements(d, cmp_equal, std::numeric_limits<double>::max());
					gsl_vector* isinf_d_indices = find(isinf_d);
					gsl_vector_view Q_upto_q = gsl_vector_subvector(Q, 0, q + 1);
					gsl_vector_memcpy(&Q_upto_q.vector, isinf_d_indices);
					gsl_vector_free(isinf_d_indices);
					gsl_vector_free(isinf_d);
					break;
				}
				
				// V=find(D==minD);
				gsl_vector_free(V);
				gsl_vector* d_equal_min_d = compare_elements(d, cmp_equal, min_d);
				V = find(d_equal_min_d);
				gsl_vector_free(d_equal_min_d);
			}
		}
		
		// DP=zeros(n,1);
		gsl_vector* dp = gsl_vector_calloc(m->size1);
		
		// for w=Q(1:n-1);
		for (int Q_index = 0; Q_index < m->size1 - 1; Q_index++) {
			int w = (int)gsl_vector_get(Q, Q_index);
			
			// BC(w)=BC(w)+DP(w)
			double bcw = gsl_vector_get(betweenness, w);
			double dpw = gsl_vector_get(dp, w);
			gsl_vector_set(betweenness, w, bcw + dpw);
			
			// for v=find(P(w,:))
			gsl_vector_view p_row = gsl_matrix_row(p, w);
			gsl_vector* found_p_row = find(&p_row.vector);
			for (int p_index = 0; found_p_row != NULL && p_index < found_p_row->size; p_index++) {
				int v = (int)gsl_vector_get(found_p_row, p_index);
				
				// DP(v)=DP(v)+(1+DP(w)).*NP(v)./NP(w);
				double dpv = gsl_vector_get(dp, v);
				double npv = gsl_vector_get(np, v);
				double npw = gsl_vector_get(np, w);
				gsl_vector_set(dp, v, dpv + (1 + dpw) * npv / npw);
			}
			if (found_p_row != NULL) {
				gsl_vector_free(found_p_row);
			}
		}
		
		gsl_vector_free(dp);
		gsl_vector_free(V);
		gsl_matrix_free(copy_m);
		gsl_vector_free(Q);
		gsl_matrix_free(p);
		gsl_vector_free(s);
		gsl_vector_free(np);
		gsl_vector_free(d);
	}
	return betweenness;
}
