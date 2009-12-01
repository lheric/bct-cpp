#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <limits>

/*
 * Computes node betweenness centrality for a binary graph.  Results are
 * returned in a vector where each element is the betweenness centrality of the
 * corresponding node.
 */
gsl_vector* bct::betweenness_bin(const gsl_matrix* m) {
	gsl_vector* betweenness = gsl_vector_alloc(m->size1);
	node_and_edge_betweenness_bin(m, betweenness, NULL);
	return betweenness;
}

/*
 * Computes edge betweenness centrality for a binary graph.  Results are
 * returned in a matrix where each element is the betweenness centrality of the
 * corresponding edge.
 */
gsl_matrix* bct::edge_betweenness_bin(const gsl_matrix* m) {
	gsl_matrix* edge_betweenness = gsl_matrix_alloc(m->size1, m->size2);
	node_and_edge_betweenness_bin(m, NULL, edge_betweenness);
	return edge_betweenness;
}

/*
 * Computes node and edge betweenness centrality for a binary graph.  Results
 * are stored in a vector (node betweenness) and a matrix (edge betweenness)
 * that are provided by the caller.
 */
void bct::node_and_edge_betweenness_bin(const gsl_matrix* m, gsl_vector* node_betweenness, gsl_matrix* edge_betweenness) {
	if (safe_mode) check_status(m, BINARY, "node_and_edge_betweenness_bin");
	if (m->size1 != m->size2) {
		throw size_exception();
	}	
	bool free_node_betweenness = false;
	bool free_edge_betweenness = false;
	
	// BC=zeros(n,1);
	if (node_betweenness == NULL) {
		free_node_betweenness = true;
		node_betweenness = gsl_vector_calloc(m->size1);
	} else {
		gsl_vector_set_zero(node_betweenness);
	}
	
	// EBC=zeros(n);
	if (edge_betweenness == NULL) {
		free_edge_betweenness = true;
		edge_betweenness = gsl_matrix_calloc(m->size1, m->size2);
	} else {
		gsl_matrix_set_zero(edge_betweenness);
	}
	
	// for u=1:n
	for (int u = 0; u < m->size1; u++) {
		
		// D=false(1,n); D(u) = 1;
		gsl_vector* d = gsl_vector_calloc(m->size1);
		gsl_vector_set(d, u, 1.0);
		
		// NP=zeros(1,n); NP(u)=1;
		gsl_vector* np = gsl_vector_calloc(m->size1);
		gsl_vector_set(np, u, 1.0);
		
		// P=false(n);
		gsl_matrix* p = gsl_matrix_calloc(m->size1, m->size2);
		
		// Q=zeros(1,n); q=n;
		gsl_vector* Q = gsl_vector_calloc(m->size1);
		int q = m->size1 - 1;
		
		// Gu=G;
		gsl_matrix* copy_m = copy(m);
		
		// V=u;
		gsl_vector* V = gsl_vector_alloc(1);
		gsl_vector_set(V, 0, u);
		
		// while V
		while (V != NULL) {
			
			// Gu(:,V)=0;
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
				
				// W=find(Gu(v,:));
				gsl_vector_view copy_m_row = gsl_matrix_row(copy_m, v);
				gsl_vector* W = find(&copy_m_row.vector);
				
				// for w=W
				for (int W_index = 0; W != NULL && W_index < W->size; W_index++) {
					int w = (int)gsl_vector_get(W, W_index);
					
					// if D(w)
					if (fp_nonzero(gsl_vector_get(d, w))) {
						
						// NP(w)=NP(w)+NP(v);
						double npw = gsl_vector_get(np, w);
						double npv = gsl_vector_get(np, v);
						gsl_vector_set(np, w, npw + npv);
						
						// P(w,v)=1;
						gsl_matrix_set(p, w, v, 1.0);
					
					// else
					} else {
						
						// D(w)=1;
						gsl_vector_set(d, w, 1.0);
						
						// NP(w)=NP(v);
						gsl_vector_set(np, w, gsl_vector_get(np, v));

						// P(w,v)=1;
						gsl_matrix_set(p, w, v, 1.0);
					}
				}
				if (W != NULL) {
					gsl_vector_free(W);
				}
			}
			
			// V=find(any(Gu(V,:),1));
			gsl_vector* column_indices = sequence(0, m->size2 - 1);
			gsl_matrix* copy_m_rows = index(copy_m, V, column_indices);
			gsl_vector* any_copy_m_rows = any(copy_m_rows);
			gsl_vector_free(V);
			V = find(any_copy_m_rows);
			gsl_vector_free(any_copy_m_rows);
			gsl_matrix_free(copy_m_rows);
			gsl_vector_free(column_indices);
		}
		
		// if ~all(D)
		if (all(d) == 0) {
			
			// Q(1:q)=find(~D);
			gsl_vector* not_d = logical_not(d);
			gsl_vector* not_d_indices = find(not_d);
			gsl_vector_view Q_upto_q = gsl_vector_subvector(Q, 0, q + 1);
			gsl_vector_memcpy(&Q_upto_q.vector, not_d_indices);
			gsl_vector_free(not_d_indices);
			gsl_vector_free(not_d);
		}
		
		// DP=zeros(n,1);
		gsl_vector* dp = gsl_vector_calloc(m->size1);
		
		// for w=Q(1:n-1);
		for (int Q_index = 0; Q_index < m->size1 - 1; Q_index++) {
			int w = (int)gsl_vector_get(Q, Q_index);
			
			// BC(w)=BC(w)+DP(w)
			double bcw = gsl_vector_get(node_betweenness, w);
			double dpw = gsl_vector_get(dp, w);
			gsl_vector_set(node_betweenness, w, bcw + dpw);
			
			// for v=find(P(w,:))
			gsl_vector_view p_row = gsl_matrix_row(p, w);
			gsl_vector* found_p_row = find(&p_row.vector);
			for (int p_index = 0; found_p_row != NULL && p_index < found_p_row->size; p_index++) {
				int v = (int)gsl_vector_get(found_p_row, p_index);
				
				// DPvw=(1+DP(w)).*NP(v)./NP(w);
				double npv = gsl_vector_get(np, v);
				double npw = gsl_vector_get(np, w);
				double dpvw = (1 + dpw) * npv / npw;
				
				// DP(v)=DP(v)+DPvw;
				double dpv = gsl_vector_get(dp, v);
				gsl_vector_set(dp, v, dpv + dpvw);
				
				// EBC(v,w)=EBC(v,w)+DPvw;
				double ebcvw = gsl_matrix_get(edge_betweenness, v, w);
				gsl_matrix_set(edge_betweenness, v, w, ebcvw + dpvw);
			}
			if (found_p_row != NULL) {
				gsl_vector_free(found_p_row);
			}
		}
		
		gsl_vector_free(dp);
		gsl_matrix_free(copy_m);
		gsl_vector_free(Q);
		gsl_matrix_free(p);
		gsl_vector_free(np);
		gsl_vector_free(d);
	}
	if (free_node_betweenness) {
		gsl_vector_free(node_betweenness);
	}
	if (free_edge_betweenness) {
		gsl_matrix_free(edge_betweenness);
	}
}
