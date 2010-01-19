#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

/*
 * Returns a 'latticized' graph R, with equivalent degree sequence to the original 
 * weighted undirected graph G, and with preserved connectedness 
 * (hence the input graph must be connected).
 *
 * Each edge is rewired (on average) ITER times. The strength distributions 
 * are not preserved for weighted graphs.
 *
 * Rewiring algorithm: Maslov and Sneppen (2002) Science 296:910
 * Latticizing algorithm: Sporns and Zwi (2004); Neuroinformatics 2:145
 */

gsl_matrix* bct::latmio_und_connected(const gsl_matrix* m, const int iters) {
	//It seems D could also be loaded from a file. How to consider that?
	gsl_matrix* R = copy(m);
	int n = m->size1;
	gsl_matrix* D = zeros(n);
	gsl_vector* v1 = sequence(1, 1, n-1); 
	gsl_vector* v2 = sequence(n-1, -1, 1); 
	gsl_matrix* u1 = concatenate_columns(v1, v2); 
	gsl_vector* v3 = min(u1); 
	gsl_vector* u = concatenate(0.0, v3); 
	gsl_vector_free(v1);
	gsl_vector_free(v2);
	gsl_vector_free(v3);
	gsl_matrix_free(u1);
	int upper_lim = (int)ceil((n-1)/2);
	for(int v = 0;v <= upper_lim;v++) {
		gsl_vector* index1 = sequence((v+1), 1, (n-1)); 
		gsl_vector* index2 = sequence(0, 1, v); 
		gsl_vector* u_seg1 = ordinal_index(u, index1); 
		gsl_vector* u_seg2 = ordinal_index(u, index2); 
		gsl_vector* seg_splice = concatenate(u_seg1, u_seg2); 
		gsl_vector_view bottom_row = gsl_matrix_row(D, (n-1-v));
		gsl_vector_memcpy(&bottom_row.vector, seg_splice);
		gsl_vector_view top_row = gsl_matrix_row(D, v);
		gsl_vector* rev_seg_splice = reverse(seg_splice); 
		gsl_vector_memcpy(&top_row.vector, rev_seg_splice);
		gsl_vector_free(index1);
		gsl_vector_free(index2);
		gsl_vector_free(u_seg1);
		gsl_vector_free(u_seg2);
		gsl_vector_free(seg_splice);
		gsl_vector_free(rev_seg_splice);
		
	}
	gsl_vector_free(u);
	
	//[i j]=find(tril(R));
	gsl_matrix* tril_R = tril(R);
	gsl_matrix* R_ij = find_ij(tril_R); 
	gsl_matrix_free(tril_R);
	gsl_vector_view i = gsl_matrix_column(R_ij, 0); 
	gsl_vector_view j = gsl_matrix_column(R_ij, 1); 
	int K = R_ij->size1;
	int tot_iters = iters;
	tot_iters *= K;
	int rewire = 0;
	srand(time(0));
	for(int iter = 1;iter <= tot_iters;iter++) {
		while(1) {
			rewire = 1;
			int a,b,c,d;
			int e1, e2;
			while(1) {
				e1 = ceil((K-1) * (((double)rand())/((double)RAND_MAX)));
				e2 = ceil((K-1) * (((double)rand())/((double)RAND_MAX)));
				while(e2 == e1) {
					e2 = ceil((K-1) * (((double)rand())/((double)RAND_MAX)));
				}
				a = gsl_vector_get(&i.vector, e1);
				b = gsl_vector_get(&j.vector, e1);
				c = gsl_vector_get(&i.vector, e2);
				d = gsl_vector_get(&j.vector, e2);
				//if all(a~=[c d]) && all(b~=[c d]);
				if(a != c && a != d && b!=c && b!=d) {  //all four vertices must be different
					break;
				}
			}
			double val = ((double)rand())/((double)RAND_MAX);
			//flip edge c-d with 50% probability to explore all potential rewirings
			if(val > 0.5) {
	            gsl_vector_set(&i.vector, e2, d);
	            gsl_vector_set(&j.vector, e2, c);
	            c = gsl_vector_get(&i.vector, e2); 
	            d = gsl_vector_get(&j.vector, e2);
	        }
	        
	        //if ~(R(a,d) || R(c,b))
	        //rewiring condition
	        if(!(gsl_matrix_get(R, a, d) || gsl_matrix_get(R, c, b))) {
	        	//if (D(a,b)+D(c,d))>=(D(a,d)+D(c,b))
	        	double val1 = gsl_matrix_get(D, a, b);
	        	double val2 = gsl_matrix_get(D, c, d);
	        	double val3 = gsl_matrix_get(D, a, d);
	        	double val4 = gsl_matrix_get(D, c, b);
	        	if((val1 + val2) >= (val3 + val4)) { //lattice condition
	        		//if ~(R(a,c) || R(b,d))
	        		if(!(gsl_matrix_get(R, a, c) || gsl_matrix_get(R, b, d))) { //connectedness condition
	        			//P=R([a d],:);
	        			gsl_vector* rows = gsl_vector_alloc(2);
	        			gsl_vector_set(rows, 0, a);
	        			gsl_vector_set(rows, 1, d);
	        			gsl_vector* columns = sequence(0, n-1);
	        			gsl_matrix* P = ordinal_index(m, rows, columns); 
	        			//P(1,b)=0; P(2,c)=0;
	        			gsl_matrix_set(P, 0, b, 0.0);
	        			gsl_matrix_set(P, 1, c, 0.0);
	        			gsl_matrix* PN = copy(P); 
	        			//PN(:,d)=1; PN(:,a)=1;
	        			gsl_vector_view dcol = gsl_matrix_column(PN, d);
	        			gsl_vector_view acol = gsl_matrix_column(PN, a);
	        			gsl_vector_set_all(&dcol.vector, 1.0);
	        			gsl_vector_set_all(&acol.vector, 1.0);
	        			
	        			while(1) {
	        				//P(1,:)=any(R(P(1,:)~=0,:),1);
	        				gsl_vector_view row = gsl_matrix_row(P, 0); 
	        				gsl_vector* row_ind = compare_elements(&row.vector, fp_not_equal, 0.0);  
	        				gsl_vector_free(columns);
	        				columns = sequence(0, n-1); 
	        				gsl_matrix* R_indxd = log_ord_index(R, row_ind, columns); 
	        				gsl_vector* any_row;
	        				if(R_indxd == NULL) { //Refer comments above the defintion of log_ord_index
	        					gsl_matrix* any_row_mat = zeros(1, columns->size);
	        					any_row = to_vector(any_row_mat);
	        					gsl_matrix_free(any_row_mat);
	        				}
	        				else {
		        				any_row = any(R_indxd, 1); 
		        			}
	        				gsl_vector_memcpy(&row.vector, any_row);
	        				gsl_vector_free(row_ind);
	        				if(R_indxd != NULL) {
		        				gsl_matrix_free(R_indxd);
		        			}
	        				gsl_vector_free(any_row);
	        				//P(2,:)=any(R(P(2,:)~=0,:),1);	        				
	        				row = gsl_matrix_row(P, 1);
	        				row_ind = compare_elements(&row.vector, fp_not_equal, 0.0); 
	        				R_indxd = log_ord_index(R, row_ind, columns); 
	        				if(R_indxd == NULL) { //Refer comments above the defintion of log_ord_index
	        					gsl_matrix* any_row_mat = zeros(1, columns->size);
	        					any_row = to_vector(any_row_mat);
	        					gsl_matrix_free(any_row_mat);
	        				}
	        				else {
		        				any_row = any(R_indxd, 1); 
		        			}
	        				gsl_vector_memcpy(&row.vector, any_row);
	        				//P=P.*(~PN);
	        				gsl_matrix* not_PN = logical_not(PN); 
	        				gsl_matrix_mul_elements(P, not_PN);
	        				//if ~all(any(P,2))
	        				gsl_vector* any_P = any(P, 2); 
	        				int val = all(any_P);
	        				if(!val) {
	        					rewire = 0;
	        					break;
	        				}
							//elseif any(any(P(:,[b c])))	        					
		    				else {
		    					gsl_vector_free(rows);
		    					gsl_vector_free(columns);
		    					rows = sequence(0, 1); //P has only 2 rows
		    					columns = gsl_vector_alloc(2);
		    					gsl_vector_set(columns, 0, b);
		    					gsl_vector_set(columns, 1, c);
		    					gsl_matrix* P_indxd = ordinal_index(P, rows, columns); 
		    					gsl_vector* P_any = any(P_indxd);
		    					val = any(P_any);
		    					if(val) {
		    						break;
		    					}
		    					gsl_vector_free(P_any);
		    					gsl_matrix_free(P_indxd);
		    				}
		    				gsl_matrix_add(PN, P);
		    				gsl_vector_free(row_ind);
		    				if(R_indxd != NULL) {
			    				gsl_matrix_free(R_indxd);
			    			}
		    				gsl_vector_free(any_row);
		    				gsl_matrix_free(not_PN);
	        			}
	        			gsl_vector_free(rows);
						gsl_vector_free(columns);
						gsl_matrix_free(P);
						gsl_matrix_free(PN);
	        		}
	        		
		    		if(rewire) {
		    			//R(a,d)=R(a,b); R(a,b)=0;
		    			gsl_matrix_set(R, a, d, gsl_matrix_get(R, a, b));
		    			gsl_matrix_set(R, a, b, 0.0);
			    		//R(d,a)=R(b,a); R(b,a)=0;
			    		gsl_matrix_set(R, d, a, gsl_matrix_get(R, b, a));
		    			gsl_matrix_set(R, b, a, 0.0);
	                	//R(c,b)=R(c,d); R(c,d)=0;
	                	gsl_matrix_set(R, c, b, gsl_matrix_get(R, c, d));
		    			gsl_matrix_set(R, c, d, 0.0);
	                	//R(b,c)=R(d,c); R(d,c)=0;
	                	gsl_matrix_set(R, b, c, gsl_matrix_get(R, d, c));
		    			gsl_matrix_set(R, d, c, 0.0);
		    			gsl_vector_set(&j.vector, e1, d);
		    			gsl_vector_set(&j.vector, e2, b);
		    			break;
		    		}
	        	}
			}
		}	 
	}
	gsl_matrix_free(R_ij);
	gsl_matrix_free(D);
	return R;
}      	
	        					
	        				
	        	
	            
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
		
	
				










