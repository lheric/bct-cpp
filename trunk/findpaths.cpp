#include "bct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_histogram.h>
#include <cstdio>
#include <math.h>

/* Original comments: Paths are sequences of linked nodes, that never visit a single
 * node more than once. This function finds all paths that start at a set of 
 * source vertices, up to a specified length. 
 * Warning: very memory-intensive.
 * See e.g. Sporns (2002). Contributor: OS.
 */
 
/* NOTE: Presently, the findpaths function returns only the 'allpaths' variable */

/* IMPORTANT: Treat very carefully the following index variables:
 * path_len
 * They are different from the way they are treated in matlab (wherever these variables 
 * are used in matlab, they should be decremented by 1 in the corresponding C++ code).
 */
 
/* NOTE: The nodes in this program will range from 0 to N-1, unlike in matlab where they
 * range from 1 to N. Also, just so space will be saved, the corresponding indices for
 * path lengths will range from 0 to path_len_max-1 though the lengths themselves range
 * from 1 to path_len_max.
 */

gsl_matrix* bct::findpaths(const gsl_matrix* m, const gsl_vector* sources_input, int path_len_max, int savepaths, gsl_matrix** ret_pq, long int* ret_tpath, gsl_vector* ret_plq, int* ret_qstop, gsl_matrix* ret_util) {
	
	if (m->size1 != m->size2) {
		throw size_exception();
	}
	
	
	gsl_matrix* CIJ = binary(m);
	gsl_vector* sources = copy(sources_input);
	
	int N = CIJ->size1;
	
	long int tpath = 0;
	
	gsl_matrix* paths = NULL;
	gsl_matrix* util = gsl_matrix_calloc(N, path_len_max);
	
	//all_pathlens array shall be indexed by path_len. Each slot in all_pathlens contains a NxN gsl matrix
	gsl_matrix* node_to_node_len[path_len_max];
	for(int i=0;i < path_len_max;i++) {
		node_to_node_len[i] = gsl_matrix_calloc(N, N);
	}
	
	//Initially 'paths' will hold only all paths of length 1
	int path_len=1;
	//The following 2 variables shall be used for the 'unique' values identification logic
	int unique_nodes_cnt = 0;
	gsl_vector* paths_node_visited = gsl_vector_calloc(N);
	
	gsl_vector_add_constant(sources, -1.0); //decrement by 1 so nodes start from 0
	
	int num_sources = sources->size;
	
	for(int j=0;j < N;j++) {
		for(int i=0;i < sources->size;i++) {
			int source_node = gsl_vector_get(sources, i); 
			if(gsl_matrix_get(CIJ, source_node, j) == 1) {
				gsl_matrix* path = gsl_matrix_alloc(path_len+1, 1);
				gsl_matrix_set(path, 0, 0, source_node);
				gsl_matrix_set(path, path_len, 0, j);
				gsl_matrix* concat_paths = concatenate_rows(paths, path);
				gsl_matrix_free(path);
				if(paths != NULL) {
					gsl_matrix_free(paths); //Hopefully, this will be freed after originally "owned" by 'concat_paths'. See code below.
				}
				paths = concat_paths; 
				//perform the following logic after 'paths' is updated in the big loop below
				int visited = gsl_vector_get(paths_node_visited, j);
				if(visited == 0) {
					unique_nodes_cnt += 1;
					gsl_vector_set(paths_node_visited, j, 1);
				}
			}
		}
	}
	
	//allocate histogram, set ranges, increment, access and update 'util'
	gsl_histogram * hist = gsl_histogram_alloc (N);
	gsl_histogram_set_ranges_uniform(hist, 0, N); //nodes fall in the range [0,N). "["=closed
	//Space between 2 bins is computed as (max-min)/(no. of bins). Hence, max=N so spacing between bins is 1.
	//IMPORTANT: '0' is the first node and NOT '1'
	
	for(int i=0;i < paths->size1;i++) {
		for(int j=0;j < paths->size2;j++) {
			gsl_histogram_increment(hist, gsl_matrix_get(paths, i, j));
		}
	}
	
	for(int i=0;i < util->size1;i++) {
		int hist_value = gsl_histogram_get(hist, i); //0th bin contains the histogram value for node '1'
		gsl_matrix_set(util, i, path_len-1, hist_value); //index (path_len-1) will hold the value for path_len
	}
	
	for(int i=0;i < paths->size2;i++) {
		int from_node = gsl_matrix_get(paths, 0, i);
		int to_node = gsl_matrix_get(paths, path_len, i);
		int count = gsl_matrix_get(node_to_node_len[path_len-1], from_node, to_node);
		count += 1;
		gsl_matrix_set(node_to_node_len[path_len-1], from_node, to_node, count);
	}
	
	gsl_matrix* all_paths;
	if(savepaths == 1) {
		all_paths = gsl_matrix_alloc(paths->size1, paths->size2);
		gsl_matrix_memcpy(all_paths, paths);
	}
	else {
		all_paths = NULL;
	}
	
	//The "big loop"
	for(path_len=2; path_len <= path_len_max;path_len++) {
		long int len_npaths;
		if(N < 40)
			len_npaths = N*path_len*10;
		else //10^(q+1)
			len_npaths = pow(10, path_len);
		
		if(len_npaths > 1900000)
			len_npaths = 1900000; //gsl_matrix_alloc can't seem to allocate (5*2000000)
		gsl_matrix* npaths = gsl_matrix_alloc(path_len+1, len_npaths); //a path of path_len will have (path_len+1) nodes
		//Don't change to 'alloc' as a gsl_matrix_isnull is run on it later
		
		//Collect the unique nodes of the 'destination row' of 'paths'
		gsl_vector* end_nodes = gsl_vector_alloc(unique_nodes_cnt);
		for(int i=0, end_node_index=0;i < N && end_node_index < unique_nodes_cnt;i++) {
			int visited = gsl_vector_get(paths_node_visited, i);
			if(visited == 1) {
				gsl_vector_set(end_nodes, end_node_index++, i);
			}
		}
		
		//paths_node_visited is used, so reset
		gsl_vector_set_zero(paths_node_visited);
		unique_nodes_cnt = 0;
		
		int npathscnt = 0;
		
		for(int i = 0;i < end_nodes->size;i++) {
			int end_node = gsl_vector_get(end_nodes, i);
			gsl_vector_view end_nodes_row = gsl_matrix_row(paths, path_len-1);
			gsl_vector* end_nodes_ind = compare_elements(&end_nodes_row.vector, cmp_equal, end_node);
			gsl_vector* end_nodes_col = find(end_nodes_ind);
			
			//find the nodes connected directly to the end node
			gsl_vector_view next_node_row = gsl_matrix_row(CIJ, end_node);
			gsl_vector* next_node_ind = compare_elements(&next_node_row.vector, cmp_equal, 1.0);
			gsl_vector* next_node_col = find(next_node_ind);
			
			if(next_node_col != NULL) { //'find' method returns NULL if nothing was found
				for(int nn=0;nn < next_node_col->size;nn++) {
					int next_node = gsl_vector_get(next_node_col, nn); //the column indices are also the nodes, in this case
					
					//pb_temp = pb(sum(j==pths(2:q,pb),1)==0);
					gsl_vector* row_indices = sequence(0, path_len-1); //check for NULL return from sequence?
					gsl_matrix* paths_term_at_endnode_all = index(paths, row_indices, end_nodes_col);
					gsl_matrix_view paths_term_at_endnode = gsl_matrix_submatrix(paths_term_at_endnode_all, 1, 0, path_len-1,\
																				 paths_term_at_endnode_all->size2);
					gsl_matrix* temp1 = compare_elements(&paths_term_at_endnode.matrix, cmp_equal, next_node);
					gsl_vector* temp2 = sum(temp1, 1);
					gsl_vector* no_prev_visits_ind = compare_elements(temp2, cmp_equal, 0.0);
					gsl_vector* no_prev_visits_col = find(no_prev_visits_ind);
					
					if(no_prev_visits_col != NULL) {
						//npths(:,npthscnt+1:npthscnt+length(pb_temp)) = [pths(:,pb_temp)' ones(length(pb_temp),1)*j]';
						// Improved version: [pths(:,pb_temp); ones(1,length(pb_temp))*j]
						row_indices = sequence(0, path_len-1);
						gsl_matrix* paths_sub = index(paths_term_at_endnode_all, row_indices, no_prev_visits_col);
					
						gsl_matrix* new_node_matrix = yens(1, no_prev_visits_col->size, (double)next_node);
						gsl_matrix* new_paths = concatenate_columns(paths_sub, new_node_matrix);
						if(npathscnt + no_prev_visits_col->size >= len_npaths) {
							len_npaths += (no_prev_visits_col->size * 4);
							gsl_matrix* npaths_ext = gsl_matrix_alloc(path_len+1, len_npaths);
							gsl_matrix_view npaths_ext_view = gsl_matrix_submatrix(npaths_ext, 0, 0, npaths->size1, npaths->size2);
							gsl_matrix_memcpy(&npaths_ext_view.matrix, npaths);
							gsl_matrix_free(npaths);				
							npaths = npaths_ext;			
						}
						gsl_matrix_view npaths_sub = gsl_matrix_submatrix(npaths, 0, npathscnt, path_len+1, no_prev_visits_col->size);
						gsl_matrix_memcpy(&npaths_sub.matrix, new_paths);
						npathscnt += no_prev_visits_col->size;
					
						gsl_vector_view sources_to_nextnode = gsl_matrix_row(paths_sub, 0);
						//allocate histogram, set ranges, increment, access and update 'util'
						gsl_histogram * hist = gsl_histogram_alloc (N);
						gsl_histogram_set_ranges_uniform(hist, 0, N); //nodes fall in the range (0,N]
						//Space between 2 bins is computed as (max-min)/(no. of bins). Hence, max=N+1 so spacing between bins is 1.
						//IMPORTANT: '0' is the first node and NOT '1'
						for(int pos=0;pos < (&sources_to_nextnode.vector)->size;pos++) {
							gsl_histogram_increment(hist, gsl_vector_get(&sources_to_nextnode.vector, pos));
						}
					
						//Pq(1:N,j,q) = Pq(1:N,j,q)+(hist(pths(1,pb_temp),1:N))';
						for(int source_node=0;source_node < N;source_node++) {
							int hist_count = gsl_histogram_get(hist, source_node);
							int curr_count = gsl_matrix_get(node_to_node_len[path_len-1], source_node, next_node);
							gsl_matrix_set(node_to_node_len[path_len-1], source_node, next_node, curr_count+hist_count);
						}
						gsl_matrix_free(paths_sub);
						gsl_matrix_free(new_node_matrix);
						gsl_matrix_free(new_paths);
						gsl_vector_free(no_prev_visits_col);
						gsl_histogram_free(hist);
					}
					
					gsl_vector_free(row_indices);
					gsl_vector_free(temp2);
					gsl_matrix_free(paths_term_at_endnode_all);
					gsl_matrix_free(temp1);
				}
			}
			
			gsl_vector_free(end_nodes_ind);
			if(end_nodes_col != NULL) {
				gsl_vector_free(end_nodes_col);
			}
			gsl_vector_free(next_node_ind);
			if(next_node_col != NULL) {
				gsl_vector_free(next_node_col);			
			}
		}
		
	    //all_paths is updated from npaths here
	    if (savepaths==1) {
		    gsl_matrix* yens_row = yens(1, all_paths->size2, -1.0);
	    	gsl_matrix* concat_all_paths = concatenate_columns(all_paths, yens_row);
	    	gsl_matrix_view npaths_sub = gsl_matrix_submatrix(npaths, 0, 0, path_len+1, npathscnt);
	    	gsl_matrix* concat_all_paths_npths = concatenate_rows(concat_all_paths, &npaths_sub.matrix);
	    	all_paths = concat_all_paths_npths;	  
	    }
	    
	    //util(1:N,q) = util(1:N,q) + hist(reshape(npths,1,size(npths,1)*size(npths,2)),1:N)' - diag(Pq(:,:,q));
	    gsl_histogram * hist = gsl_histogram_alloc (N);
		gsl_histogram_set_ranges_uniform(hist, 0, N); 
		for(int i=0;i < npaths->size1;i++) {
			for(int j=0;j < npathscnt;j++) {
				gsl_histogram_increment(hist, gsl_matrix_get(npaths, i, j));
			}
		}
		
		for(int node=0;node < util->size1;node++) {
			int hist_value = gsl_histogram_get(hist, node);
			int cycle_count = gsl_matrix_get(node_to_node_len[path_len-1], node, node); //Correct for cycles
			gsl_matrix_set(util, node, path_len-1, (hist_value - cycle_count)); 
		}
		
		//paths is updated from npaths here. After this npaths won't be needed
		int PATHS_UPDATED = 0;
	    if(npathscnt > 0) {
			gsl_matrix_free(paths);
			gsl_matrix* paths_temp = gsl_matrix_alloc(npaths->size1, npathscnt);
			int paths_column_index = 0;
			for(int i=0;i < npathscnt;i++) {
				int first_node = gsl_matrix_get(npaths, 0, i);
				int last_node = gsl_matrix_get(npaths, path_len, i);
				if(first_node != last_node) {
					gsl_vector_view no_cycles_path = gsl_matrix_column(npaths, i);
					gsl_vector_view paths_column = gsl_matrix_column(paths_temp, paths_column_index++);
					gsl_vector_memcpy(&paths_column.vector, &no_cycles_path.vector);
					PATHS_UPDATED = 1;
					int visited = gsl_vector_get(paths_node_visited, last_node);
					if(visited == 0) {
						unique_nodes_cnt += 1;
						gsl_vector_set(paths_node_visited, last_node, 1);
					}
				}
			}
			if(PATHS_UPDATED) {
				paths = gsl_matrix_alloc(paths_temp->size1, paths_column_index);
				gsl_matrix_view paths_copy = gsl_matrix_submatrix(paths_temp, 0, 0, paths->size1, paths->size2);
				gsl_matrix_memcpy(paths, &paths_copy.matrix);
				gsl_matrix_free(paths_temp);
			}
			gsl_matrix_free(npaths);
		}
		else {
			gsl_matrix_free(paths);
			paths = NULL;
		}
		
		//stop finding more paths here. Compute qstop and plq
		if(!PATHS_UPDATED) {
			//tpath = sum(sum(sum(Pq)));
			gsl_matrix* first_sum = NULL;
			for(int i=0;i < path_len;i++) {
				gsl_vector* temp_sum = sum(node_to_node_len[i]);
				gsl_matrix* concat_sum = concatenate_columns(first_sum, temp_sum);
				if(first_sum != NULL) {
					gsl_matrix_free(first_sum);
				}
				gsl_vector_free(temp_sum);
				first_sum = concat_sum;
			}
			gsl_vector* second_sum = sum(first_sum, 2);
			tpath = sum(second_sum);
			gsl_matrix_free(first_sum);
			
			//plq = reshape(sum(sum(Pq)),1,qmax);
			gsl_vector* plq = second_sum; //no need of a reshape, the way second_sum is computed is exactly what is wanted as 'plq'
			
			break; //end of the program, return values after this
		}
				
	} //end of the "big loop"
	
	//End of the program. Compute qstop and plq
	int qstop;
	if(path_len > path_len_max) {
		qstop = path_len_max;
	}
	else {
		qstop = path_len;
	}
	
	//tpath = sum(sum(sum(Pq)));
	gsl_matrix* first_sum = NULL;
	for(int i=0;i < path_len_max;i++) {
		gsl_vector* temp_sum = sum(node_to_node_len[i]);
		gsl_matrix* concat_sum = concatenate_columns(first_sum, temp_sum);
		if(first_sum != NULL) {
			gsl_matrix_free(first_sum);
		}
		gsl_vector_free(temp_sum);
		first_sum = concat_sum;
	}
	gsl_vector* second_sum = sum(first_sum, 2);
	tpath = sum(second_sum);
	if(first_sum != NULL) {
		gsl_matrix_free(first_sum);
	}
	
	//plq = reshape(sum(sum(Pq)),1,qmax);
	gsl_vector* plq = second_sum; //no need of a reshape, the way second_sum is computed is exactly what is wanted as 'plq'
	
	gsl_matrix_add_constant(all_paths, 1.0); //Revert back to matlab's node range: from (0,N-1) to (1,N)
	
	if(paths != NULL) {
		gsl_matrix_free(paths);
	}
	
	//Assign the local return values to the corresponding pointer parameters and free the local variables
	if(plq != NULL) {
		if(ret_plq == NULL) {
			ret_plq = gsl_vector_alloc(plq->size);
		}
		gsl_vector_memcpy(ret_plq, plq);
		gsl_vector_free(plq);
	}
	
	if(node_to_node_len != NULL) {
		if(ret_pq == NULL) {
			ret_pq = new gsl_matrix* [path_len_max];
			for(int i=0;i < path_len_max;i++) {
				ret_pq[i] = gsl_matrix_alloc(N, N);
			}
		}
		for(int i=0;i < path_len_max;i++) {	
			gsl_matrix_memcpy(ret_pq[i],node_to_node_len[i]);
		}
		for(int i=0;i < path_len_max;i++) {	
			gsl_matrix_free(node_to_node_len[i]);
		}
	}
	
	if(util != NULL) {
		if(ret_util == NULL) {
			ret_util = gsl_matrix_alloc(util->size1, util->size2);
		}
		gsl_matrix_memcpy(ret_util, util);
		gsl_matrix_free(util);
	}	
	
	if(ret_tpath == NULL) {
		ret_tpath = new long int;
	}
	*ret_tpath = tpath;
	
	if(ret_qstop == NULL) {
		ret_qstop = new int;
	}
	*ret_qstop = qstop;
	
	//This shall be considered the 'main' return value. NOTE that this will be NULL if 'savepaths' = 0
	return all_paths;
		
} //THE END

