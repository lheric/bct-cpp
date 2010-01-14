#include <cstdlib>
#include <ctime>

namespace matlab {
	typedef int (*cmp_fn)(void*, void*);
	
	/*
	 * Stable quicksort.
	 */
	template<typename T> void quicksort(T* array, size_t count, cmp_fn compare) {
		quicksort_internal(NULL, array, count, compare);
	}

	/*
	 * Stable indirect quicksort.
	 */
	template<typename T> void quicksort_index(size_t* indices, const T* array, size_t count, cmp_fn compare) {
		T array_copy[count];
		for (int i = 0; i < count; i++) {
			array_copy[i] = array[i];
			indices[i] = i;
		}
		quicksort_internal(indices, array_copy, count, compare);
	}

	template<typename T> void quicksort_internal(size_t* indices, T* array, size_t count, cmp_fn compare) {
		
		// Empty or one-element array is already sorted
		if (count <= 1) {
			return;
		}
		
		// Initialize less/greater arrays
		T* less = new T[count - 1];
		T* greater = new T[count - 1];
		size_t* less_indices;
		size_t* greater_indices;
		if (indices != NULL) {
			less_indices = new size_t[count - 1];
			greater_indices = new size_t[count - 1];
		}
		int n_less = 0;
		int n_greater = 0;
		
		// Choose pivot randomly
		// Make sure it is the last among equal elements
		std::srand(std::time(NULL));
		int i_pivot = std::rand() % count;
		T pivot = array[i_pivot];
		for (int i = i_pivot + 1; i < count; i++) {
			T value = array[i];
			if (compare(&value, &pivot) == 0) {
				i_pivot = i;
				pivot = value;
			}
		}
		int pivot_index;
		if (indices != NULL) {
			pivot_index = indices[i_pivot];
		}
		
		// Separate into values less/greater than pivot
		// Elements equal to pivot are considered less
		for (int i = 0; i < count; i++) {
			if (i == i_pivot) {
				continue;
			}
			T value = array[i]; 
			if (compare(&value, &pivot) <= 0) {
				less[n_less] = value;
				if (indices != NULL) {
					less_indices[n_less] = indices[i];
				}
				n_less++;
			} else {
				greater[n_greater] = value;
				if (indices != NULL) {
					greater_indices[n_greater] = indices[i];
				}
				n_greater++;
			}
		}
		
		// Store values back into original array
		for (int i = 0; i < n_less; i++) {
			array[i] = less[i];
			if (indices != NULL) {
				indices[i] = less_indices[i];
			}
		}
		array[n_less] = pivot;
		if (indices != NULL) {
			indices[n_less] = pivot_index;
		}
		for (int i = 0; i < n_greater; i++) {
			array[n_less + i + 1] = greater[i];
			if (indices != NULL) {
				indices[n_less + i + 1] = greater_indices[i];
			}
		}
		delete[] less;
		delete[] greater;
		if (indices != NULL) {
			delete[] less_indices;
			delete[] greater_indices;
		}
		
		// Sort less/greater portions separately
		if (indices == NULL) {
			quicksort_internal(NULL, array, n_less, compare);
			quicksort_internal(NULL, &array[n_less + 1], n_greater, compare);
		} else {
			quicksort_internal(indices, array, n_less, compare);
			quicksort_internal(&indices[n_less + 1], &array[n_less + 1], n_greater, compare);
		}
	}
}
