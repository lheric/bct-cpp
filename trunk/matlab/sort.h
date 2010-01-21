#ifndef SORT_H
#define SORT_H

#include <algorithm>
#include <vector>

namespace matlab {
	template<class T> void stable_sort(T*, size_t);
	template<class T> void stable_sort(T*, size_t, bool (*)(T, T));
	template<class T> void stable_sort_index(size_t*, const T*, size_t);
	template<class T> void stable_sort_index(size_t*, const T*, size_t, bool (*)(T, T));
}

template <class T1, class T2> struct transparent_pair {
	T1 first;
	T2 second;
	transparent_pair(const T1& x, const T2& y) : first(x), second(y) { }
	operator T1() const { return first; }
};

/*
 * Performs an in-place stable sort.
 */
template<class T> void matlab::stable_sort(T* array, size_t count) {
	std::vector<T> v;
	v.assign(array, array + count);
	std::stable_sort(v.begin(), v.end());
	std::copy(v.begin(), v.end(), array);
}

/*
 * Performs an in-place stable sort using a custom comparison function.
 */
template<class T> void matlab::stable_sort(T* array, size_t count, bool (*compare)(T, T)) {
	std::vector<T> v;
	v.assign(array, array + count);
	std::stable_sort(v.begin(), v.end(), compare);
	std::copy(v.begin(), v.end(), array);
}

/*
 * Performs an indirect stable sort.
 */
template<class T> void matlab::stable_sort_index(size_t* indices, const T* array, size_t count) {
	std::vector<transparent_pair<T, size_t> > v;
	for (size_t i = 0; i < count; i++) {
		v.push_back(transparent_pair<T, size_t>(array[i], i));
	}
	std::stable_sort(v.begin(), v.end());
	for (size_t i = 0; i < count; i++) {
		indices[i] = v[i].second;
	}
}

/*
 * Performs an indirect stable sort using a custom comparison function.
 */
template<class T> void matlab::stable_sort_index(size_t* indices, const T* array, size_t count, bool (*compare)(T, T)) {
	std::vector<transparent_pair<T, size_t> > v;
	for (size_t i = 0; i < count; i++) {
		v.push_back(transparent_pair<T, size_t>(array[i], i));
	}
	std::stable_sort(v.begin(), v.end(), compare);
	for (size_t i = 0; i < count; i++) {
		indices[i] = v[i].second;
	}
}

#endif
