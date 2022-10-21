#pragma once

#include <cmath>

// Math utils
using T = double;
using RawArray = std::unique_ptr<T[]>;

RawArray make_raw_array(size_t size) {
	return RawArray(new T[size * sizeof(T)]);
}

RawArray make_raw_array(size_t size, T defaultValue) {
	RawArray arr(new T[size * sizeof(T)]);

	for (size_t k = 0; k < size; ++k) arr[k] = defaultValue;

	return arr;
}

void print_array(T* const arr, size_t size) {
	std::cout << "{ ";
	for (size_t k = 0; k < size - 1; ++k) std::cout << arr[k] << ", ";
	std::cout << arr[size - 1] << " }\n";
}


constexpr T PI = 3.14159265358979323846;


template<typename Type>
constexpr Type sqr(Type value) { return value * value; } // screw you C++, I want my sqr()


constexpr T norm_of_difference_L2(T* const a, T* const b, size_t size) {
	T res(0);
	for (size_t i = 0; i < size; ++i) res += sqr(a[i] - b[i]);
	return std::sqrt(res);
};
