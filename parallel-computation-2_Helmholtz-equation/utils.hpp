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

void print_array(T* arr, size_t size) {
	std::cout << "{ ";
	for (size_t k = 0; k < size - 1; ++k) std::cout << arr[k] << ", ";
	std::cout << arr[size - 1] << " }\n";
}

constexpr T PI = 3.14159265358979323846;

template<typename T>
constexpr T sqr(T value) { return value * value; } // screw you C++, I want my sqr()


template<typename T>
constexpr T norm_of_difference_L2(T* a, T* b, size_t size) {
	T res(0);
	for (size_t i = 0; i < size; ++i) res += sqr(a[i] - b[i]);
	return std::sqrt(res);
};

// Formating utils
template< typename T >
void print_system(const T* bot,
	const T* lower,
	const T* diag,
	const T* upper,
	const T* top,
	const T* right_part, size_t size)
{
	const int step = sqrt(size);
	const int params[2] = { 8, 3 };
	std::cout << "---------------------------------------------------- \n";
	std::cout << "|";
	std::cout << std::setw(params[0]) << std::setprecision(params[1]) << diag[0] << " ";
	std::cout << std::setw(params[0]) << std::setprecision(params[1]) << upper[0] << " ";
	for (int i = 2; i < size; ++i)
		if (i == step)
		{
			std::cout << std::setw(params[0]) << std::setprecision(params[1]) << top[0] << " ";
		}
		else
		{
			std::cout << std::setw(params[0]) << std::setprecision(params[1]) << T(0) << " ";
		}
	std::cout << "| " << std::setw(params[0]) << std::setprecision(params[1]) << right_part[0];
	std::cout << "|\n";
	for (int i = 1; i < size - 1; ++i)
	{
		std::cout << "|";
		for (int j = 0; j < i - 1; ++j)
			if (j == i - step)
			{
				std::cout << std::setw(params[0]) << std::setprecision(params[1]) << bot[j] << " ";
			}
			else
			{
				std::cout << std::setw(params[0]) << std::setprecision(params[1]) << T(0) << " ";
			}
		std::cout << std::setw(params[0]) << std::setprecision(params[1]) << lower[i - 1] << " ";
		std::cout << std::setw(params[0]) << std::setprecision(params[1]) << diag[i] << " ";
		std::cout << std::setw(params[0]) << std::setprecision(params[1]) << upper[i] << " ";
		for (int j = i + 2; j < size; ++j)
			if (j == step + i)
			{
				std::cout << std::setw(params[0]) << std::setprecision(params[1]) << top[i] << " ";
			}
			else
			{
				std::cout << std::setw(params[0]) << std::setprecision(params[1]) << T(0) << " ";
			}
		std::cout << "| " << std::setw(params[0]) << std::setprecision(params[1]) << right_part[i];
		std::cout << "|\n";
	}
	std::cout << "|";
	for (int i = 0; i < size - 2; ++i)
		if (i == size - 1 - step)
		{
			std::cout << std::setw(params[0]) << std::setprecision(params[1]) << bot[i] << " ";
		}
		else
		{
			std::cout << std::setw(params[0]) << std::setprecision(params[1]) << T(0) << " ";
		}
	std::cout << std::setw(params[0]) << std::setprecision(params[1]) << lower[size - 2] << " ";
	std::cout << std::setw(params[0]) << std::setprecision(params[1]) << diag[size - 1] << " ";
	std::cout << "| " << std::setw(params[0]) << std::setprecision(params[1]) << right_part[size - 1];
	std::cout << "|\n";
	std::cout << "---------------------------------------------------- \n";
}