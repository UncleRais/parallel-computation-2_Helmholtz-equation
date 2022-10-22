#pragma once

#include "utils.hpp"


inline T norm_of_difference_L2(T* const a, T* const b, size_t size) {
	T res(0);
	for (size_t i = 0; i < size; ++i) res += sqr(a[i] - b[i]);
	return std::sqrt(res);
};


// Iterative Jacobi method for linear system solution
// A * x = f, 
// (L + U) * x_k + D * x_k+1 = f, 
// L + U = A - D,
// A = L + D + U, L - upper triangular, D - diagonal, U - lower triangular
// x_k+1 = D^{-1} (f - (A - D) x_k)
// x_k+1 = y + C x_k,
// y = D^{-1} * f, C = D^{-1} * (D - A)
// Ñonvergence condition: diagonal dominance
UniquePtrArray jacobi_seq(T alpha, T* const f, size_t num_var, size_t internalN, T epsilon) {

	const T inverseAlpha = T(1) / alpha;

	auto x = make_raw_array(num_var, 0); // initial guess in zero-vector
	auto x0 = make_raw_array(num_var);

	auto x_ptr = x.get();
	auto x0_ptr = x0.get();

	do {
		// x0 = x
		std::swap(x_ptr, x0_ptr);

		// Rows i = 0
		x_ptr[0] = inverseAlpha * (
			f[0]
			+ x0_ptr[1]
			+ x0_ptr[internalN]
			);

		// Rows i = 1, internalN - 1
		for (size_t i = 1; i < internalN; ++i)
			x_ptr[i] = inverseAlpha * (
				f[i]
				+ x0_ptr[i - 1]
				+ x0_ptr[i + 1]
				+ x0_ptr[i + internalN]
				);

		// Rows i = internalN, N - interanlN
		for (size_t i = internalN; i < num_var - internalN; ++i)
			x_ptr[i] = inverseAlpha * (
				f[i]
				+ x0_ptr[i - 1]
				+ x0_ptr[i + 1]
				+ x0_ptr[i + internalN]
				+ x0_ptr[i - internalN]
				);

		// Rows i = N - interanlN + 1, N - 1
		for (size_t i = num_var - internalN; i < num_var - 1; ++i)
			x_ptr[i] = inverseAlpha * (
				f[i]
				+ x0_ptr[i - 1]
				+ x0_ptr[i + 1]
				+ x0_ptr[i - internalN]
				);

		// Rows i = N
		x_ptr[num_var - 1] = inverseAlpha * (
			f[num_var - 1]
			+ x0_ptr[num_var - 2]
			+ x0_ptr[num_var - 1 - internalN]
			);

		// Account for 0's on upper/lower diagonals
		for (size_t i = internalN; i <= num_var - internalN; i += internalN) {
			x_ptr[i - 1] -= inverseAlpha * x0_ptr[i];
			x_ptr[i] -= inverseAlpha * x0_ptr[i - 1];
		}

		// Uncomment to debug
		/*std::cout << "x0 = \n";
		print_array(x0.get(), num_var);

		std::cout << "x = \n";
		print_array(x.get(), num_var);

		std::cout << "err = " << norm_of_difference_L2(x.get(), x0.get(), num_var) << "\n";*/

	} while (norm_of_difference_L2(x_ptr, x0_ptr, num_var) > epsilon);

	if (x.get() == x_ptr) return x;
	else return x0;
}