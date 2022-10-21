#pragma once

#include <omp.h>

#include "utils.hpp"


inline T norm_of_difference_L2_par(T* const a, T* const b, int size) {
	T res(0);

	#pragma omp parallel for schedule(guided) firstprivate(a, b) reduction(+:res)
	for (int i = 0; i < size; ++i) res += sqr(a[i] - b[i]);

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
UniquePtrArray jacobi_par(T alpha, T* const f, int num_var, int internalN, T epsilon) {
	const T inverseAlpha = T(1) / alpha;

	auto x = make_raw_array(num_var, 0); // initial guess in zero-vector
	auto x0 = make_raw_array(num_var);

	auto x_ptr = x.get();
	auto x0_ptr = x0.get();

	#pragma omp parallel
	while (true) {
		// x0 = x
		memcpy(x0.get(), x.get(), num_var * sizeof(T));

		// Rows i = 0
		x[0] = inverseAlpha * (
			f[0]
			+ x0[1]
			+ x0[internalN]
			);

		//firstprivate(num_var, internalN, f, inverseAlpha, x_ptr, x0_ptr)

		//#pragma omp parallel default(none) firstprivate(num_var, internalN, f, inverseAlpha, x_ptr, x0_ptr)
		//{
		//	// Rows i = 1, internalN - 1
		//	#pragma omp for nowait schedule(static)
		//	for (int i = 1; i < internalN; ++i)
		//		x_ptr[i] = inverseAlpha * (
		//			f[i]
		//			+ x0_ptr[i - 1]
		//			+ x0_ptr[i + 1]
		//			+ x0_ptr[i + internalN]
		//			);

		//	// Rows i = internalN, N - interanlN
		//	#pragma omp for nowait schedule(static)
		//	for (int i = internalN; i < num_var - internalN; ++i)
		//		x_ptr[i] = inverseAlpha * (
		//			f[i]
		//			+ x0_ptr[i - 1]
		//			+ x0_ptr[i + 1]
		//			+ x0_ptr[i + internalN]
		//			+ x0_ptr[i - internalN]
		//			);

		//	// Rows i = N - interanlN + 1, N - 1
		//	#pragma omp for nowait schedule(static)
		//	for (int i = num_var - internalN; i < num_var - 1; ++i)
		//		x_ptr[i] = inverseAlpha * (
		//			f[i]
		//			+ x0_ptr[i - 1]
		//			+ x0_ptr[i + 1]
		//			+ x0_ptr[i - internalN]
		//			);
		//}

		// Rows i = 1, internalN - 1
		#pragma omp parallel for schedule(static) firstprivate(num_var, internalN, f, inverseAlpha, x_ptr, x0_ptr)
		for (int i = 1; i < internalN; ++i)
			x_ptr[i] = inverseAlpha * (
				f[i]
				+ x0_ptr[i - 1]
				+ x0_ptr[i + 1]
				+ x0_ptr[i + internalN]
				);

		// Rows i = internalN, N - interanlN
		#pragma omp parallel for schedule(static) firstprivate(num_var, internalN, f, inverseAlpha, x_ptr, x0_ptr)
		for (int i = internalN; i < num_var - internalN; ++i)
			x_ptr[i] = inverseAlpha * (
				f[i]
				+ x0_ptr[i - 1]
				+ x0_ptr[i + 1]
				+ x0_ptr[i + internalN]
				+ x0_ptr[i - internalN]
				);

		// Rows i = N - interanlN + 1, N - 1
		#pragma omp parallel for schedule(static) firstprivate(num_var, internalN, f, inverseAlpha, x_ptr, x0_ptr)
		for (int i = num_var - internalN; i < num_var - 1; ++i)
			x_ptr[i] = inverseAlpha * (
				f[i]
				+ x0_ptr[i - 1]
				+ x0_ptr[i + 1]
				+ x0_ptr[i - internalN]
				);

		// Rows i = N
		x[num_var - 1] = inverseAlpha * (
			f[num_var - 1]
			+ x0[num_var - 2]
			+ x0[num_var - 1 - internalN]
			);

		// Account for 0's on upper/lower diagonals
		for (int i = internalN; i <= num_var - internalN; i += internalN) {
			x[i - 1] -= inverseAlpha * x0[i];
			x[i] -= inverseAlpha * x0[i - 1];
		}

	}/// while (norm_of_difference_L2_par(x.get(), x0.get(), num_var) > epsilon);

	return x;
}