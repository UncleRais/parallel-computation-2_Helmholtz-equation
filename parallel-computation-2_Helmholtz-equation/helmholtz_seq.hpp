#pragma once

#include <functional>
#include <memory>

#include "utils.hpp"


// SLAE format
// [0] top_diag   -> -1
// [1] upper_diag -> ((i+1) % internalN) ? -1 : 0 
// [2] diag       -> alpha
// [3] lower_diag -> ((i+1) % internalN) ? -1 : 0 
// [4] bot_diag   -> -1
// [5] right_part -> beta * f_ij
T get_SLAE_alpha(T k, T L, size_t N) {
	const T h = L / (N - 1);

	return T(4) + sqr(h) * sqr(k);
}
RawArray get_SLAE_right_part(std::function<T(T, T)> f, T L, size_t N) {
	const T h = L / (N - 1);

	const size_t internalN = N - 2;
	const size_t num_var = sqr(internalN);

	auto res = make_raw_array(num_var);

	const T beta = sqr(h);

	for (size_t j = 1; j <= internalN; ++j)
		for (size_t i = 1; i <= internalN; ++i)
			res[(j - 1) * internalN + (i - 1)] = beta * f(i * h, (static_cast<T>(internalN) - j + 1) * h);

	return res;
}


// Iterative Jacobi method for linear system solution
// A * x = f, 
// (L + U) * x_k + D * x_k+1 = f, 
// L + U = A - D,
// A = L + D + U, L - upper triangular, D - diagonal, U - lower triangular
// x_k+1 = D^{-1} (f - (A - D) x_k)
// x_k+1 = y + C x_k,
// y = D^{-1} * f, C = D^{-1} * (D - A)
// Ñonvergence condition: diagonal dominance
RawArray Jacobi(T alpha, T* const f, size_t num_var, size_t internalN, T epsilon) {
	T normC = 0.9999;

	const T inverseAlpha = T(1) / alpha;

	auto x = make_raw_array(num_var, 0); // initial guess in zero-vector
	auto x0 = make_raw_array(num_var);

	///std::unique_ptr<bool[]> ibool(new bool[num_var * sizeof(bool)]);

	do {
		// x0 = x
		memcpy(x0.get(), x.get(), num_var * sizeof(T));

		// Rows i = 0
		x[0] = inverseAlpha * (
			f[0]
			+ x0[1]
			+ x0[internalN]
			);

		// Rows i = 1, internalN - 1
		for (size_t i = 1; i < internalN; ++i)
			x[i] = inverseAlpha * (
				f[i]
				+ x0[i - 1]
				+ x0[i + 1]
				+ x0[i + internalN]
				);

		// Rows i = internalN, N - interanlN
		for (size_t i = internalN; i < num_var - internalN; ++i)
			x[i] = inverseAlpha * (
				f[i]
				+ x0[i - 1]
				+ x0[i + 1]
				+ x0[i + internalN]
				+ x0[i - internalN]
				);

		// Rows i = N - interanlN + 1, N - 1
		for (size_t i = num_var - internalN; i < num_var - 1; ++i)
			x[i] = inverseAlpha * (
				f[i]
				+ x0[i - 1]
				+ x0[i + 1]
				+ x0[i - internalN]
				);

		// Rows i = N
		x[num_var - 1] = inverseAlpha * (
			f[num_var - 1]
			+ x0[num_var - 2]
			+ x0[num_var - 1 - internalN]
			);

		// Account for 0's on upper/lower diagonals
		for (size_t i = internalN; i <= num_var - internalN; i += internalN) {
			x[i - 1] -= inverseAlpha * x0[i];
			x[i] -= inverseAlpha * x0[i - 1];
		}

		// Uncomment to debug
		/*std::cout << "x0 = \n";
		print_array(x0.get(), num_var);

		std::cout << "x = \n";
		print_array(x.get(), num_var);

		std::cout << "err = " << norm_of_difference_L2(x.get(), x0.get(), num_var) << "\n";*/

	} while (norm_of_difference_L2(x.get(), x0.get(), num_var) > (1 - normC) / normC * epsilon);
	
	return x;
}


// Solves equation of given form:
// -u_xx - u_yy + k^2 u(x, y) = f(x, y)
// k - wave_number
// f - right_part
RawArray Helmholtz_solve(T k, std::function<T(T, T)> f, T L, size_t N, T epsilon = 1e-5) {
	auto alpha = get_SLAE_alpha(k, L, N);
	auto right_part = get_SLAE_right_part(f, L, N);
	const size_t num_var = sqr(N - 2);
	const size_t internalN = N - 2;

	RawArray sol = Jacobi(alpha, right_part.get(), num_var, internalN, epsilon);

	return sol;
}