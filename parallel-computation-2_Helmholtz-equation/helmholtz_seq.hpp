#pragma once

#include <functional>

#include "jacobi_seq.hpp"
#include "jacobi_par.hpp"


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

UniquePtrArray get_SLAE_right_part(std::function<T(T, T)> f, T L, size_t N) {
	const T h = L / (N - 1);

	const size_t internalN = N - 2;
	const size_t num_var = sqr(internalN);

	auto res = make_raw_array(num_var);

	const T beta = sqr(h);

	for (size_t j = 1; j <= internalN; ++j)
		for (size_t i = 1; i <= internalN; ++i)
			res[(j - 1) * internalN + (i - 1)] = beta * f(i * h, h * (internalN - j + 1));

	return res;
}


enum class Method {
	SEQUENTIAL_JACOBI,
	PARALLEL_JACOBI,
	SEQUENTIAL_SEIDEL,
	PARALLEL_SEIDEL,
};

// Solves equation of given form:
// -u_xx - u_yy + k^2 u(x, y) = f(x, y)
// k - wave_number
// f - right_part
UniquePtrArray Helmholtz_solve(T k, std::function<T(T, T)> f, T L, size_t N, Method method, T epsilon = 1e-5) {
	auto alpha = get_SLAE_alpha(k, L, N);
	auto right_part = get_SLAE_right_part(f, L, N);
	const size_t num_var = sqr(N - 2);
	const size_t internalN = N - 2;

	switch (method) {
	case Method::SEQUENTIAL_JACOBI:
		return jacobi_seq(alpha, right_part.get(), num_var, internalN, epsilon);
	case Method::PARALLEL_JACOBI:
		return jacobi_par(alpha, right_part.get(), num_var, internalN, epsilon);
	case Method::SEQUENTIAL_SEIDEL:
		return jacobi_seq(alpha, right_part.get(), num_var, internalN, epsilon);
	case Method::PARALLEL_SEIDEL:
		return jacobi_seq(alpha, right_part.get(), num_var, internalN, epsilon);
	default:
		exit_with_error("Unknown method");
	}
}