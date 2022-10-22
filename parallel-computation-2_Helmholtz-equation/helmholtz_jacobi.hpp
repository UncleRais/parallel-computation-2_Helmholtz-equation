#pragma once

#include <functional>
#include <omp.h>

#include "utils.hpp"

#define SEQUENTIAL_MODE 0
#define PARALLEL_MODE 1


// # Solves Helholtz equation #
// -u_xx - u_yy + k^2 u(x, y) = f(x, y)
// k - wave_number
// f - right_part
// using
// # Iterative Jacobi method for linear system solution #
// A * x = f, 
// (L + U) * x_k + D * x_k+1 = f, 
// L + U = A - D,
// A = L + D + U, L - upper triangular, D - diagonal, U - lower triangular
// x_k+1 = D^{-1} (f - (A - D) x_k)
// x_k+1 = y + C x_k,
// y = D^{-1} * f, C = D^{-1} * (D - A)
// Ñonvergence condition: diagonal dominance
UniquePtrArray helholtz_jacobi(T k, std::function<T(T, T)> f, T L, size_t N, T epsilon,
	std::function<T(T)> BC_left, std::function<T(T)> BC_right, 
	std::function<T(T)> BC_bot, std::function<T(T)> BC_top,
	bool parallel
) {
	T h = T(1) / (N - 1);

	const T alpha = T(4) + sqr(h) * sqr(k);
	const T inverseAlpha = T(1) / alpha;
	const T beta = sqr(h);

	const int internalN = N - 2;
	const int num_var = sqr(N);

	auto x = make_raw_array(num_var, 0); // initial guess in zero-vector
	auto x0 = make_raw_array(num_var);

	auto x_ptr = x.get();
	auto x0_ptr = x0.get();

	// Tabulate boundaries
	// Left
	for (int i = 0; i < N; ++i) {
		x_ptr[i * N + 0] = BC_left(i * h);
		x0_ptr[i * N + 0] = BC_left(i * h);
	}
	// Right
	for (int i = 0; i < N; ++i) {
		x_ptr[i * N + N - 1] = BC_right(i * h);
		x0_ptr[i * N + N - 1] = BC_right(i * h);
	}
	// Top
	for (int j = 1; j < N - 1; ++j) {
		x_ptr[0 * N + j] = BC_top(j * h);
		x0_ptr[0 * N + j] = BC_top( j * h);
	}
	// Bottom
	for (int j = 1; j < N - 1; ++j) {
		x_ptr[(N - 1) * N + j] = BC_bot(j * h);
		x0_ptr[(N - 1) * N + j] = BC_bot(j * h);
	}

	T norm_diff_2 = 0;

	// Jacobi method
	do {
		// x0 = x
		std::swap(x_ptr, x0_ptr);

		// Internal loop
		#pragma omp parallel for if(parallel)
		for (int i = 1; i < N - 1; ++i) {
			for (int j = 1; j < N - 1; ++j) {
				const int indexIJ = i * N + j;

				x_ptr[indexIJ] = inverseAlpha * (
					x0_ptr[indexIJ - N] + // up
					x0_ptr[indexIJ + N] + // down
					x0_ptr[indexIJ - 1] + // left
					x0_ptr[indexIJ + 1] + // right
					beta * f(j * h, i * h)
					);
			}
		}

		// Stop condition
		norm_diff_2 = 0;
		#pragma omp parallel for if(parallel) reduction(+:norm_diff_2)
		for (size_t i = 0; i < num_var; ++i) norm_diff_2 += sqr(x_ptr[i] - x0_ptr[i]);

	} while (sqrt(norm_diff_2) > epsilon);
	
	// Return vector from last iteration
	if (x.get() == x_ptr) return x;
	else return x0;
}