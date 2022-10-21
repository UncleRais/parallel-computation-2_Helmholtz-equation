#include <iostream>
#include <iomanip>

#include "static_timer.hpp"
#include "helmholtz_seq.hpp"


// # Config #
// Wave number
const T k = 10;

// Area size
const T L = 1;

// Grid size
const size_t N = 102;
const size_t internalN = N - 2;
const size_t internalPointCount = sqr(internalN);

// Precision
const T precision = 1e-9;

// Right part
T right_part(T x, T y) { return(2 * sin(PI * y) + k * k * (1 - x) * x * sin(PI * y) + PI * PI * (1 - x) * x * sin(PI * y)); }

// Precise solution (for analythical purposes)
T analythical_solution(T x, T y) { return (1 - x) * x * sin(PI * y); }

RawArray get_exact_solution() {
	RawArray sol_exact(new T[internalPointCount * sizeof(T)]);

	T h = T(1) / (N - 1);

	for (size_t j = 1; j < internalN + 1; ++j)
		for (size_t i = 1; i < internalN + 1; ++i)
			sol_exact[(j - 1) * internalN + (i - 1)] = analythical_solution(
				i * h,
				(static_cast<double>(internalN) - j + 1) * h
			);

	return sol_exact;
}


int main(int argc, char** argv) {
	// Numeric, Jacobi
	StaticTimer::start();
	RawArray sol_numeric = Helmholtz_solve(k, right_part, T(1), N, precision);
	std::cout << "Jacobi done in " << StaticTimer::end() << " sec\n";

	// Exact
	RawArray sol_exact = get_exact_solution();

	T L2_norm_of_difference(0);
	T L2_norm_of_exact_sol (0);
	for (size_t i = 0; i < internalPointCount; ++i) {
		L2_norm_of_difference += sqr(sol_numeric[i] - sol_exact[i]);
		L2_norm_of_exact_sol += sqr(sol_exact[i]);
	}


	std::cout << "Error L2_Norm = " << sqrt(L2_norm_of_difference / L2_norm_of_exact_sol) << "\n";
}

