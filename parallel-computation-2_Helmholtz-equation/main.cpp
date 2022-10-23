#include <iostream>
#include <iomanip>

#include "static_timer.hpp"
#include "helmholtz_jacobi.hpp"
#include "helmholtz_seidel.hpp"
#include "table.hpp"

#define SEQUENTIAL_MODE 1


/// Uncomment to check for memory leaks
//#define _DEBUG
//#define _CRTDBG_MAP_ALLOC
//#include <stdlib.h>
//#include <crtdbg.h>


// # Config #
// Area size
const T L = 1;

// Grid size
const size_t N = 2002;
const size_t internalN = N - 2;

const size_t pointCount = sqr(N);
const size_t internalPointCount = sqr(internalN);

// Wave number
const T c1 = 10;
const T h = L / (N - 1);

const T k = sqrt(c1) / h;

// Precision
const T precision = 1e-6;

// Right part
T right_part(T x, T y) { return(2 * sin(PI * y) + k * k * (1 - x) * x * sin(PI * y) + PI * PI * (1 - x) * x * sin(PI * y)); }

// Boundaries
T zero_boundary(T value) { return 0; }

// Precise solution (for analythical purposes)
T analythical_solution(T x, T y) { return (1 - x) * x * sin(PI * y); }

UniquePtrArray get_exact_solution() {
	auto sol_exact = make_raw_array(pointCount);

	T h = T(1) / (N - 1);

	for (size_t i = 0; i < N; ++i)
		for (size_t j = 0; j < N; ++j)
			sol_exact[i * N + j] = analythical_solution(j * h, i * h);

	return sol_exact;
}

T get_relative_error_L2(T* const sol_numeric, T* const sol_exact) {
	T L2_norm_of_difference(0);
	T L2_norm_of_exact_sol(0);

	for (size_t i = 0; i < pointCount; ++i) {
		L2_norm_of_difference += sqr(sol_numeric[i] - sol_exact[i]);
		L2_norm_of_exact_sol += sqr(sol_exact[i]);
	}

	return sqrt(L2_norm_of_difference / L2_norm_of_exact_sol);
}

void print_sol(T* const sol) {
	for (size_t i = 0; i < N; ++i) {
		std::cout << "[";
		for (size_t j = 0; j < N; ++j)
			std::cout << std::setw(15) << sol[i * N + j];
		std::cout << "]\n";
	}
}


int main(int argc, char** argv) {
	// Set up OpenMP
	/*const int THREAD_CAP = 4;

	const int MAX_THREADS = omp_get_max_threads();
	const int NUM_THREADS = std::min(THREAD_CAP, MAX_THREADS);

	omp_set_num_threads(NUM_THREADS);*/

	std::cout
		<< "N = " << N << "\n"
		<< "k^2 h^2 = " << c1 << "\n"
		<< "precision = " << precision << "\n\n";
	//	<< "Using " << NUM_THREADS << " threads out of " << MAX_THREADS << " available\n\n";

	table_add_1("Method");
	table_add_2("Time (sec)");
	table_add_3("Rel. error");
	table_add_4("Speedup");
	table_add_5("Threads");
	table_hline();

	double jacobiSeqTime = 0;
	double jacobiParTime = 0;
	double SeidelSeqTime = 0;
	double SeidelParTime = 0;

	double MAXTHREADS = 8;

	auto solution_exact = get_exact_solution();

	// 1) Sequential Jacobi
	{
		table_add_1("Jacobi");

		StaticTimer::start();
		auto solution = helholtz_jacobi(
			k, right_part, L, N, precision,
			zero_boundary, zero_boundary, zero_boundary, zero_boundary,
			SEQUENTIAL_MODE
		);
		jacobiSeqTime = StaticTimer::end();

		table_add_2(jacobiSeqTime);
		table_add_3(get_relative_error_L2(solution.get(), solution_exact.get()));
		table_add_4(1);
		table_add_5(1);
	}

	// 2) Parallel Jacobi
	for (int PARALLEL_MODE = 2; PARALLEL_MODE <= MAXTHREADS; PARALLEL_MODE++)
	{
		table_add_1("Parallel Jacobi");

		StaticTimer::start();
		auto solution = helholtz_jacobi(
			k, right_part, L, N, precision,
			zero_boundary, zero_boundary, zero_boundary, zero_boundary,
			PARALLEL_MODE
		);
		jacobiParTime = StaticTimer::end();

		table_add_2(jacobiParTime);
		table_add_3(get_relative_error_L2(solution.get(), solution_exact.get()));
		table_add_4(jacobiSeqTime / jacobiParTime);
		table_add_5(PARALLEL_MODE);
	}

	// 3)Sequential Seidel
	{
		table_add_1("Seidel");

		StaticTimer::start();
		auto solution = helholtz_seidel(
			k, right_part, L, N, precision,
			zero_boundary, zero_boundary, zero_boundary, zero_boundary,
			SEQUENTIAL_MODE
		);
		SeidelSeqTime = StaticTimer::end();

		table_add_2(SeidelSeqTime);
		table_add_3(get_relative_error_L2(solution.get(), solution_exact.get()));
		table_add_4(SeidelSeqTime / SeidelSeqTime);
		table_add_5(SEQUENTIAL_MODE);
	}

	// 2) Parallel Seidel
	for (int PARALLEL_MODE = 2; PARALLEL_MODE <= MAXTHREADS; PARALLEL_MODE++)
	{
		table_add_1("Parallel Seidel");

		StaticTimer::start();
		auto solution = helholtz_seidel(
			k, right_part, L, N, precision,
			zero_boundary, zero_boundary, zero_boundary, zero_boundary,
			PARALLEL_MODE
		);
		SeidelParTime = StaticTimer::end();

		table_add_2(SeidelParTime);
		table_add_3(get_relative_error_L2(solution.get(), solution_exact.get()));
		table_add_4(SeidelSeqTime / SeidelParTime);
		table_add_5(PARALLEL_MODE);
	}

	/// Uncomment to check for memory leaks
	// solution_exact.reset();
	/*if (_CrtDumpMemoryLeaks()) { std::cout << "\nMemory leaks!\n"; }
	else { std::cout << "\n No leaks\n"; }*/
}

