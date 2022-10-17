#define _USE_MATH_DEFINES
#define M_PI
#include <iostream>
#include <functional>
#include <cmath>
#include <iomanip>

template< typename T >
constexpr T sqr(T value) { return value * value; } // screw you C++, I want my sqr()

template< typename T >
constexpr T norm_minus(T* a, T* b, size_t size) 
{ 
	T res = T(0);
	for (size_t i = 0; i < size; ++i) res += sqr(a[i] - b[i]);
	return(sqrt(res));
};


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
		std::cout << std::setw(params[0]) << std::setprecision(params[1]) << lower[i-1] << " ";
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

// Assemble system:
// A * x = b
// Returns five pointers to diagonals of A and a pointer to b 
template<typename T>
std::vector< T* > assemble_system(T k, std::function<T(T, T)> f, T L, size_t N)
{
	const T h = L / (N - 1);

	const size_t internalN = N - 2;
	const size_t num_var = sqr(internalN);


	// Init SLAE matrix/rigthpart
	T* top_diag = new T[(num_var - internalN) * sizeof(T)];
	T* upper_diag = new T[(num_var - 1) * sizeof(T)];
	T* diag = new T[(num_var) * sizeof(T)];
	T* lower_diag = new T[(num_var - 1) * sizeof(T)];
	T* bot_diag = new T[(num_var - internalN) * sizeof(T)];

	T* right_part = new T[num_var];

	// Fill SLAE matrix/rigthpart
	const T diag_coef = 4 + sqr(h) * sqr(k);
	const T right_coef = sqr(h);

	for (int i = 0; i < num_var - internalN; ++i) top_diag[i] = -1;
	for (int i = 0; i < num_var - 1; ++i)		  upper_diag[i] = -1;
	for (int i = 0; i < num_var; ++i)			  diag[i] = diag_coef;
	for (int i = 0; i < num_var - 1; ++i)		  lower_diag[i] = -1;
	for (int i = 0; i < num_var - internalN; ++i) bot_diag[i] = -1;

	for (size_t j = 1; j <= internalN ; ++j)
		for (size_t i = 1; i <= internalN ; ++i)
			right_part[(j - 1) * internalN + (i - 1)] = right_coef * f(i * h, (internalN - j + 1) * h);

	// Add 0's from boundary conditions
	for (int i = internalN - 1; i < num_var - 1; i += internalN) upper_diag[i] = 0;
	for (int i = internalN - 1; i < num_var - 1; i += internalN) lower_diag[i] = 0;

	std::vector<T*> res = { bot_diag, lower_diag, diag, upper_diag, top_diag, right_part };

	return(res);

}


// Iterative Jacobi method for linear system solution
// A * x = f, 
// (L + U) * x_k + D * x_k+1 = f, 
// L + U = A - D,
// A = L + D + U, L - upper triangular, D - diagonal, U - lower triangular
// x_k+1 = D^{-1} (f - (A - D) x_k)
// x_k+1 = y + C x_k,
// y = D^{-1} * f, C = D^{-1} * (D - A)
// Сonvergence condition: diagonal dominance
template < typename T >
T* Jacobi(std::vector<T*>& system, size_t num_var,  T epsilon) 
{
	const size_t len = num_var - size_t(sqrt(num_var));
	const size_t step = size_t(sqrt(num_var)) - 2;
	T* res = new T[num_var * sizeof(T)];

	T temp = 1 / system[2][0];

	res[0] = system[5][0] * temp;
	system[3][0] = -system[3][0] * temp;
	system[4][0] = -system[4][0] * temp;
	system[5][0] = system[5][0] * temp;

	for (size_t i = 1; i < num_var - 1; ++i)
	{
		res[i] = system[5][i] * temp; // initial vector
		if(i - 2 - step >= 0 && i - 2 - step < num_var) system[0][i  - 2 - step] = -system[0][i - 2 - step] * temp; //bot

		system[1][i - 1] = -system[1][i - 1] * temp; //lower
		system[3][i] = -system[3][i] * temp; // upper

		if(i <= len - 1) system[4][i] = -system[4][i] * temp; // top
		// y = D^{-1} * f
		system[5][i] = system[5][i] * temp; //right_part
	}

	res[num_var - 1] = system[5][num_var - 1] * temp;
	system[0][len - 1] = -system[0][len - 1] * temp;
	system[1][num_var - 2] = -system[1][num_var - 2] * temp;
	system[5][num_var - 1] = system[5][num_var - 1] * temp;

	//Print your system if you want 
	//print_system(system[0], system[1], system[2], system[3], system[4], system[5], num_var);

	T normC = 4 * system[0][0]; /// Dima, plz, calculate this value, i am fed up of it

	size_t counter = 0;
	T* prev = new T[num_var * sizeof(T)];
	for (size_t i = 0; i < num_var; ++i) prev[i] = T(0); 

	std::cout << "Iter estimate: " << log((1 - normC) * epsilon / norm_minus(res, prev, num_var)) / log(normC) << "\n";
	if (normC <= 1)
	{
		while (norm_minus(res, prev, num_var) >= (1 - normC) * (epsilon) / normC)
		{
			for (size_t i = 0; i < num_var; ++i) prev[i] = res[i];

			res[0] = system[5][0] + system[3][0] * prev[1] + system[4][0] * prev[step + 2] ;
			for (size_t i = 1; i < num_var - 1; ++i)
			{
				res[i] = system[5][i] + system[1][i - 1] * prev[i - 1] + system[3][i] * prev[i + 1]; // y + (lower + upper) * x_k
				if (i - 2 - step >= 0 && i - 2 - step < num_var) res[i] += system[0][i - 2 - step] * prev[i - 2 - step]; //bot
				if (i <= len - 1) res[i] += system[4][i] * prev[i + step + 2]; // top
	
			};
			res[num_var - 1] = system[5][num_var - 1] + system[0][len - 1] * prev[num_var - 3 - step] + system[1][num_var - 2] * prev[num_var - 2];
			++counter;
		}
	}
	std::cout << "Iterations: " << counter << "\n";

	return res;
}

// Solves equation of given form:
// -u_xx - u_yy + k^2 u(x, y) = f(x, y)
// k - wave_number
// f - right_part
template<typename T>
T* Helmholtz_solve(T wave_number, std::function<T(T, T)> right_part, T length, size_t num_points, T epsilon = 1e-5)
{
	std::vector< T* > system(6);
	system = assemble_system(wave_number, right_part, length, num_points);

	const size_t num_var = sqr(num_points - 2);
	//Print your system if you want 
	print_system(system[0], system[1], system[2], system[3], system[4], system[5], num_var);

	T* sol = Jacobi(system, num_var, epsilon);

	for (size_t i = 0; i < 6; ++i)
		delete system[i];

	return(sol);
}


typedef double T;

static T k = 10;

T right_part(T x, T y)
{
	return(2 * sin(M_PI * y) + k * k * (1 - x) * x * sin(M_PI * y) + M_PI * M_PI * (1 - x) * x * sin(M_PI * y));
}

T solution(T x, T y)
{
	return((1 - x) * x * sin(M_PI * y));
}

int main(int argc, char** argv)
{
	std::function<T (T, T)> f = right_part;
	size_t N = 5;
	size_t num_var = N - 2;
	size_t internalN = sqr(num_var);
	T* sol = Helmholtz_solve(k, f, T(1), N, 1e-5);

	T* exact = new T[internalN * sizeof(T)];
	T h = T(1) / (N - 1);
	for (size_t j = 1; j < num_var + 1; ++j)
		for (size_t i = 1; i < num_var + 1; ++i) 
			exact[(j - 1) * num_var + (i - 1)] = solution(i * h, (num_var - j + 1) * h);


	T LtwoNorm = T(0), LtwoExact = T(0);
	for (size_t i = 0; i < internalN; ++i)
	{
		LtwoNorm += sqr(sol[i] - exact[i]);
		LtwoExact += sqr(exact[i]);
	}

	const int params[2] = { 5, 8 };
	/*std::cout << "Jacobi solution : ";
	for (size_t i = 0; i < internalN; ++i)
	{
		std::cout << std::setw(params[0]) << std::setprecision(params[1]) << sol[i] << " ";
	}
	std::cout << "\n\n";
	std::cout << "Exact solution : ";
	for (size_t i = 0; i < internalN; ++i)
	{
		std::cout << std::setw(params[0]) << std::setprecision(params[1]) << exact[i] << " ";
	}
	std::cout << "\n\n";*/
	std::cout << "Error L2_Norm = " << sqrt(LtwoNorm / LtwoExact) << "\n";
		
	

}

