# parallel-computation-2_Helmholtz-equation

Contains serial and parallel (openMP) implementations of following algorithms:

* Helholtz equation solver (internal iterations using Jacobi method)
* Helholtz equation solver (internal iterations using Seidel method with black-red ordering)

Note that present implementations are intended for academic purposes, as such they are not meant to be used in any sort of high-performance production code.

## Compilation

* Recommended compiler: Intel C++ Compiler
* Requires C++17 support
* Requires openMP implementation

## Usage

Helmholtz equation:<br>
$$-u_{xx} - u_{yy} + k^2 u(x, y) = f(x, y), \quad (x, y) \in D = [0, L_x]\times[0, L_y], $$
$$u(0, y)  = \phi(y),\quad u(x, 0) = \xi(x),\quad  u(L_x, y) = \theta(y),\quad  u(x, L_y) = \eta(x),$$
where<br>
$k$ - wave_number, $f$ - right_part.<br>

Adjust grid size $N$, area sizes $L_x$ and $L_y$, wave number $k$, precision and boundary conditions in "main.cpp" to configure testing parameters. Parallel implementations assume $(N - 2)$ to be a multiple of num_threads.
