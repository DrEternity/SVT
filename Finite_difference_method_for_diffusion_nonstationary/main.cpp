#include "inmost.h"
#include <cstdlib>
#include <iostream>
#include <chrono>

using namespace std;
using namespace INMOST;

int main(int argc, char *argv[])
{
	// INMOST is generally used for MPI-parallel computations
	// For some users errors may occur if there's no Solver::Initialize()
	Solver::Initialize(&argc, &argv);
	double delta_time = atof(argv[1]);

	// Create a sparse matrix (CSR form) and load it
	Sparse::Matrix A;
	string loc = ""; // put your location here or pass in argv
	A.Load(loc + "./matrix.mtx");
        cout << "Размер системы: " << A.Size() << endl;
	// Create a (sparse) vector and load it
	Sparse::Vector rhs;
	rhs.Load(loc + "./rhs.mtx");
	// Create solution vector
	// Initialize as rhs so they have equal size
	Sparse::Vector sol = rhs;
        
	// Create linear solver and set parameters
	Solver S(Solver::INNER_ILU2);
	S.SetParameter("drop_tolerance", "0");
	S.SetParameter("absolute_tolerance", "1e-9");
	S.SetParameter("relative_tolerance", "1e-15");
	// Set matrix AND compute preconditioner (!)
	
	auto begin = std::chrono::steady_clock::now(); 
	S.SetMatrix(A);
	auto end = std::chrono::steady_clock::now();
	auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
	std::cout << "ILU(2): " << elapsed_ms.count() << " ms\n";

        // Run preconditioned BiCGStab
	begin = std::chrono::steady_clock::now();
	bool solved = S.Solve(rhs, sol);
	end = std::chrono::steady_clock::now();
	elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
	std::cout << "BiCGStab: " << elapsed_ms.count() << " ms\n";

	if(!solved){
	    cout << "Linear solver failed: " << S.GetReason() << endl;
	}
    sol.Save("./solution.mtx");
	cout << "Number of iterations: " << S.Iterations() << endl;
	cout << "Residual: " << S.Residual() << endl;
	Solver::Finalize();
	return 0;
}