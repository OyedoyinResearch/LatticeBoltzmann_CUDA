#ifndef D2Q9_HEADER
#define D2Q9_HEADER

#include "BoltzmannSolver.hpp"

/// @brief 2D Lattice boltzmann solver with the D2Q9 implementation
class D2Q9 : public BoltzmannSolver
{
protected:

public:
	/// <summary>
	/// Base constructor
	/// </summary>
	D2Q9() : BoltzmannSolver(9, 1, 1, 1) {}
	/// <summary>
	/// Construct the D2Q9 solver with the specified sizes
	/// </summary>
	/// <param name="nx"></param>
	/// <param name="ny"></param>
	D2Q9(int nx, int ny) : BoltzmannSolver(9, nx, ny, 1) {}
	/// @brief Setup 2D solver properties
	void InitialiseSolver() override;
};


#endif // D2Q9_HEADER




