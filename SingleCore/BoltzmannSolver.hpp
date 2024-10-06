#ifndef LBM_HEADER
#define LBM_HEADER

#include <cmath>
#include <cstdio>
#include <malloc.h>
#include "vector3d.hpp"

enum CollisionMode
{
	PureBGK,
	BGK_LES,
	MRT
};

/// @brief Base class for the lattice boltzmann solver
class BoltzmannSolver
{
public:
	/// <summary>
	///
	/// </summary>
	CollisionMode collisionMode;
	/// @brief Cell size along the horizontal/longitudinal axis
	int NX = 10;
	/// @brief Cell size along the vertical axis
	int NY = 10;
	/// @brief Cell size along the depth/lateral axis
	int NZ = 1;
	/// @brief Number of cell neighbours to keep track of
	int Q = 9;
	/// <summary>
	/// Lattice Reynolds number
	/// </summary>
	double ReynoldsNumber = 1000.0f;

	/// <summary>
	/// Lattice Rayleigh Number
	/// </summary>
	double RayleighNumber = 1000.0f;

	/// <summary>
	/// Lattice Prandtl Number
	/// </summary>
	double PrandtlNumber = 0.71f;

	/// <summary>
	/// Default density in the lattice
	/// </summary>
	double AmbientDensity = 2.7f;

	/// <summary>
	/// Ambient temperature on lattice initialisation [lattice units]
	/// </summary>
	double AmbientTemperature = 0.5f;

	/// <summary>
	/// Maximum temperature within the lattice [lattice units]
	/// </summary>
	double MaximumTemperature = 1.0f;

	/// <summary>
	/// [m/s2] * Thermal diffusion constant
	/// </summary>
	double gBeta = 0.01f;

	/// <summary>
	/// Cs to compute eddy viscosity with
	/// </summary>
	double SmagorinskyConstant = 0.16f;

	// Density field collision parameters
	double LatticeViscosity = 0.0f;
	double BaseDensityRelaxationTime = 0.0f;
	double TurbulentDensityRelaxationTime = 0.0f;
	double EffectiveDensityRelaxationTime = 0.0f;

	// Temperature field collision parameters
	double BaseTemperatureRelaxationTime = 0.0f;
	double EffectiveTemperatureRelaxationTime = 0.0f;

	// Density distribution fields
	double *density_field{};
	double *density_field_equilibrum{};
	double *density_field_new{};

	// Temperature distribution fields
	double *temperature_field{};
	double* heatsource_field{};

	// Temporary field for streaming
	double *field_storage{};

	// Macroscopic properties
	double *U{};
	double *V{};
	double *W{};
	double *Density{};
	double *Temperature{};
	/// <summary>
	/// Magnitude of the strain rate |S|
	/// </summary>
	double *Sigma{};
	// Solver specifics
	double *Weights{};
	int *CX{};
	int *CY{};
	int *CZ{};
	int *InverseDirections{};
	double *eX{};
	double *eY{};
	double *eZ{};
	double *X{};
	double *Y{};
	double *Z{};
	double* InletVelocities{};

	double ConvergenceError = 0.0f;

	double *FieldSetInput;
	double *FieldSetOutput;
	double *M;
	double *MInverse;
	double *S;

	/// @brief Get the position of a matrix value in the flattened array
	/// @param i
	/// @param j
	/// @return
	int GetMatrixIndex(int i, int j)
	{
		return j + (i * Q);
	}

	void TransformVector(double *matrix, double *vector, double *&result)
	{
		for (int i = 0; i < Q; i++)
		{
			result[i] = 0.0f;
			for (int j = 0; j < Q; j++)
			{
				result[i] += matrix[GetMatrixIndex(i, j)] * vector[j];
			}
		}
	}

public:
	BoltzmannSolver(int q, int nx, int ny, int nz);
	~BoltzmannSolver();
	void Cleanup();
	virtual void InitialiseSolver();
	void ComputeFieldStream(double *&field);
	void ComputeDensityFieldCollision();
	void ComputeMacroscopicProperites();
	void ComputeTemperatureFieldCollision();
	void ComputeMacroscopicTemperature();

	void ComputeTemperatureBoundaryConditions();
	void ComputeDensityBoundaryConditions(bool full_bounce);

	void ComputeAtmosphereDensityBC();
	void ComputeAtmosphereTemperatureBC();

#pragma region GET FUNCTIONS

	/// @brief Get the position of the specified point in the flattened 1D array
	/// int Row-Major order mode
	/// @param x
	/// @param y
	/// @param z
	/// @return
	int GetScalarIndex(int i, int j)
	{
		int np1 = NX + 1;
		return j + (i * np1);
	}

	/// @brief Get the position of the specified point in the flattened field array
	/// @param x
	/// @param y
	/// @param z
	/// @param k
	/// @return
	int GetFieldIndex(int i, int j, int q)
	{
		int ds = GetScalarIndex(i, j);
		return q + (ds * 9);
	}

#pragma endregion
};

#endif // !LBM_HEADER