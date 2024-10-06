#include "BoltzmannSolver.hpp"

/// <summary>
/// Base constructor for the solver.  Configure the required parameters and assign the needed memory
/// </summary>
/// <param name="q"></param>
/// <param name="nx"></param>
/// <param name="ny"></param>
/// <param name="nz"></param>
BoltzmannSolver::BoltzmannSolver(int q, int nx, int ny, int nz)
{
	// Set input values
	Q = q;
	NX = nx;
	NY = ny;
	NZ = nz;

	// Allocate memory for the arrays
	int np1 = NX + 1;
	int mp1 = NY + 1;
	int scalar_size = np1 * mp1;
	int field_size = np1 * mp1 * Q;

	// Scalars
	U = (double*)malloc(scalar_size * sizeof(double));
	V = (double*)malloc(scalar_size * sizeof(double));
	W = (double*)malloc(scalar_size * sizeof(double));
	Density = (double*)malloc(scalar_size * sizeof(double));
	Temperature = (double*)malloc(scalar_size * sizeof(double));
	Sigma = (double*)malloc(scalar_size * sizeof(double));

	X = (double*)malloc(np1 * sizeof(double));
	Y = (double*)malloc(mp1 * sizeof(double));
	Z = (double*)malloc(NZ * sizeof(double));
	InletVelocities = (double*)malloc(mp1 * sizeof(double));

	Weights = (double*)malloc(Q * sizeof(double));
	CX = (int*)malloc(Q * sizeof(int));
	CY = (int*)malloc(Q * sizeof(int));
	CZ = (int*)malloc(Q * sizeof(int));
	eX = (double*)malloc(Q * sizeof(double));
	eY = (double*)malloc(Q * sizeof(double));
	eZ = (double*)malloc(Q * sizeof(double));
	InverseDirections = (int*)malloc(Q * sizeof(int));

	int matrix_size = Q * Q;
	FieldSetInput = (double*)malloc(Q * sizeof(double));
	FieldSetOutput = (double*)malloc(Q * sizeof(double));
	S = (double*)malloc(Q * sizeof(double));
	M = (double*)malloc(matrix_size * sizeof(double));
	MInverse = (double*)malloc(matrix_size * sizeof(double));

	// Density Fields
	density_field = (double*)malloc(field_size * sizeof(double));
	density_field_new = (double*)malloc(field_size * sizeof(double));
	density_field_equilibrum = (double*)malloc(field_size * sizeof(double));
	// Temperature Fields
	temperature_field = (double*)malloc(field_size * sizeof(double));
	field_storage = (double*)malloc(field_size * sizeof(double));
	heatsource_field = (double*)malloc(field_size * sizeof(double));

	// Clean the assigned arrays
	for (int x = 0; x <= NX; x++)
	{
		for (int y = 0; y <= NY; y++)
		{
			for (int z = 0; z < NZ; z++)
			{
				// Set initial scalar values
				int ds = GetScalarIndex(x, y);
				U[ds] = 0.0f;
				V[ds] = 0.0f;
				W[ds] = 0.0f;
				Density[ds] = 0.0f;
				Temperature[ds] = 0.0f;
				Sigma[ds] = 0.0f;

				// Set initial field values
				for (int k = 0; k < Q; k++)
				{
					int dk = k + ds * Q;
					density_field[dk] = 0.0f;
					temperature_field[dk] = 0.0f;
					density_field_equilibrum[dk] = 0.0f;
					field_storage[dk] = 0.0f;
					density_field_new[dk] = 0.0f;
					heatsource_field[dk] = 0.0f;
				}
			}
		}
	}
}

/// @brief Destructor,
BoltzmannSolver::~BoltzmannSolver()
{
	// Cleanup();
}

/// @brief deallocate array memories
void BoltzmannSolver::Cleanup()
{
	free(density_field);
	free(temperature_field);
	free(density_field_equilibrum);
	free(density_field_new);
	free(field_storage);

	free(U);
	free(V);
	free(W);
	free(Density);
	free(Temperature);
	free(Sigma);

	free(Weights);
	free(CX);
	free(CY);
	free(CZ);
	free(eX);
	free(eY);
	free(eZ);
	free(X);
	free(Y);
	free(Z);
}

/// @brief Initialise the distribution field to the equilibrum field. must be called after the velocity arrays have been configured
void BoltzmannSolver::InitialiseSolver()
{
	// Halfway grid
	float dx = 1.0f / NX;
	float dy = 1.0f / NY;
	float dz = 1.0f / NZ;
	X[0] = 0.0f;
	Y[0] = 0.0f;
	Z[0] = 0.0f;
	for (int i = 1; i <= NX; i++)
	{
		X[i] = X[i - 1] + dx;
	}

	for (int i = 1; i <= NY; i++)
	{
		Y[i] = Y[i - 1] + dy;
	}

	// Set the field properties
	for (int i = 0; i <= NX; i++)
	{
		for (int j = 0; j <= NY; j++)
		{
			// Initialise scalar values
			int dc = GetScalarIndex(i, j);
			Density[dc] = AmbientDensity;
			Temperature[dc] = AmbientTemperature;

			// Initialise field values
			for (int a = 0; a < Q; a++)
			{
				int df = GetFieldIndex(i, j, a);
				double ginit = Weights[a] * AmbientTemperature;
				temperature_field[df] = ginit;

				double finit = Weights[a] * AmbientDensity;
				density_field[df] = finit;
				density_field_equilibrum[df] = finit;
			}
		}
	}
}

#pragma region INNER SOLVER FUNCTIONS

/// @brief Streams the field values along the specified velocity indices
/// It actually works!! prepped on 22-09-2024 12:40am
void BoltzmannSolver::ComputeFieldStream(double*& field)
{
	// loop over all interior voxels
	for (int x = 0; x <= NY; x++)
	{
		for (int y = 0; y <= NX; y++)
		{
			for (int a = 0; a < Q; a++)
			{
				// Negative because I am checking the values at the previous location i.e t - dt
				int iS = x - CX[a];
				int jS = y - CY[a];
				int neighbour_field_index = GetFieldIndex(iS, jS, a);
				int current_field_index = GetFieldIndex(x, y, a);

				if (
					iS < 0 || iS > NX ||
					jS < 0 || jS > NY)
				{
					// Copy the elements at the boundaries and leave unchanged
					field_storage[current_field_index] = field[current_field_index];
				}
				else
				{
					// Push out the elements in the interior cells (current = value along previous direction)
					// e.g value at current 1 axis = value at previous 1 position
					field_storage[current_field_index] = field[neighbour_field_index];
				}
			}
		}
	}

	// Copy the elements from the temporary array back into the functional array
	for (int x = 0; x <= NY; x++)
	{
		for (int y = 0; y <= NX; y++)
		{
			for (int a = 0; a < Q; a++)
			{
				int current_field_index = GetFieldIndex(x, y, a);
				field[current_field_index] = field_storage[current_field_index];
			}
		}
	}
}

/// @brief STEP:X Compute the new distribution function by relaxing the streamed field towards the equilibrium field
void BoltzmannSolver::ComputeDensityFieldCollision()
{
	// loop over all interior voxels
	for (int i = 0; i <= NX; i++)
	{
		for (int j = 0; j <= NY; j++)
		{
			// Get working index
			int ds = GetScalarIndex(i, j);
			double temperature = Temperature[ds];
			double density = Density[ds];

			for (int a = 0; a < Q; a++)
			{
				// 1.0 Calculate Equilibrium distribution field
				double edotu =
					eX[a] * U[ds] +
					eY[a] * V[ds] +
					eZ[a] * W[ds];
				double udotu =
					U[ds] * U[ds] +
					V[ds] * V[ds] +
					W[ds] * W[ds];
				double feq = density * Weights[a] *
					(1.0f + (3.0f * edotu) + (4.50f * edotu * edotu) - 1.50f * udotu);

				// 2.0 Calculate buoyancy force
				double wk = Weights[a];
				double force = 3.0f * wk * gBeta * (temperature - AmbientTemperature) * eY[a] * density;
				if (i == 0 || i == NX)
				{
					force = 0.0;
				}
				if (j == 0 || j == NY)
				{
					force = 0.0;
				}
				int df = GetFieldIndex(i, j, a);

				if (collisionMode == CollisionMode::BGK_LES)
				{
					// Calculate the turbulent viscosity relaxation time [eqn 22]
					double delta_x_sqr = 1.0f;
					double Q_bar = Sigma[ds];
					double Cs_sqr = SmagorinskyConstant * SmagorinskyConstant;
					double tau_sqr = BaseDensityRelaxationTime * BaseDensityRelaxationTime;
					double tau_RHS = tau_sqr + 18.0f * Cs_sqr * delta_x_sqr * Q_bar;
					TurbulentDensityRelaxationTime = 0.5f * (sqrt(tau_RHS) - BaseDensityRelaxationTime);

					// Calculate the effective relaxation time
					EffectiveDensityRelaxationTime = BaseDensityRelaxationTime + TurbulentDensityRelaxationTime;
					// Note that tau_eff now varies from point to point in the domain, and is
					// larger for large strain rates. If the strain rate is zero, tau_eff = 0 and we
					// revert back to the original (laminar) LBM scheme where tau_eff = tau.
				}
				else if (collisionMode == CollisionMode::MRT)
				{

				}
				else
				{
					EffectiveDensityRelaxationTime = BaseDensityRelaxationTime;
				}


				double ux = U[ds];
				double uy = V[ds];
				double p = Density[ds];

				double m0 = p;
				double m1 = p - (3.0f * p * ((ux * ux) + (uy * uy)));
				double m2 = (9.0f * p * ux * ux * uy * uy) - (3.0f * p * ((ux * ux) + (uy * uy))) + p;
				double m3 = p * ux;
				double m4 = (3.0f * p * ux * ux * ux) - (p * ux);
				double m5 = p * uy;
				double m6 = (3.0f * p * uy * uy * uy) - (p * uy);
				double m7 = p * ((ux * ux) - (uy * uy));
				double m8 = p * (ux * uy);

				// 3.0 Implement the BGK collision scheme with Eddy viscosity computed
				// using the Smagorinsky model
				// BGK Collision
				double fc = density_field[df];
				double omega = 1.0f / EffectiveDensityRelaxationTime;
				double fi = (omega * feq) + ((1.0f - omega) * fc) + force;
				// Store
				density_field_equilibrum[df] = feq;
				density_field[df] = fi;
			}
		}
	}
}

/// @brief STEP:X Calculate the new velocity components and density from the streamed field
void BoltzmannSolver::ComputeMacroscopicProperites()
{
	// loop over all interior voxels
	for (int i = 0; i <= NX; i++)
	{
		for (int j = 0; j <= NY; j++)
		{
			// Get working index
			int dc = GetScalarIndex(i, j);

			// Compute convergence error with the L2 norm
			double temp1 = 0;
			double temp2 = 0;
			for (int a = 0; a < Q; a++)
			{
				int dk = GetFieldIndex(i, j, a);
				double fi = density_field_new[dk];
				double fo = density_field[dk];
				temp1 += pow(fi - fo, 2.0f);
				temp2 += pow(fi, 2.0f);
			}

			ConvergenceError = pow(temp1 / temp2, 0.5f);

			// 1.0 Push the density field into the error compute container
			for (int a = 0; a < Q; a++)
			{
				int df = GetFieldIndex(i, j, a);
				density_field_new[df] = density_field[df];
			}

			// 2.0 Calculate the interior density values
			double summation_fi = 0.0f;
			for (int a = 0; a < Q; a++)
			{
				int df = GetFieldIndex(i, j, a);
				double fi_temp = density_field[df];
				summation_fi += fi_temp;
			}

			Density[dc] = summation_fi;

			// 3.0 Calculate the interior velocity components
			double summation_fi_x = 0.0f;
			double summation_fi_y = 0.0f;
			double summation_fi_z = 0.0f;
			for (int a = 0; a < Q; a++)
			{
				double eix = eX[a];
				double eiy = eY[a];
				double eiz = eZ[a];
				int df = GetFieldIndex(i, j, a);
				double fi_temp = density_field[df];

				summation_fi_x += (fi_temp * eix);
				summation_fi_y += (fi_temp * eiy);
				summation_fi_z += (fi_temp * eiz);
			}

			// Compute the velocity from the momentum sum
			U[dc] = summation_fi_x / summation_fi;
			V[dc] = summation_fi_y / summation_fi;
			W[dc] = summation_fi_z / summation_fi;

			// 4.0 Calculate the strain rate components
			double sum_xx = 0.0, sum_xy = 0.0, sum_xz = 0.0;
			double sum_yx = 0.0, sum_yy = 0.0, sum_yz = 0.0;
			double sum_zx = 0.0, sum_zy = 0.0, sum_zz = 0.0;
			for (int a = 1; a < Q; a++)
			{
				int df = GetFieldIndex(i, j, a);
				double eix = eX[a];
				double eiy = eY[a];
				double eiz = eZ[a];

				double fi = density_field[df];
				double feq = density_field_equilibrum[df];
				double fi_feq = fi - feq;
				sum_xx += (fi_feq * eix * eix);
				sum_xy += (fi_feq * eix * eiy);
				sum_xz += (fi_feq * eix * eiz);
				sum_yx += (fi_feq * eiy * eix);
				sum_yy += (fi_feq * eiy * eiy);
				sum_yz += (fi_feq * eiy * eiz);
				sum_zx += (fi_feq * eiz * eix);
				sum_zy += (fi_feq * eiz * eiy);
				sum_zz += (fi_feq * eiz * eiz);
			}

			double Sigma_Inner =
				pow(sum_xx, 2) + pow(sum_xy, 2) + pow(sum_xz, 2) +
				pow(sum_yx, 2) + pow(sum_yy, 2) + pow(sum_yz, 2) +
				pow(sum_zx, 2) + pow(sum_zy, 2) + pow(sum_zz, 2);
			Sigma[dc] = pow(Sigma_Inner, 0.5);
		}
	}
}

/// @brief STEP:X Compute the new temperature distribution function by relaxing the streamed field towards the equilibrum field
void BoltzmannSolver::ComputeTemperatureFieldCollision()
{
	// loop over all interior voxels
	for (int i = 0; i <= NX; i++)
	{
		for (int j = 0; j <= NY; j++)
		{
			// Get working index
			int dc = GetScalarIndex(i, j);
			double temperature = Temperature[dc];

			for (int a = 0; a < Q; a++)
			{
				// Get field index
				int df = GetFieldIndex(i, j, a);
				double Wi = Weights[a];

				// 1.0 Calculate Equilibrium distribution field
				double edotu = (eX[a] * U[dc]) + (eY[a] * V[dc]) + (eZ[a] * W[dc]);
				double geq = temperature * Wi * (1.0f + 3.0f * edotu);

				// Leave space for more complex relaxation implementation
				EffectiveTemperatureRelaxationTime = BaseTemperatureRelaxationTime;

				// 2.0 Implement the BGK collision scheme
				double Qi = Wi * heatsource_field[df];
				double gc = temperature_field[df];
				double omega_t = 1.0f / EffectiveTemperatureRelaxationTime;
				double gi = (gc * (1.0f - omega_t)) + (omega_t * geq);
				// Store
				temperature_field[df] = gi + Qi;
			}
		}
	}
}

/// @brief STEP:X Calculate the new temperature from the streamed field
void BoltzmannSolver::ComputeMacroscopicTemperature()
{
	// loop over all interior voxels
	for (int i = 0; i <= NX; i++)
	{
		for (int j = 0; j <= NY; j++)
		{
			// Get working index
			int ds = GetScalarIndex(i, j);

			double summation_gi = 0.0f;
			for (int a = 0; a < Q; a++)
			{
				int df = GetFieldIndex(i, j, a);
				double gi_temp = temperature_field[df];
				summation_gi += gi_temp;
			}

			Temperature[ds] = summation_gi;
		}
	}
}

#pragma endregion

#pragma region BOUNDARY CONDITIONS

/// @brief
void BoltzmannSolver::ComputeTemperatureBoundaryConditions()
{
	int i, j, k;
	// Left boundary, T = tw (usually tw = 1.0 for this boundary condition)
	for (j = 0; j <= NY; j++)
	{
		int dk1 = GetFieldIndex(0, j, 1);
		int dk5 = GetFieldIndex(0, j, 5);
		int dk8 = GetFieldIndex(0, j, 8);

		int dk3 = GetFieldIndex(0, j, 3);
		int dk7 = GetFieldIndex(0, j, 7);
		int dk6 = GetFieldIndex(0, j, 6);

		temperature_field[dk1] = MaximumTemperature * (Weights[1] + Weights[3]) - temperature_field[dk3];
		temperature_field[dk5] = MaximumTemperature * (Weights[5] + Weights[7]) - temperature_field[dk7];
		temperature_field[dk8] = MaximumTemperature * (Weights[8] + Weights[6]) - temperature_field[dk6];
	}

	// Right boundary, T = 0.0 (reflective condition with zero temperature)
	for (j = 0; j <= NY; j++)
	{
		int dk6 = GetFieldIndex(NX, j, 6);
		int dk3 = GetFieldIndex(NX, j, 3);
		int dk7 = GetFieldIndex(NX, j, 7);

		int dk8 = GetFieldIndex(NX, j, 8);
		int dk1 = GetFieldIndex(NX, j, 1);
		int dk5 = GetFieldIndex(NX, j, 5);

		temperature_field[dk6] = -temperature_field[dk8];
		temperature_field[dk3] = -temperature_field[dk1];
		temperature_field[dk7] = -temperature_field[dk5];
	}

	// Top boundary, adiabatic (no temperature flux, temperature gradient = 0)
	for (i = 1; i <= NX - 1; i++)
	{
		for (k = 0; k < Q; k++)
		{
			int dkl = GetFieldIndex(i, NY, k);
			int dkr = GetFieldIndex(i, NY - 1, k);
			temperature_field[dkl] = temperature_field[dkr];
		}
	}

	// Bottom boundary, T = tw (fixed temperature, tw is typically set to 1.0)
	for (i = 1; i <= NX - 1; i++)
	{
		for (k = 0; k < Q; k++)
		{
			int dkl = GetFieldIndex(i, 0, k);
			int dkr = GetFieldIndex(i, 1, k);
			temperature_field[dkl] = temperature_field[dkr];
		}
	}
}

/// @brief
void BoltzmannSolver::ComputeDensityBoundaryConditions(bool full_bounce)
{
	int i, j;
	for (j = 0; j <= NY; j++)
	{
		int dk1 = GetFieldIndex(0, j, 1);
		int dk5 = GetFieldIndex(0, j, 5);
		int dk8 = GetFieldIndex(0, j, 8);

		int dk3 = GetFieldIndex(0, j, 3);
		int dk7 = GetFieldIndex(0, j, 7);
		int dk6 = GetFieldIndex(0, j, 6);

		density_field[dk1] = density_field[dk3];
		density_field[dk5] = density_field[dk7];
		density_field[dk8] = density_field[dk6];
	}

	for (i = 0; i <= NX; i++)
	{
		if (full_bounce)
		{
			int dk4 = GetFieldIndex(i, NY, 4);
			int dk8 = GetFieldIndex(i, NY, 8);
			int dk7 = GetFieldIndex(i, NY, 7);
			//
			int dk2 = GetFieldIndex(i, NY, 2);
			int dk6 = GetFieldIndex(i, NY, 6);
			int dk5 = GetFieldIndex(i, NY, 5);
			//
			density_field[dk4] = density_field[dk2];
			density_field[dk8] = density_field[dk6];
			density_field[dk7] = density_field[dk5];
		}
		else
		{
			int NY_1 = NY;
			int dc0 = GetFieldIndex(i, NY_1, 0);
			int dc1 = GetFieldIndex(i, NY_1, 1);
			int dc2 = GetFieldIndex(i, NY_1, 2);
			int dc3 = GetFieldIndex(i, NY_1, 3);
			int dc5 = GetFieldIndex(i, NY_1, 5);
			int dc6 = GetFieldIndex(i, NY_1, 6);

			double f0 = density_field[dc0];
			double f1 = density_field[dc1];
			double f2 = density_field[dc2];
			double f3 = density_field[dc3];
			double f5 = density_field[dc5];
			double f6 = density_field[dc6];

			double rhoN = f0 + f1 + f3 + (2.0f * (f2 + f6 + f5));

			int dk4 = GetFieldIndex(i, NY_1, 4);
			int dk8 = GetFieldIndex(i, NY_1, 8);
			int dk7 = GetFieldIndex(i, NY_1, 7);

			int dk2 = GetFieldIndex(i, NY_1, 2);
			int dk6 = GetFieldIndex(i, NY_1, 6);
			int dk5 = GetFieldIndex(i, NY_1, 5);

			density_field[dk4] = f2;
			density_field[dk8] = f6 + (rhoN * 0.1f / 6.0f);
			density_field[dk7] = f5 - (rhoN * 0.1f / 6.0f);
		}
	}

	for (i = 0; i <= NX; i++)
	{
		int dk2 = GetFieldIndex(i, 0, 2);
		int dk5 = GetFieldIndex(i, 0, 5);
		int dk6 = GetFieldIndex(i, 0, 6);

		int dk4 = GetFieldIndex(i, 0, 4);
		int dk7 = GetFieldIndex(i, 0, 7);
		int dk8 = GetFieldIndex(i, 0, 8);

		density_field[dk2] = density_field[dk4];
		density_field[dk5] = density_field[dk7];
		density_field[dk6] = density_field[dk8];
	}

	for (j = 0; j <= NY; j++)
	{
		int dk3 = GetFieldIndex(NX, j, 3);
		int dk7 = GetFieldIndex(NX, j, 7);
		int dk6 = GetFieldIndex(NX, j, 6);

		int dk1 = GetFieldIndex(NX, j, 1);
		int dk5 = GetFieldIndex(NX, j, 5);
		int dk8 = GetFieldIndex(NX, j, 8);

		density_field[dk3] = density_field[dk1];
		density_field[dk7] = density_field[dk5];
		density_field[dk6] = density_field[dk8];
	}
}

#pragma region ATMOSPHERE EXTENSION

/// <summary>
/// Compute the boundary conditions for the atmospheric model
/// </summary>
void BoltzmannSolver::ComputeAtmosphereDensityBC()
{
	// Left Wall => inlet with computed atmosphere boundary layer velocity 
			// Zho-He Boundary with a specified velocity
	int i, j;
	for (j = 0; j <= NY; j++)
	{
		int dk1 = GetFieldIndex(0, j, 1);
		int dk5 = GetFieldIndex(0, j, 5);
		int dk8 = GetFieldIndex(0, j, 8);

		double f0 = density_field[GetFieldIndex(0, j, 0)];
		double f2 = density_field[GetFieldIndex(0, j, 2)];
		double f4 = density_field[GetFieldIndex(0, j, 4)];
		double f3 = density_field[GetFieldIndex(0, j, 3)];
		double f6 = density_field[GetFieldIndex(0, j, 6)];
		double f7 = density_field[GetFieldIndex(0, j, 7)];

		double uw = InletVelocities[j];
		double rho_w = (1.0f / (1.0f - uw)) * (f0 + f2 + f4 + (2.0f * (f3 + f6 + f7)));
		double f1 = f3 + (0.6666667f * rho_w * uw);
		double f5 = f7 - (0.5f * (f2 - f4)) + (0.16666667f * rho_w * uw);
		double f8 = f6 + (0.5f * (f2 - f4)) + (0.16666667f * rho_w * uw);

		density_field[dk1] = f1;
		density_field[dk5] = f5;
		density_field[dk8] = f8;
	}

	// Top Wall => Free-Slip boundary condition
	for (i = 0; i <= NX; i++)
	{
		int dk4 = GetFieldIndex(i, NY, 4);
		int dk8 = GetFieldIndex(i, NY, 8);
		int dk7 = GetFieldIndex(i, NY, 7);
		//           
		int dk2 = GetFieldIndex(i, NY, 2);
		int dk6 = GetFieldIndex(i, NY, 6);
		int dk5 = GetFieldIndex(i, NY, 5);
		//           
		density_field[dk4] = density_field[dk2];
		density_field[dk8] = density_field[dk5];
		density_field[dk7] = density_field[dk6];
	}

	// Bottom wall => No-Slip boundary condition
	for (i = 0; i <= NX; i++)
	{
		int dk2 = GetFieldIndex(i, 0, 2);
		int dk5 = GetFieldIndex(i, 0, 5);
		int dk6 = GetFieldIndex(i, 0, 6);

		int dk4 = GetFieldIndex(i, 0, 4);
		int dk7 = GetFieldIndex(i, 0, 7);
		int dk8 = GetFieldIndex(i, 0, 8);

		density_field[dk2] = density_field[dk4];
		density_field[dk5] = density_field[dk7];
		density_field[dk6] = density_field[dk8];
	}
	// Right Wall (Outlet => density)
	// Zho-He Boundary with a specified pressure (density)
	for (j = 0; j <= NY; j++)
	{
		int dk3 = GetFieldIndex(NX, j, 3);
		int dk7 = GetFieldIndex(NX, j, 7);
		int dk6 = GetFieldIndex(NX, j, 6);

		double f1 = density_field[GetFieldIndex(NX, j, 1)];
		double f5 = density_field[GetFieldIndex(NX, j, 5)];
		double f8 = density_field[GetFieldIndex(NX, j, 8)];
		double f0 = density_field[GetFieldIndex(NX, j, 0)];
		double f2 = density_field[GetFieldIndex(NX, j, 2)];
		double f4 = density_field[GetFieldIndex(NX, j, 4)];

		double rho_outlet = AmbientDensity;
		double ux = -1.0f + (f0 + f2 + f4 + (2.0f * (f1 + f5 + f8))) / rho_outlet;
		double f3 = f1 - (0.666667f * rho_outlet * ux);
		double f7 = f5 + (0.5f * (f2 - f4)) - (0.1666667 * rho_outlet * ux);
		double f6 = f8 - (0.5f * (f2 - f4)) - (0.1666667 * rho_outlet * ux);

		density_field[dk3] = f3;
		density_field[dk7] = f7;
		density_field[dk6] = f6;
	}
}

/// <summary>
/// Compute the temperature boundary conditions for the atmospheric model
/// </summary>
void BoltzmannSolver::ComputeAtmosphereTemperatureBC()
{
	int i, j, k;
	// Left boundary, T = ambient temperature
	for (j = 0; j <= NY; j++)
	{
		int dk1 = GetFieldIndex(0, j, 1);
		int dk5 = GetFieldIndex(0, j, 5);
		int dk8 = GetFieldIndex(0, j, 8);

		int dk3 = GetFieldIndex(0, j, 3);
		int dk7 = GetFieldIndex(0, j, 7);
		int dk6 = GetFieldIndex(0, j, 6);

		temperature_field[dk1] = AmbientTemperature * (Weights[1] + Weights[3]) - temperature_field[dk3];
		temperature_field[dk5] = AmbientTemperature * (Weights[5] + Weights[7]) - temperature_field[dk7];
		temperature_field[dk8] = AmbientTemperature * (Weights[8] + Weights[6]) - temperature_field[dk6];
	}

	// Right boundary, T = ambient temperature
	for (j = 0; j <= NY; j++)
	{
		int dk6 = GetFieldIndex(NX, j, 6);
		int dk3 = GetFieldIndex(NX, j, 3);
		int dk7 = GetFieldIndex(NX, j, 7);

		int dk8 = GetFieldIndex(NX, j, 8);
		int dk1 = GetFieldIndex(NX, j, 1);
		int dk5 = GetFieldIndex(NX, j, 5);

		temperature_field[dk6] = AmbientTemperature * (Weights[6] + Weights[8]) - temperature_field[dk8];
		temperature_field[dk3] = AmbientTemperature * (Weights[3] + Weights[1]) - temperature_field[dk1];
		temperature_field[dk7] = AmbientTemperature * (Weights[7] + Weights[5]) - temperature_field[dk5];
	}

	// Top boundary, adiabatic (no temperature flux, temperature gradient = 0)
	for (i = 1; i <= NX - 1; i++)
	{
		for (k = 0; k < Q; k++)
		{
			int dkl = GetFieldIndex(i, NY, k);
			int dkr = GetFieldIndex(i, NY - 1, k);
			temperature_field[dkl] = temperature_field[dkr];
		}
	}

	// Bottom boundary, T = tw (fixed temperature, tw is typically set to 1.0)
	for (i = 1; i <= NX - 1; i++)
	{
		for (k = 0; k < Q; k++)
		{
			int dkl = GetFieldIndex(i, 0, k);
			int dkr = GetFieldIndex(i, 1, k);
			temperature_field[dkl] = temperature_field[dkr];
		}
	}
}

#pragma endregion

#pragma endregion
