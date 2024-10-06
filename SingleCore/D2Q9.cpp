#include "D2Q9.hpp"


/// @brief Initialise the weights and streaming velocity vectors
void D2Q9::InitialiseSolver()
{
	// Structure

	//      6        2        5
	//               ^
	//               |
	//               |
	//      3 <----- 0 -----> 1
	//               |
	//               |
	//               v
	//      7        4        8

	// Set inner properties
	Q = 9;

	// Configure the weights
	const float w0 = 4.0f / 9.0f;
	const float wa = 1.0f / 9.0f;
	const float wd = 1.0f / 36.0f;

	// Central Point
	Weights[0] = w0;
	// Axial weights
	Weights[1] = wa;
	Weights[2] = wa;
	Weights[3] = wa;
	Weights[4] = wa;
	// Diagonal weights
	Weights[5] = wd;
	Weights[6] = wd;
	Weights[7] = wd;
	Weights[8] = wd;

	// Configure the streaming velocity directions
	const float c_float = 1.0f;
	float cmx[9] = { 0, c_float, 0, -c_float, 0, c_float, -c_float, -c_float, c_float };
	float cmy[9] = { 0, 0, c_float, 0, -c_float, c_float, c_float, -c_float, -c_float };
	for (int i = 0; i < Q; i++)
	{
		CX[i] = (int)cmx[i];
		CY[i] = (int)cmy[i];
		CZ[i] = 0;

		eX[i] = cmx[i];
		eY[i] = cmy[i];
		eZ[i] = 0.0f;
	}

	// Define opposite directions (useful for implementing bounce back)
	InverseDirections[0] = 0;
	InverseDirections[1] = 3;
	InverseDirections[2] = 4;
	InverseDirections[3] = 1;
	InverseDirections[4] = 2;
	InverseDirections[5] = 7;
	InverseDirections[6] = 8;
	InverseDirections[7] = 5;
	InverseDirections[8] = 6;

	double m_vals[81] = 
	{ 
		0, 0, 0, 0, 0, 1, -1, 1, -1,
		0, 1, -1, 1, -1, 0, 0, 0, 0,
		0, 0, -2, 0, 2, 1, 1, -1, -1,
		0, 0, 1, 0, -1, 1, 1, -1, -1,
		0, -2, 0, 2, 0, 1, -1, -1, 1,
		0, 1, 0, -1, 0, 1, -1, -1, 1,
		4, -2, -2, -2, -2, 1, 1, 1, 1,
		-4, -1, -1, -1, -1, 2, 2, 2, 2,
		1, 1, 1, 1, 1, 1, 1, 1, 1
	};
	M = m_vals;

	double f0_0 = 0.0f;
	double f1_4 = 1.0f/4.0f;
	double f1_6 = 1.0f/6.0f;
	double f1_9 = 1.0/9.0f;
	double f1_12 = 1.0f/12.0f;
	double f1_18 = 1.0f/18.0f;
	double f1_36 = 1.0f/36.0f;

	double m_inv[81] =
	{
		f1_9, f1_18, f1_36, f1_6, f1_12, -f1_6, -f1_12, f0_0, -f1_4,
		f1_9, f1_18, f1_36, -f1_6, -f1_12, -f1_6, -f1_12, f0_0, f1_4,
		f1_9, f1_18, f1_36, -f1_6, -f1_12, f1_6, f1_12, f0_0, -f1_4,
		f1_9, f1_18, f1_36, f1_6, f1_12, f1_6, f1_12, f0_0, f1_4,
		f1_9, -f1_36, -f1_18, f0_0, f0_0, -f1_6, f1_6, -f1_4, f0_0,
		f1_9, -f1_36, -f1_18, -f1_6, f1_6, f0_0, f0_0, f1_4, f0_0,
		f1_9, -f1_36, -f1_18, f0_0, f0_0, f1_6, -f1_6, -f1_4, f0_0,
		f1_9, -f1_36, -f1_18, f1_6, -f1_6, f0_0, f0_0, f1_4, f0_0,
		f1_9, -f1_9, f1_9, f0_0, f0_0, f0_0, f0_0, f0_0, f0_0,
	};
	MInverse = m_inv;

	double se = 1.64f;
	double sv = 1.0f/BaseDensityRelaxationTime;
	double sE = 1.54f;
	double sq = (8.0f * (2.0f - sv))/(8.0f - sv);
	double s_vector[9] = {0.0f, se, sE, 0, sq, 0.0, sq, sv, sv};
	S = s_vector;

	// Initiliase the distribution fields 
	BoltzmannSolver::InitialiseSolver();
}
