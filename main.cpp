#include <cmath>
#include <cstdlib>
#include "External/glew/include/GL/glew.h"
#include "External/freeglut/include/GL/freeglut.h"
#include "SingleCore/D2Q9.hpp"

// Lattice Parameters
int LatticeSize = 100;
/// <summary>
/// Test Re
/// </summary>
float ReynoldsNumber = 1000.0f;

// Boundary Conditions
float LidVelocity = 0.20f;
/// <summary>
/// Default density [lattice units]
/// </summary>
float Rho0 = 16.0f;

D2Q9 Solver{};

int NSTEPS = 3000;
int Counter = 0;

const int width = 800;
const int height = 600;

/// <summary>
/// Draw fluid stuff
/// </summary>
void FixedUpdate()
{
	Solver.ComputeDensityFieldCollision();
	Solver.ComputeFieldStream(Solver.density_field);
	Solver.ComputeDensityBoundaryConditions(false);
	Solver.ComputeMacroscopicProperites();
	printf("Step: %d => Error: %g\n", Counter, Solver.ConvergenceError);
	Counter++;

	// Request redisplay
	glutPostRedisplay();
}

void DrawGizmos()
{
	// Clear the buffer
	glClear(GL_COLOR_BUFFER_BIT); 

	float cellSize = width / Solver.NX;

	// Loop through the grid
	for (int x = 0; x < LatticeSize; x++)
	{
		for (int y = 0; y < LatticeSize; y++)
		{
			int dc = Solver.GetScalarIndex(x, y);

			// Compute the origin of the current cell
			float originX = x / 2.0f;
			float originY = y / 2.0f;

			// Get the temperature/density and velocity components
			float temperature = Solver.Temperature[dc];
			float u = Solver.U[dc];
			float v = Solver.V[dc];

			float scale = sqrt(u * u + v * v) / LidVelocity;
			scale *= 100.0f;

			// Set color for velocity vector (magenta)
			glColor3f(1.0f, 0.0f, 1.0f);

			// Draw velocity vector as a line (from origin)
			glBegin(GL_LINES);
			glVertex2f(originX, originY); // Start of the line
			glVertex2f(originX + u * scale, originY + v * scale); // End of the line
			glEnd();
		}
	}

	// Swap the buffers to display the result
	glutSwapBuffers(); 
}

void InitialiseDisplay()
{
	// Initialize GLUT
	int argc = 1;
	char* argv[1] = { (char*)"Display" };
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(width, height);
	glutCreateWindow("OpenGL Lattice Visualization");

	// Initialize GLEW
	glewInit();

	// Set up the drawing function and initialization
	glutDisplayFunc(DrawGizmos);

	// Set the background color
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	// Set up 2D orthographic projection
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0, LatticeSize / 2.0f, 0, LatticeSize / 2.0f);  

	glutIdleFunc(FixedUpdate);

	// Enter the GLUT main loop
	glutMainLoop();
}

int main()
{
	printf("Simulating Lid Driven Cavity Case with the LBM \n");
	printf("Domain size: %ux%u\n", LatticeSize, LatticeSize);
	printf("Re: %g\n", ReynoldsNumber);
	printf("\n");

#pragma region VOID START

	int np1 = LatticeSize + 1;
	double LatticeViscosity = LidVelocity * ((double)np1 / ReynoldsNumber);
	double BaseDensityRelaxationTime = 3.0f * LatticeViscosity + 0.5f;

	Solver = D2Q9(LatticeSize, LatticeSize);
	Solver.collisionMode = CollisionMode::BGK_LES;
	Solver.AmbientTemperature = 0.5f;
	Solver.MaximumTemperature = 1.0f;

	// Set Collision Parameters
	Solver.LatticeViscosity = LatticeViscosity;
	Solver.BaseDensityRelaxationTime = BaseDensityRelaxationTime;
	Solver.BaseTemperatureRelaxationTime = 0.5f;

	// Initialise the fluid distribution field arrays
	Solver.InitialiseSolver();

	InitialiseDisplay();

#pragma endregion

	return 0;
}
