/**
 * \file linear_elastic_csm.cpp
 * \brief CSM program executable
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

#include <stdio.h>
#include "../../quasi_1d_euler/inner_prod_vector.hpp"
#include "../lecsm.hpp"
#include "../1D_mesh_tools.hpp"
using namespace std;

// =====================================================================

int main() {

	// Declare the solver
	int nnp = 5;
	LECSM csm(nnp);

	// Define material properties
	double E = 100000000;		// Young's Modulus (Rubber) (Pa)
	double t = 0.03;				// Thickness of the beam elements (meters)
	double w = 1;					  // Width of the domain (meters)
	double h = 1;						// Height of the domain (meters)
	csm.set_material(E, t, w, h);

	// Create problem mesh
	double length = 10.0;
  InnerProdVector x_coord(nnp, 0.0);
  InnerProdVector y_coord(nnp, 0.0);
  InnerProdVector area(nnp, 0.0);
  for (int i = 0; i < nnp; i++) {
    // evenly spaced nodes along the x
    x_coord(i) = i*length/(nnp-1);
    // parabolic nozzle wall for y coords
    y_coord(i) = 0.01*(length-x_coord(i))*x_coord(i);
  }
  csm.GenerateMesh(x_coord, y_coord);

  // determine the nodal structural boundary conditions
  InnerProdVector BCtype(3*nnp, 0.0);
  InnerProdVector BCval(3*nnp, 0.0);
  for (int i=0; i<nnp; i++) {
    BCtype(3*i) = 0;
    BCtype(3*i+1) = -1;
    BCtype(3*i+2) = -1;
    BCval(3*i) = 0;
    BCval(3*i+1) = 0;
    BCval(3*i+2) = 0;
  }
  BCtype(0) = 0;
  BCtype(1) = 0;
  BCtype(2) = -1;
  BCtype(3*nnp-3) = 0;
  BCtype(3*nnp-2) = 0;
  BCtype(3*nnp-1) = -1;
  csm.SetBoundaryConds(BCtype, BCval);

  // Specify pressure
  InnerProdVector press(nnp, 20.0);
  csm.set_press(press);

	// Call FEA solver
	csm.Solve();

	csm.CalcResidual();
	InnerProdVector & res = csm.get_res();

	printf("Printing residual for inspection:\n");
	for(int i=0; i<3*nnp; i++) {
		printf("|  %f  |\n", res(i));
	}

	return 0;
}