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
	int nnp = 51;
	LECSM csm(nnp);

	// Define material properties
	double E = 69000000000; // Young's Modulus (Rubber) (Pa)
	double t = 0.1;				  // Thickness of the beam elements (meters)
	double w = 1;					  // Width of the beam elements (meters)
	double h = 1;						// Maximum nozzle height (meters)
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
    y_coord(i) = 0.1*(length-x_coord(i))*x_coord(i);
  }
  csm.GenerateMesh(x_coord, y_coord);

// =====================================================================
// BOUNDARY CONDITIONS
// =====================================================================

#if 0
  // ~~~~~ CANTILEVER BEAM ~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  InnerProdVector BCtype(3*nnp, 0.0);
  InnerProdVector BCval(3*nnp, 0.0);
  for (int i=0; i<nnp; i++) {
    BCtype(3*i) = -1;
    BCtype(3*i+1) = -1;
    BCtype(3*i+2) = -1;
    BCval(3*i) = 0;
    BCval(3*i+1) = 0;
    BCval(3*i+2) = 0;
  }
  BCtype(0) = 0;
  BCtype(1) = 0;
  BCtype(2) = 0;
  BCtype(3*nnp-2) = 1;
  BCval(3*nnp-2) = -1000;
#else
  // ~~~~~ PARABOLIC NOZZLE ~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
  BCtype(2) = 0;
  BCtype(3*nnp-3) = 0;
  BCtype(3*nnp-2) = 0;
  BCtype(3*nnp-1) = 0;
#endif
  csm.SetBoundaryConds(BCtype, BCval);

  // ~~~~~ NODAL PRESSURES ~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if 0
  InnerProdVector press(nnp, 50.0);
#else
  InnerProdVector press(nnp, 0.0);
  double maxP = 10000000;
  for (int i=0; i<nnp; i++)
    press(i) = (maxP*(x_coord(i)-length)*x_coord(i)/25)+maxP;
#endif
  csm.set_press(press);

// =====================================================================
// VALIDATING THE PARTIAL DERIVATIVES
// =====================================================================

#if 1
  //csm.InspectMesh();

  InnerProdVector wrkU(3*nnp,1.0), wrkP(nnp,1.0);
  InnerProdVector outU(nnp,0.0), outP(3*nnp,0.0);
  InnerProdVector vU(nnp,0.0), vP(3*nnp,0.0);
  double delta = 1.e-3;

  // Perform (dS/dp) product with the built-in routine
  csm.Calc_dSdp_Product(wrkP, vP);
  
  // Calculate (dS/dp) with finite differencing
  csm.CalcResidual();
  InnerProdVector res1 = csm.get_res();
  vector< vector<double> > dSdp(3*nnp, vector<double>(nnp));
  for (int i=0; i<nnp; i++) {
    // preturb the pressure on i-th node
    press(i) += delta;
    csm.set_press(press);
    csm.CalcResidual();
    InnerProdVector res2 = csm.get_res();
    for (int j=0; j<3*nnp; j++) {
      dSdp[j][i] = (res2(j) - res1(j))/delta;
    }
    // reset the perturbation before testing next node
    press(i) -= delta;
  }

  // Multiply (dS/dP)*wrkP
  for (int i=0; i<3*nnp; i++) {
    //printf("|");
    for (int j=0; j<nnp; j++) {
      outP(i) += dSdp[i][j] * wrkP(j);
      //printf(" %f ", dSdp[i][j]);
    }
    //printf("|\n");
  }

  // Validate by seeing if the difference is zero
  // printf("(dS/dp) product:\n");
  // printf("---------------------\n");
  InnerProdVector diffP(3*nnp, 0.0);
  for (int i=0; i<3*nnp; i++) {
    diffP(i) = vP(i) - outP(i);
    // printf("%i th difference: %f\n", i, diffP(i));
  }
  printf("normalized L2 norm of the dS/dp product difference: %E\n",
         diffP.Norm2()/vP.Norm2());

  // Perform (dA/du) product with the built-in routine
  csm.Calc_dAdu_Product(wrkU, vU);

  // Calculate (dA/du) with finite differencing
  csm.CalcCoordsAndArea();
  InnerProdVector area1 = csm.get_area();
  vector< vector<double> > dAdu(nnp, vector<double>(3*nnp));
  InnerProdVector u(3*nnp, 0.0);
  for (int i=0; i<3*nnp; i++) {
    // perturb displacement in the i-th degree of freedom
    u(i) += delta;
    csm.set_u(u);
    csm.CalcCoordsAndArea();
    InnerProdVector area2 = csm.get_area();
    for (int j=0; j<nnp; j++) {
      dAdu[j][i] = (area2(j) - area1(j))/delta;
    }
    // undo the displacement before testing the next degree of freedom
    u(i) = 0;
  }

  // Multiply (dA/du)*wrkU
  for (int i=0; i<nnp; i++) {
    for (int j=0; j<3*nnp; j++) {
      outU(i) += dAdu[i][j] * wrkU(j);
    }
  }

  // Validate by seeing if the difference is zero
  // printf("(dA/du) product:\n");
  // printf("---------------------\n");
  InnerProdVector diffU(nnp, 0.0);
  for (int i=0; i<nnp; i++) {
    diffU(i) = vU(i) - outU(i);
    // printf("%i th difference: %f\n", i, diffU(i));
  }
  printf("L2 norm of the dA/du product difference: %E\n", diffU.Norm2());
#else
  // Call FEA solver
  csm.Solve();
  InnerProdVector u = csm.get_u();
#endif

// =====================================================================
// INSPECT THE RESIDUAL
// =====================================================================  

#if 0
  csm.CalcResidual();
  InnerProdVector & res = csm.get_res();

  printf("Printing residual for inspection:\n");
  for(int i=0; i<3*nnp; i++) {
    printf("|  %f  |\n", res(i));
  }
#endif

// =====================================================================
// GRID REFINEMENT TEST
// =====================================================================

#if 0
  // Declare the solver
  int nnp_fine = 41;
  LECSM csm_fine(nnp_fine);

  csm_fine.set_material(E, t, w, h);

  // Create problem mesh
  InnerProdVector x_fine(nnp, 0.0);
  InnerProdVector y_fine(nnp, 0.0);
  for (int i = 0; i < nnp_fine; i++) {
    // evenly spaced nodes along the x
    x_fine(i) = i*length/(nnp_fine-1);
    // parabolic nozzle wall for y coords
    y_fine(i) = 0.1*(length-x_fine(i))*x_fine(i);
  }
  csm_fine.GenerateMesh(x_fine, y_fine);

  InnerProdVector BCtype_fine(3*nnp_fine, 0.0);
  InnerProdVector BCval_fine(3*nnp_fine, 0.0);
  for (int i=0; i<nnp_fine; i++) {
    BCtype_fine(3*i) = 0;
    BCtype_fine(3*i+1) = -1;
    BCtype_fine(3*i+2) = -1;
    BCval_fine(3*i) = 0;
    BCval_fine(3*i+1) = 0;
    BCval_fine(3*i+2) = 0;
  }
  BCtype_fine(0) = 0;
  BCtype_fine(1) = 0;
  BCtype_fine(2) = 0;
  BCtype_fine(3*nnp_fine-3) = 0;
  BCtype_fine(3*nnp_fine-2) = 0;
  BCtype_fine(3*nnp_fine-1) = 0;
  csm_fine.SetBoundaryConds(BCtype_fine, BCval_fine);

  InnerProdVector press_fine(nnp_fine, 0.0);
  for (int i=0; i<nnp; i++)
    press_fine(i) = (maxP*(x_fine(i)-length)*x_fine(i)/25)+maxP;
  csm_fine.set_press(press_fine);

  csm_fine.Solve();
  InnerProdVector u_fine = csm_fine.get_u();
  double error = u_fine.Norm2() - u.Norm2();
  printf("Coarse grid deformation L2 norm: %E\n", u.Norm2());
  printf("Fine grid deformation L2 norm: %E\n", u_fine.Norm2());
  printf("Grid size error: %E\n", error);
#endif

	return 0;
}
