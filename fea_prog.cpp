/**
 * \file fea_prog.cpp
 * \FEA program driver function
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

#include <stdio.h>
#include "./fea_prog.hpp"
#include "./1D_mesh_tools.hpp"
#include "./setup_eq.hpp"
#include "./matrix_tools.hpp"
#include "./output_tools.hpp"
using namespace std;

// =====================================================================

void FEA(Mesh nozzle, double* props, double* P)
{
  // Core Finite Element Analaysis procedure.
  //
  // Inputs:
  //    mesh              - pMeshMdl mesh object
  //    geom              - pGeomMdl geom object
  //    D                 - material properties matrix
  //    ntyp              - vector indicating the type of DoF for each
  //                        component at the node
  //    FG                - vector of nodal point forces or nodal
  //                        prescribed essential boundary conditions
  //    f_init            - body forces

  // Process the problem mesh.
  int nnp = nozzle.nnp;
  int nel = nozzle.nel;

  // Initialize the problem equation.
  printf("Setting up global equation mapping...\n");
  vector< vector< vector<double> > > id(3, vector< vector<double> >(nnp, vector<double>(2)));
  vector<double> G;
  vector<double> F;
  int ndof, ndog;
  setup_eq(nozzle, id, G, F, ndof, ndog);
  printf("DONE\n");
  printf("Allocating global stiffness matrix...\n");
  vector< vector<double> > K(ndof, vector<double>(ndof));
  for (int a = 0; a < ndof; a++)
  {
    for (int b = 0; b < ndof; b++)
      {K[a][b]=0;} // Zero out the global stiffness matrix
  }
  printf("DONE\n");

  // Recover material/geometric properties
  double E = props[0];
  double w = props[1];
  double t = props[2];

  // Loop over all elements, assuming that each face is an element.
  printf("Starting element iteration...\n");
  Element elem;
  for (int i=0; i<nel; i++)
  {
    // Get information about the element.
    elem = nozzle.allElems[i];
    int nen = elem.nen;
    int nee = nen*3;
    printf("Constructing the element stiffness matrix and force vector...\n");
    vector< vector< vector<double> > > lm(3, vector< vector<double> >(nen, vector<double>(2)));
    vector< vector<double> > KE(nee, vector<double>(nee));
    vector<double> FE(nee);
    double locP = P[i];
    elem.GetStiff(E, w, t, locP, id, lm, KE, FE);
    printf("DONE\n");
    
    // Assemble the element contributions into the global matrices.
    printf("Assembling the element contributions into global matrices...\n");
    elem.assemble(lm, KE, FE, G, K, F);
    printf("DONE\n");

    // Clear vectors
    lm.clear();
    KE.clear();
    FE.clear();
  }// finish looping over elements
  printf("Element iteration complete!\n");

  // Print matrices for inspection.
  printf("    Global Stiffness Matrix K:\n");
  printMatrix(K, ndof, ndof);
  printf("    Global Force Vector F:\n");
  for (int n = 0; n < ndof; n++)
    {printf("|  %f  |\n", F[n]);}

  // Solve the global Kd = F system.
  vector<double> disp(ndof);
  int maxIt = 100000;
  printf("Starting CG Solver...\n");
  CGSolve(K, ndof, ndof, F, ndof, maxIt, disp);
  printf("DONE\n");
  printf("    Global Displacement Vector (disp):\n");
  for (int n = 0; n < ndof; n++)
    {printf("|  %f  |\n", disp[n]);}

  // Print out the node displacements.
  vector< vector<double> > nodeDisp(nnp, vector<double>(2));
  printf("Outputting displacements...\n");
  output_disp(nnp, G, id, disp, nodeDisp);

  printf("SUCCESS: Finite Element Analysis complete!\n");
}