/**
 * \file fea_prog.cpp
 * \FEA program driver function
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

#include <stdio.h>
#include "./lecsm.hpp"
#include "./matrix_tools.hpp"
#include "./output_tools.hpp"
using namespace std;

// =====================================================================

void LECSM::GetStiff()
{
  // Initialize the problem equation
  int nnp = geometry.nnp;
  vector< vector< vector<int> > > gm(3, vector< vector<int> >(nnp, vector<int>(2)));
  geometry.SetupEq(G, F, gm);

  // Pre-allocate out the global stiffness matrix
  for (int a = 0; a < ndof; a++)
  {
    for (int b = 0; b < ndof; b++)
      {K[a][b]=0;}
  }

  // Loop over all elements, assuming that each face is an element.
  Element elem;
  int nel = geometry.nel;
  for (int i=0; i<nel; i++)
  {
    // Get information about the element.
    elem = geometry.allElems[i];
    int nen = elem.nen;
    vector< vector<double> > KE(nen*3, vector<double>(nen*3));
    vector<double> FE(nen*3);
    vector< vector< vector<int> > > lm(3, vector< vector<int> >(nen, vector<int>(2)));
    double locP = elem.pressure;
    elem.GetElemStiff(E, w, t, locP, gm, lm, KE, FE);
    
    // Assemble the element contributions into the global matrices.
    elem.Assemble(KE, FE, lm, G, F, K);

    // Clean-up
    KE.clear();
    FE.clear();
    lm.clear();
  }
}

void LECSM::CalcStiffProd(InnerProdVector& u_csm, InnerProdVector& v_csm)
{

  // Strip u_csm of fixed nodes
  Node node;
  int nnp = geometry.nnp;
  int ndof = geometry.ndof;
  vector<double> u_dof(ndof);
  int p = 0;
  for (int i=0; i<nnp; i++) {
    for (int j=0; j<3; j++) {
      if (gm[j][i][0] == 1) {
        u_dof[p] = u_csm(i*j);
        p++;
      }
    }
  }

  // Perform the (dS/du)*u_dof operation
  vector<double> v_dof(ndof);
}

void LECSM::Solve()
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

  // Get the global LHS stiffness matrix
  GetStiff();

  // Recover mesh properties
  int nnp = geometry.nnp;
  int ndof = geometry.ndof;

  // Solve the global Kd = F system.
  vector<double> disp(ndof);
  int maxIt = 100000;
  CGSolve(K, ndof, ndof, F, ndof, maxIt, disp);
  for (int n = 0; n < ndof; n++)
    {printf("|  %f  |\n", disp[n]);}

  // Print out the node displacements.
  vector< vector<double> > nodeDisp(nnp, vector<double>(2));
  printf("Outputting displacements...\n");
  output_disp(nnp, G, gm, disp, nodeDisp);
}