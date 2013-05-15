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

void LECSM::CalcStiffProd(InnerProdVector& u_csm, InnerProdVector& v_csm)
{
  // Initialize the problem equation
  int nnp = geometry.nnp;
  vector< vector< vector<int> > > gm(3, vector< vector<int> >(nnp, vector<int>(2)));
  vector<double> G;
  vector<double> F;
  geometry.SetupEq(G, F, gm);

  // Zero out the global stiffness matrix
  int ndof = geometry.ndof;
  int ndog = geometry.ndog;
  vector< vector<double> > K(ndof, vector<double>(ndof));
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

  }// finish looping over elements

  // Clean-up
  G.clear();
  F.clear();
  gm.clear();
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

  // Initialize the problem equation
  int nnp = geometry.nnp;
  vector< vector< vector<int> > > gm(3, vector< vector<int> >(nnp, vector<int>(2)));
  vector<double> G;
  vector<double> F;
  geometry.SetupEq(G, F, gm);

  // Zero out the global stiffness matrix
  int ndof = geometry.ndof;
  int ndog = geometry.ndog;
  vector< vector<double> > K(ndof, vector<double>(ndof));
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

  }// finish looping over elements

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