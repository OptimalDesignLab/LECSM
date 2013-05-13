/**
 * \file fea_prog.cpp
 * \FEA program driver function
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

#include <stdio.h>
#include "./fea_prog.h"
using namespace std;

// =====================================================================

void FEA(Mesh nozzle, double props, double P,
         vector< vector<double> > D,
         vector< vector<int> > ntyp,
         vector< vector<double> > FG,
         vector<double> f_init,
         double P)
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
  vector< vector< vector<double> > > id(nsd, vector< vector<double> >(nnp, vector<double>(2)));
  vector<double> G;
  vector<double> F;
  int ndof, ndog;
  setup_eq(nozzle, FG, id, G, F, ndof, ndog)
  printf("Allocating global stiffness matrix...\n");
  vector< vector<double> > K(ndof, vector<double>(ndof));
  for (int a = 0; a < ndof; a++)
  {
    for (int b = 0; b < ndof; b++)
      {K[a][b]=0;} // Zero out the global stiffness matrix
  }
  printf("DONE\n");
  printf("Governing equation initialized.\n");

  // Loop over all elements, assuming that each face is an element.
  printf("Starting element iteration...\n");
  Element elem;
  for (int i=0; i<nel; i++)
  {
    // Get the element stiffness matrix.
    elem = nozzle.allElems[i];
    printf("Element Number: %i\n", elem.id);
    int nee = nen*2;
    vector< vector<double> > KE(nee, vector<double>(nee));
    vector<double> FE(nee);
    elem.GetStiff(E, w, t, P, KE, FE);

    // Assemble the element contributions into the global matrices.
    printf("Assembling the element contributions into global matrices...\n");
    assemble(nsd, nen, lm, KE, FE, G, K, F);
    printf("DONE\n");

    // Clear vectors
    ien.clear();
    weights.clear();
    intPts.clear();
    nodeCoords.clear();
    lm.clear();
    KE.clear();
    FE.clear();
    Se.clear();
    iterEnd = FMDB_PartEntIter_GetNext(iter, elem);
  }// finish looping over elements
  FMDB_PartEntIter_Del(iter);
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
  vector<double> nodeDisp(2);
  printf("Outputting displacements...\n");
  output_disp(nnp, nsd, G, id, disp, nodeDisp);

  // Calculate and print out stresses/fluxes.
  printf("Calculating stresses at each element...\n");
  FMDB_PartEntIter_Init(part, FMDB_FACE, FMDB_ALLTOPO, iter);
  iterEnd = FMDB_PartEntIter_GetNext(iter, elem);
  while(!iterEnd)
  {
    // Get information about the element.
    elemID = FMDB_Ent_ID(elem);
    printf("Element Number: %i\n", elemID);
    int nint, nen, nee;
    Elem_CheckType(elem, nint, nen);
    nee = nen*nsd;
    vector<int> ien(nen);
    vector<double> weights(nint);
    vector< vector<double> > intPts(nint, vector<double>(nsd));
    vector< vector<double> > nodeCoords(nen, vector<double>(nsd));
    Elem_GetInfo(elem, nsd, nen, nint, nee, ien, weights, intPts, nodeCoords);
    // Calculate the element stress vector.
    vector<double> SIG(3);
    vector< vector< vector<double> > > lm(nsd, vector< vector<double> >(nen, vector<double>(2)));
    out_flux(nsd, nen, nint, ien, id, S[elemID], disp, G, lm, SIG);

    // Print out the results.
    printf("    Stresses at element %d:\n", elemID);
    for (int n = 0; n < 3; n++)
      {printf("|  %f  |\n", SIG[n]);}

    iterEnd = FMDB_PartEntIter_GetNext(iter, elem);
    ien.clear();
    weights.clear();
    intPts.clear();
    lm.clear();
    nodeCoords.clear();
    SIG.clear();
  }// finish looping over elements
  FMDB_PartEntIter_Del(iter);
  printf("DONE\n");
  printf("SUCCESS: Finite Element Analysis complete!\n");
}