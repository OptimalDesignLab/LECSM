/**
 * \file fea_prog.cpp
 * \FEA program driver function
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

#include <stdio.h>
#include <math.h>
#include "./lecsm.hpp"
#include "./matrix_tools.hpp"
#include "./output_tools.hpp"
using namespace std;

// =====================================================================

void LECSM::InitGlobalVecs(vector<double>& G, vector<double>& F,
                           vector< vector<double> >& K)
{
  Node nd;
  int nnp = geometry.nnp;
  for (int b = 0; b < nnp; b++)
  {
    nd = geometry.allNodes[b];
    for (int i = 0; i < 3; i++)
    {
      if (nd.type[i]==1)        // DoF - possible nodal load
      {
        F.push_back(nd.forceBC[i]);     // store the nodal load
      }
      else                      // prescribed BC
      {
        if (nd.type[i]==2)      // DoG - non-zero BC
        {
          G.push_back(nd.dispBC[i]);   // store the essential BC
        }
      }
    }
  }
  int ndof = geometry.ndof;
  for (int a = 0; a < ndof; a++)
  {
    for (int b = 0; b < ndof; b++)
      {K[a][b]=0;}
  }
}

void LECSM::GetStiff(vector< vector< vector<int> > >& gm,
                     vector<double>& G, vector<double>& F,
                     vector< vector<double> >& K)
{
  // Loop over all elements, assuming that each face is an element.
  Element elem;
  int nel = geometry.nel;
  for (int i=0; i<nel; i++)
  {
    // Get information about the element.
    elem = geometry.allElems[i];
    int nen = elem.nen;

    // Calculate the stiffness matrix
    vector< vector< vector<int> > > lm(3, vector< vector<int> >(nen, vector<int>(2)));
    vector< vector<double> > KE(nen*3, vector<double>(nen*3));
    vector<double> FE(nen*3);
    elem.GetElemStiff(E, w, t, P, gm, lm, KE, FE);
    
    // Assemble the element contributions into the global matrices.
    elem.Assemble(KE, FE, lm, G, F, K);
  }
}

void LECSM::Calc_dSdu_Product(InnerProdVector& u_csm, InnerProdVector& v_csm)
{
  // Map out the global equation numbers
  int nnp = geometry.nnp;
  vector< vector< vector<int> > > gm(3, vector< vector<int> >(nnp, vector<int>(2)));
  geometry.SetupEq(gm);

  // Calculate the stiffness matrix (dS/du)
  int ndof = geometry.ndof;
  int ndog = geometry.ndog;
  vector<double> G(ndog), F(ndof);
  vector< vector<double> > K(ndof, vector<double>(ndof));
  InitGlobalVecs(G, F, K);
  GetStiff(gm, G, F, K);
  G.clear();
  F.clear();

  // Strip u_csm of fixed nodes
  vector<double> u_dof(ndof);
  int p = 0;
  for (int i=0; i<nnp; i++) {
    for (int j=0; j<3; j++) {
      if (gm[j][i][0] == 1) {
        u_dof[p] = u_csm(i*3+j);
        p++;
      }
    }
  }

  // Perform the (dS/du)*u_dof operation
  vector<double> v_dof(ndof);
  matrixVecMult(K, ndof, ndof, u_dof, ndof, v_dof);
  u_dof.clear();

  // Insert results into v_csm
  int a, b;
  for (int i=0; i<nnp; i++) {
    for (int j=0; j<3; j++) {
      a = (i*3)+j;
      if (gm[j][i][0] == 0) {       // node dof is fixed (zero)
        v_csm(a) = 0;
      }
      else if (gm[j][i][0] == 1) {  // node dof is free
        b = gm[j][i][1];
        v_csm(a) = v_dof[b];
      }
      else if (gm[j][i][0] == 2) {  // node dof is fixed (prescribed)

      }
    }
  }

  // Clean-up
  gm.clear();
  K.clear();
  v_dof.clear();

}

void LECSM::Calc_dAdu_Product(InnerProdVector& u_csm, InnerProdVector& wrk)
{
  // Calculate dA/du
  int nnp = geometry.nnp;
  int p = 0;
  vector< vector<double> > dAdu(nnp, vector<double>(nnp*3));
  for (int i=0; i<nnp; i++) {
    dAdu[i][i*p] = 0;
    dAdu[i][i*p+1] = -w;
    dAdu[i][i*p+2] = 0;
    wrk(i) = 0;     // initialize the output vector
  }

  // Calculate wrk = (dA/du)*u_csm
  for (int i=0; i<nnp; i++) {
    for (int j=0; j<nnp*3; j++) {
      wrk(i) += dAdu[i][j] * u_csm(j);
    }
  }
}

void LECSM::Calc_dSdp_Product(InnerProdVector& wrk, InnerProdVector& u_cfd)
{
  // Initialize the global derivative matrix
  int nnp = geometry.nnp;
  vector< vector<double> > dSdu(nnp*3, vector<double>(nnp));
  for (int i=0; i<nnp*3; i++) {
    for (int j=0; j<nnp; j++) {
      dSdu[i][j] = 0;
    }
  }

  // Loop over elements, calculating dS/dp at the element level before
  // adding the contributions into u_cfd
  int nel = geometry.nel;
  Element elem;
  Node nodeL, nodeR;
  int idL, idR;
  double len, c, s, dFxdp, dFydp;
  vector< vector<double> > dSdu_ele(6, vector<double>(2));
  for (int i=0; i<nel; i++) {
    // Initialize element parameters
    elem = geometry.allElems[i];
    c = elem.cosine;
    s = elem.sine;
    len = elem.length;

    // Calculate the element node contribution
    dFxdp = -w*len*c/4;
    dFydp = -w*len*s/4;

    // Add the element node contributions to the global derivative
    nodeL = elem.adjNodes[0];
    nodeR = elem.adjNodes[1];
    idL = nodeL.id;
    idR = nodeR.id;

    dSdu[3*idL][idL] += dFxdp;
    dSdu[3*idL+1][idL] += dFydp;
    // dSdu[3*idL+2][idL] = 0 // Moment term (pressure independent)
    dSdu[3*idL][idR] += dFxdp;
    dSdu[3*idL+1][idR] += dFydp;
    // dSdu[3*idL+2][idR] = 0 // Moment term (pressure independent)
    dSdu[3*idR][idL] += dFxdp;
    dSdu[3*idR+1][idR] += dFydp;
    // dSdu[3*idR+2][idR] = 0 // Moment term (pressure independent)
    dSdu[3*idR][idL] += dFxdp;
    dSdu[3*idR+1][idL] += dFydp;
    // dSdu[3*idR+2][idL] = 0 // Moment term (pressure independent)
  }

  // Initialize the output vector
  for (int i=0; i<nnp*3; i++) {
    u_cfd(i) = 0;
  }

  // Calculate (dS/du)*wrk
  for (int i=0; i<nnp*3; i++) {
    for (int j=0; j<nnp; j++) {
      u_cfd(i) += dSdu[i][j] * wrk(j);
    }
  }

  // Clean-up
  dSdu.clear();
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

  int nnp = geometry.nnp;
  vector< vector< vector<int> > > gm(3, vector< vector<int> >(nnp, vector<int>(2)));
  geometry.SetupEq(gm);

  int ndof = geometry.ndof;
  int ndog = geometry.ndog;
  vector<double> G(ndog), F(ndof);
  vector< vector<double> > K(ndof, vector<double>(ndof));
  InitGlobalVecs(G, F, K);

  // Get the global LHS stiffness matrix
  GetStiff(gm, G, F, K);

  // Solve the global Kd = F system
  vector<double> disp(ndof);
  int maxIt = 100000;
  CGSolve(K, ndof, ndof, F, ndof, maxIt, disp);

  // Print out the node displacements.
  vector< vector<double> > nodeDisp(nnp, vector<double>(2));
  printf("Outputting displacements...\n");
  output_disp(nnp, G, gm, disp, nodeDisp);
}