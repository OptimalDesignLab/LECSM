/**
 * \file fea_prog.cpp
 * \brief FEA program driver function
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

void LECSM::GenerateMesh(const InnerProdVector & x, const InnerProdVector & y)
{
  // Create the mesh nodes
  Node node;
  vector<Node> nodes;
  for (int i=0; i<nnp_; i++) {
    double c[2] = {x(i), y(i)};
    node.CreateNode(i,c)
    nodes.push_back(node);
  }

  // Create 2D beam elements from the nodess
  Element elem;
  vector<Element> elems;
  Node nodeL;
  Node nodeR;
  vector<Node> elemNodes(2);
  int nel = nnp_ - 1;
  for (int i=0; i<nel; i++) {
    nodeL = nodes[i];
    nodeR = nodes[i+1];
    elemNodes[0] = nodeL;
    elemNodes[1] = nodeR;
    elem.CreateElem(j, elemNodes);
    elems.push_back(elem);
  }

  // Create the problem mesh
  geom_.CreateMesh(elems);

  // Clean-up
  nodes.clear();
  elems.clear();
}

// =====================================================================

void LECSM::SetBoundarConds(const InnerProdVector & BCtype, 
                            const InnerProdVector & BCval)
{
  // Loop over all mesh nodes
  Node node;
  double type[3], val[3];
  for (int i=0; i<nnp_; i++) {
    node = geom_.allNodes[i];
    for (int j=0; j<3; j++) {
      type[i] = BCtype(3*i+j);  // extract BC type for the node
      val[i] = BCval(3*i+j);  // extract BC value for the node
    }
    // define the BC
    node.DefineBCs(type, val);
  }
}

// =====================================================================

void LECSM::InitGlobalVecs(vector<double>& G, vector<double>& F,
                           vector< vector<double> >& K)
{
  Node nd;
  int nnp = geom_.nnp;
  for (int b = 0; b < nnp; b++)
  {
    nd = geom_.allNodes[b];
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
  int ndof = geom_.ndof;
  for (int a = 0; a < ndof; a++)
  {
    for (int b = 0; b < ndof; b++)
      {K[a][b]=0;}
  }
}

// =====================================================================

void LECSM::GetStiff(vector< vector< vector<int> > >& gm,
                     vector<double>& G, vector<double>& F,
                     vector< vector<double> >& K)
{
  // Loop over all elements, assuming that each face is an element.
  Element elem;
  int nel = geom_.nel;
  Node nodeL, nodeR;
  for (int i=0; i<nel; i++)
  {
    // Get information about the element.
    elem = geom_.allElems[i];
    int nen = elem.nen;

    // Calculate the stiffness matrix
    vector< vector< vector<int> > > lm(3, vector< vector<int> >(nen, vector<int>(2)));
    vector< vector<double> > KE(nen*3, vector<double>(nen*3));
    vector<double> FE(nen*3);
    nodeL = elem.adjNodes[0];
    nodeR = elem.adjNodes[1];
    vector<double> locP(2);
    locP[0] = P_(nodeL.id);
    locP[1] = P_(nodeR.id);
    elem.GetElemStiff(E, w, t, locP, gm, lm, KE, FE);
    
    // Assemble the element contributions into the global matrices.
    elem.Assemble(KE, FE, lm, G, F, K);
  }
}

// =====================================================================

void LECSM::Calc_dSdu_Product(InnerProdVector& u_csm, InnerProdVector& v_csm)
{
  // Map out the global equation numbers
  int nnp = geom_.nnp;
  vector< vector< vector<int> > > gm(3, vector< vector<int> >(nnp, vector<int>(2)));
  geom_.SetupEq(gm);

  // Calculate the stiffness matrix (dS/du)
  int ndof = geom_.ndof;
  int ndog = geom_.ndog;
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

// =====================================================================

void LECSM::Calc_dAdu_Product(InnerProdVector& u_csm, InnerProdVector& wrk)
{
  // Calculate dA/du
  int nnp = geom_.nnp;
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

// =====================================================================

void LECSM::Calc_dSdp_Product(InnerProdVector& wrk, InnerProdVector& u_cfd)
{
  // Initialize the global derivative matrix
  int nnp = geom_.nnp;
  vector< vector<double> > dSdu(nnp*3, vector<double>(nnp));
  for (int i=0; i<nnp*3; i++) {
    for (int j=0; j<nnp; j++) {
      dSdu[i][j] = 0;
    }
  }

  // Loop over elements, calculating dS/dp at the element level before
  // adding the contributions into u_cfd
  int nel = geom_.nel;
  Element elem;
  Node nodeL, nodeR;
  int idL, idR;
  double len, c, s, dFxdp, dFydp;
  vector< vector<double> > dSdu_ele(6, vector<double>(2));
  for (int i=0; i<nel; i++) {
    // Initialize element parameters
    elem = geom_.allElems[i];
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

// =====================================================================

void LECSM::CalcStateVars()
{
  int nnp = geom_.nnp;
  double y, h;
  for (int i=0; i<nnp; i++) {
    xCoords(i) = u_(3*i);
    y = u_(3*i+1);
    h = 2*(0.5 - y);
    area(i) = w*h;
  }
}

// =====================================================================

void LECSM::CalcResidual()
{
  // Update the mesh from the latest displacements
  UpdateMesh(u_);

  // Generate the global equation number mapping
  int nnp = geom_.nnp;
  vector< vector< vector<int> > > gm(3, vector< vector<int> >(nnp, vector<int>(2)));
  geom_.SetupEq(gm);

  // Initiate global vectors used in the solver
  int ndof = geom_.ndof;
  int ndog = geom_.ndog;
  vector<double> G(ndog), F(ndof);
  vector< vector<double> > K(ndof, vector<double>(ndof));
  InitGlobalVecs(G, F, K);

  // Calculate the stiffness matrix and the forcing vector
  GetStiff(gm, G, F, K);

  int p;
  vector<double> u_dof(ndof);
  for (int i=0; i<nnp; i++) {
    for (int j=0; j<3; j++) {
      if (gm[j][i][0] != 1) {
        p = gm[j][i][1];
        u_dof[p] = u_(3*i+j);
      }
    }
  }

  // Calculate the K*u product
  vector<double> v_dof(ndof);
  matrixVecMult(K, ndof, ndof, u_dof, ndof, v_dof);
  u_dof.clear();

  // Form the Ku-f residual for free nodes
  vector<double> res_dof(ndof);
  for (int i=0; i<ndof; i++)
    res_dof[i] = v_dof[i] - F[i];

  // Assemble the whole residual
  for (int i=0; i<nnp; i++) {
    for (int j=0; j<3; j++) {
      if (gm[j][i][0] == 1)   // node is free
        res_(3*i+j) = res[gm[j][i][1]];
      else  // node is fixed
        res_(3*i+j) = 0;
    }
  }

}

// =====================================================================

void LECSM::Solve()
{
  // Generate the global equation number mapping
  int nnp = geom_.nnp;
  vector< vector< vector<int> > > gm(3, vector< vector<int> >(nnp, vector<int>(2)));
  geom_.SetupEq(gm);

  // Initiate global vectors used in the solver
  int ndof = geom_.ndof;
  int ndog = geom_.ndog;
  vector<double> G(ndog), F(ndof);
  vector< vector<double> > K(ndof, vector<double>(ndof));
  InitGlobalVecs(G, F, K);

  // Calculate the stiffness matrix and the forcing vector
  GetStiff(gm, G, F, K);

  // Solve the global Kd = F system
  vector<double> disp(ndof);
  int maxIt = 100000;
  CGSolve(K, ndof, ndof, F, ndof, maxIt, disp);

  // Assemble the nodal displacements
  vector< vector<double> > du(nnp, vector<double>(2));
  printf("Outputting displacements...\n");
  output_disp(nnp, G, gm, disp, du);
}