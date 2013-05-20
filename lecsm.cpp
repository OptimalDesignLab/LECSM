/**
 * \file lecsm.cpp
 * \brief Linear Elastic CSM Solver (2D beam)
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

#include <stdio.h>
#include <math.h>
#include "./lecsm.hpp"
#include "./matrix_tools.hpp"
#include "./output_tools.hpp"

using namespace std;

void LECSM::set_u(const InnerProdVector & u_csm)
{ 
  u_ = u_csm;
  geom_.Update(u_);
  for (int i=0; i<nnp_; i++) {
    Node node = geom_.allNodes[i];
    xCoords_(i) = node.coords[0];
    yCoords_(i) = node.coords[1];
  }
}

// =====================================================================

void LECSM::GenerateMesh(const InnerProdVector & x, const InnerProdVector & y)
{
  // Save coordinates for later
  xCoords_ = x;
  yCoords_ = y;

  // Create the mesh nodes
  vector<Node> nodes;
  for (int i=0; i<nnp_; i++) {
    double c[2] = {x(i), y(i)};
    Node node(i, c);
    nodes.push_back(node);
  }

  // Create 2D beam elements from the nodess
  int nel = nnp_ - 1;
  vector<Element> elems(nel);
  vector<Node> elemNodes(2);
  for (int i=0; i<nel; i++) {
    Node nodeL = nodes[i];
    Node nodeR = nodes[i+1];
    elemNodes[0] = nodeL;
    elemNodes[1] = nodeR;
    Element elem(i, elemNodes);
    elems[i] = (elem);
  }

  // Create the problem mesh
  geom_.CreateMesh(elems, nodes);

  // Clean-up
  nodes.clear();
  elems.clear();
}

// =====================================================================

void LECSM::SetBoundaryConds(const InnerProdVector & BCtype, 
                             const InnerProdVector & BCval)
{
  // Loop over all mesh nodes
  int type[3];
  double val[3];
  for (int i=0; i<nnp_; i++) {;
    for (int j=0; j<3; j++) {
      type[j] = BCtype(3*i+j);  // extract BC type for the node
      val[j] = BCval(3*i+j);  // extract BC value for the node
    }
    // define the BC
    geom_.allNodes[i].DefineBCs(type, val);
  }

  // Cascade the nodal changes into their elements
  for (int i=0; i<geom_.nel; i++) {
    geom_.allElems[i].adjNodes[0] = geom_.allNodes[i];
    geom_.allElems[i].adjNodes[1] = geom_.allNodes[i+1];
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
    elem.GetElemStiff(E_, w_, t_, locP, gm, lm, KE, FE);
    
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
  int nnp = geom_.nnp;
  vector< vector<double> > dAdu(nnp, vector<double>(3*nnp));
  for (int i=0; i<nnp; i++) {
    Node node = geom_.allNodes[i];
    dAdu[i][3*i] = 0;
    if (node.type[1] == 1)
      dAdu[i][3*i+1] = -2*w_;
    else
      dAdu[i][3*i+1] = 0;
    dAdu[i][3*i+2] = 0;
  }

  for (int i=0; i<nnp; i++) {
    for (int j=0; j<3*nnp; j++) {
      wrk(i) += dAdu[i][j] * u_csm(j);
    }
  }
}

// =====================================================================

void LECSM::Calc_dSdp_Product(InnerProdVector& wrk, InnerProdVector& u_cfd)
{
  // Initialize the global derivative matrix
  int nnp = geom_.nnp;
  vector< vector<double> > dSdp(nnp*3, vector<double>(nnp));
  for (int i=0; i<nnp*3; i++) {
    for (int j=0; j<nnp; j++) {
      dSdp[i][j] = 0;
    }
  }

  // Loop over elements, calculating dS/dp at the element level before
  // adding the contributions into u_cfd
  int nel = geom_.nel;
  Element elem;
  Node nodeL, nodeR;
  int idE, idL, idR, type[3];
  double x1, x2, y1, y2, len, c, s, dFxdp, dFydp;
  vector< vector<double> > dSdp_elem(nnp*3, vector<double>(nnp));
  for (int i=0; i<nel; i++) {
    // Initialize element parameters
    elem = geom_.allElems[i];
    idE = elem.id;
    nodeL = elem.adjNodes[0];
    nodeR = elem.adjNodes[1];
    idL = nodeL.id;
    idR = nodeR.id;
    vector< vector<double> > dSdp_elem(6, vector<double>(2));

    // Calculate element length and orientation
    x1 = nodeL.coords[0];
    x2 = nodeR.coords[0];
    y1 = nodeL.coords[1];
    y2 = nodeR.coords[1];
    len = sqrt(pow(x2-x1,2)+pow(y2-y1,2));
    c = (x2 - x1)/len;
    s = (y2 - y1)/len;

    // Calculate the element node contribution
    dFxdp = -w_*len*s/4;
    dFydp = w_*len*c/4;
    dSdp_elem[0][0] = dFxdp;
    dSdp_elem[1][0] = dFydp;
    dSdp_elem[2][0] = 0; // Moment term (pressure independent)
    dSdp_elem[0][1] = dFxdp;
    dSdp_elem[1][1] = dFydp;
    dSdp_elem[2][1] = 0; // Moment term (pressure independent)

    for (int k=0; k<3; k++) {
      if (nodeL.type[k] == 1) {
        dSdp[3*idL+k][idE] += dSdp_elem[k][0];
        dSdp[3*idL+k][idE+1] += dSdp_elem[k][1];
      }
      if (nodeR.type[k] == 1) {
        dSdp[3*idR+k][idE] += dSdp_elem[k][0];
        dSdp[3*idR+k][idE+1] += dSdp_elem[k][1];
      }
    }

    dSdp_elem.clear();
  }

  // Calculate (dS/du)*wrk
  for (int i=0; i<nnp*3; i++) {
    for (int j=0; j<nnp; j++) {
      u_cfd(i) += dSdp[i][j] * wrk(j); // perform the multiplication
    }
  }

  // Clean-up
  dSdp.clear();
}

// =====================================================================

void LECSM::CalcArea()
{
  int nnp = geom_.nnp;
  double y, meshH;
  for (int i=0; i<nnp; i++) {
    y = yCoords_(i);
    meshH = 2*(0.5*h_ - y);
    area_(i) = w_*meshH;
  }
}

// =====================================================================

void LECSM::CalcResidual()
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
  G.clear();

  int p;
  vector<double> u_dof(ndof);
  for (int i=0; i<nnp; i++) {
    for (int j=0; j<3; j++) {
      if (gm[j][i][0] == 1) {
        p = gm[j][i][1];
        u_dof[p] = u_(3*i+j);
      }
    }
  }

  // Calculate the K*u product
  vector<double> v_dof(ndof);
  matrixVecMult(K, ndof, ndof, u_dof, ndof, v_dof);
  u_dof.clear();
  K.clear();

  // Form the Ku-f residual for free nodes
  vector<double> res_dof(ndof);
  for (int i=0; i<ndof; i++)
    res_dof[i] = v_dof[i] - F[i];
  v_dof.clear();
  F.clear();

  // Assemble the whole residual
  for (int i=0; i<nnp; i++) {
    for (int j=0; j<3; j++) {
      if (gm[j][i][0] == 1)   // node is free
        res_(3*i+j) = res_dof[gm[j][i][1]];
      else  // node is fixed
        res_(3*i+j) = 0;
    }
  }
  res_dof.clear();
}

// =====================================================================

void LECSM::Solve()
{
#if 1
  geom_.InspectNodes();
#endif

  // Generate the global equation number mapping
  int nnp = geom_.nnp;
  vector< vector< vector<int> > > gm(3, vector< vector<int> >(nnp, vector<int>(2)));
  geom_.SetupEq(gm);
#if 1
  printf("Printing global mapping for inspection:\n");
  for (int i=0; i<nnp; i++) {
    for (int j=0; j<3; j++) {
      printf("    gm[%i][%i][%i] = %i\n", j, i, 0, gm[j][i][0]);
      printf("    gm[%i][%i][%i] = %d\n", j, i, 1, gm[j][i][1]);
    }
  }
#endif

  // Initiate global vectors used in the solver
  int ndof = geom_.ndof;
  int ndog = geom_.ndog;
  vector<double> G(ndog), F(ndof);
  vector< vector<double> > K(ndof, vector<double>(ndof));
  InitGlobalVecs(G, F, K);

  // Calculate the stiffness matrix and the forcing vector
  GetStiff(gm, G, F, K);
#if 1
  printf("Printing the stiffness matrix for inspection:\n");
  printMatrix(K, ndof, ndof);
  printf("Printing the forcing vector for inspection:\n");
  for (int i=0; i<ndof; i++) {
    printf("|   %f   |\n", F[i]);
  }
#endif

  // Solve the global Kd = F system
  vector<double> disp(ndof);
  int maxIt = 100000;
  CGSolve(K, ndof, ndof, F, ndof, maxIt, disp);

  // Assemble the nodal displacements
  vector< vector<double> > du(nnp, vector<double>(2));
  printf("Outputting displacements...\n");
  output_disp(nnp, G, gm, disp, du);
}