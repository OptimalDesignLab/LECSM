/**
 * \file lecsm.cpp
 * \brief Linear Elastic CSM Solver (2D beam)
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

#include <stdio.h>
#include <math.h>
#include <fstream>
#include <string>
#include "./lecsm.hpp"
#include "./matrix_tools.hpp"

using namespace std;

void LECSM::set_u(const InnerProdVector & u_csm)
{ 
  u_ = u_csm;
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

void LECSM::ResetCoords()
{
  int nnp = geom_.nnp;
  for (int i=0; i<nnp; i++) {
    Node node = geom_.allNodes[i];
    xCoords_(i) = node.coords[0];
    yCoords_(i) = node.coords[1];
  }
}

// =====================================================================

void LECSM::UpdateMesh()
{
  int nnp = geom_.nnp;
  for (int i=0; i<nnp; i++) {
    geom_.allNodes[i].coords[0] = xCoords_(i);
    geom_.allNodes[i].coords[1] = yCoords_(i);
  }
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

void LECSM::InitGlobalVecs(vector<double>& G, vector<double>& F)
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

void LECSM::Precondition(InnerProdVector& in, InnerProdVector& out)
{
  // Map out the global equation numbers
  int nnp = geom_.nnp;
  vector< vector< vector<int> > > gm(3, vector< vector<int> >(nnp, vector<int>(2)));
  geom_.SetupEq(gm);

  // Calculate the stiffness matrix (dS/du)
  int ndof = geom_.ndof;
  int ndog = geom_.ndog;
  vector<double> G(ndog), F(ndof);
  vector< vector<double> > K(ndof, vector<double>(ndof, 0.0));
  InitGlobalVecs(G, F);
  GetStiff(gm, G, F, K);
  G.clear();
  F.clear();

  for (int i=0; i<nnp; i++) {
    Node node = geom_.allNodes[i];
    for (int j=0; j<3; j++) {
      if (node.type[j] == 1) {
        int p = gm[j][node.id][1];
        out(3*i+j) = in(3*i+j)/K[p][p];
      }
      else
        out(3*i+j) = in(3*i+j);
    }
  }
  K.clear();

}

// =====================================================================

void LECSM::Calc_dSdu_Product(const InnerProdVector& in, InnerProdVector& out)
{
  // Generate the global equation number mapping
  int nnp = geom_.nnp;
  vector< vector< vector<int> > > gm(3, vector< vector<int> >(nnp, vector<int>(2)));
  geom_.SetupEq(gm);

  // Initiate global vectors used in the solver
  int ndof = geom_.ndof;
  int ndog = geom_.ndog;
  vector<double> G(ndog), F(ndof);
  vector< vector<double> > K(ndof, vector<double>(ndof, 0.0));
  InitGlobalVecs(G, F);

  // Calculate the stiffness matrix and the forcing vector
  GetStiff(gm, G, F, K);
  G.clear();
  F.clear();

  int p;
  vector<double> u_dof(ndof);
  for (int i=0; i<nnp; i++) {
    for (int j=0; j<3; j++) {
      if (gm[j][i][0] == 1) {
        p = gm[j][i][1];
        u_dof[p] = in(3*i+j);
      }
    }
  }

  // Calculate the K*u product
  vector<double> v_dof(ndof);
  matrixVecMult(K, ndof, ndof, u_dof, ndof, v_dof);
  u_dof.clear();
  K.clear();

  // Assemble the whole product
  for (int i=0; i<nnp; i++) {
    for (int j=0; j<3; j++) {
      if (gm[j][i][0] == 1)   // node is free
        out(3*i+j) = v_dof[gm[j][i][1]];
      else  // node is fixed
        out(3*i+j) = 0;
    }
  }

  // Clean-up
  gm.clear();
  K.clear();
  v_dof.clear();
}

// =====================================================================

void LECSM::Calc_dAdu_Product(InnerProdVector& in, InnerProdVector& out)
{
  int nnp = geom_.nnp;
  vector< vector<double> > dAdu(nnp, vector<double>(3*nnp, 0.0));
  for (int i=0; i<nnp; i++) {
    Node node = geom_.allNodes[i];
    dAdu[i][3*i] = 0;
    out(i) = 0;
    if (node.type[1] == 1)
      dAdu[i][3*i+1] = -2*w_;
    else
      dAdu[i][3*i+1] = 0;
    dAdu[i][3*i+2] = 0;
  }

  for (int i=0; i<nnp; i++) {
    for (int j=0; j<3*nnp; j++) {
      out(i) += dAdu[i][j] * in(j);
    }
  }

  dAdu.clear();
}

// =====================================================================

void LECSM::CalcTrans_dAdu_Product(InnerProdVector& in, InnerProdVector& out)
{
  int nnp = geom_.nnp;
  vector< vector<double> > dAdu(nnp, vector<double>(3*nnp, 0.0));
  for (int i=0; i<nnp; i++) {
    Node node = geom_.allNodes[i];
    dAdu[i][3*i] = 0;
    if (node.type[1] == 1)
      dAdu[i][3*i+1] = -2*w_;
    else
      dAdu[i][3*i+1] = 0;
    dAdu[i][3*i+2] = 0;
  }

  for (int i=0; i<3*nnp; i++) {
    out(i) = 0;
    for (int j=0; j<nnp; j++) {
      out(i) += dAdu[j][i] * in(j);
    }
  }

  dAdu.clear();
}

// =====================================================================

void LECSM::Calc_dudA_Product(InnerProdVector& in, InnerProdVector& out)
{
  int nnp = geom_.nnp;
  vector< vector<double> > dudA(3*nnp, vector<double>(nnp, 0.0));
  for (int i=0; i<nnp; i++) {
    Node node = geom_.allNodes[i];
    if (node.type[1] == 1)
      dudA[3*i+1][i] = -1/(2*w_);
  }
  for (int i=0; i<3*nnp; i++) {
    for (int j=0; j<nnp; j++)
      out(i) += dudA[i][j] * in(j);
  }
  dudA.clear();
}

// =====================================================================

void LECSM::CalcTrans_dudA_Product(InnerProdVector& in, InnerProdVector& out)
{
  int nnp = geom_.nnp;
  vector< vector<double> > dudA(3*nnp, vector<double>(nnp, 0.0));
  for (int i=0; i<nnp; i++) {
    Node node = geom_.allNodes[i];
    if (node.type[1] == 1)
      dudA[3*i+1][i] = -1/(2*w_);
  }
  for (int i=0; i<nnp; i++) {
    for (int j=0; j<3*nnp; j++)
      out(i) += dudA[j][i] * in(j);
  }
  dudA.clear();
}
// =====================================================================

void LECSM::Calc_dSdp_Product(InnerProdVector& in, InnerProdVector& out)
{
  // Initialize the global derivative matrix and zero out the resultant
  int nnp = geom_.nnp;
  vector< vector<double> > dSdp(nnp*3, vector<double>(nnp, 0.0));

  // Loop over elements, calculating dS/dp at the element level before
  // adding the contributions into u_cfd
  int nel = geom_.nel;
  Element elem;
  Node nodeL, nodeR;
  int idE, idL, idR, type[3];
  double x1, x2, y1, y2, len, c, s, dFxdp, dFydp;
  for (int i=0; i<nel; i++) {
    // Initialize element parameters
    elem = geom_.allElems[i];
    idE = elem.id;
    nodeL = elem.adjNodes[0];
    nodeR = elem.adjNodes[1];
    idL = nodeL.id;
    idR = nodeR.id;
    vector< vector<double> > dSdp_elem(6, vector<double>(2));

#if 0
  double q1 = -P[0]*w;
  double q2 = -P[1]*w;
  double f1 = (length/6)*((2*q1)+q2);
  double f2 = (length/6)*(q1+(2*q2));
  FE[0] = (-sine*f1) + nodeL.forceBC[0];
  FE[1] = (cosine*f1) + nodeL.forceBC[1];
  FE[2] = nodeL.forceBC[2];
  FE[3] = (-sine*f2) + nodeR.forceBC[0];
  FE[4] = (cosine*f2) + nodeR.forceBC[1];
  FE[5] = nodeR.forceBC[2];
#endif
  
    // Calculate element length and orientation
    x1 = nodeL.coords[0];
    x2 = nodeR.coords[0];
    y1 = nodeL.coords[1];
    y2 = nodeR.coords[1];
    len = sqrt(pow(x2-x1,2)+pow(y2-y1,2));
    c = (x2 - x1)/len;
    s = (y2 - y1)/len;

    // Calculate the element node contribution
    dSdp_elem[0][0] = len*w_*s/3;
    dSdp_elem[1][0] = -len*w_*c/3;
    dSdp_elem[2][0] = 0; // Moment term (pressure independent)
    dSdp_elem[3][0] = len*w_*s/6;
    dSdp_elem[4][0] = -len*w_*c/6;
    dSdp_elem[5][0] = 0;
    dSdp_elem[0][1] = len*w_*s/6;
    dSdp_elem[1][1] = -len*w_*c/6;
    dSdp_elem[2][1] = 0; // Moment term (pressure independent)
    dSdp_elem[3][1] = len*w_*s/3;
    dSdp_elem[4][1] = -len*w_*c/3;
    dSdp_elem[5][1] = 0;

    // subtract element terms, because -f is on the left side
    for (int k=0; k<3; k++) {
      if (nodeL.type[k] == 1) {
        dSdp[3*idL+k][idE] -= dSdp_elem[k][0];
        dSdp[3*idL+k][idE+1] -= dSdp_elem[k][1];
      }
      if (nodeR.type[k] == 1) {
        dSdp[3*idR+k][idE] -= dSdp_elem[3+k][0];
        dSdp[3*idR+k][idE+1] -= dSdp_elem[3+k][1];
      }
    }

    dSdp_elem.clear();
  }

  // Calculate (dS/du)*wrk
  for (int i=0; i<nnp*3; i++) {
    out(i) = 0;
    for (int k=0; k<nnp; k++) {
      out(i) += dSdp[i][k] * in(k); // perform the multiplication
    }
  }

  // Clean-up
  dSdp.clear();
}

// =====================================================================

void LECSM::CalcTrans_dSdp_Product(InnerProdVector& in, InnerProdVector& out)
{
  // Initialize the global derivative matrix and zero out the resultant
  int nnp = geom_.nnp;
  vector< vector<double> > dSdp(nnp*3, vector<double>(nnp, 0.0));

  // Loop over elements, calculating dS/dp at the element level before
  // adding the contributions into u_cfd
  int nel = geom_.nel;
  Element elem;
  Node nodeL, nodeR;
  int idE, idL, idR, type[3];
  double x1, x2, y1, y2, len, c, s, dFxdp, dFydp;
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
    dSdp_elem[0][0] = len*w_*s/3;     // dF1x/dP1
    dSdp_elem[1][0] = -len*w_*c/3;    // dF1y/dP1
    dSdp_elem[2][0] = 0; // Moment term (pressure independent)
    dSdp_elem[3][0] = len*w_*s/6;     //
    dSdp_elem[4][0] = -len*w_*c/6;
    dSdp_elem[5][0] = 0; // Moment term (pressure independent)
    dSdp_elem[0][1] = len*w_*s/6;
    dSdp_elem[1][1] = -len*w_*c/6;
    dSdp_elem[2][1] = 0; // Moment term (pressure independent)
    dSdp_elem[3][1] = len*w_*s/3;
    dSdp_elem[4][1] = -len*w_*c/3;
    dSdp_elem[5][1] = 0; // Moment term (pressure independent)

    for (int k=0; k<3; k++) {
      if (nodeL.type[k] == 1) {
        dSdp[3*idL+k][idL] -= dSdp_elem[k][0];
        dSdp[3*idL+k][idR] -= dSdp_elem[k][1];
      }
      if (nodeR.type[k] == 1) {
        dSdp[3*idR+k][idL] -= dSdp_elem[3+k][0];
        dSdp[3*idR+k][idR] -= dSdp_elem[3+k][1];
      }
    }
    dSdp_elem.clear();
  }

  // Calculate [(dS/du)^T]*u_csm
  for (int i=0; i<nnp; i++) {
    out(i) = 0.0;
    for (int k=0; k<3*nnp; k++) {
      out(i) += dSdp[k][i] * in(k); // perform the multiplication
    }
  }

  // Clean-up
  dSdp.clear();
}
// =====================================================================

void LECSM::CalcCoordsAndArea()
{
  int nnp = geom_.nnp;
  for (int i=0; i<nnp_; i++) {
    Node node = geom_.allNodes[i];
    xCoords_(i) = node.coords[0];
    if (node.type[0] == 1)
      xCoords_(i) += u_(3*i);
    yCoords_(i) = node.coords[1];
    if (node.type[1] == 1)
      yCoords_(i) += u_(3*i+1);
    area_(i) = w_*(h_ - 2.0*yCoords_(i));
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
  vector< vector<double> > K(ndof, vector<double>(ndof, 0.0));
  InitGlobalVecs(G, F);

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

int LECSM::SolveFor(InnerProdVector & rhs)
{
  // Generate the global equation number mapping
  int nnp = geom_.nnp;
  vector< vector< vector<int> > > gm(3, vector< vector<int> >(nnp, vector<int>(2)));
  geom_.SetupEq(gm);

  // Make sure the RHS vector has zeros for all fixed degrees of freedom
  for (int i=0; i<nnp; i++) {
    Node node = geom_.allNodes[i];
    for (int j=0; j<3; j++) {
      if (node.type[j] == 1)
        rhs(3*i+j) = 0;
    }
  }

  kona::MatrixVectorProduct<InnerProdVector>* 
      mat_vec = new StiffnessVectorProduct(this);
  kona::Preconditioner<InnerProdVector>*
      precond = new ApproxStiff(this);

  string filename = "lecsm_krylov.dat";
  ofstream fout(filename.c_str());
  int max_iter = 40;
  double tol = 1e-6;
  int precond_calls = 0;
  u_ = 0.0;
  kona::FGMRES(max_iter, tol, rhs, u_, *mat_vec, *precond,
               precond_calls, fout);  
  return precond_calls;
  
#if 0
  // Generate the global equation number mapping
  int nnp = geom_.nnp;
  vector< vector< vector<int> > > gm(3, vector< vector<int> >(nnp, vector<int>(2)));
  geom_.SetupEq(gm);

  // Initiate global vectors used in the solver
  int ndof = geom_.ndof;
  int ndog = geom_.ndog;
  vector<double> G(ndog), F(ndof);
  vector< vector<double> > K(ndof, vector<double>(ndof, 0.0));
  InitGlobalVecs(G, F);

  // Calculate the stiffness matrix and the forcing vector
  GetStiff(gm, G, F, K);
  F.clear();

  // extract values from the RHS matrix according to degree of freedoms
  int p;
  vector<double> F_new(ndof);
  for (int i=0; i<nnp; i++) {
    for (int j=0; j<3; j++) {
      if (gm[j][i][0] == 1) {
        p = gm[j][i][1];
        F_new[p] = rhs(3*i+j);
      }
    }
  }

  // Solve the global Kd = F system
  vector<double> disp(ndof);
  int maxIt = 1000;
  int iter = CGSolve(K, ndof, ndof, F_new, ndof, maxIt, disp);
  printf("LECSM: Solver converged in %i iterations!\n", iter);

  // Assemble the nodal displacements
  for (int A = 0; A < nnp; A++)
  {
    for (int i = 0; i < 3; i++)
    {
      int t = gm[i][A][0];
      double P = gm[i][A][1];
      if (t == 1)             // dof
        {u_(3*A+i) = disp[P];}
      else
        {
          if (t == 2)          // dog
            {u_(3*A+i) = G[P];}
          else
            {u_(3*A+i) = 0;}
        }
    }
  }
#endif
}

// =====================================================================

void LECSM::Solve()
{
#if 0
  geom_.InspectElements();
#endif

  // Generate the global equation number mapping
  int nnp = geom_.nnp;
  vector< vector< vector<int> > > gm(3, vector< vector<int> >(nnp, vector<int>(2)));
  geom_.SetupEq(gm);
#if 0
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
  vector< vector<double> > K(ndof, vector<double>(ndof, 0.0));
  InitGlobalVecs(G, F);

  // Calculate the stiffness matrix and the forcing vector
  GetStiff(gm, G, F, K);
#if 0
  printf("Printing the stiffness matrix for inspection:\n");
  printMatrix(K, ndof, ndof);
  printf("Printing the forcing vector for inspection:\n");
  for (int i=0; i<ndof; i++) {
    printf("|   %f   |\n", F[i]);
  }
#endif

  // Solve the global Kd = F system
  vector<double> disp(ndof);
  int maxIt = 100;
  int iter = CGSolve(K, ndof, ndof, F, ndof, maxIt, disp);
  printf("LECSM: Solver converged in %i iterations!\n", iter);

  // Assemble the nodal displacements
  printf("Directions:\n");
  printf("  0 - x-axis\n");
  printf("  1 - y-axis\n");
  printf("  2 - rotation about z-axis\n");
  for (int A = 0; A < nnp; A++)
  {
    for (int i = 0; i < 3; i++)
    {
      int t = gm[i][A][0];
      double P = gm[i][A][1];
      if (t == 1)             // dof
        {u_(3*A+i) = disp[P];}
      else
        {
          if (t == 2)          // dog
            {u_(3*A+i) = G[P];}
          else
            {u_(3*A+i) = 0;}
        }
      printf("    Node %d displaced %f in direction %d\n", A, u_(3*A+i), i);
    }
  }
}

// ======================================================================

void StiffnessVectorProduct::operator()(const InnerProdVector & u, 
                                        InnerProdVector & v) { 
  solver_->Calc_dSdu_Product(u, v);
}
