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
#include <boost/property_tree/ptree.hpp>

using namespace std;

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

void LECSM::set_coords(const InnerProdVector & x, const InnerProdVector & y)
{
  xCoords_ = x;
  yCoords_ = y;
}

// =====================================================================

void LECSM::UpdateMesh()
{
  int nnp = geom_.nnp;
  int nel = geom_.nel;
  for (int i=0; i<nnp; i++) {
    geom_.allNodes[i].coords[0] = xCoords_(i);
    geom_.allNodes[i].coords[1] = yCoords_(i);
  }
  for (int i=0; i<nel; i++) {
    geom_.allElems[i].adjNodes[0] = geom_.allNodes[i];
    geom_.allElems[i].adjNodes[1] = geom_.allNodes[i+1];
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
    elem.Assemble<double>(KE, FE, lm, G, F, K);
  }
}

// =====================================================================

template <typename type>
void LECSM::GetStiff(const vector<type>& x, const vector<type>& y,
                     vector< vector< vector<int> > >& gm,
                     vector<type>& G, vector<type>& F,
                     vector< vector<type> >& K)
{
  // Loop over all elements, assuming that each face is an element.
  Element elem;
  int nel = geom_.nel;
  Node nodeL, nodeR;
  type tE = static_cast<type>(E_);
  type tw = static_cast<type>(w_);
  type tt = static_cast<type>(t_);
  for (int i=0; i<nel; i++)
  {
    // Get information about the element.
    elem = geom_.allElems[i];
    int nen = elem.nen;

    // Calculate the stiffness matrix
    vector< vector< vector<int> > >
        lm(3, vector< vector<int> >(nen, vector<int>(2)));
    vector< vector<type> > KE(nen*3, vector<type>(nen*3));
    vector<type> FE(nen*3);
    nodeL = elem.adjNodes[0];
    nodeR = elem.adjNodes[1];
    vector<type> locP(2);
    locP[0] = static_cast<type>(P_(nodeL.id));
    locP[1] = static_cast<type>(P_(nodeR.id));
    type x1 = x[nodeL.id];
    type x2 = x[nodeR.id];
    type y1 = y[nodeL.id];
    type y2 = y[nodeR.id];
    elem.GetElemStiff<type>(x1, x2, y1, y2, tE, tw, tt, locP, gm, lm, KE, FE);
    
    // Assemble the element contributions into the global matrices.
    elem.Assemble<type>(KE, FE, lm, G, F, K);
  }
}

// explicit instantiations
template void LECSM::GetStiff<double>(const vector<double>& x, const vector<double>& y,
                                      vector< vector< vector<int> > >& gm,
                                      vector<double>& G, vector<double>& F,
                                      vector< vector<double> >& K);
template void LECSM::GetStiff<complex<double> >(
    const vector<complex<double> >& x, const vector<complex<double> >& y,
    vector< vector< vector<int> > >& gm,
    vector<complex<double> >& G, vector<complex<double> >& F,
    vector< vector<complex<double> > >& K);

// =====================================================================

void LECSM::Precondition(InnerProdVector& in, InnerProdVector& out)
{
  kona::Preconditioner<InnerProdVector>*
          precond = new GaussSeidelPrecond(*this, this->geom_);
  (*precond)(in, out);
  delete precond;
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
        out(3*i+j) = 0.0;
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
    if (node.type[1] == 1)
      dAdu[i][3*i+1] = -2*w_;
    else
      dAdu[i][3*i+1] = 0;
    dAdu[i][3*i+2] = 0;
  }

  for (int i=0; i<nnp; i++) {
    out(i) = 0;
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

void LECSM::Calc_dydA_Product(InnerProdVector& in, InnerProdVector& out)
{
  int nnp = geom_.nnp;
  for (int i=0; i<nnp; i++) {
    Node node = geom_.allNodes[i];
    if (node.type[1] == 1)
      out(i) = -1.0*in(i)/(2.0*w_);
    else
      out(i) = 0.0;
  }
}

// =====================================================================

void LECSM::CalcFD_dSdy_Product(InnerProdVector& in, InnerProdVector& out)
{
  // Calculate the Residual at the unperturbed state
  CalcResidual();
  InnerProdVector res0 = get_res();

  // Perturb the y coordinates
  double eps = kona::CalcEpsilon(yCoords_.Norm2(), in.Norm2());
  InnerProdVector save_y = yCoords_;
  for (int i = 0; i < nnp_; i++)
    yCoords_(i) += eps*in(i);

  // Re-evaluate residual
  UpdateMesh();
  CalcResidual();
  out = get_res();
  out -= res0;
  out /= eps;

  // Reset coordinates
  yCoords_ = save_y;
  UpdateMesh();
  
#if 0
  // NOTE: Cannot find the paper on easier finite differencing
  // ~~~ FIX THIS LATER ~~~
  ResetCoords();
  InnerProdVector save_y = yCoords_;
  double delta = 1.e-6;
  InnerProdVector disp(3*nnp_, 1.0);
  set_u(disp);

  CalcResidual();
  InnerProdVector res0 = get_res();

  // Calculate derivative via finite differencing
  vector< vector<double> > dSdy(3*nnp_, vector<double>(nnp_, 0.0));
  for (int i=0; i < nnp_; i++) {
    yCoords_(i) += delta;
    UpdateMesh();
    CalcResidual();
    InnerProdVector res1 = get_res();
    for (int j=0; j < 3*nnp_; j++) {
      dSdy[j][i] = (res1(j) - res0(j))/delta;
    }
    yCoords_(i) = save_y(i);
    UpdateMesh();
  }

  // Perform the multiplication
  for (int i=0; i < 3*nnp_; i++) {
    out(i) = 0;
    for (int j=0; j < nnp_; j++)
      out(i) += dSdy[i][j] * in(j);
  }

  // Clean-up
  dSdy.clear();
#endif
}

// =====================================================================

void LECSM::CalcCmplx_dSdy_Product(InnerProdVector& in, InnerProdVector& out)
{
  // complex Perturb the y coordinates
  double eps = 1e-40;
  vector<complex<double> > x_cmplx(nnp_), y_cmplx(nnp_), res_cmplx(3*nnp_);
  for (int i = 0; i < nnp_; i++) {
    x_cmplx[i] = complex<double>(geom_.allNodes[i].coords[0], 0.0);
    y_cmplx[i] = complex<double>(geom_.allNodes[i].coords[1], eps*in(i));
  }
  
  // Calculate the Residual at the unperturbed state
  CalcResidual<complex<double> >(x_cmplx, y_cmplx, res_cmplx);

  // compute the product
  for (int i = 0; i < 3*nnp_; i++)
    out(i) = imag(res_cmplx[i])/eps;

  x_cmplx.clear();
  y_cmplx.clear();
  res_cmplx.clear();
}

// =====================================================================

void LECSM::CalcTransFD_dSdy_Product(InnerProdVector& in, InnerProdVector& out)
{ 
  // Set up parameters necessary for FD
  ResetCoords();
  InnerProdVector save_y = yCoords_;
  double delta = 1.e-7;

  CalcResidual();
  InnerProdVector res0 = get_res();

  // Calculate derivative via finite differencing
  vector< vector<double> > dSdy(3*nnp_, vector<double>(nnp_, 0.0));
  for (int i=0; i < nnp_; i++) {
    double eps = std::max(delta*yCoords_(i), kona::kEpsilon);
    yCoords_(i) += eps;
    UpdateMesh();
    CalcResidual();
    InnerProdVector res1 = get_res();
    for (int j=0; j < 3*nnp_; j++) {
      dSdy[j][i] = (res1(j) - res0(j))/eps;
    }
    yCoords_(i) = save_y(i);
    UpdateMesh();
  }

  // Perform the transpose multiplication
  for (int i=0; i < nnp_; i++) {
    out(i) = 0;
    for (int j=0; j < 3*nnp_; j++)
      out(i) += dSdy[j][i] * in(j);
  }

  // Clean-up
  dSdy.clear();  
}

// =====================================================================

void LECSM::CalcTransCmplx_dSdy_Product(InnerProdVector& in,
                                        InnerProdVector& out)
{
  double eps = 1e-40;
  vector<complex<double> > x_cmplx(nnp_), y_cmplx(nnp_), res_cmplx(3*nnp_);
  for (int i = 0; i < nnp_; i++) {
    x_cmplx[i] = complex<double>(geom_.allNodes[i].coords[0], 0.0);
    y_cmplx[i] = complex<double>(geom_.allNodes[i].coords[1], 0.0);
  }
  
  // Calculate matrix via complex step
  vector< vector<double> > dSdy(3*nnp_, vector<double>(nnp_, 0.0));
  for (int i=0; i < nnp_; i++) {
    y_cmplx[i] = complex<double>(real(y_cmplx[i]), eps);
    CalcResidual<complex<double> >(x_cmplx, y_cmplx, res_cmplx);
    for (int j=0; j < 3*nnp_; j++)
      dSdy[j][i] = imag(res_cmplx[j])/eps;
    y_cmplx[i] = complex<double>(real(y_cmplx[i]), 0.0);
  }

  // Perform the transpose multiplication
  for (int i=0; i < nnp_; i++) {
    out(i) = 0;
    for (int j=0; j < 3*nnp_; j++)
      out(i) += dSdy[j][i] * in(j);
  }

  // Clean-up
  dSdy.clear();
  x_cmplx.clear();
  y_cmplx.clear();
  res_cmplx.clear();
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
        res_(3*i+j) = 0.0;
    }
  }
  res_dof.clear();
}

// =====================================================================

template <typename type>
void LECSM::CalcResidual(const vector<type>& x, const vector<type>& y,
                         vector<type>& res)
{
  // Generate the global equation number mapping
  int nnp = geom_.nnp;
  vector< vector< vector<int> > >
      gm(3, vector< vector<int> >(nnp, vector<int>(2)));
  geom_.SetupEq(gm);

  // Initiate global vectors used in the solver
  int ndof = geom_.ndof;
  int ndog = geom_.ndog;
  vector<type> G(ndog), F(ndof);
  vector< vector<type> > K(ndof, vector<type>(ndof, 0.0));
  Node nd;
  for (int b = 0; b < nnp; b++) {
    nd = geom_.allNodes[b];
    for (int i = 0; i < 3; i++) {
      if (nd.type[i]==1)        // DoF - possible nodal load
        F.push_back(static_cast<type>(nd.forceBC[i])); // store the nodal load
      else                       // prescribed BC
        if (nd.type[i]==2)      // DoG - non-zero BC
          G.push_back(static_cast<type>(nd.dispBC[i]));// store the essential BC
    }
  }
  //InitGlobalVecs(G, F);

  // Calculate the stiffness matrix and the forcing vector
  GetStiff(x, y, gm, G, F, K);
  G.clear();

  int p;
  vector<type> u_dof(ndof);
  for (int i=0; i<nnp; i++) {
    for (int j=0; j<3; j++) {
      if (gm[j][i][0] == 1) {
        p = gm[j][i][1];
        u_dof[p] = static_cast<type>(u_(3*i+j));
      }
    }
  }

  // Calculate the K*u product
  vector<type> v_dof(ndof);
  for (int i = 0; i < ndof; i++) {
    v_dof[i] = 0;
    for (int j = 0; j < ndof; j++)
      v_dof[i] += K[i][j]*u_dof[j];
  }  
  //matrixVecMult(K, ndof, ndof, u_dof, ndof, v_dof);
  u_dof.clear();
  K.clear();

  // Form the Ku-f residual for free nodes
  vector<type> res_dof(ndof);
  for (int i=0; i<ndof; i++)
    res_dof[i] = v_dof[i] - F[i];
  v_dof.clear();
  F.clear();

  // Assemble the whole residual
  for (int i=0; i<nnp; i++) {
    for (int j=0; j<3; j++) {
      if (gm[j][i][0] == 1)   // node is free
        res[3*i+j] = res_dof[gm[j][i][1]];
      else  // node is fixed
        res[3*i+j] = static_cast<type>(0.0);
    }
  }
  res_dof.clear();
}

// explicit instantiations
template void LECSM::CalcResidual<double>(const vector<double>& x, const vector<double>& y,
                                          vector<double>& res);
template void LECSM::CalcResidual<complex<double> >(
    const vector<complex<double> >& x, const vector<complex<double> >& y,
    vector<complex<double> >& res);

// =====================================================================

int LECSM::SolveFor(InnerProdVector & rhs, const int & max_iter,
                    const double & tol)
{
  // Make sure the RHS vector has zeros for all fixed degrees of freedom  
  int nnp = geom_.nnp;
  for (int i=0; i<nnp; i++) {
    Node node = geom_.allNodes[i];
    for (int j=0; j<3; j++) {
      if (node.type[j] != 1)
        rhs(3*i+j) = 0.0;
    }
  }

  kona::MatrixVectorProduct<InnerProdVector>* 
      mat_vec = new StiffnessVectorProduct(*this, this->geom_);
  kona::Preconditioner<InnerProdVector>*
      precond = new GaussSeidelPrecond(*this, this->geom_);

  string filename = "lecsm_krylov.dat";
  ofstream fout(filename.c_str());
  int precond_calls = 0;
  u_ = 0.0;
  kona::FGMRES<InnerProdVector>(max_iter, tol, rhs, u_, *mat_vec, *precond,
             precond_calls, fout);
//  kona::MINRESSolver<InnerProdVector> solver;
//  solver.SubspaceSize(max_iter);
//  using boost::property_tree::ptree;
//  ptree input_params, output_params;
//  input_params.put<double>("tol", tol);
//  input_params.put<bool>("check", true);
//  solver.Solve(input_params, rhs, u_, *mat_vec, *precond, output_params, fout);
//  precond_calls = output_params.get<int>("iters");
  return precond_calls;
}

// =====================================================================

void LECSM::Solve(bool info)
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
  int maxIt = 10000;
  int iter = CGSolve(K, ndof, ndof, F, ndof, maxIt, disp);
  if (info) printf("LECSM: Solver converged in %i iterations!\n", iter);

#if 0
  // Assemble the nodal displacements
  printf("Directions:\n");
  printf("  0 - x-axis\n");
  printf("  1 - y-axis\n");
  printf("  2 - rotation about z-axis\n");
#endif
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
            {u_(3*A+i) = 0.0;}
        }
      //printf("    Node %d displaced %f in direction %d\n", A, u_(3*A+i), i);
    }
  }
}

// ======================================================================

void LECSM::StiffDiagProduct(const InnerProdVector & in,
                             InnerProdVector & out)
{
   // Set up the global adjacency mapping
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
   vector<double> Kdiag(ndof, 0.0);
   for (int i=0; i<ndof; i++)
      Kdiag[i] = K[i][i];
   K.clear();

   // Perform the diagonal multiplication
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
   for (int i=0; i<ndof; i++)
      v_dof[i] = Kdiag[i]*u_dof[i];
   u_dof.clear();
   Kdiag.clear();

   // Assemble the whole product
   for (int i=0; i<nnp; i++) {
     for (int j=0; j<3; j++) {
       if (gm[j][i][0] == 1)   // node is free
         out(3*i+j) = v_dof[gm[j][i][1]];
       else  // node is fixed
         out(3*i+j) = 0.0;
     }
   }

   // Clean-up
   gm.clear();
   v_dof.clear();
}

// ======================================================================

StiffnessVectorProduct::StiffnessVectorProduct(LECSM& solver, Mesh& geom) {
  // Generate the global equation number mapping
  nnp_ = geom.nnp;
  gm_.resize(3, vector< vector<int> >(nnp_, vector<int>(2)));
  geom.SetupEq(gm_);
  // Build the stiffness matrix
  ndof_ = geom.ndof;
  int ndog = geom.ndog;
  vector<double> G(ndog), F(ndof_);
  K_.resize(ndof_, vector<double>(ndof_, 0.0));
  u_dof_.resize(ndof_, 0.0);
  v_dof_.resize(ndof_, 0.0);
  solver.InitGlobalVecs(G, F);
  solver.GetStiff(gm_, G, F, K_);  
} 

// ======================================================================

void StiffnessVectorProduct::operator()(const InnerProdVector & u, 
                                        InnerProdVector & v) {
  // move u into reduced vector
  for (int i = 0; i < nnp_; i++) {
    for (int j = 0; j < 3; j++) {
      if (gm_[j][i][0] == 1)
        u_dof_[gm_[j][i][1]] = u(3*i+j);
    }
  }  
  // Perform the multiplication.
  for (int i = 0; i < ndof_; i++) {
    v_dof_[i] = 0.0;
    for (int j = 0; j < ndof_; j++)
      v_dof_[i] += K_[i][j]*u_dof_[j];
  }
  // move v_dof_ into the full vector
  for (int i = 0; i < nnp_; i++) {
    for (int j = 0; j < 3; j++) {
      if (gm_[j][i][0] == 1) // node is free
        v(3*i+j) = v_dof_[gm_[j][i][1]];
      else // node is fixed
        v(3*i+j) = 0.0;
    }
  }
}

// ======================================================================

DiagonalPrecond::DiagonalPrecond(LECSM& solver, Mesh& geom) {
  // Generate the global equation number mapping
  nnp_ = geom.nnp;
  gm_.resize(3, vector< vector<int> >(nnp_, vector<int>(2)));
  geom.SetupEq(gm_);
  // Build the stiffness matrix
  ndof_ = geom.ndof;
  int ndog = geom.ndog;
  vector<double> G(ndog), F(ndof_);
  vector< vector<double> > K(ndof_, vector<double>(ndof_, 0.0));
  Kdiag_.resize(ndof_, 0.0);
  u_dof_.resize(ndof_, 0.0);
  v_dof_.resize(ndof_, 0.0);
  solver.InitGlobalVecs(G, F);
  solver.GetStiff(gm_, G, F, K);
  G.clear();
  F.clear();

  // Calculate the stiffness matrix and the forcing vector
  for (int i = 0; i < ndof_; i++) {
    if (K[i][i] < 1e-16)
      Kdiag_[i] = 1.0;
    else
      Kdiag_[i] = 1.0/K[i][i];
  }
  K.clear();
}

// ======================================================================

void DiagonalPrecond::operator()(InnerProdVector & u, 
                             InnerProdVector & v) { 
  //solver_->StiffDiagProduct(u, v);
  // move u into reduced vector
  for (int i = 0; i < nnp_; i++) {
    for (int j = 0; j < 3; j++) {
      if (gm_[j][i][0] == 1)
        u_dof_[gm_[j][i][1]] = u(3*i+j);
    }
  }  
  // Perform the multiplication with the diagonal matrix.
  for (int i = 0; i < ndof_; i++)
    v_dof_[i] = Kdiag_[i]*u_dof_[i];
  // move v_dof_ into the full vector
  for (int i = 0; i < nnp_; i++) {
    for (int j = 0; j < 3; j++) {
      if (gm_[j][i][0] == 1) // node is free
        v(3*i+j) = v_dof_[gm_[j][i][1]];
      else // node is fixed
        v(3*i+j) = 0.0;
    }
  }
}

// ======================================================================

GaussSeidelPrecond::GaussSeidelPrecond(LECSM& solver, Mesh& geom) {
  // Generate the global equation number mapping
  nnp_ = geom.nnp;
  gm_.resize(3, vector< vector<int> >(nnp_, vector<int>(2)));
  geom.SetupEq(gm_);
  // Build the stiffness matrix
  ndof_ = geom.ndof;
  int ndog = geom.ndog;
  vector<double> G(ndog), F(ndof_);
  K_.resize(ndof_, vector<double>(ndof_, 0.0));
  u_dof_.resize(ndof_, 0.0);
  v_dof_.resize(ndof_, 0.0);
  solver.InitGlobalVecs(G, F);
  solver.GetStiff(gm_, G, F, K_);
  max_iters_ = 10;
}

// ======================================================================

void GaussSeidelPrecond::operator()(InnerProdVector & u,
                                    InnerProdVector & v) {
    
  // move RHS into reduced vector
  for (int i = 0; i < nnp_; i++) {
    for (int j = 0; j < 3; j++) {
      if (gm_[j][i][0] == 1)
        u_dof_[gm_[j][i][1]] = u(3*i+j);
    }
  }

  // Gauss Seidel iterations
  for (int k = 0; k < max_iters_; k++) {
    // iterations over DOFs
    for (int i = 0; i < ndof_; i++) {
      v_dof_[i] = u_dof_[i];
      for (int j = 0; j < ndof_; j++) {
        if (j != i) v_dof_[i] -= K_[i][j]*v_dof_[j];
      }
      v_dof_[i] *= 1.0/K_[i][i];
    }
  }

  // move the solution into the output vector
  v = 0.0;
  for (int i = 0; i < nnp_; i++) {
    for (int j = 0; j < 3; j++) {
      if (gm_[j][i][0] == 1) // node is free
        v(3*i+j) = v_dof_[gm_[j][i][1]];
      else // node is fixed
        v(3*i+j) = 0.0;
    }
  }
}