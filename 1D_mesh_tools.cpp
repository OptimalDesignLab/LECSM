/**
 * \file 1D_mesh_tools.cpp
 * \brief 1D mesh library
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

#include <stdio.h>
#include <math.h>
#include "./1D_mesh_tools.hpp"
#include "./matrix_tools.hpp"
using namespace std;

// =====================================================================

void Node::CreateNode(int num, double* c) {
  id = num;
  for (int i=0; i<3; i++) {
    coords[i] = c[i];
    type[i] = 1;            // node is initially free
    dispBC[i] = NULL;       // no displacement BCs defined
    forceBC[i] = 0;         // force and moment BCs are zero
  }             
}

// =====================================================================

void Node::DefineBCs(int* BCtype, double* BCval) {
  for (int i=0; i<3; i++) {    // loop over the three dimensions of BCs
    if (BCtype[i] == 0) {      // displacement BC
      dispBC[i] = BCval[i];
      if (BCval[i] == 0) {     // zero BC (fixed node)
        type[i] = 0;
      }
      else {                   // non-zero BC (prescribed displacement)
        type[i] = 2;
      }
    }
    else if (BCtype[i] == 1) { // force or moment BC
      forceBC[i] = BCval[i];
    }
  }
}

// =====================================================================

void Element::CreateElem(int num, vector<Node> nodes) {
  id = num;
  nen = nodes.size();
  adjNodes = nodes;
  Node nodeL = adjNodes[0];
  Node nodeR = adjNodes[1];
  double x1 = nodeL.coords[0];
  double x2 = nodeR.coords[0];
  double y1 = nodeL.coords[1];
  double y2 = nodeR.coords[1];
  length = sqrt(pow(x2-x1,2)+pow(y2-y1,2));
  cosine = (x2 - x1)/length;
  sine = (y2 - y1)/length;
}

// =====================================================================

void Element::GetElemStiff(double E, double w, double t, vector<double>& P,
                           vector< vector< vector<int> > >& gm,
                           vector< vector< vector<int> > >& lm,
                           vector< vector<double> >& KE, vector<double>& FE) {

  // Recover the left and right nodes of the element
  Node nodeL = adjNodes[0];
  Node nodeR = adjNodes[1];
  int idL = nodeL.id;
  int idR = nodeR.id;

  // Generate the local equation mapping
  lm[0][0][0] = gm[0][idL][0];
  lm[0][0][1] = gm[0][idL][1];
  lm[0][1][0] = gm[0][idR][0];
  lm[0][1][1] = gm[0][idR][1];
  lm[1][0][0] = gm[1][idL][0];
  lm[1][0][1] = gm[1][idL][1];
  lm[1][1][0] = gm[1][idR][0];
  lm[1][1][1] = gm[1][idR][1];
  lm[2][0][0] = gm[2][idL][0];
  lm[2][0][1] = gm[2][idL][1];
  lm[2][1][0] = gm[2][idR][0];
  lm[2][1][1] = gm[2][idR][1];

  // Calculate the local element stiffness matrix
  double A = w*t; // cross section area of the element
  double I = w*pow(t,3)/12; // area moment of inertia of the beam element x-section
  vector< vector<double> > KEloc(nen*3, vector<double>(nen*3));
  KEloc[0][0] = A*E/length;
  KEloc[0][1] = 0;
  KEloc[0][2] = 0;
  KEloc[0][3] = -A*E/length;
  KEloc[0][4] = 0;
  KEloc[0][5] = 0;
  KEloc[1][0] = 0;
  KEloc[1][1] = 12*E*I/pow(length,3);
  KEloc[1][2] = 6*E*I/pow(length,2);
  KEloc[1][3] = 0;
  KEloc[1][4] = -12*E*I/pow(length,3);
  KEloc[1][5] = 6*E*I/pow(length,2);
  KEloc[2][0] = 0;
  KEloc[2][1] = 6*E*I/pow(length,2);
  KEloc[2][2] = 4*E*I/length;
  KEloc[2][3] = 0;
  KEloc[2][4] = -6*E*I/pow(length,2);
  KEloc[2][5] = 2*E*I/length;
  KEloc[3][0] = -A*E/length;
  KEloc[3][1] = 0;
  KEloc[3][2] = 0;
  KEloc[3][3] = A*E/length;
  KEloc[3][4] = 0;
  KEloc[3][5] = 0;
  KEloc[4][0] = 0;
  KEloc[4][1] = -12*E*I/pow(length,3);
  KEloc[4][2] = -6*E*I/pow(length,2);
  KEloc[4][3] = 0;
  KEloc[4][4] = 12*E*I/pow(length,3);
  KEloc[4][5] = -6*E*I/pow(length,2);
  KEloc[5][0] = 0;
  KEloc[5][1] = 6*E*I/pow(length,2);
  KEloc[5][2] = 2*E*I/length;
  KEloc[5][3] = 0;
  KEloc[5][4] = -6*E*I/pow(length,2);
  KEloc[5][5] = 4*E*I/length;

  // Create the local to global transformation matrix
  vector< vector<double> > T(nen*3, vector<double>(nen*3));
  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      T[i][j] = 0;
    }
  }
  T[0][0] = cosine;
  T[0][1] = sine;
  T[1][0] = -sine;
  T[1][1] = cosine;
  T[2][2] = 1;
  T[3][3] = cosine;
  T[3][4] = sine;
  T[4][3] = -sine;
  T[4][4] = cosine;
  T[5][5] = 1;

  // Calculate the global element stiffness matrix
  vector< vector<double> > Tt(nen*3, vector<double>(nen*3));
  matrixTranspose(T, nen*3, nen*3, Tt);
  vector< vector<double> > KT(nen*3, vector<double>(nen*3));
  matrixMult(KEloc, nen*3, nen*3, T, nen*3, nen*3, KT);  
  matrixMult(Tt, nen*3, nen*3, KT, nen*3, nen*3, KE);

  // Create the element forcing vector due to pressure
  // Add in the nodal force contributions
  double p1 = P[0];
  double p2 = P[1];
  double Pave = 0.5*(p1+p2);
  double f_hat = -Pave*length*w/2;
  double fy = f_hat*cosine;
  double fx = f_hat*sine;
  FE[0] = fx + nodeL.forceBC[0];
  FE[1] = fy + nodeL.forceBC[1];
  FE[2] = nodeL.forceBC[2];
  FE[3] = fx + nodeR.forceBC[0];
  FE[4] = fy + nodeR.forceBC[1];
  FE[5] = nodeR.forceBC[2];

  // Clean-up
  T.clear();
  Tt.clear();
  KEloc.clear();
}

// =====================================================================

void Element::Assemble(vector< vector<double> >& KE, vector<double>& FE,
                       vector< vector< vector<int> > >& lm,
                       vector<double>& G, vector<double>& F,
                       vector< vector<double> >& K)
{
  int p = 0;
  for (int a = 0; a < nen; a++)
  {
    for (int i = 0; i < 3; i++)
    {
      if (lm[i][a][0] == 1) // dof
      {
        int P = lm[i][a][1];
        F[P] = F[P] + FE[p];
        int q = 0;
        for (int b = 0; b < nen; b++)
        {
          for (int j = 0; j < 3; j++)
          {
            int Q = lm[j][b][1];
            if (lm[j][b][0] == 1) // dof
              {K[P][Q] = K[P][Q] + KE[p][q];}
            else if (lm[j][b][0] == 2) // dog
              {F[P] = F[P] - G[Q]*KE[p][q];}
            q++;
          } // end j loop over nsd
        } // end b loop over nen (columns)
      }
      p++;
    } // end i loop over nsd
  } // end a loop over nen (rows)
}

// =====================================================================

void Mesh::CreateMesh(vector<Element>& elems) {
  allElems = elems;
  nel = allElems.size();
  nnp = 0;
  Element elem = allElems.at(0);
  Node nodeL = elem.adjNodes[0];
  Node nodeR = elem.adjNodes[1];
  allNodes.push_back(nodeL);
  allNodes.push_back(nodeR);
  for (int i=1; i<nel; i++) {
    elem = allElems.at(i);
    allNodes.push_back(elem.adjNodes[1]);    
  }
  nnp = allNodes.size();
}

// =====================================================================

void Mesh::SetupEq(vector< vector< vector<int> > >& gm)
{
  ndof = 0;
  ndog = 0;
  Node nd;
  for (int b = 0; b < nnp; b++)
  {
    nd = allNodes[b];
    int a = nd.id;
    for (int i = 0; i < 3; i++)
    {
      if (nd.type[i]==1)        // DoF - possible nodal load
      {
        gm[i][a][0] = 1;        // store the node type
        gm[i][a][1] = ndof;     // store the equation number
        ndof++;
      }
      else                      // prescribed BC
      {
        if (nd.type[i]==2)      // DoG - non-zero BC
        {
          gm[i][a][0] = 2;
          gm[i][a][1] = ndog;
          ndog++;
        }
        else                    // zero BC (node.type[i] == 0)
        {
          gm[i][a][0] = 0;
          gm[i][a][1] = 0;
        }
      }
    }
  }
}

// =====================================================================

void Mesh::Update(const InnerProdVector& u_csm)
{
  // Loop over mesh nodes and update their coordinates
  Node node;
  for (int i=0; i<nnp; i++) {
    node = allNodes[i];
    node.coords[0] += u_csm(3*i);
    node.coords[1] += u_csm(3*i+1);
  }
}