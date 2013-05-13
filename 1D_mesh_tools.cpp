/**
 * \file 1D_mesh_tools.cpp
 * \a 1D mesh library
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

#include <stdio.h>
#include <math.h>
#include "./1D_mesh_tools.hpp"
#include "./matrix_tools.hpp"
using namespace std;

// =====================================================================

void Node::CreateNode() {
  id = -1;
  type = -1;
}

void Node::CreateNode(int id_, int type_, double coords_) {
  id = id_;
  type = 1; // node is free to move by default
  coords = coords_;
  dispBC = NULL;
  forceBC = {0,0,0};
}

void Node::DefineBCs(int BCtype, double BCval) {
  if (BCtype == 0) { // displacement BC
    if (BCval == NULL) { // zero BC
      type = 0;
      dispBC = {0,0,0};
    }
    else { // non-zero BC
      type = 2;
      dispBC = BCval;
    }
  }
  else { // nodal force BC
    forceBC = BCval;
  }
}

// =====================================================================

void Element::CreateElem() {
  id = -1;
  nen = 0;
}

void Element::CreateElem(int id_, vector<Node> nodes) {
  id = id_;
  nen = nodes.size();
  adjNodes = nodes;
}

// =====================================================================

void Mesh::CreateMesh() {
  nnp = 0;
  nel = 0;
}

void Mesh::CreateMesh(vector<Element> elems) {
  allElems = elems;
  nel = allElems.size();
  nnp = 0;
  Element elem = allElems.at(0);
  Node nodeL = elem.adjNodes[0];
  Node nodeR = elem.adjNodes[1];
  vector<Node> allNodes = {nodeL, nodeR};
  for (int i=1; i<nel; i++) {
    elem = allElems.at(i);
    allNodes.push_back(elem.adjNodes[1]);    
  }
  nnp = allNodes.size();
}

// =====================================================================

void Element::GetStiff(double E, double w, double t, double P,
                       vector< vector< vector<double> > > id,
                       vector< vector< vector<double> > > lm, 
                       vector< vector<double> >& KE,
                       vector<double>& FE) {
  // Recover the left and right nodes of the element
  Node nodeL = adjNodes[0];
  Node nodeR = adjNodes[1];
  // Generate the local equation mapping
  lm[0][0][0] = id[0][nodeL.id][0];
  lm[0][0][1] = id[0][nodeL.id][1];
  lm[0][1][0] = id[0][nodeR.id][0];
  lm[0][1][1] = id[0][nodeR.id][1];
  lm[1][0][0] = id[1][nodeL.id][0];
  lm[1][0][1] = id[1][nodeL.id][1];
  lm[1][1][0] = id[1][nodeR.id][0];
  lm[1][1][1] = id[1][nodeR.id][1];
  lm[2][0][0] = id[1][nodeL.id][0];
  lm[2][0][1] = id[1][nodeL.id][1];
  lm[2][1][0] = id[1][nodeR.id][0];
  lm[2][1][1] = id[1][nodeR.id][1];
  // Calculate the local element stiffness matrix
  double x1 = nodeL.coords[0];
  double y1 = nodeL.coords[1];
  double x2 = nodeR.coords[0];
  double y2 = nodeR.coords[1];
  double len = sqrt(pow(x2-x1,2) + pow(y2-y1,2));
  double c = (x2 - x1)/len;
  double s = (y2 - y1)/len;
  double A = w*t; // cross section area of the element
  double I = w*pow(h,3)/12 // area moment of inertia of the beam element x-section
  vector< vector<double> > KEloc(6, vector<double>(6));
  KEloc[0][0] = A*E/len;
  KEloc[0][1] = 0;
  KEloc[0][2] = 0;
  KEloc[0][3] = -A*E/len;
  KEloc[0][4] = 0;
  KEloc[0][5] = 0;
  KEloc[1][0] = 0;
  KEloc[1][1] = 12*E*I/pow(len,3);
  KEloc[1][2] = 6*E*I/pow(len,2);
  KEloc[1][3] = 0;
  KEloc[1][4] = -12*E*I/pow(len,3);
  KEloc[1][5] = 6*E*I/pow(len,2);
  KEloc[2][0] = 0;
  KEloc[2][1] = 6*E*I/pow(len,2);
  KEloc[2][2] = 4*E*I/len;
  KEloc[2][3] = 0;
  KEloc[2][4] = -6*E*I/pow(len,2);
  KEloc[2][5] = 2*E*I/len;
  KEloc[3][0] = -A*E/len;
  KEloc[3][1] = 0;
  KEloc[3][2] = 0;
  KEloc[3][3] = A*E/len;
  KEloc[3][4] = 0;
  KEloc[3][5] = 0;
  KEloc[4][0] = 0;
  KEloc[4][1] = -12*E*I/pow(len,3);
  KEloc[4][2] = -6*E*I/pow(len,2);
  KEloc[4][3] = 0;
  KEloc[4][4] = 12*E*I/pow(len,3);
  KEloc[4][5] = -6*E*I/pow(len,2);
  KEloc[5][0] = 0;
  KEloc[5][1] = 6*E*I/pow(len,2);
  KEloc[5][2] = 2*E*I/len;
  KEloc[5][3] = 0;
  KEloc[5][4] = -6*E*I/pow(len,2);
  KEloc[5][5] = 4*E*I/len;
  // Create the local to global transformation matrix
  vector< vector<double> > T(6, vector<double>(6))
  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      T[i][j] = 0;
    }
  }
  T[0][0] = c;
  T[0][1] = s;
  T[1][0] = -s;
  T[1][1] = c;
  T[2][2] = 1;
  T[3][3] = c;
  T[3][4] = s;
  T[4][3] = -s;
  T[4][4] = c;
  T[5][5] = 1;
  // Calculate the global element stiffness matrix
  vector< vector<double> > Tt(6, vector<double>(6))
  matrixTranspose(T, 6, 6, Tt);
  vector< vector<double> > KT(6, vector<double>(6))
  matrixMult(KEloc, 6, 6, T, 6, 6, KT);
  matrixMult(Tt, 6, 6, KT, 6, 6, KE);
  // Create the element forcing vector due to pressure
  // Add in the nodal force contributions
  double f_hat = -P[id]*len*w/2;
  double fy = f_hat*c;
  double fx = f_hat*s;
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

void Element::assemble(vector< vector< vector<double> > > lm,
                       vector< vector<double> > KE,
                       vector<double> FE, vector<double>& G,
                       vector< vector<double> >& K,
                       vector<double>& F)
{
  // This procedure assemples the element contributions into
  // the global matrices.
  //
  // Inputs:
  //    nsd             - number of space dimensions
  //    nen             - number of element nodes
  //    lm              - element location vector, identifying the
  //                      equation numbers and the location of
  //                      essential BCs for the element
  //    KE              - element stiffness matrix
  //    FE              - element force vector
  //    G               - global prescribed displacement vector
  // Outputs:
  //    K               - global stiffness matrix
  //    F               - global force matrix
  //
  printf("    Testing lm values:\n");
  for (int c = 0;c<2;c++)
  {
    for (int a = 0;a<nen;a++)
    { 
      for (int b = 0;b<nsd;b++)
      { 
        printf("      lm[%i][%i][%i]=%f \n",b,a,c,lm[b][a][c]);
      }
    }
    printf("      --------------------\n");
  }
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