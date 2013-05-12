/**
 * \file 1D_mesh_tools.cpp
 * \a 1D mesh library
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

#include <stdio.h>
#include <math.h>
#include "./1D_mesh_tools.hpp"
using namespace std;

// =====================================================================

Node::CreateNode() {
  id = -1;
  type = -1;
}

Node::CreateNode(int id_, int type_, double coords_) {
  id = id_;
  type = 1; // node is free to move by default
  coords = coords_;
  dispBC = NULL;
  forceBC = {0,0};
}

Node::DefineBCs(int BCtype, double BCval) {
  if (BCtype == 0) { // displacement BC
    if (BCval == NULL) { // zero BC
      type = 0;
      dispBC = {0,0};
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

Element::CreateElem() {
  id = -1;
  nen = 0;
}

Element::CreateElem(int id_, vector<Node> nodes) {
  id = id_;
  nen = nodes.size();
  adjNodes = nodes;
}

// =====================================================================

Mesh::CreateMesh() {
  nnp = 0;
  nel = 0;
}

Mesh::CreateMesh(vector<Element> elems) {
  allElems = elems;
  nel = allElems.size();
  // count number of mesh nodes, assuming 1-D mesh
  nnp = 0;
  Element elem;
  for (int i=0; i<nel; i++) {
    elem = allElems.at(i);
    nnp += elem.nen - 1;
  }
  nnp++;
}

// =====================================================================

void Element::GetStiff(double E, double w, double t, double P, 
                       vector< vector<double> >& KE,
                       vector<double>& FE) {
  nodeL = adjNodes[0];
  nodeR = adjNodes[1];
  double x1 = nodeL.coords[0];
  double y1 = nodeL.coords[1];
  double x2 = nodeR.coords[0];
  double y2 = nodeR.coords[1];
  double len = sqrt(pow(x2-x1,2) + pow(y2-y1,2));
  double c = (x2 - x1)/len;
  double s = (y2 - y1)/len;
  // Create the element stiffness matrix
  double A = w*t;
  double blk[2][2];
  blk[0][0] = E*A*pow(c,2)/len;
  blk[0][1] = E*A*c*s/len;
  blk[1][0] = E*A*c*s/len;
  blk[1][1] = E*A*pow(s,2)/len;
  for (int i=0; i<4; i++) {
    if (i < 2) {
      KE[i][0] = blk[i][0];
      KE[i][1] = blk[i][1];
      KE[i][2] = -blk[i][0];
      KE[i][3] = -blk[i][1];
    }
    else {
      KE[i][0] = -blk[i][0];
      KE[i][1] = -blk[i][1];
      KE[i][2] = blk[i][0];
      KE[i][3] = blk[i][1];
    }
  }
  // Create the element forcing vector due to pressure
  double f_hat = -P[id]*len*w/2;
  double fy = f_hat*c;
  double fx = f_hat*s;
  FE[0] = fx;
  FE[1] = fy;
  FE[2] = fx;
  FE[3] = fy;
  // Add in nodal force contributions
  FE[0] += nodeL.forceBC[0];
  FE[1] += nodeL.forceBC[1];
  FE[2] += nodeR.forceBC[2];
  FE[3] += nodeR.forceBC[3];
}