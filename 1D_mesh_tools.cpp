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
  // Calculate the element stiffness matrix
  double x1 = nodeL.coords[0];
  double y1 = nodeL.coords[1];
  double x2 = nodeR.coords[0];
  double y2 = nodeR.coords[1];
  double len = sqrt(pow(x2-x1,2) + pow(y2-y1,2));
  double c = (x2 - x1)/len;
  double s = (y2 - y1)/len;
  double A = w*t;     // cross section area of the element
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