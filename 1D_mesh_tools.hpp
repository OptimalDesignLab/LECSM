/**
 * \file 1D_mesh_tools.hpp
 * \a 1D mesh library
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

#pragma once
#include <vector>
using namespace std;

// =====================================================================

class Node {
public:
	int id;
	int type[3]; // 0 for zero BC, 1 for free node, 2 for prescribed displacements
	double coords[2];
	double dispBC[3]; // 3 degrees of freedom per node (axial, trans, rotational)
	double forceBC[3]; // forcing in two directions plus moment about the third

	void CreateNode(int num, double* c);
	void DefineBCs(int* BCtype, double* BCval);
};

// =====================================================================

class Element {
public:
	int id;
	int nen;
	vector<Node> adjNodes;

	void CreateElem(int num, vector<Node> nodes);

	void GetStiff(double E, double w, double t, double P,
                  vector< vector< vector<double> > > gm,
                  vector< vector< vector<double> > >& lm, 
                  vector< vector<double> >& KE,
                  vector<double>& FE);

	void assemble(vector< vector< vector<double> > > lm,
                  vector< vector<double> > KE,
                  vector<double> FE, vector<double>& G,
                  vector< vector<double> >& K,
                  vector<double>& F);
};

// =====================================================================

class Mesh {
public:
	int nnp;
	int nel;
	vector<Element> allElems;
	vector<Node> allNodes;

	void CreateMesh(vector<Element> elems);
};