/**
 * \file 1D_mesh_tools.hpp
 * \a 1D mesh library
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

using namespace std;

// =====================================================================

class Node {
public:
	int id;
	int type[2];
	double coords[2];
	double dispBC[2];
	double forceBC[2];

	CreateNode();
	CreateNode(int id_, int type_, double coords_);
	DefineBCs(int BCtype, double BCval);
};

// =====================================================================

class Element {
public:
	int id;
	int nen;
	vector<Node> adjNodes;

	CreateElem();

	CreateElem(int id_, vector<Node> nodes);

	void GetStiff(double E, double w, double t, double P, 
								vector< vector<double> >& KE,
				  			vector<double>& FE);
};

// =====================================================================

class Mesh {
public:
	int nnp;
	int nel;
	vector<Element> allElems;
	vector<Node> allNodes;

	CreateMesh();
	CreateMesh(vector<Element> elems);
};