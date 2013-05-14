/**
 * \file linear_elastic_csm.cpp
 * \CSM program executable
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

#include <stdio.h>
#include "./fea_prog.hpp"
#include "./1D_mesh_tools.hpp"
using namespace std;

// =====================================================================

int main() {
	// Define material properties
	double E = 10000000; 		// Pascals (Rubber)
	double w = 2;						// meters
	double t = 0.03; 				// meters
	double props[3] = {E, w, t};

	// Create problem mesh
	Node node;
	vector<Node> nodes;
	double x[6] = {0.0, 0.5, 1.2, 1.8, 2.5, 3.0};
	double y[6] = {0.0, 0.35, 0.5, 0.5, 0.35, 0.0};
	for (int i=0; i<6; i++) {
		double c[2] = {x[i],y[i]};
		node.CreateNode(i, c);
		if ((i==0)||(i==5)) {
			int BCtype[3] = {0,0,NULL};
			double BCval[3] = {0,0,NULL};
			node.DefineBCs(BCtype, BCval);
		}
		nodes.push_back(node);
	}
	Element elem;
	vector<Element> elems;
	Node nodeL;
	Node nodeR;
	vector<Node> elemNodes(2);
	for (int j=0; j<5; j++) {
		nodeL = nodes[j];
		nodeR = nodes[j+1];
		elemNodes[0] = nodeL;
		elemNodes[1] = nodeR;
		elem.CreateElem(j, elemNodes);
		elems.push_back(elem);
	}
	Mesh nozzle;
	nozzle.CreateMesh(elems);
	nodes.clear();
	elems.clear();

	// Create sample pressure vector
	double P[5] = {20,20,20,20,20}; // Pascals

	// Call FEA solver
	FEA(nozzle, props, P);

	return 0;
}