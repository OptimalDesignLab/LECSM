/**
 * \file linear_elastic_csm.cpp
 * \CSM program executable
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

#include <stdio.h>
#include "../quasi_1d_euler/inner_prod_vector.hpp"
#include "./lecsm.hpp"
#include "./1D_mesh_tools.hpp"
using namespace std;

// =====================================================================

int main() {

	// Declare the solver
	LECSM csm;

	// Define material properties
	csm.E = 100000000;		// Young's Modulus (Rubber) (Pa)
	csm.w = 2;					// Width of the geometry (meters)
	csm.t = 0.03;				// Thickness of the beam elements (meters)

	// Create problem mesh
	Node node;
	vector<Node> nodes;
	double x[6] = {0.0, 0.5, 1.2, 1.8, 2.5, 3.0};
	double y[6] = {0.0, 0.35, 0.5, 0.5, 0.35, 0.0};
	for (int i=0; i<6; i++) {
		double c[2] = {x[i],y[i]};
		node.CreateNode(i, c);
		csm.P.push_back(20); 	// define nodal pressures
		if ((i==0)||(i==5)) {
			int BCtype[3] = {0,0,1};
			double BCval[3] = {0,0,0};
			node.DefineBCs(BCtype, BCval);
		}
		else {
			int BCtype[3] = {0,1,1};
			double BCval[3] = {0,0,0};
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
	csm.geometry = nozzle;

	// Call FEA solver
	csm.Solve();

	return 0;
}