/**
 * \file fea_prog.hpp
 * \FEA program header file
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

#pragma once
#include <vector>
#include "./1D_mesh_tools.hpp"
#include "../quasi_1d_euler/inner_prod_vector.hpp"
using namespace std;

// =====================================================================

class LECSM {
public:
	double E, w, t;			// Material properties
	Mesh geometry;		  // Problem geometry and BCs
	vector< vector<double> > K; // Global LHS stiffness matrix
	vector<double> F, G; // Global RHS and prescribed BC vector

	void GetStiff();

	void CalcStiffProd(InnerProdVector& u_csm, InnerProdVector& v_csm);

	void Solve();
};