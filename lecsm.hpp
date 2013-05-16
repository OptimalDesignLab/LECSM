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
	double h;						// Distance to axis of symmetry
	vector<double> P;		// Nodal pressures
	Mesh geometry;		  // Problem geometry and BCs

	void InitGlobalVecs(vector<double>& G, vector<double>& F,
											vector< vector<double> >& K);

	void GetStiff(vector< vector< vector<int> > >& gm,
								vector<double>& G, vector<double>& F,
                vector< vector<double> >& K);

	void Calc_dSdu_Product(InnerProdVector& u_csm, InnerProdVector& v_csm);

	void Calc_dAdu_Product(InnerProdVector& u_csm, InnerProdVector& wrk);

	void Calc_dSdp_Product(InnerProdVector& wrk, InnerProdVector& u_cfd);

	void Solve();
};