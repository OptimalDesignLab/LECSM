/**
 * \file fea_prog.hpp
 * \FEA program header file
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

#pragma once
#include <vector>
#include "./1D_mesh_tools.hpp"
using namespace std;

// =====================================================================

class LECSM {
public:
	double E, w, t;			// Material properties
	Mesh geometry;		  // Problem geometry and BCs

	void Solve();
};