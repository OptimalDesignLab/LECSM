/**
 * \file linear_elastic_csm.cpp
 * \main CSM program
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
	double E = 10000; // Pascals
	double w = 2;			// meters
	double t = 0.03; 	// meters
	double props = {E, w, t};

	// Create problem mesh
	
	// Call FEA solver
}



