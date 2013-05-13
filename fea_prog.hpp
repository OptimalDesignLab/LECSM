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

void FEA(Mesh nozzle, double* props, double P,
         vector< vector<double> > FG);