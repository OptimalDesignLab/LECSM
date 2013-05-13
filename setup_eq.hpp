/**
 * \file setup_eq.cpp
 * \setup equations header file
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

#pragma once
#include <vector>
#include "./1D_mesh_tools.hpp"
using namespace std;

void setup_eq(Mesh nozzle, vector< vector<double> >FG,
             vector< vector< vector<double> > >& id,
             vector<double>& G, vector<double>& F,
             int& ndof,      int& ndog);