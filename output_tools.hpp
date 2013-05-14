/**
 * \file output_tools.hpp
 * \output tools header file
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

#pragma once
#include <vector>
using namespace std; 

void output_disp(int nnp, vector<double> G,
                 vector< vector< vector<double> > > id,
                 vector<double> disp, vector< vector<double> >& nodeDisp);