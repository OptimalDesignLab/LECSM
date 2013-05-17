/**
 * \file output_tools.hpp
 * \brief output tools header file
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

#pragma once
#include <vector>
using namespace std; 

/*!
 * \brief outputs the nodal displacements of a mesh
 * \param[in] nnp - number of mesh nodes
 * \param[in] G - vector of prescribed nodal displacements
 * \param[in] gm - global equation number mapping
 * \param[in] disp - CSM solution vector
 * \param[out] nodeDisp - matrix of all nodal displacements
 */
void output_disp(int nnp, vector<double>& G,
                 vector< vector< vector<int> > >& gm,
                 vector<double>& disp, vector< vector<double> >& nodeDisp);