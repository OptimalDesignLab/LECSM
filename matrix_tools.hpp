/**
 * \file matrix_tools.hpp
 * \matrix tools header file
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

#pragma once
#include <vector>
using namespace std;

void matrixMult(vector< vector<double> > A, int rowA, int colA,
                vector< vector<double> > B, int rowB, int colB,
                vector< vector<double> >& AB);

void matrixVecMult(vector< vector<double> > A, int rowA, int colA,
                   vector<double> B, int rowB,
                   vector<double>& AB);

void matrixTranspose(vector< vector<double> > A, int rowA, int colA,
                     vector< vector<double> >& Atrans);

void printMatrix(vector< vector<double> > A, int rowA, int colA);

void CGSolve(vector< vector<double> > K, int rowK, int colK,
             vector<double> F, int rowF, int maxIt, vector<double>& Disp);