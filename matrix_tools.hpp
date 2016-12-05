/**
 * \file matrix_tools.hpp
 * \brief matrix tools header file
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

#pragma once
#include <vector>
using namespace std;

/*!
 * \brief performs a matrix-matrix multiplication
 * \param[in] A - left matrix
 * \param[in] rowA - row size of A
 * \param[in] colA - column size of A
 * \param[in] B - right matrix
 * \param[in] rowB - row size of B
 * \param[in] colB - column size of B
 * \param[out] AB - resultant matrix
 */
void matrixMult(vector< vector<double> >& A, int rowA, int colA,
                vector< vector<double> >& B, int rowB, int colB,
                vector< vector<double> >& AB);

/*!
 * \brief performs a matrix-vector multiplication
 * \param[in] A - matrix
 * \param[in] rowA - row size of A
 * \param[in] colA - column size of A
 * \param[in] B - vector
 * \param[in] rowB - row size of B
 * \param[out] AB - resultant vector
 */
void matrixVecMult(vector< vector<double> >& A, int rowA, int colA,
                   vector<double>& B, int rowB,
                   vector<double>& AB);

/*!
 * \brief transposes the given matrix
 * \param[in] A - matrix
 * \param[in] rowA - row size of A
 * \param[in] colA - column size of A
 * \param[out] Atrans - transpose of A
 */
void matrixTranspose(vector< vector<double> >& A, int rowA, int colA,
                     vector< vector<double> >& Atrans);

/*!
 * \brief prints the given matrix
 * \param[in] A - matrix
 * \param[in] rowA - row size of A
 * \param[in] colA - column size of A
 */
void printMatrix(vector< vector<double> >& A, int rowA, int colA);

/*!
 * \brief solves a K*disp=F system with conjugate gradient iterations
 * \param[in] K - LHS matrix
 * \param[in] rowK - row size of K
 * \param[in] colK - column size of K
 * \param[in] F - RHS vector
 * \param[in] rowF - row size of F
 * \param[in] maxIt - maximum number of iterations
 * \param[out] Disp - solution vector
 */
int CGSolve(vector< vector<double> >& K, int rowK, int colK,
            vector<double>& F, int rowF, int maxIt, vector<double>& Disp, 
            bool info=false);