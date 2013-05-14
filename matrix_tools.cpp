/**
 * \file matrix_tools.cpp
 * \a set of matrix operation tools
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

#include <stdio.h>
#include <stdlib.h>
#include "./matrix_tools.hpp"
using namespace std;

// ======================================================================

void matrixMult(vector< vector<double> > A, int rowA, int colA,
                vector< vector<double> > B, int rowB, int colB,
                vector< vector<double> >& AB)
{
  // Performs a matrix multiplication for the given matrices
  // A and B.
  //
  // Sanity check on matrix dimensions.
  if (colA != rowB)
  {
    printf("ERROR: Matrix inner dimensions don't match.\
      Multiplication failed!\n");
    exit(EXIT_FAILURE);
  }
  // Pre-allocate the resultant matrix.
  for (int i = 0; i < rowA; i++)
  {
    for (int j = 0; j < colB; j++)
      {AB[i][j] = 0;}
  }
  // Perform the multiplication.
  for (int i = 0; i < rowA; i++)
  {
    for (int j = 0; j < colB; j++)
    {
      for (int k = 0; k < colA; k++)
      {
        AB[i][j] = AB[i][j] + A[i][k]*B[k][j];
      }
    }
  }
}

// ======================================================================

void matrixVecMult(vector< vector<double> > A, int rowA, int colA,
                   vector<double> B, int rowB,
                   vector<double>& AB)
{
  // Performs a matrix-vector multiplication for the given matrix
  // A and given vector B.

  // Sanity check on matrix dimensions.
  if (colA != rowB)
  {
    printf("ERROR: Matrix inner dimensions don't match.\
      Multiplication failed!\n");
    exit(EXIT_FAILURE);
  }
  // Pre-allocate the resultant vector.
  for (int i = 0; i < rowB; i++)
    {AB[i] = 0;}
  // Perform the multiplication.
  for (int i = 0; i < rowA; i++)
  {
    for (int j = 0; j < colA; j++)
    {
      AB[i] = AB[i] + A[i][j]*B[i];
    }
  }
}

// ======================================================================

void matrixTranspose(vector< vector<double> > A, int rowA, int colA,
                     vector< vector<double> >& Atrans)
{
  for (int i = 0; i < rowA; i++)
  {
    for (int j = 0; j < colA; j++)
    {
      Atrans[j][i] = A[i][j];
    }
  }
}

// ======================================================================

void printMatrix(vector< vector<double> > A, int rowA, int colA)
{
  for (int r = 0; r < rowA; r++)
  {
    printf("|  ");
    for (int c = 0; c < colA; c++)
    {
      double val = A[r][c];
      printf("%f  ", val);
      if (c == colA-1)
        {printf("|\n");}
    }
  }
}

// ======================================================================

void CGSolve(vector< vector<double> > K, int rowK, int colK,
             vector<double> F, int rowF, int maxIt, vector<double>& Disp)
{
  // Iterative Conjugate Gradient solver for an Ax=b system.
  //
  // Inputs:
  //    K         - left hand side matrix
  //    rowK      - vertical size of A
  //    colK      - horizontal size of A
  //    F         - right hand side vector
  //    rowF      - vertical size of b
  //    maxIt     - number of maximum iterations
  // Outputs:
  //    Disp      - approximate solution vector
  
  // Sanity check for system dimensions
  if (colK != rowF)
  {
    printf("ERROR: Matrix dimensions don't match. Cannot solve!\n");
    exit(EXIT_FAILURE);
  }
  double tol = 0;
  for(int i = 0; i < rowF ;i++)
  {
    Disp[i] = 0;
  }
  double r[rowF];
  double p[rowF];
  double rtr = 0;
  for(int i = 0; i < rowF; i++)
  {
   r[i] = F[i];
   p[i] = r[i];
   rtr = rtr + r[i]*r[i];
  }
  double alpha=0;
  double Ap[rowF];
  double ptAp=0;
  int iter = 0;
  for(int k = 0; k < maxIt; k++)
  {
    iter++;
    for(int i = 0; i<rowF ;i++)
    {Ap[i] = 0;}
    ptAp = 0;
    for(int i = 0; i < rowF; i++)
    {
      for(int j = 0; j< rowF; j++)
      {
        Ap[i] = Ap[i] + K[i][j]*p[j];
      }
    ptAp = ptAp + p[i]*Ap[i];
    }
    alpha = rtr/(ptAp);
    double rtr_new = 0;
    for(int i = 0; i < rowF; i++)
    {
      Disp[i] = Disp[i]+alpha*p[i];
      r[i] = r[i] - alpha*Ap[i];
      rtr_new = rtr_new + r[i]*r[i];
    }
    if(rtr_new < 1e-6)
    {
      tol = rtr_new;
      break;
    }
    for(int i = 0; i < rowF; i++)
    {
      p[i] = r[i] + rtr_new*p[i]/rtr;
    }
    rtr = rtr_new;
    tol = rtr_new;
  }
}