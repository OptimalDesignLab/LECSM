/**
 * \file elem_calcs.cpp
 * \functions that perform element-level calculations
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

#include <stdio.h>
#include "./elem_calcs.h"
using namespace std;

// ======================================================================


void ele_stiff(pMeshEnt elem, int nsd, int nint, int nen, int nee,
               vector< vector<double> > D, vector<double> f_init,
               vector< vector<double> > hOn, vector<double> h,
               vector< vector<double> > intPts, vector<double> weights,
               vector< vector<double> > nodeCoords, vector<int> ien,
               vector< vector< vector<double> > > id,
               vector< vector< vector<double> > >& lm,
               vector< vector<double> >& KE, vector<double>& FE,
               vector< vector< vector<double> > >& Se)
{
  // Main routine for calculating the element matrices.
  //
  // Inputs:
  //    elem           - element to be operated on
  //    nsd            - number of space dimensions
  //    D              - material properties matrix
  //    id             - destination matrix which, for each node,
  //                     stores the type of DoG and the equation number
  // Outputs:
  //    lm (nsd,nen,2) - element location vector identifying the
  //                     equation numbers and location of prescribed
  //                     essential BC for the element
  //    KE (nee,nee)   - element stiffness matrix
  //    FE (nee)       - element load vector

  // Set-up the location vector which will keep track 
  // of the mapping between the local DoF numbers and 
  // the global equation or the specified essential 
  // boundary condition numbers.
  printf("Setting up the location vector...\n");
  for (int a = 0; a < nen; a++)
  {
    int A = ien[a];
    for (int i = 0; i < nsd; i++)
    {
      lm[i][a][0] = id[i][A][0];
      lm[i][a][1] = id[i][A][1];
    }
  }
  printf("DONE\n");
  // Zero out the stiffness matrix and the load vector.
  printf("Zeroing out the stiffness matrix and the load vector...\n");
  for (int p = 0; p < nee; p++)
  {
    FE[p] = 0;
    for (int q = 0; q < nee; q++)
    {
      KE[p][q] = 0;
    }
  }
  printf("DONE\n");
  // Loop over the integration points and add into the
  // appropriate contributions.
  double eta, eps, detJ;
    printf("Looping over integration points...\n");
  for (int i = 0; i < nint; i++)
  {
    printf("   Integration point %d...\n", i);
    printf("       - Weight = %d\n", weights[i]);
    // Calculate element shape functions at the integration point.
    eps = intPts[i][0];
    eta = intPts[i][1];
    vector<double> N(nen);
    vector<double> dNdx(nen);
    vector<double> dNdy(nen);
    Elem_GenShp(elem, eps, eta, nsd, nen, nint, nodeCoords,
                detJ, N, dNdx, dNdy);
    printf("       - detJ = %d\n", detJ);
    // Construct the B matrix (strain-displacement relationships
    vector< vector<double> > B(3, vector<double>(nsd*nen));
    vector< vector<double> > Btrans(nsd*nen, vector<double>(3));
    vector< vector<double> > S(3, vector<double>(nsd*nen));
    vector< vector<double> > BS(nsd*nen, vector<double>(nsd*nen));
    for (int j = 0; j < nen; j++)
    {
      B[0][2*j]   = dNdx[j];
      B[0][2*j+1] = 0;
      B[1][2*j]   = 0;
      B[1][2*j+1] = dNdy[j];
      B[2][2*j]   = dNdy[j];
      B[2][2*j+1] = dNdx[j];
    }
    N.clear();
    dNdx.clear();
    dNdy.clear();
    printf("    B (strain-disp relationships):\n");
    printMatrix(B, 3, nsd*nen);
    // Construct the transpose of B
    matrixTranspose(B, 3, nsd*nen, Btrans);
    printf("    Transpose of B:\n");
    printMatrix(Btrans, nsd*nen, 3);
    // Construct S = D*B
    matrixMult(D, 3, 3, B, 3, nsd*nen, S);  
    printf("    S = D*B:\n");
    printMatrix(S, 3, nsd*nen);
    Se[i] = S;
    // Calculate BS = Btrans*D*B = Btrans*S
    matrixMult(Btrans, nsd*nen, 3, S, 3, nsd*nen, BS);
    printf("    BS = Btrans*D*B = Btrans*S:\n");
    printMatrix(BS, nsd*nen, nsd*nen);
    // Add contributions to the element stiffness matrix.
    for (int a = 0; a < nee; a++)
    {
      for (int b = 0; b < nee; b++)
      {
        KE[a][b] = KE[a][b] + BS[a][b]*weights[i]/detJ;
      }
    }
    B.clear();
    Btrans.clear();
    S.clear();
    BS.clear();
    // printf("    Stiffness contributions at integration point %d\n", nint);
    // printMatrix(KE, nee, nee);
    // Calculate the body force contributions.
    for (int a = 0; a < nen; a++)
    {
      for (int j = 0; j < nsd; j++)
      {
        int p = nsd*a + i;
        for (int i = 0; i < nint; i++)
        {
          eps = intPts[i][0];
          eta = intPts[i][1];
          vector<double> N(nen);
          vector<double> dNdx(nen);
          vector<double> dNdy(nen);
          Elem_GenShp(elem, eps, eta, nsd, nen, nint, nodeCoords,
                      detJ, N, dNdx, dNdy);
          FE[p] = FE[p] + N[a]*weights[i]*detJ*f_init[j];
          N.clear();
          dNdx.clear();
          dNdy.clear();
        }
      }
    }
    // Calculate the edge load contributions.
    for (int a = 0; a < nen; a++)
    {
      for (int n = 0; n < nsd; n++)
      {
        int p = nsd*a + n; 
        int elemID = ien[a];
        if (hOn[elemID][n] == 1)
        {
          // Calculate 1-D integration points.
          int d;
          if ((nen == 3)||(nen == 4))      // linear element
          {
            d = 2;
            vector< vector<double> > int1D(2, vector<double>(2));
            vector<double> w1D(2);
            w1D[0] = 1;
            w1D[1] = 1;
            if (n == 0)      // x-direction
            {
              int1D[0][0] =  0;
              int1D[0][1] = -0.577350269189626;
              int1D[1][0] =  0;
              int1D[1][1] =  0.577350269189626;
            }
            else if (n == 1) // y -direction
            {
              int1D[0][1] =  0;
              int1D[0][0] = -0.577350269189626;
              int1D[1][1] =  0;
              int1D[1][0] =  0.577350269189626;
            }
            // Add the contributions.
            for (int i = 0; i < d; i++)
            {
              eps = int1D[i][0];
              eta = int1D[i][1];
              vector<double> N(nen);
              vector<double> dNdx(nen);
              vector<double> dNdy(nen);
              Elem_GenShp(elem, eps, eta, nsd, nen, nint, nodeCoords,
                          detJ, N, dNdx, dNdy);
              FE[p] = FE[p] + N[a]*w1D[i]*detJ*h[n];
              N.clear();
              dNdx.clear();
              dNdy.clear();
            }
          }
          else if ((nen == 6)||(nen == 8)) // quadratic element
          {
            d = 3;
            vector< vector<double> > int1D(3, vector<double>(2));
            vector<double> w1D(3);
            w1D[0] = 5/9;
            w1D[1] = 8/9;
            w1D[2] = 5/9;
            if (n == 0)      // x-direction
            {
              int1D[0][0] =  1;
              int1D[0][1] =  -0.774596669241483;
              int1D[1][0] =  1;
              int1D[1][1] =  0;
              int1D[2][0] =  1;
              int1D[2][1] =  0.774596669241483;
            }
            else if (n == 1) // y-direction
            {
              int1D[0][1] =  1;
              int1D[0][0] = -0.774596669241483;
              int1D[1][1] =  1;
              int1D[1][0] =  0;
              int1D[2][1] =  1;
              int1D[2][0] =  0.774596669241483;
            }
            // Add the contributions.
            for (int i = 0; i < d; i++)
            {
              eps = int1D[i][0];
              eta = int1D[i][1];
              vector<double> N(nen);
              vector<double> dNdx(nen);
              vector<double> dNdy(nen);
              Elem_GenShp(elem, eps, eta, nsd, nen, nint, nodeCoords,
                          detJ, N, dNdx, dNdy);
              FE[p] = FE[p] + N[a]*w1D[i]*detJ*h[n];
              N.clear();
              dNdx.clear();
              dNdy.clear();
            }
          }
        }
      }
    }
  } // end loop over integration points
  printf("DONE\n");
} // end of ele_stiff
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void assemble(int nsd,           int nen,
              vector< vector< vector<double> > > lm,
              vector< vector<double> > KE,
              vector<double> FE, vector<double>& G,
              vector< vector<double> >& K,
              vector<double>& F)
{
  // This procedure assemples the element contributions into
  // the global matrices.
  //
  // Inputs:
  //    nsd             - number of space dimensions
  //    nen             - number of element nodes
  //    lm              - element location vector, identifying the
  //                      equation numbers and the location of
  //                      essential BCs for the element
  //    KE              - element stiffness matrix
  //    FE              - element force vector
  //    G               - global prescribed displacement vector
  // Outputs:
  //    K               - global stiffness matrix
  //    F               - global force matrix
  //
  printf("    Testing lm values:\n");
  for (int c = 0;c<2;c++)
  {
    for (int a = 0;a<nen;a++)
    { 
      for (int b = 0;b<nsd;b++)
      { 
        printf("      lm[%i][%i][%i]=%f \n",b,a,c,lm[b][a][c]);
      }
    }
    printf("      --------------------\n");
  }
  int p = 0;
  for (int a = 0; a < nen; a++)
  {
    for (int i = 0; i < nsd; i++)
    {
      if (lm[i][a][0] == 1) // dof
      {
        int P = lm[i][a][1];
        F[P] = F[P] + FE[p];
        int q = 0;
        for (int b = 0; b < nen; b++)
        {
          for (int j = 0; j < nsd; j++)
          {
            int Q = lm[j][b][1];
            if (lm[j][b][0] == 1) // dof
              {K[P][Q] = K[P][Q] + KE[p][q];}
            else if (lm[j][b][0] == 2) // dog
              {F[P] = F[P] - G[Q]*KE[p][q];}
            q++;
          } // end j loop over nsd
        } // end b loop over nen (columns)
      }
      p++;
    } // end i loop over nsd
  } // end a loop over nen (rows)
}