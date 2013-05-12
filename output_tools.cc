/**
 * \file output_tools.cpp
 * \a set of tools to output FEA analysis results
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

#include <stdio.h>
#include "./output_tools.h"
using namespace std;

// ======================================================================

void output_disp(int nnp, int nsd, vector<double> G,
                 vector< vector< vector<double> > > id,
                 vector<double> disp, vector<double>& nodeDisp)
{
  // This procedure prints the displacements at nodes.
  //
  // Inputs:
  //    nnp               - number of mesh nodes
  //    nsd               - number of space dimensions
  //    G                 - global prescribed displacement vector
  //    id                - destination matrix which, for each node,
  //                        stores the type of DoG and the equation number
  //    disp              - global displacement vector
  // Outputs:
  //    nodeDisp          - displacement at a particular node
  //
   for (int A = 0; A < nnp; A++)
   {
      for (int i = 0; i < nsd; i++)
      {
         int t = id[i][A][0];
         double P = id[i][A][1];
         if (t == 1)             // dof
            {nodeDisp[i] = disp[P];}
         else
         {
            if (t == 2)          // dog
               {nodeDisp[i] = G[P];}
            else
               {nodeDisp[i] = 0;}
         }
         printf("    Node %d displaced %f in direction %d\n", A, nodeDisp[i], i);
      }
   }
}

// ======================================================================

void out_flux(int nsd, int nen, int nint, vector<int> ien,
              vector< vector< vector<double> > > id,
              vector< vector< vector<double> > > Se,
              vector<double> disp, vector<double> G,
              vector< vector< vector<double> > >& lm,
              vector<double>& SIG)
{
  // This procedure calculates the stresses for each element.
  //
  // Inputs:
  //    nsd               - number of space dimensions
  //    nint              - number of integration points
  //    nen               - number of element nodes
  //    id                - destination matrix which, for each node,
  //                        stores the type of DoG and the equation number
  //    lm                - element location matrix, identifying the
  //                        equation numbers and the location of
  //                        essential BCs for the element
  //    Se                - element stress displacement matrix
  //    disp              - global displacement vector
  //    G                 - global prescribed displacement vector
  // Outputs:
  //    SIG               - element stress vector
  //
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
  vector<double> DE(nsd*nen);
  int p = 0;
  printf("Recovering the element displacement vector...\n");
  // Recover the element displacement vector out of the global one.
  for (int a = 0; a < nen; a++)
  {
    for (int i = 0; i < nsd; i++)
    {
      int P = lm[i][a][1];
      if (lm[i][a][0] == 1)         // dof
        {DE[p] = disp[P];}
      else
      {
        if (lm[i][a][0] == 2)       // dog
          {DE[p] = G[P];}
        else
          {DE[p] = 0;}
      }
      p++;
    }
  }
  printf("DONE\n");
  // Initialize the element stress vector.
  SIG[0] = 0;
  SIG[1] = 0;
  SIG[2] = 0;
  // Loop over integration points, calculate the stresses and add
  // the contributions to the element stress vector.
  printf("Calculating stresses...\n");
  for (int i = 0; i < nint; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      for (int k = 0; k < nsd*nen; k++)
        {SIG[j] = SIG[j] + Se[i][j][k]*DE[k];}
    }
  }
  printf("DONE\n");
}