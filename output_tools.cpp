/**
 * \file output_tools.cpp
 * \a set of tools to output FEA analysis results
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

#include <stdio.h>
#include "./output_tools.hpp"
using namespace std;

// ======================================================================

void output_disp(int nnp, vector<double> G,
                 vector< vector< vector<double> > > id,
                 vector<double> disp, vector< vector<double> >& nodeDisp)
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
  //    nodeDisp          - nodal displacement vector
  printf("Directions:\n");
  printf("  0 - x-axis\n");
  printf("  1 - y-axis\n");
  printf("  2 - rotation about z-axis\n");
  vector<double> nd(2);
  for (int A = 0; A < nnp; A++)
  {
    for (int i = 0; i < 3; i++)
    {
      int t = id[i][A][0];
      double P = id[i][A][1];
      if (t == 1)             // dof
        {nd[i] = disp[P];}
      else
        {
          if (t == 2)          // dog
            {nd[i] = G[P];}
          else
            {nd[i] = 0;}
        }
      printf("    Node %d displaced %f in direction %d\n", A, nd[i], i);
      nodeDisp[A][i] = nd[i];
    }
  }
}