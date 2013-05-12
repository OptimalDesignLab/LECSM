/**
 * \file setup_eq.cc
 * \determine the local to global equation mapping for FEA
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

#include <stdio.h>
#include "./setup_eq.h"
using namespace std;

// ======================================================================

void setup_eq(int nnp, int nel, int nsd,
             vector< vector<int> > ntyp, vector< vector<double> >FG,
             vector< vector< vector<double> > >& id,
             vector<double>& G, vector<double>& F,
             int& ndof,      int& ndog)
{
  // This procedure is responsible for specifying the
  // equation numbers at each of the possible nodal
  // degrees of freedom and identifying the essential
  // boundary conditions.
  //
  // Inputs:
  //    nnp           - number of mesh nodes
  //    nel           - number of mesh elements
  //    nsd           - number of space dimensions
  //    ntyp          - vector indicating the type of DoF for each
  //                    component at the node
  //    FG            - vector of nodal point forces or nodal
  //                    prescribed essential boundary conditions
  // Outputs:
  //    id            - destination matrix which, for each node,
  //                    stores the type of DoG and the equation number
  //    K(ndof,ndof)  - global stiffness matrix
  //    F(ndof)       - global force vector
  //    G(ndog)       - global prescribed displacement vector
  //    ndof          - number of DoF without BCs
  //    ndog          - number of DoF with non-zero essential BCs
  //
  printf("Allocating global force and prescribed displacement vectors...\n");
  // Loop over the nodes setting the equation or prescribed
  // essential BC number and load the global nodal loads or
  // prescribed essential BCs.
  ndof = 0;
  ndog = 0;
  for (int a = 0; a < nnp; a++)
  {
    for (int i = 0; i < nsd; i++)
    {
      if (ntyp[a][i]==1)       // DoF - possible nodal load
      {
        id[i][a][0] = 1;        // store the node type
        id[i][a][1] = ndof;     // store the equation number
        ndof++;
        F.push_back(FG[a][i]);     // store the nodal load
      }
      else                      // prescribed BC
      {
        if (ntyp[a][i]==2)     // DoG - non-zero BC
        {
          id[i][a][0] = 2;
          id[i][a][1] = ndog;
          ndog++;
          G.push_back(FG[a][i]);   // store the essential BC
        }
        else                    // zero BC (nodeTyle[i] == 0)
        {
          id[i][a][0] = 0;
          id[i][a][1] = 0;
        }
      }
    } // end loop over possible dof at node
  } // end loop over nodes
  printf("DONE\n");
  printf("    Testing id values:\n");
  for (int c = 0; c < 2; c++)
  {
    for (int a = 0; a < nnp; a++)
    { 
      for (int b = 0; b < nsd; b++)
      { 
        printf("      id[%i][%i][%i]=%f \n",b,a,c,id[b][a][c]);
      }
    }
    printf("      --------------------\n");
  }
  printf("DONE\n");
} // end of setup_eq