/**
 * \file setup_eq.cpp
 * \determine the local to global equation mapping for FEA
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

#include <stdio.h>
#include "./setup_eq.hpp"
using namespace std;

// ======================================================================

void setup_eq(Mesh nozzle, vector< vector<double> >FG,
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
  Node nd;
  for (int b = 0; a < nnp; a++)
  {
    nd = mesh.allNodes[b];
    a = nd.id;
    for (int i = 0; i < 3; i++)
    {
      if (nd.type[i]==1)       // DoF - possible nodal load
      {
        id[i][a][0] = 1;        // store the node type
        id[i][a][1] = ndof;     // store the equation number
        ndof++;
        F.push_back(nd.forceBC[i]);     // store the nodal load
      }
      else                      // prescribed BC
      {
        if (nd.type[i]==2)     // DoG - non-zero BC
        {
          id[i][a][0] = 2;
          id[i][a][1] = ndog;
          ndog++;
          G.push_back(nd.dispBC[i]);   // store the essential BC
        }
        else                    // zero BC (nodeTyle[i] == 0)
        {
          id[i][a][0] = 0;
          id[i][a][1] = 0;
        }
      }
    } // end loop over possible dof at node
  } // end loop over nodes
  // printf("DONE\n");
  // printf("    Testing id values:\n");
  // for (int c = 0; c < 2; c++)
  // {
  //   for (int a = 0; a < nnp; a++)
  //   { 
  //     for (int b = 0; b < 3; b++)
  //     { 
  //       printf("      id[%i][%i][%i]=%f \n",b,a,c,id[b][a][c]);
  //     }
  //   }
  //   printf("      --------------------\n");
  // }
  // printf("DONE\n");
}