/**
 * \file linear_elastic_csm.cpp
 * \main CSM program
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

#include <stdio.h>
#include "./linear_elastic_csm.hpp"
using namespace std;

// =====================================================================

int main(int argc, char* argv[]) {
   // Create the mesh for the example problem.
   printf("Creating test problem...\n");
   printf("Generating mesh...\n");
   int nnp, nel, nsd;
   int type = 3; // linear quad mesh (large)
   mesh_problem(type, (pGeomMdl)NULL, mesh, nnp, nel, nsd);
   printf("SUCCESS\n");
   // Define the problem specifics.
   vector< vector<double> > D(3, vector<double>(3));
   vector< vector<int> > ntyp(8, vector<int>(2));
   vector< vector<double> > FG(8, vector<double>(2));
   vector< vector<double> > hOn(8, vector<double>(2));
   vector<double> f_init(2);
   vector<double> h(2);
   printf("Preparing problem variables...\n");
   define_problem(type, D, ntyp, FG, f_init, h, hOn);
   printf("SUCCESS\n");
   // Call the FEA Program.
   printf("Starting FEA Program...\n");
   FEA(mesh, (pGeomMdl)NULL, nsd, D, ntyp, FG, f_init, hOn, h);
}



