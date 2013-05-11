#include <stdio.h>
#include <mpi.h>
#include "FMDB.h"
using namespace std;
//#include "TopoModel.h"

int loadMesh(pMeshMdl& mesh, pGeomMdl geom, const char* meshFileName) {
   FMDB_Mesh_Create(geom, mesh);

   if (FMDB_Mesh_LoadFromFile(mesh, meshFileName, SCUTIL_CommSize() - 1)) {
      fprintf(stderr,"ERROR [%d]: FMDB_Mesh_LoadFromFile failed to load %s ... exiting\n", SCUTIL_CommRank(), meshFileName);
      FMDB_Mesh_Del(mesh);
      SCUTIL_Finalize();
      return SCUtil_FAILURE;      
   }
/*
   int isValid = 0;
   FMDB_Mesh_Verify(mesh, &isValid);
   if ( 0 == isValid ) {
      fprintf(stderr,"ERROR [%d]: Mesh is invalid ... exiting\n", SCUTIL_CommRank());
      FMDB_Mesh_Del(mesh);
      SCUTIL_Finalize();
      return SCUtil_FAILURE;
   }
*/
   return SCUtil_SUCCESS;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void define_problem(int type,
                    vector< vector<double> >& D,
                    vector< vector<int> >& ntyp,
                    vector< vector<double> >& FG,
                    vector<double>& f_init,
                    vector<double>& h,
                    vector< vector<double> >& hOn)
{
  // Creates the example problem.
  //
  // Inputs:
  //    geom          - pGeomMdl geometry object for the domain
  // Outputs:
  //    mesh          - pMeshMdl mesh object created for the problem
  //    nsd           - number of space dimensions (1D, 2D or 3D)
  //    D             - material properties matrix
  //    ntyp          - vector indicating the type of DoF for each
  //                    component at the node
  //    FG            - vector of nodal point forces or nodal
  //                    prescribed essential boundary conditions
  //    f_init        - distributed forces on the geometry
  //
  double E = 10000000;                   // Young's Modulus
  double v = 0.33;                       // Poisson's Ratio
  double lambda = v*E/((1+v)*(1-2*v));
  double mu = E/(1+v);
  // Construct the material properties matrix
  printf("Constructing material props matrix...\n");
  D[0][0] = lambda + 2*mu;
  D[0][1] = lambda;
  D[0][2] = 0;
  D[1][0] = lambda;
  D[1][1] = lambda + 2*mu;
  D[1][2] = 0;
  D[2][0] = 0;
  D[2][1] = 0;
  D[2][2] = mu;
  printf("DONE\n");
  if (type == 0)        // linear quadrilateral
  {
    // Define types of DoF at each node
    // 0 for zero BC (fixed point), 1 for DoF, 2 for non-zero BC
    printf("Defining types of DoF at each node...\n");
    ntyp[0][0] = 0;
    ntyp[0][1] = 0;
    ntyp[1][0] = 1;
    ntyp[1][1] = 1;
    ntyp[2][0] = 1;
    ntyp[2][1] = 1;
    ntyp[3][0] = 1;
    ntyp[3][0] = 1;
    printf("DONE\n");
    // Define FG for the problem
    printf("Prescribing nodal loads...\n");
    FG[0][0] = 0;
    FG[0][1] = 0;
    FG[1][0] = 0;
    FG[1][1] = 0;
    FG[2][0] = 0;
    FG[2][1] = 0;
    FG[3][0] = 0;
    FG[3][1] = 0; 
    printf("DONE\n");
    // Define body forces, if any
    printf("Defining body forces...\n");
    f_init[0] = 0;
    f_init[1] = 0;
    printf("DONE\n");
    // Define edge loads, if any
    printf("Defining edge loads...\n");
    h[0] = 0;
    h[1] = 0;
    hOn[0][0] = 0;
    hOn[0][1] = 0;
    hOn[1][0] = 0;
    hOn[1][1] = 0;
    hOn[2][0] = 1;
    hOn[2][1] = 0;
    hOn[3][0] = 1;
    hOn[3][1] = 0;
    printf("DONE\n");
  }
  else if (type == 1)        // quadratic quadrilateral
  {
    // Define types of DoF at each node
    // 0 for zero BC (fixed point), 1 for DoF, 2 for non-zero BC
    printf("Defining types of DoF at each node...\n");
    ntyp[0][0] = 0;
    ntyp[0][1] = 0;
    ntyp[1][0] = 1;
    ntyp[1][1] = 1;
    ntyp[2][0] = 1;
    ntyp[2][1] = 1;
    ntyp[3][0] = 0;
    ntyp[3][1] = 0;
    ntyp[4][0] = 1;
    ntyp[4][1] = 1;
    ntyp[5][0] = 1;
    ntyp[5][1] = 1;
    ntyp[6][0] = 1;
    ntyp[6][1] = 1;
    ntyp[7][0] = 0;
    ntyp[7][1] = 0;
    printf("DONE\n");
    // Define FG for the problem
    printf("Prescribing nodal loads...\n");
    FG[0][0] = 0;
    FG[0][1] = 0;
    FG[1][0] = 0;
    FG[1][1] = 0;
    FG[2][0] = 0;
    FG[2][1] = 10000;
    FG[3][0] = 0;
    FG[3][1] = 0; 
    FG[4][0] = 0;
    FG[4][1] = 0;
    FG[5][0] = 0;
    FG[5][1] = 0;
    FG[6][0] = 0;
    FG[6][1] = 0;
    FG[7][0] = 0;
    FG[7][1] = 0; 
    printf("DONE\n");
    // Define body forces, if any
    printf("Defining body forces...\n");
    f_init[0] = 0;
    f_init[1] = 0;
    printf("DONE\n");
    // Define edge loads, if any
    printf("Defining edge loads...\n");
    h[0] = 0;
    h[1] = 0;
    hOn[0][0] = 0;
    hOn[0][1] = 0;
    hOn[1][0] = 0;
    hOn[1][1] = 0;
    hOn[2][0] = 0;
    hOn[2][1] = 0;
    hOn[3][0] = 0;
    hOn[3][1] = 0;
    hOn[0][0] = 0;
    hOn[0][1] = 0;
    hOn[1][0] = 0;
    hOn[1][1] = 0;
    hOn[2][0] = 0;
    hOn[2][1] = 0;
    hOn[3][0] = 0;
    hOn[3][1] = 0;
    printf("DONE\n");
  }
  else if (type == 2)        // linear triangular
  {
    // Define types of DoF at each node
    // 0 for zero BC (fixed point), 1 for DoF, 2 for non-zero BC
    printf("Defining types of DoF at each node...\n");
    ntyp[0][0] = 0;
    ntyp[0][1] = 0;
    ntyp[1][0] = 1;
    ntyp[1][1] = 1;
    ntyp[2][0] = 0;
    ntyp[2][1] = 0;
    ntyp[3][0] = 1;
    ntyp[3][1] = 1;
    printf("DONE\n");
    // Define FG for the problem
    printf("Prescribing nodal loads...\n");
    FG[0][0] = 0;
    FG[0][1] = 0;
    FG[1][0] = 0;
    FG[1][1] = 10000;
    FG[2][0] = 0;
    FG[2][1] = 0;
    FG[3][0] = 0;
    FG[3][1] = 0;
    printf("DONE\n");
    // Define body forces, if any
    printf("Defining body forces...\n");
    f_init[0] = 0;
    f_init[1] = 0;
    printf("DONE\n");
    // Define edge loads, if any
    printf("Defining edge loads...\n");
    h[0] = 0;
    h[1] = 0;
    hOn[0][0] = 0;
    hOn[0][1] = 0;
    hOn[1][0] = 0;
    hOn[1][1] = 0;
    hOn[2][0] = 0;
    hOn[2][1] = 0;
    hOn[3][0] = 0;
    hOn[3][1] = 0;
    printf("DONE\n");
  }
  else if (type == 3)        // quadratic triangular
  {
    // Define types of DoF at each node
    // 0 for zero BC (fixed point), 1 for DoF, 2 for non-zero BC
    printf("Defining types of DoF at each node...\n");
    ntyp[0][0] = 0;
    ntyp[0][1] = 0;
    ntyp[1][0] = 1;
    ntyp[1][1] = 1;
    ntyp[2][0] = 0;
    ntyp[2][1] = 0;
    ntyp[3][0] = 1;
    ntyp[3][1] = 1;
    ntyp[4][0] = 1;
    ntyp[4][1] = 1;
    ntyp[5][0] = 0;
    ntyp[5][1] = 0;
    printf("DONE\n");
    // Define FG for the problem
    printf("Prescribing nodal loads...\n");
    FG[0][0] = 0;
    FG[0][1] = 0;
    FG[1][0] = 0;
    FG[1][1] = 10000;
    FG[2][0] = 0;
    FG[2][1] = 0;
    FG[3][0] = 0;
    FG[3][1] = 0;
    FG[4][0] = 0;
    FG[4][1] = 0;
    FG[5][0] = 0;
    FG[5][1] = 0;
    printf("DONE\n");
    // Define body forces, if any
    printf("Defining body forces...\n");
    f_init[0] = 0;
    f_init[1] = 0;
    printf("DONE\n");
    // Define edge loads, if any
    printf("Defining edge loads...\n");
    h[0] = 0;
    h[1] = 0;
    hOn[0][0] = 0;
    hOn[0][1] = 0;
    hOn[1][0] = 0;
    hOn[1][1] = 0;
    hOn[2][0] = 0;
    hOn[2][1] = 0;
    hOn[3][0] = 0;
    hOn[3][1] = 0;
    hOn[4][0] = 0;
    hOn[4][1] = 0;
    hOn[5][0] = 0;
    hOn[5][1] = 0;
    printf("DONE\n");
  }
  else if (type == 4)     // larger linear quadratic problem
  {
    // Define types of DoF at each node
    // 0 for zero BC (fixed point), 1 for DoF, 2 for non-zero BC
    printf("Defining types of DoF at each node...\n");
    ntyp[0][0] = 0;
    ntyp[0][1] = 0;
    ntyp[1][0] = 1;
    ntyp[1][1] = 1;
    ntyp[2][0] = 1;
    ntyp[2][1] = 1;
    ntyp[3][0] = 2;
    ntyp[3][1] = 1;
    ntyp[4][0] = 0;
    ntyp[4][1] = 0;
    ntyp[5][0] = 1;
    ntyp[5][1] = 1;
    ntyp[6][0] = 1;
    ntyp[6][1] = 1;
    ntyp[7][0] = 2;
    ntyp[7][1] = 1;
    printf("DONE\n");
    // Define FG for the problem
    printf("Prescribing nodal loads...\n");
    FG[0][0] = 0;
    FG[0][1] = 0;
    FG[1][0] = 0;
    FG[1][1] = 0;
    FG[2][0] = 0;
    FG[2][1] = 0;
    FG[3][0] = 1.5;
    FG[3][1] = 0; 
    FG[4][0] = 0;
    FG[4][1] = 0;
    FG[5][0] = 0;
    FG[5][1] = 0;
    FG[6][0] = 0;
    FG[6][1] = 0;
    FG[7][0] = 1.5;
    FG[7][1] = 0;
    printf("DONE\n");
    // Define body forces, if any
    printf("Defining body forces...\n");
    f_init[0] = 0;
    f_init[1] = 0;
    printf("DONE\n");
    // Define edge loads, if any
    printf("Defining edge loads...\n");
    h[0] = 0;
    h[1] = -1000;
    hOn[0][0] = 0;
    hOn[0][1] = 1;
    hOn[1][0] = 0;
    hOn[1][1] = 1;
    hOn[2][0] = 0;
    hOn[2][1] = 1;
    hOn[3][0] = 0;
    hOn[3][1] = 1;
    hOn[4][0] = 0;
    hOn[4][1] = 0;
    hOn[5][0] = 0;
    hOn[5][1] = 0;
    hOn[6][0] = 0;
    hOn[6][1] = 0;
    hOn[7][0] = 0;
    hOn[7][1] = 0;
    printf("DONE\n");
  }
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void mesh_problem(int type, pGeomMdl geom, pMeshMdl& mesh,
                  int& nnp, int& nel, int& nsd)
{
  FMDB_Mesh_Create(geom, mesh);
  pPart part;
  FMDB_Mesh_GetPart(mesh, 0, part);
  nsd = 2;
  if (type == 0)    // linear quadrilateral mesh
  {
    nnp = 4;
    nel = 1;
    printf("Creating vertices...\n");
    double params[] = {0, 0, 0};
    double coords1[] = {-1, -1, 0};
    double coords2[] = {1, -1, 0};
    double coords3[] = {1, 1, 0};
    double coords4[] = {-1, 1, 0};
    pMeshEnt* new_vtx = new pMeshEnt[4];
    FMDB_Vtx_Create(part, (pGeomEnt)NULL, coords1, params, new_vtx[0]);
    FMDB_Ent_SetID(new_vtx[0], 0);
    FMDB_Vtx_Create(part, (pGeomEnt)NULL, coords2, params, new_vtx[1]);
    FMDB_Ent_SetID(new_vtx[1], 1);
    FMDB_Vtx_Create(part, (pGeomEnt)NULL, coords3, params, new_vtx[2]);
    FMDB_Ent_SetID(new_vtx[2], 2);
    FMDB_Vtx_Create(part, (pGeomEnt)NULL, coords4, params, new_vtx[3]);
    FMDB_Ent_SetID(new_vtx[3], 3);
    printf("DONE\n");
    printf("Creating elements/faces...\n");
    pMeshEnt new_face[1];
    FMDB_Face_Create(part, (pGeomEnt)NULL, FMDB_QUAD, new_vtx, 0, new_face[0]);
    FMDB_Ent_SetID(new_face[0], 0);
    printf("MESH COMPLETE\n");
  }
  else if (type == 1)    // quadratic quadrilateral mesh
  {
    nnp = 8;
    nel = 1;
    printf("Creating vertices...\n");
    double params[] = {0, 0, 0};
    double coords1[] = {-1, -1, 0};
    double coords2[] = {1, -1, 0};
    double coords3[] = {1, 1, 0};
    double coords4[] = {-1, 1, 0};
    pMeshEnt* new_vtx = new pMeshEnt[4];
    FMDB_Vtx_Create(part, (pGeomEnt)NULL, coords1, params, new_vtx[0]);
    FMDB_Ent_SetID(new_vtx[0], 0);
    FMDB_Vtx_Create(part, (pGeomEnt)NULL, coords2, params, new_vtx[1]);
    FMDB_Ent_SetID(new_vtx[1], 1);
    FMDB_Vtx_Create(part, (pGeomEnt)NULL, coords3, params, new_vtx[2]);
    FMDB_Ent_SetID(new_vtx[2], 2);
    FMDB_Vtx_Create(part, (pGeomEnt)NULL, coords4, params, new_vtx[3]);
    FMDB_Ent_SetID(new_vtx[3], 3);
    printf("DONE\n");
    printf("Creating elements/faces...\n");
    pMeshEnt new_face[1];
    FMDB_Face_Create(part, (pGeomEnt)NULL, FMDB_QUAD, new_vtx, 0, new_face[0]);
    FMDB_Ent_SetID(new_face[0], 0);
    printf("DONE\n");
    printf("Creating edge nodes...\n");
    pPartEntIter iter;
    FMDB_PartEntIter_Init(part, FMDB_EDGE, FMDB_ALLTOPO, iter);
    pMeshEnt edge;
    int nodeID = 4;
    while (SCUtil_SUCCESS == FMDB_PartEntIter_GetNext(iter, edge))
    {
      vector<pMeshEnt> edgeVtx;
      FMDB_Ent_GetAdj(edge, FMDB_VERTEX, 1, edgeVtx);
      double coordsA[3];
      double coordsB[3];
      FMDB_Vtx_GetCoord(edgeVtx[0], coordsA);
      FMDB_Vtx_GetCoord(edgeVtx[1], coordsB);
      double x = (coordsA[0]+coordsB[0])/2;
      double y = (coordsA[1]+coordsB[1])/2;
      double coordsC[] = {x, y, 0};
      FMDB_Ent_SetID(edge, nodeID);
      pNode node;
      FMDB_Edge_CreateNode(edge, coordsC, node);
      nodeID++;
    }
    FMDB_PartEntIter_Del(iter);
    printf("DONE\n");
  }
  else if (type == 2)    // linear triangular mesh
  {
    nnp = 3;
    nel = 2;
    printf("Creating vertices...\n");
    double params[] = {0, 0, 0};
    double coords1[] = {0, 0, 0};
    double coords2[] = {1, 1, 0};
    double coords3[] = {0, 1, 0};
    double coords4[] = {1, 0, 0};
    pMeshEnt* new_vtx = new pMeshEnt[4];
    FMDB_Vtx_Create(part, (pGeomEnt)NULL, coords1, params, new_vtx[0]);
    FMDB_Ent_SetID(new_vtx[0], 0);
    FMDB_Vtx_Create(part, (pGeomEnt)NULL, coords2, params, new_vtx[1]);
    FMDB_Ent_SetID(new_vtx[1], 1);
    FMDB_Vtx_Create(part, (pGeomEnt)NULL, coords3, params, new_vtx[2]);
    FMDB_Ent_SetID(new_vtx[2], 2);
    FMDB_Vtx_Create(part, (pGeomEnt)NULL, coords4, params, new_vtx[3]);
    FMDB_Ent_SetID(new_vtx[3], 3);
    printf("DONE\n");
    printf("Creating elements/faces...\n");
    pMeshEnt* new_face1 = new pMeshEnt[3];
    new_face1[0] = new_vtx[0];
    new_face1[1] = new_vtx[1];
    new_face1[2] = new_vtx[2];
    pMeshEnt* new_face2 = new pMeshEnt[3];
    new_face2[0] = new_vtx[0];
    new_face2[1] = new_vtx[1];
    new_face2[2] = new_vtx[3];
    pMeshEnt* new_face = new pMeshEnt[2];
    FMDB_Face_Create(part, (pGeomEnt)NULL, FMDB_TRI, new_face1, 0, new_face[0]);
    FMDB_Ent_SetID(new_face[0], 0);
    FMDB_Face_Create(part, (pGeomEnt)NULL, FMDB_TRI, new_face2, 0, new_face[1]);
    FMDB_Ent_SetID(new_face[1], 1);
    printf("DONE\n");
  }
  else if (type == 3)    // quadratic triangular mesh
  {
    nnp = 6;
    nel = 1;
    printf("Creating vertices...\n");
    double params[] = {0, 0, 0};
    double coords1[] = {0, 0, 0};
    double coords2[] = {1, 0, 0};
    double coords3[] = {0, 1, 0};
    pMeshEnt* new_vtx = new pMeshEnt[3];
    FMDB_Vtx_Create(part, (pGeomEnt)NULL, coords1, params, new_vtx[0]);
    FMDB_Ent_SetID(new_vtx[0], 0);
    FMDB_Vtx_Create(part, (pGeomEnt)NULL, coords2, params, new_vtx[1]);
    FMDB_Ent_SetID(new_vtx[1], 1);
    FMDB_Vtx_Create(part, (pGeomEnt)NULL, coords3, params, new_vtx[2]);
    FMDB_Ent_SetID(new_vtx[2], 2);
    printf("DONE\n");
    printf("Creating elements/faces...\n");
    pMeshEnt new_face[1];
    FMDB_Face_Create(part, (pGeomEnt)NULL, FMDB_TRI, new_vtx, 0, new_face[0]);
    FMDB_Ent_SetID(new_face[0], 0);
    printf("DONE\n");
    printf("Creating edge nodes...\n");
    pPartEntIter iter;
    FMDB_PartEntIter_Init(part, FMDB_EDGE, FMDB_ALLTOPO, iter);
    pMeshEnt edge;
    int nodeID = 3;
    while (SCUtil_SUCCESS == FMDB_PartEntIter_GetNext(iter, edge))
    {
      vector<pMeshEnt> edgeVtx;
      FMDB_Ent_GetAdj(edge, FMDB_VERTEX, 1, edgeVtx);
      double coordsA[3];
      double coordsB[3];
      FMDB_Vtx_GetCoord(edgeVtx[0], coordsA);
      FMDB_Vtx_GetCoord(edgeVtx[1], coordsB);
      double x = (coordsA[0]+coordsB[0])/2;
      double y = (coordsA[1]+coordsB[1])/2;
      double coordsC[] = {x, y, 0};
      FMDB_Ent_SetID(edge, nodeID);
      pNode node;
      FMDB_Edge_CreateNode(edge, coordsC, node);
      nodeID++;
    }
    FMDB_PartEntIter_Del(iter);
    printf("DONE\n");
  }
  else if (type == 4)     // large linear quad problem
  {
    nnp = 8;
    nel = 3;
    // Calculate the coordinates. The coordinate system is centered
    // at node 5, bottom-left of the geometry. Create the vertices.
    printf("Creating vertices...\n");
    double params[] = {0, 0, 0};
    pMeshEnt* new_vtx = new pMeshEnt[8];
    double coords1[] = {0, 4, 0};
    FMDB_Vtx_Create(part, (pGeomEnt)NULL, coords1, params, new_vtx[0]);
    FMDB_Ent_SetID(new_vtx[0], 0);
    double coords2[] = {4, 4, 0};
    FMDB_Vtx_Create(part, (pGeomEnt)NULL, coords2, params, new_vtx[1]);
    FMDB_Ent_SetID(new_vtx[1], 1);
    double coords3[] = {8, 4, 0};
    FMDB_Vtx_Create(part, (pGeomEnt)NULL, coords3, params, new_vtx[2]);
    FMDB_Ent_SetID(new_vtx[2], 2);
    double coords4[] = {12, 4, 0};
    FMDB_Vtx_Create(part, (pGeomEnt)NULL, coords4, params, new_vtx[3]);
    FMDB_Ent_SetID(new_vtx[3], 3);
    double coords5[] = {0, 0, 0};
    FMDB_Vtx_Create(part, (pGeomEnt)NULL, coords5, params, new_vtx[4]);
    FMDB_Ent_SetID(new_vtx[4], 4);
    double coords6[] = {4, 0, 0};
    FMDB_Vtx_Create(part, (pGeomEnt)NULL, coords6, params, new_vtx[5]);
    FMDB_Ent_SetID(new_vtx[5], 5);
    double coords7[] = {8, 0, 0};
    FMDB_Vtx_Create(part, (pGeomEnt)NULL, coords7, params, new_vtx[6]);
    FMDB_Ent_SetID(new_vtx[6], 6);
    double coords8[] = {12, 0, 0};
    FMDB_Vtx_Create(part, (pGeomEnt)NULL, coords8, params, new_vtx[7]);
    FMDB_Ent_SetID(new_vtx[7], 7);
    printf("DONE\n");
    // Manually create the 4 elements/faces.
    printf("Creating elements/faces...\n");
    pMeshEnt elems[3];
    pMeshEnt face1vtx[4] = {new_vtx[0],new_vtx[4],new_vtx[5],new_vtx[1]};
    FMDB_Face_Create(part, (pGeomEnt)NULL, FMDB_QUAD, face1vtx, 0, elems[0]);
    FMDB_Ent_SetID(elems[0],0);
    printf("Element 1 DONE\n");
    pMeshEnt face2vtx[4] = {new_vtx[1],new_vtx[5],new_vtx[6],new_vtx[2]};
    FMDB_Face_Create(part, (pGeomEnt)NULL, FMDB_QUAD, face2vtx, 0, elems[1]);
    FMDB_Ent_SetID(elems[1],1);
    printf("Element 2 DONE\n");
    pMeshEnt face3vtx[4] = {new_vtx[2],new_vtx[6],new_vtx[7],new_vtx[3]};
    FMDB_Face_Create(part, (pGeomEnt)NULL, FMDB_QUAD, face3vtx, 0, elems[2]);
    FMDB_Ent_SetID(elems[2],2);
    printf("Element 3 DONE\n");
  }
  printf("MESH COMPLETE\n");
  // Print mesh statistics.
  int numEdge;
  FMDB_Part_GetNumEnt(part, FMDB_EDGE, FMDB_ALLTOPO, &numEdge);
  int numVtx;
  FMDB_Part_GetNumEnt(part, FMDB_VERTEX, FMDB_ALLTOPO, &numVtx);
  int numFace;
  FMDB_Part_GetNumEnt(part, FMDB_FACE, FMDB_ALLTOPO, &numFace);
  printf("    Mesh Properties:\n");
  printf("      - Number of edges: %d\n",numEdge);
  printf("      - Number of nodes: %d\n",numVtx);
  printf("      - Number of faces: %d\n",numFace);
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void process_mesh(pMeshMdl mesh, pGeomMdl geom, int nsd,
                  pPart& part, int& nnp, int& nel)
{
  // Analyzes a mesh in order to determine the number of nodes
  // and the number of elements. Faces are assumed to be elements
  // in 2-D meshes. Regions are assumed to be elements in 3-D.
  FMDB_Mesh_GetPart(mesh, 0, part);
  if (nsd == 2)         // 2-D mesh
  {
    FMDB_Part_GetNumEnt(part, FMDB_FACE, FMDB_ALLTOPO, &nel);
    FMDB_Part_GetNumEnt(part, FMDB_VERTEX, FMDB_ALLTOPO, &nnp);
    pPartEntIter iter;
    FMDB_PartEntIter_Init(part, FMDB_EDGE, FMDB_ALLTOPO, iter);
    pMeshEnt edge;
    int iterEnd = FMDB_PartEntIter_GetNext(iter, edge);
    while (!iterEnd)
    {
      int numNode;
      FMDB_Edge_GetNumNode(edge, &numNode);
      if (numNode > 0)
        {nnp = nnp + numNode;}
      iterEnd = FMDB_PartEntIter_GetNext(iter, edge);
    }
    FMDB_PartEntIter_Del(iter);
  }
  else
  {
    printf("ERROR: Unsupported mesh dimension!\n");
    exit(EXIT_FAILURE);
  }
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
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
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
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
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void matrixVecMult(vector< vector<double> > A, int rowA, int colA,
                   vector<double> B, int rowB,
                   vector<double>& AB)
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
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void matrixTranspose(vector< vector<double> > A, int rowA, int colA,
                     vector< vector<double> >& Atrans)
{
  // Transposes the given matrix A.
  //
  // Pre-allocate the transpose.
  for (int i = 0; i < rowA; i++)
  {
    for (int j = 0; j < colA; j++)
    {
      Atrans[j][i] = A[i][j];
    }
  }
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void printMatrix(vector< vector<double> > A, int rowA, int colA)
{
  // Given a matrix A, displays the matrix values in a
  // formatted way.
  //
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
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void directSolve(vector< vector<double> > Ainv, int rowA, int colA,
                 vector<double> b,              int rowB,
                 vector<double>& x)
{
  // Performs a direct solution to the equation Ax = B, given
  // the right hand side vector b and the inverse of left hand side
  // matrix A.
  //
  // Sanity check on matrix dimensions.
  if (colA != rowB)
  {
    printf("ERROR: Matrix inner dimensions don't match.\
      Direct solve failed!\n");
    exit(EXIT_FAILURE);
  }
  // Pre-allocate the solution vector.
  for (int i = 0; i < rowB; i++)
  {
    x[i] = 0;
  }
  // Carry out the direct-solve.
  for (int i = 0; i < rowB; i++)
  {
    for (int j = 0; j < colA; j++)
    {
      x[i] += Ainv[i][j]*b[j];
    }
  }
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void CGSolve(vector< vector<double> > K, int rowK, int colK,
             vector<double> F, int rowF, int maxIt, vector<double>& Disp)
{
  // Iterative Conjugate Gradient solver for an Ax=b system.
  //
  // Inputs:
  //    A         - left hand side matrix
  //    rowA      - vertical size of A
  //    colA      - horizontal size of A
  //    b         - right hand side vector
  //    rowb      - vertical size of b
  //    maxIt     - number of maximum iterations
  // Outputs:
  //    x         - approximate solution vector
  //
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
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void Elem_CheckType(pMeshEnt elem, int& nint, int& nen)
{
  // Check whether a given element is triangular or square, and
  // determine the order of the edges.
  //
  // Input:
  //    elem        - pMeshEnt element object
  // Outputs:
  //    type        - element type (0=tri, 1=rect, 2=unknown)
  //    quad        - edge order (0=linear, 1=quadratic, 2=higher)
  //
  int numEdge;
  FMDB_Ent_GetNumAdj(elem, FMDB_EDGE, &numEdge);
  printf("    Number of edges: %d.\n",numEdge);
  vector<pMeshEnt> adjEdges;
  FMDB_Ent_GetAdj(elem, FMDB_EDGE, 1, adjEdges);
  int type, quad;
  if (numEdge == 3)
  {
    type = 0;
  }
  else if (numEdge == 4)
  {
    type = 1;
  }
  else
    {type = 2;} // polygonal element (NOT SUPPORTED)
  quad = 0;     // assume linear element at first
  for (int i = 0; i < adjEdges.size(); i++)
  {
    int numNode;
    FMDB_Edge_GetNumNode(adjEdges[i], &numNode);
    if (numNode > 0)
    {
      quad = 1;
      break;
    }
  }
  // Figure out the number of integration points and element nodes.
  if (type == 0) // tri
    {
      printf("    Type: Triangular\n");
      if (quad == 0) // linear
      {
        printf("    Order: Linear\n");
        nint = 1;
        nen = 3;
      }
      else if (quad == 1) //quadratic
      {
        printf("    Order: Quadratic\n");
        nint = 3;
        nen = 6;
      }
      else
      {
        printf("ERROR: Unsupported high-order element!\n");
        exit(EXIT_FAILURE);
      }
    }
    else if (type == 1) // quad
    {
      printf("    Type: Quadrilateral\n");
      if (quad == 0) // linear
      {
        printf("    Order: Linear\n");
        nint = 4;
        nen = 4;
      }
      else if (quad == 1) //quadratic
      {
        printf("    Order: Quadratic\n");
        nint = 9;
        nen = 8;
      }
      else
      {
        printf("ERROR: Unsupported high-order element!\n");
        exit(EXIT_FAILURE);
      }
    }
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void Elem_GetInfo(pMeshEnt elem, int nsd, int nen,
                  int nint, int nee, vector<int>& ien,
                  vector<double>& weights,
                  vector< vector<double> >& intPts,
                  vector< vector<double> >& nodeCoords)
{
  // Retreive element and nodal coordinate information to be used 
  // in calculating the element stiffness matrices.
  //
  // Inputs:
  //    elem            - pMeshEnt object for the element
  //    nsd             - number of space dimensions
  // Outputs:
  //    nint            - number of integration points
  //    nen             - number of nodes for the element
  //    nee             - size of the element stiffness matrix
  //    ien(nen)        - vector of global node IDs for nodes
  //                      that define the element
  //    weights(nint)   - Gaussian integration weights
  //    intPts(nint,2)  - local epsilon and eta coordinates for each
  //                      integration point
  //
  printf("Calculating integration points and weights...\n");
  if (nen == 3)          // lineartriangular element
  {
    // M=1 Gaussian quadrature for right-isosceles triangle domain
    // Local coordinate system centered at the 90-degree corner
    // Single integration point at the centroid
    weights[0] = 0.5;
    intPts[0][0] = 0.33333333333333;
    intPts[0][1] = 0.33333333333333;
  }
  else if (nen == 6)   // quadratic triangular element
  {
    // M=3 Gaussian quadrature for right-isosceles triangle domain
    // Local coordinate system centered at the 90-degree corner
    // Integration points start bottom left and go clockwise
    weights[0] = 0.16666666666667; 
    intPts[0][0] = 0.16666666666667;
    intPts[0][1] = 0.16666666666667;
    weights[2] = 0.66666666666666;
    intPts[2][0] = 0.16666666666667;
    intPts[2][1] = 0.66666666666666;;
    weights[1] = 0.66666666666666;;
    intPts[1][0] = 0.66666666666666;;
    intPts[1][1] = 0.16666666666667;
  }
  else if (nen == 4)     // linear quadrilateral element
  {   
    // M=2 Gaussian quadrature for square domain
    // Local coordinate system centered at the center of square
    for (int i = 0; i < nint; i++)
      {weights[i]=1;}           // all weights are 1
    double val = 0.57735027;    // 1/sqrt(3)
    intPts[0][0] = -val; 
    intPts[0][1] = -val;
    intPts[1][0] =  val;
    intPts[1][1] = -val;
    intPts[2][0] =  val;
    intPts[2][1] =  val;
    intPts[3][0] = -val;
    intPts[3][1] =  val;
  }
  else if (nen = 8)    // quadratic quadrilateral element
  {
    // M=3 Gaussian quadrature for square domain
    // Local coordinate system centered at the center of square
    // Integration points start from top-left and go by rows
    double val = 0.77459667; // sqrt(3/5)
    weights[3] = 0.30864197530864; // top-left
    intPts[3][0] = -val;
    intPts[3][1] = val;
    weights[6] = 0.49382716049383; // top-center
    intPts[6][0] = 0;
    intPts[6][1] = val;
    weights[2] = 0.30864197530864; // top-right
    intPts[2][0] = val;
    intPts[2][1] = val;
    weights[7] = 0.49382716049383; // center-left
    intPts[7][0] = -val;
    intPts[7][1] = 0;
    weights[8] = 0.79012345679012; // center-center
    intPts[8][0] = 0;
    intPts[8][1] = 0;
    weights[5] = 0.49382716049383; // center-right
    intPts[5][0] = val;
    intPts[5][1] = 0;
    weights[0] = 0.30864197530864; // bottom-left
    intPts[0][0] = -val;
    intPts[0][1] = -val;
    weights[4] = 0.49382716049383; // bottom-center
    intPts[4][0] = 0;
    intPts[4][1] = -val;
    weights[1] = 0.30864197530864; // bottom-right
    intPts[1][0] = val;
    intPts[1][1] = -val;

  }
  else                  // higher-order edges
  {
    printf("WARNING: Higher-order element. Elem_GetInfo() \
      couldn't be called by ele_stiff!\n");
    exit(EXIT_FAILURE);
  }
  printf("DONE\n");
  nee = nen * nint;       // size of the element stiffness matrix
  // Build a list of element node coordinates and create ien
  printf("Retrieving element node coordinates and form ien vector...\n");
  double coords[3];
  int localID = 0;
  vector<pMeshEnt> adjVtx;
  FMDB_Ent_GetAdj(elem, FMDB_VERTEX, 0, adjVtx);
  printf("Going through vertices...\n");
  for (int n = 0; n < adjVtx.size(); n++)
  {
    FMDB_Vtx_GetCoord(adjVtx[n], coords);
    nodeCoords[localID][0] = coords[0];
    nodeCoords[localID][1] = coords[1];
    ien[localID] = FMDB_Ent_ID(adjVtx[n]);
    localID++;
  }
  if ((nen == 6)||(nen == 8))
  {
    vector<pMeshEnt> adjEdges;
    FMDB_Ent_GetAdj(elem, FMDB_EDGE, 0, adjEdges);
    printf("Going through edge nodes...\n");
    for (int k = 0; k < adjEdges.size(); k++)
    {
      pNode node;
      FMDB_Edge_GetNode(adjEdges[k], 0, node);
      FMDB_Node_GetCoord(node, coords);
      nodeCoords[localID][0] = coords[0];
      nodeCoords[localID][1] = coords[1];
      ien[localID] = FMDB_Ent_ID(adjEdges[k]);
      localID++;
    }
  }
  printf("DONE\n");
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void Elem_GenShp(pMeshEnt elem,       double eps,      double eta,
                 int nsd,             int nen,         int nint,
                 vector< vector<double> > nodeCoords, double& detJ,
                 vector<double>& N,   vector<double>& dNdx,
                 vector<double>& dNdy)
{
  // Generates the shape function and its derivatives at a given
  // integration point on the element.
  //
  // Input:
  //    elem        - pMeshEnt element object
  //    eps         - local epsilon coordinate for the integration point
  //    eta         - local eta coordinate for the integration point
  //    nsd         - number of space dimensions
  //    nen         - number of element nodes
  //    nint        - number of integration points
  //    nodeCoords  - coordinate vector for element nodes
  // Output:
  //    detJ        - determinant of the Jacobian
  //    N           - shape function
  //    dNdx        - partial derivative of the shape function with
  //                  respect to global x
  //    dNdy        - partial derivative of the shape function with
  //                  respect to global y
  //
  // printf("Calculating shape functions and their derivatives...\n");
  double dN_dEps[nen], dN_dEta[nen];
  // 1-D Linear Element: (NOT YET SUPPORTED - DO NOT USE)
  if (nen == 2)
  {
    // Define the shape functions.
    N[0] = 0.5*(1-eps);
    N[1] = 0.5*(eps+1);
    // Calculate the local derivatives of the shape functions.
    dN_dEps[0] = -0.5;
    dN_dEps[1] =  0.5;
  }
  // 2-D Linear Triangle Element:
  else if (nen == 3)
  {
    // Define the shape functions.
    N[0] = 1-eps-eta;
    N[1] = eps;
    N[2] = eta;
    // Calculate the local derivatives of the shape functions.
    dN_dEps[0] = -1;
    dN_dEta[0] = -1;
    dN_dEps[1] =  1;
    dN_dEta[1] =  0;
    dN_dEps[2] =  0;
    dN_dEps[2] =  1;
  }
  // 2-D Quadratic Triangle Element:
  else if (nen == 6)
  {
    // Define the shape functions.
    N[0] = 2*(1-eps-eta)*(0.5-eps-eta);
    N[1] = 2*eps*(eps-0.5);
    N[2] = 2*eta*(eta-0.5);
    N[3] = 4*(1-eps-eta)*eps;
    N[4] = 4*eps*eta;
    N[5] = 4*(1-eps-eta)*eta;
    // Calculate the local derivatives of the shape functions.
    dN_dEps[0] = -2*(0.5 - eps - eta) - 2*(1 - eps - eta);
    dN_dEta[0] = -2*(0.5 - eps - eta) - 2*(1 - eps - eta);
    dN_dEps[1] =  2*(-0.5 + eps) + 2*eps;
    dN_dEta[1] =  0;
    dN_dEps[2] =  0;
    dN_dEta[2] =  2*(-0.5 + eta) + 2*eta;
    dN_dEps[3] =  4 - 8*eps - 4*eta;
    dN_dEta[3] = -4*eps;
    dN_dEps[4] =  4*eta;
    dN_dEta[4] =  4*eps;
    dN_dEps[5] = -4*eta;
    dN_dEta[5] =  4*(1 - eps - eta) - 4*eta;
  }
  // 2-D Linear Quadrilateral Element:
  else if (nen == 4)
  {
    // Define shape functions.
    N[0] = 0.25*(1-eps)*(1-eta);
    N[1] = 0.25*(1+eps)*(1-eta);
    N[2] = 0.25*(1+eps)*(1+eta);
    N[3] = 0.25*(1-eps)*(1+eta);
    // Calculate the local derivatives of the shape functions.
    dN_dEps[0] = -0.25*(1 - eta);
    dN_dEta[0] = -0.25*(1 - eps);
    dN_dEps[1] =  0.25*(1 - eta);
    dN_dEta[1] = -0.25*(1 + eps);
    dN_dEps[2] =  0.25*(1 + eta);
    dN_dEta[2] =  0.25*(1 + eps);
    dN_dEps[3] = -0.25*(1 + eta);
    dN_dEta[3] =  0.25*(1 - eps);

  }
  // 2-D Quadratic Quadrilateral Element:
  else if (nen == 8)
  {
    // Define the shape functions.
    N[0] = 0.25*(1-eps)*(1-eta)*(-eps-eta-1);
    N[1] = 0.25*(1+eps)*(1-eta)*(eps-eta-1);
    N[2] = 0.25*(1+eps)*(1+eta)*(eps+eta-1);
    N[3] = 0.25*(1-eps)*(1+eta)*(-eps+eta-1);
    N[4] = 0.5*(1-eps*eps)*(1-eta);
    N[5] = 0.5*(1+eps)*(1-eta*eta);
    N[6] = 0.5*(1-eps*eps)*(1+eta);
    N[7] = 0.5*(1-eps)*(1-eta*eta);
    // Calculate the local derivatives of the shape functions.
    dN_dEps[0] = -0.25*(1 - eps)*(1 - eta) - 0.25*(1 - eta)*(-1 - eps - eta);
    dN_dEta[0] = -0.25*(1 - eps)*(1 - eta) - 0.25*(1 - eps)*(-1 - eps - eta);
    dN_dEps[1] =  0.25*(1 + eps)*(1 - eta) + 0.25*(1 - eta)*(-1 + eps - eta);
    dN_dEta[1] = -0.25*(1 + eps)*(1 - eta) - 0.25*(1 + eps)*(-1 + eps - eta);
    dN_dEps[2] =  0.25*(1 + eps)*(1 + eta) + 0.25*(1 + eta)*(-1 + eps + eta);
    dN_dEta[2] =  0.25*(1 + eps)*(1 + eta) + 0.25*(1 + eps)*(-1 + eps + eta);
    dN_dEps[3] = -0.25*(1 - eps)*(1 + eta) - 0.25*(1 + eta)*(-1 - eps + eta);
    dN_dEta[3] =  0.25*(1 - eps)*(1 + eta) + 0.25*(1 - eps)*(-1 - eps + eta);
    dN_dEps[4] = -eps*(1 - eta);
    dN_dEta[4] = -0.5*(1 - eps*eps);
    dN_dEps[5] =  0.5*(1 - eta*eta);
    dN_dEta[5] = -eta*(1 + eps);
    dN_dEps[6] = -eps*(1 + eta);
    dN_dEta[6] =  0.5*(1 - eps*eps);
    dN_dEps[7] = -0.5*(1 - eta*eta);
    dN_dEta[7] = -eta*(1 - eps);
  }
  // Other element types: (NOT YET SUPPORTED)
  else
  {
    printf("ERROR: Unsupported element type in Elem_GenShp!\n");
    exit(EXIT_FAILURE);
  }
  // printf("Shape functions...DONE\n");
  // Collect local derivatives into a matrix.
  vector< vector<double> > dN_local(nsd, vector<double>(nen));
  for (int i = 0; i < nsd; i++)
  {
    for (int j = 0; j < nen; j++)
    {
      if (i == 0)
        {dN_local[i][j] = dN_dEps[j];}
      else if (i = 1)
        {dN_local[i][j] = dN_dEps[j];}
    }
  }
  // printf("Local derivatives...DONE\n");
  // printf("    dN_local:\n");
  // printMatrix(dN_local, nsd, nen);
  // Calculate the Jacobian, its determinant and its inverse (if
  // applicable).
  vector< vector<double> > J(nsd, vector<double>(nsd));
  double dx_dEps = 0;
  double dx_dEta = 0;
  double dy_dEps = 0;
  double dy_dEta = 0;
  for (int i = 0; i < nen; i++)
  {
    dx_dEps = dx_dEps + dN_dEps[i]*nodeCoords[i][0];
    dx_dEta = dx_dEta + dN_dEta[i]*nodeCoords[i][0];
    dy_dEps = dy_dEps + dN_dEps[i]*nodeCoords[i][1];
    dy_dEta = dy_dEta + dN_dEta[i]*nodeCoords[i][1];
  }
  detJ = dx_dEps*dy_dEta - dx_dEta*dy_dEps;
  // Calculate the global derivatives of the shape functions.
  for (int i = 0; i < nen; i++)
  {
    dNdx[i] =  (dy_dEta*dN_dEps[i] - dy_dEps*dN_dEta[i])/detJ;
    dNdy[i] = -(dx_dEta*dN_dEps[i] - dx_dEps*dN_dEta[i])/detJ;
  }
  // printf("Global derivatives...DONE\n");
  
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
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
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
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
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
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
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void FEA(pMeshMdl mesh, pGeomMdl geom, int nsd,
         vector< vector<double> > D,
         vector< vector<int> > ntyp,
         vector< vector<double> > FG,
         vector<double> f_init,
         vector< vector<double> > hOn,
         vector<double> h)
{
  // Core Finite Element Analaysis procedure.
  //
  // Inputs:
  //    mesh              - pMeshMdl mesh object
  //    geom              - pGeomMdl geom object
  //    D                 - material properties matrix
  //    ntyp              - vector indicating the type of DoF for each
  //                        component at the node
  //    FG                - vector of nodal point forces or nodal
  //                        prescribed essential boundary conditions
  //    f_init            - distributed forces on the geometry
  //
  // Process the problem mesh.
  int nnp, nel;
  pPart part;
  printf("Processing problem mesh...\n");
  process_mesh(mesh, geom, nsd, part, nnp, nel);
  printf("DONE\n");

  // Initialize the problem equation.
  vector< vector< vector<double> > > id(nsd, vector< vector<double> >(nnp, vector<double>(2)));
  vector<double> G;
  vector<double> F;
  int ndof, ndog;
  setup_eq(nnp, nel, nsd, ntyp, FG, id, G, F, ndof, ndog);
  printf("Allocating global stiffness matrix...\n");
  vector< vector<double> > K(ndof, vector<double>(ndof));
  for (int a = 0; a < ndof; a++)
  {
    for (int b = 0; b < ndof; b++)
      {K[a][b]=0;} // Zero out the global stiffness matrix
  }
  printf("DONE\n");
  printf("Governing equation initialized.\n");

  // Initiate the global stress-displacement matrix.
  vector< vector< vector< vector<double> > > > S;

  // Loop over all elements, assuming that each face is an element.
  printf("Starting element iteration...\n");
  pPartEntIter iter;
  FMDB_PartEntIter_Init(part, FMDB_FACE, FMDB_ALLTOPO, iter);
  pMeshEnt elem;
  int iterEnd = FMDB_PartEntIter_GetNext(iter, elem);
  int elemID;
  while (!iterEnd)
  {
    // Get information about the element.
    elemID = FMDB_Ent_ID(elem);
    printf("Element Number: %i\n", elemID);
    int nint, nen, nee;
    Elem_CheckType(elem, nint, nen);
    nee = nen*nsd;
    printf("    Number of integration points: %i\n", nint);
    printf("    Number of element nodes: %i\n", nen);
    vector<int> ien(nen);
    vector<double> weights(nint);
    vector< vector<double> > intPts(nint, vector<double>(nsd));
    vector< vector<double> > nodeCoords(nen, vector<double>(nsd));
    Elem_GetInfo(elem, nsd, nen, nint, nee, ien, weights, intPts, nodeCoords);
    printf("    Element node coordinates:\n");
    for (int n = 0; n < nen; n++)
      {printf("        %f    %f\n",nodeCoords[n][0],nodeCoords[n][1]);}

    // Construct the element stiffness matrix and force vector.
    printf("Constructing the element stiffness matrix and force vector...\n");
    vector< vector< vector<double> > > lm(nsd, vector< vector<double> >(nen, vector<double>(2)));
    vector< vector<double> > KE(nee, vector<double>(nee));
    vector<double> FE(nee);
    vector< vector< vector<double> > > Se(nint, vector< vector<double> >(3, vector<double>(nsd*nen)));
    ele_stiff(elem, nsd, nint, nen, nee, D, f_init, hOn, h, 
              intPts, weights, nodeCoords, ien, id, lm, KE, FE, Se);
    // Clear out vectors
    printf("DONE\n");
    S.push_back(Se);
    printf("    Element Stiffness Matrix (KE):\n");
    printMatrix(KE, nee, nee);
    printf("    Element Force Matrix (FE):\n");
    for (int i = 0; i < FE.size(); i++)
      {printf("    |  %d  |\n", FE[i]);}

    // Assemble the element contributions into the global matrices.
    printf("Assembling the element contributions into global matrices...\n");
    assemble(nsd, nen, lm, KE, FE, G, K, F);
    printf("DONE\n");

    // Clear vectors
    ien.clear();
    weights.clear();
    intPts.clear();
    nodeCoords.clear();
    lm.clear();
    KE.clear();
    FE.clear();
    Se.clear();
    iterEnd = FMDB_PartEntIter_GetNext(iter, elem);
  }// finish looping over elements
  FMDB_PartEntIter_Del(iter);
  printf("Element iteration complete!\n");

  // Print matrices for inspection.
  printf("    Global Stiffness Matrix K:\n");
  printMatrix(K, ndof, ndof);
  printf("    Global Force Vector F:\n");
  for (int n = 0; n < ndof; n++)
    {printf("|  %f  |\n", F[n]);}

  // Solve the global Kd = F system.
  vector<double> disp(ndof);
  int maxIt = 100000;
  printf("Starting CG Solver...\n");
  CGSolve(K, ndof, ndof, F, ndof, maxIt, disp);
  printf("DONE\n");
  printf("    Global Displacement Vector (disp):\n");
  for (int n = 0; n < ndof; n++)
    {printf("|  %f  |\n", disp[n]);}

  // Print out the node displacements.
  vector<double> nodeDisp(2);
  printf("Outputting displacements...\n");
  output_disp(nnp, nsd, G, id, disp, nodeDisp);

  // Calculate and print out stresses/fluxes.
  printf("Calculating stresses at each element...\n");
  FMDB_PartEntIter_Init(part, FMDB_FACE, FMDB_ALLTOPO, iter);
  iterEnd = FMDB_PartEntIter_GetNext(iter, elem);
  while(!iterEnd)
  {
    // Get information about the element.
    elemID = FMDB_Ent_ID(elem);
    printf("Element Number: %i\n", elemID);
    int nint, nen, nee;
    Elem_CheckType(elem, nint, nen);
    nee = nen*nsd;
    vector<int> ien(nen);
    vector<double> weights(nint);
    vector< vector<double> > intPts(nint, vector<double>(nsd));
    vector< vector<double> > nodeCoords(nen, vector<double>(nsd));
    Elem_GetInfo(elem, nsd, nen, nint, nee, ien, weights, intPts, nodeCoords);
    // Calculate the element stress vector.
    vector<double> SIG(3);
    vector< vector< vector<double> > > lm(nsd, vector< vector<double> >(nen, vector<double>(2)));
    out_flux(nsd, nen, nint, ien, id, S[elemID], disp, G, lm, SIG);

    // Print out the results.
    printf("    Stresses at element %d:\n", elemID);
    for (int n = 0; n < 3; n++)
      {printf("|  %f  |\n", SIG[n]);}

    iterEnd = FMDB_PartEntIter_GetNext(iter, elem);
    ien.clear();
    weights.clear();
    intPts.clear();
    lm.clear();
    nodeCoords.clear();
    SIG.clear();
  }// finish looping over elements
  FMDB_PartEntIter_Del(iter);
  printf("DONE\n");
  printf("SUCCESS: Finite Element Analysis complete!\n");
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void printUsage(char* exe) {
   fprintf(stderr, "Usage: %s -m meshFileName [OPTIONS]\n", exe);
   fprintf(stderr, "-m specify the input mesh file name\n"
           "-g specify the input geometric model file name\n"
           "-w write the mesh to the specified file name\n");
}

int main(int argc, char* argv[]) {
   SCUTIL_Init(MPI_COMM_WORLD);

   const int rank = SCUTIL_CommRank();
   const int commSize = SCUTIL_CommSize();

   int writeMesh = 0;
   static const int MaxStringLength = 1024;
   char meshFileName[MaxStringLength] = " ";
   char geomFileName[MaxStringLength] = " ";
   char outMeshFileName[MaxStringLength] = " ";

   if (0 == rank) {
      int opt;
      while ((opt = getopt(argc, argv, "m:g:w:")) != -1) {
         switch (opt) {
            case 'm':
               strncpy(meshFileName, optarg, MaxStringLength);
               meshFileName[MaxStringLength - 1] = '\0';
               break;
            case 'g':
               strncpy(geomFileName, optarg, MaxStringLength);
               geomFileName[MaxStringLength - 1] = '\0';
               break;
            case 'w':
               writeMesh = 1;
               strncpy(outMeshFileName, optarg, MaxStringLength);
               outMeshFileName[MaxStringLength - 1] = '\0';
               break;
            default:
               printUsage(argv[0]);
               SCUTIL_Finalize();
               exit(EXIT_FAILURE);
         }
      }
   }
   // broadcast the user arguments
   MPI_Bcast(&writeMesh, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(meshFileName, MaxStringLength, MPI_CHAR, 0, MPI_COMM_WORLD);
   MPI_Bcast(outMeshFileName, MaxStringLength, MPI_CHAR, 0, MPI_COMM_WORLD);
   MPI_Bcast(geomFileName, MaxStringLength, MPI_CHAR, 0, MPI_COMM_WORLD);

   if (rank == 0) {
      printf("STATUS Inputs: meshFileName=%s geomFileName=%s outMeshFileName=%s\n",
              meshFileName, geomFileName, outMeshFileName);
      printf("STATUS # MPI Processes = %d\n", commSize);
   }

   // load geometric model
   pGeomMdl geom = NULL;

   // load mesh
   pMeshMdl mesh;
   if (strcmp(meshFileName, " ") != 0) {
      loadMesh(mesh, geom, meshFileName);
   }


   ////////////////////////////
   /// Write code below here //
   //////////////////////////// 

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


   ////////////////////////////
   /// Write code above here //
   //////////////////////////// 

   if ( writeMesh != 0 ) {
      FMDB_Mesh_WriteToFile(mesh, outMeshFileName, commSize-1);
   }   

   FMDB_Mesh_Del(mesh);

   SCUTIL_Finalize();

   return 0;
}



