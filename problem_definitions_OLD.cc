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