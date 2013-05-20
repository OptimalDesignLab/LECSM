/**
 * \file 1D_mesh_tools.hpp
 * \brief 1D mesh library
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

#pragma once
#include <vector>
#include "../quasi_1d_euler/inner_prod_vector.hpp"
using namespace std;

// =====================================================================

/*!
 * \class Node
 * \brief defines the nodes of a 2-D mesh
 */
class Node {
public:
	int id;
	int type[3]; // 0 for zero BC, 1 for free node, 2 for prescribed displacements
	double coords[2];
	double dispBC[3]; // 3 degrees of freedom per node (axial, trans, rotational)
	double forceBC[3]; // forcing in two directions plus moment about the third

  Node() {}

	/*!
   * \brief construct the node
   * \param[in] num - global node id
   * \param[in] c - vector of node coordinates
   */
  Node(int num, double* c)
  {
    id = num;
    for (int i=0; i<3; i++) {
      coords[i] = c[i];
      type[i] = 1;            // node is initially free
    }   
  }

  ~Node() {}

	/*!
   * \brief define nodal boundary conditions
   * \param[in] BCtype - type of boundary conditions
   * \param[in] BCval - value of the boundary conditions
   */
	void DefineBCs(int* BCtype, double* BCval);
};

// =====================================================================

/*!
 * \class Element
 * \brief defines a 2D first order finite element beam
 */
class Element {
public:
	int id, nen;
	vector<Node> adjNodes;

  Element() {}

	/*!
   * \brief construct the element
   * \param[in] num - global element id
   * \param[in] nodes - vector of nodes defining the element
   */
	Element(int num, vector<Node> nodes);

  ~Element() {}

	/*!
   * \brief construct the element stiffness matrix and forcing vector
   * \param[in] E - young's modulus
   * \param[in] w - element width
   * \param[in] t - element thickness
   * \param[in] P - pressures on element nodes
   * \param[in] gm - global equation number mapping
   * \param[out] lm - local equation number mapping
   * \param[out] KE - element stiffness matrix
   * \param[out] FE - element forcing vector
   */
	void GetElemStiff(double E, double w, double t, vector<double>& P,
                	  vector< vector< vector<int> > >& gm,
                	  vector< vector< vector<int> > >& lm,
                	  vector< vector<double> >& KE, vector<double>& FE);

	/*!
   * \brief assemble the element contributions into global matrix/vectors
   * \param[in] KE - element stiffness matrix
   * \param[in] FE - element forcing vector
   * \param[in] lm - local equation number mapping
   * \param[in] G - global prescribed displacements vector
   * \param[out] K - global stiffness matrix
   * \param[out] F - global forcing vector
   */
	void Assemble(vector< vector<double> >& KE, vector<double>& FE,
								vector< vector< vector<int> > >& lm,
								vector<double>& G, vector<double>& F,
                vector< vector<double> >& K);
};

// =====================================================================

/*!
 * \class Mesh
 * \brief defines a 2D finite element beam mesh
 */
class Mesh {
public:
	int nnp, nel, ndof, ndog;
	vector<Element> allElems;
	vector<Node> allNodes;

  /*!
   * \brief class constructor
   */
  Mesh() {
    nnp = 0;
    nel = 0;
  }

	/*!
   * \brief construct the mesh
   * \param[in] elems - vector of mesh elements
   */
	void CreateMesh(vector<Element>& elems, vector<Node>& nodes)
  {
    allElems = elems;
    allNodes = nodes;
    nel = allElems.size();
    nnp = allNodes.size();
  }

  /*!
   * \brief default destructor
   */
  ~Mesh() {} //

  void InspectNodes();

  void InspectElements();

	/*!
   * \brief generate the global equation number map
   * \param[out] gm - global equation number mapping
   */
	void SetupEq(vector< vector< vector<int> > >& gm);

	/*!
   * \brief update the mesh
   * \param[in] u_csm - nodal displacement vector
   */
	void Update(const InnerProdVector& u_csm);
};