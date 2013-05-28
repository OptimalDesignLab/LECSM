/**
 * \file lecsm.hpp
 * \brief LECSM solver header file
 * \author  Alp Dener <alp.dener@gmail.com>
 * \version 1.0
 */

#pragma once
#include <vector>
#include "./1D_mesh_tools.hpp"
#include "../quasi_1d_euler/inner_prod_vector.hpp"
#include <krylov.hpp>

using namespace std;

// =====================================================================

/*!
 * \class LECSM
 * \brief 2D linear elastic computational structural analysis solver
 */
class LECSM {
public:

	/*!
   * \brief class constructor
   * \param[in] nnp - number of nodes in the domain
   * \param[in] x - vector of x coordinates
   * \param[in] y - vector of y coordinates
   */
	LECSM(int nnp) :
			area_(nnp, 0.0),
			xCoords_(nnp, 0.0),
      yCoords_(nnp, 0.0),
			res_(3*nnp, 0.0),
			u_(3*nnp, 0.0),
			P_(nnp, 0.0) { nnp_ = nnp; }

	/*!
   * \brief default destructor
   */
	~LECSM() {}

  /*!
   * \brief returns a vector x coordinates for each node
   * \returns xCoords_ member value
   */
  InnerProdVector & get_x() { return xCoords_; }
  
  /*!
   * \brief returns a vector y coordinates for each node
   * \returns yCoords_ member value
   */
  InnerProdVector & get_y() { return yCoords_; }

	/*!
   * \brief returns a vector areas at each node
   * \returns area_ member value
   */
	InnerProdVector & get_area() { return area_; }

	/*!
   * \brief returns the CSM residual vector
   * \returns res_ member value
   */
	InnerProdVector & get_res() { return res_; }

  /*!
   * \brief returns the CSM displacement vector
   * \returns u_ member value
   */
  InnerProdVector & get_u() { return u_; }

	/*!
   * \brief sets the material properties for the CSM solver
   * \param[in] E - young's modulus
   * \param[in] t - 2D beam thickness
   * \param[in] w - nozzle fixed width
   * \param[in] h - nozzle max height (at y_coord = 0)
   */
	void set_material(double E, double t, double w, double h) { 
      E_ = E; t_ = t; w_ = w; h_ = h;}

	/*!
   * \brief updates the CSM mesh with displacements
   * \param[in] u_csm - vector of displacements (3*nnp)
   */
	void set_u(const InnerProdVector & u_csm) { u_ = u_csm; }

	/*!
   * \brief sets the pressure values at each node
   * \param[in] press - vector of pressure values (nnp)
   */
	void set_press(const InnerProdVector & press) { P_ = press; }

	/*!
   * \brief generates the initial problem mesh
   * \param[in] x - vector of node x coordinates
   * \param[in] y - vector of node y coordinates
   */
	void GenerateMesh(const InnerProdVector & x, const InnerProdVector & y);

  /*!
   * \brief resets the solver coordinates to the nodal coordinates of the geometry
   */
  void ResetCoords();

  /*!
   * \brief resets the solver coordinates to the nodal coordinates of the geometry
   */
  void SetCoords(const InnerProdVector & x, const InnerProdVector & y);

  /*!
   * \brief updates the geometry with the coordinates stored in the solver
   */
  void UpdateMesh();

	/*!
   * \brief set the nodal boundary conditions (displacements and forcing)
   * \param[in] BCtype - type of BC (displacement or forcing)
   * \param[in] BYval - value of BC
   */
	void SetBoundaryConds(const InnerProdVector & BCtype, const InnerProdVector & BCval);

	/*!
   * \brief initializes the global vectors used in the solver
   * \param[out] G - global prescribed displacements vector
   * \param[out] F - global prescribed nodal forcing
   */
	void InitGlobalVecs(vector<double>& G, vector<double>& F);

	/*!
   * \brief calculates the global stiffness matrix and associated vectors
   * \param[in] gm - global equation number mapping
   * \param[out] G - global prescribed displacements vector
   * \param[out] F - global prescribed nodal forcing
   * \param[out] K - global stiffness matrix
   */
	void GetStiff(vector< vector< vector<int> > >& gm,
								vector<double>& G, vector<double>& F,
                vector< vector<double> >& K);

  /*!
   * \brief preconditions the input vector with the diagonal of the system stiffness matrix
   * \param[in] in - un-preconditioned vector
   * \param[out] v_csm - preconditioned vector
   */
  void Precondition(InnerProdVector& in, InnerProdVector& out);

	/*!
   * \brief calculates the (dS/du)*vector product
   * \param[in] in - multiplied vector (3*num_nodes)
   * \param[out] out - resultant vector (3*num_nodes)
   */
  void Calc_dSdu_Product(const InnerProdVector& in, InnerProdVector& out);

	/*!
   * \brief calculates the (dA/du)*vector product
   * \param[in] in - multiplied vector (3*num_nodes)
   * \param[out] out - resultant vector (num_nodes)
   */
	void Calc_dAdu_Product(InnerProdVector& in, InnerProdVector& out);

  /*!
   * \brief calculates the [(dA/du)^T]*vector product
   * \param[in] in - multiplied vector (num_nodes)
   * \param[out] out - resultant vector (3*num_nodes)
   */
  void CalcTrans_dAdu_Product(InnerProdVector& in, InnerProdVector& out);

  /*!
   * \brief calculates the (dy/dA)*vector product
   * \param[in] in - multiplied vector (num_nodes)
   * \param[out] out - resultant vector (num_nodes)
   */
  void Calc_dydA_Product(InnerProdVector& in, InnerProdVector& out);

  /*!
   * \brief calculates the [(dy/dA)^T]*vector product
   * \param[in] in - multiplied vector (num_nodes)
   * \param[out] out - resultant vector (num_nodes)
   */
  void CalcTrans_dydA_Product(InnerProdVector& in, InnerProdVector& out);

  /*!
   * \brief product for FD derivative residual w.r.t nodal coordinates
   * \param[in] in - multiplied vector (num_nodes)
   * \param[out] out - resultant vector (num_nodes)
   */
  void CalcFD_dSdy_Product(InnerProdVector& in, InnerProdVector& out);

  /*!
   * \brief transpose product for FD derivative residual w.r.t nodal coordinates
   * \param[in] in - multiplied vector (num_nodes)
   * \param[out] out - resultant vector (num_nodes)
   */
  void CalcTransFD_dSdy_Product(InnerProdVector& in, InnerProdVector& out);

	/*!
   * \brief calculates the (dS/dp)*vector product
   * \param[in] in - multiplied vector (num_nodes)
   * \param[out] out - resultant vector (3*num_nodes)
   */
	void Calc_dSdp_Product(InnerProdVector& in, InnerProdVector& out);

  /*!
   * \brief calculates the [(dS/dp)^T]*vector product
   * \param[in] in - multiplied vector (3*num_nodes)
   * \param[out] out - resultant vector (num_nodes)
   */
  void CalcTrans_dSdp_Product(InnerProdVector& in, InnerProdVector& out);

  /*!
   * \brief calculates the displaced coordinates and nozzle area
   */
	void CalcCoordsAndArea();

	/*!
   * \brief calculates the CSM residual based on displacements in u_
   */
	void CalcResidual();

  /*!
   * \brief inspects the solver mesh
   */
  void InspectMesh() {
    geom_.InspectElements();
  }

  /*!
   * \brief solves the CSM problem for a given forcing vector using MINRES
   * \param[in] rhs - prescribed right-hand-side vector for the system
   * \param[in] max_iter - maximum number of iterations permitted
   * \param[in] tol - relative residual tolerance desired
   * \result nodal displacements are calculated and saved to u_
   * \returns number of matrix-vector products
   */
  int SolveFor(InnerProdVector & rhs, const int & max_iter = 1000,
               const double & tol = 1e-6);

	/*!
   * \brief independent solution of a CSM problem using conjugate gradient
   */
	void Solve();

private:
	InnerProdVector area_;
	InnerProdVector xCoords_;
  InnerProdVector yCoords_;
	InnerProdVector res_;
	InnerProdVector u_;
	InnerProdVector P_;
  vector< vector<double> > K_;
	int nnp_;
	double E_, w_, t_, h_;
	Mesh geom_;
};

// ======================================================================

/*!
 * \class StiffnessVectorProduct
 * \brief specialization of matrix-vector product for lecsm
 */
class StiffnessVectorProduct :
    public kona::MatrixVectorProduct<InnerProdVector> {
 public:

  /*!
   * \brief default constructor
   * \param[in] csm_solver - a linear elasticity solver (defines product)
   */
  StiffnessVectorProduct(LECSM * solver) {
    solver_ = solver; 
  } 

  ~StiffnessVectorProduct() {} ///< class destructor

  /*!
   * \brief operator that defines the Stiffness-Matrix-Vector product
   * \param[in] u - vector that is being multiplied by K
   * \param[out] v - vector that is the result of the product
   */
  void operator()(const InnerProdVector & u, InnerProdVector & v);

 private:
  LECSM * solver_; ///< used to access the Matrix-Vector product routine
};

// ======================================================================

/*!
 * \class ApproxStiff
 * \brief specialization of preconditioner for lecsm
 */
class ApproxStiff :
    public kona::Preconditioner<InnerProdVector> {
 public:

  /*!
   * \brief default constructor
   * \param[in] euler_solver - a Quasi1DEuler solver to access precond.
   */
  ApproxStiff(LECSM * solver) {
    solver_ = solver; 
  }

  ~ApproxStiff() {} ///< class destructor

  /*!
   * \brief operator that applies the approximate inverse Stiffness matrix.
   * \param[in] u - vector that is being preconditioned
   * \param[out] v - vector that is the result of the preconditioning
   */
  void operator()(InnerProdVector & u, InnerProdVector & v) {
    v = u;
  }

 private:
  LECSM * solver_; ///< used to access the Stiffness matrix
};
