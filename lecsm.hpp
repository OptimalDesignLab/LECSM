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
	void set_u(const InnerProdVector & u_csm);

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
   * \brief set the nodal boundary conditions (displacements and forcing)
   * \param[in] BCtype - type of BC (displacement or forcing)
   * \param[in] BYval - value of BC
   */
	void SetBoundaryConds(const InnerProdVector & BCtype, const InnerProdVector & BCval);

	/*!
   * \brief initializes the global vectors used in the solver
   * \param[out] G - global prescribed displacements vector
   * \param[out] F - global prescribed nodal forcing
   * \param[out] K - global stiffness matrix
   */
	void InitGlobalVecs(vector<double>& G, vector<double>& F,
											vector< vector<double> >& K);

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
   * \brief calculates the (dS/du)*u_csm product
   * \param[in] u_csm - displacement vector
   * \param[out] v_csm - product vector
   */
	void Calc_dSdu_Product(InnerProdVector& u_csm, InnerProdVector& v_csm);

	/*!
   * \brief calculates the (dA/du)*u_csm product
   * \param[in] u_csm - displacement vector
   * \param[out] v_csm - product vector
   */
	void Calc_dAdu_Product(InnerProdVector& u_csm, InnerProdVector& wrk);

	/*!
   * \brief calculates the (dS/dp)*(dp/dq)*u_cfd product
   * \param[in] wrk - (dp/dq)*u_cfd vector
   * \param[out] v_csm - product vector
   */
	void Calc_dSdp_Product(InnerProdVector& wrk, InnerProdVector& u_cfd);

	/*!
   * \brief calculates the nozzle area and stores in area_
   */
	void CalcArea();

	/*!
   * \brief calculates the CSM residual using u_ displacements
   */
	void CalcResidual();

  /*!
   * \brief inspects the solver mesh
   */
  void InspectMesh() {
    geom_.InspectElements();
  }

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
	int nnp_;
	double E_, w_, t_, h_;
	Mesh geom_;
};