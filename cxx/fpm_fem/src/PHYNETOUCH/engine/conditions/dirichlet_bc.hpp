/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */

/**
   \file dirichlet_bc.hpp
   \brief Header file for Dirichlet boundary condition.
   \author Konstantinos A. Mountris
   \date 29/07/2020
*/

#ifndef PHYNETOUCH_DIRICHLET_BC_HPP_
#define PHYNETOUCH_DIRICHLET_BC_HPP_

#include "PHYNETOUCH/engine/conditions/load_curve_basic.hpp"
#include "PHYNETOUCH/engine/conditions/load_curve_factory.hpp"

#include <IMP/Vectors>
#include <vector>
#include <memory>



namespace PNT {

/** \addtogroup Conditions \{ */

/**
 * \class DirichletBc
 * \author Konstantinos A. Mountris
 * \brief Dirichlet boundary condition. If the boundary condition is applied on a scalar field its dimension is 1.
 *        Otherwise it is the same with the dimension of the geometry.
 * /tparam DIM the dimension of the Dirichlet boundary condition.
 */
template<short DIM>
class DirichletBc {

private:

    std::vector<int> node_ids_;                     /**< The indices of the nodes where the condition is imposed */

    IMP::Vec<DIM,double> dir_;                      /**< The direction of the Dirichlet boundary condition application */

    double value_;                                  /**< The Dirichlet condition imposed value */

    std::shared_ptr<LoadCurveBasic> load_curve_;    /**< The load curve defining how the Dirichlet boundary condition is applied */


public:

    /**
     * \brief The default constructor.
     */
    DirichletBc();


    /**
     * \brief The destructor.
     */
    virtual ~DirichletBc();


    /**
     * \brief Set the indices of the nodes where the Dirichlet condition is applied.
     * \param [in] node_ids The indices of the nodes.
     * \return [void]
     */
    inline void SetNodeIds(const std::vector<int> &node_ids);


    /**
     * \brief Set the value of the Dirichlet condition.
     * \param [in] value The value of the Dirichlet condition.
     * \return [void]
     */
    inline void SetValue(double value);


    /**
     * \brief Set the application direction of the Dirichlet boundary condition.
     * \param [in] dir The application direction of the Dirichlet boundary condition application.
     * \return [void]
     */
    inline void SetDirection(const IMP::Vec<DIM,double> &dir);


    /**
     * \brief Set the Loading Curve object
     * \param [in] load_curve_type 
     * \param [in] load_time 
     * \return [void]
     */
    inline void SetLoadingCurve(LoadCurveType load_curve_type, double load_start, double load_duration);


    /**
     * \brief Define the loading increments of the dirichlet boundary condition.
     * \param [in] dt The time step of each loading increment.
     * \return [void]
     */
    inline void LoadingDefinition(double dt);


    inline void UpdateValue(int step, double dt);

    
    /**
     * \brief Apply Dirichlet boundary condition on displacement matrix.
     * \param [in] padding The nodes padding to modify the appropriate row of the displacement matrix.
     * \param [in] step The index of the integration's step.
     * \param [in] dt The time step to retrive the displament BC from the load curve.
     * \param [out] displacement The displacement matrix where the Dirichlet boundary condition will be applied.
     *  It is applied on the entries of the matrix that corresponds to the nodes indices of the boundary condition.
     * \return [void]
     */
    inline void Apply(int padding, int step, double dt, Eigen::MatrixXd &displacement) const;


    /**
     * \brief Get the indices of the nodes where the Dirichlet condition is applied.
     * \return [const std::vector<int>&] The nodes where the Dirichlet condition is applied.
     */
    inline auto & NodeIds() const { return this->node_ids_; }


    /**
     * \brief Get the value of the Dirichlet condition.
     * \return [double] The value of the Dirichlet condition.
     */
    inline double Value() const { return this->value_; }


    /**
     * \brief Get the application direction of the Dirichlet boundary condition.
     * \return [const IMP::Vec<DIM,double>&] The application direction of the Dirichlet boundary condition.
     */
    inline auto & Direction() const { return this->dir_; }


};


/** \} End of Doxygen Groups */

} // End of namespace PNT

#endif //PHYNETOUCH_DIRICHLET_BC_HPP_

#include "PHYNETOUCH/engine/conditions/dirichlet_bc.tpp"