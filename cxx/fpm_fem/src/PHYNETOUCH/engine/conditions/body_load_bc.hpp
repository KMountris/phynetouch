/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */

/**
   \file body_load_bc.hpp
   \brief Header file for body load boundary condition.
   \author Konstantinos A. Mountris
   \date 29/07/2020
*/

#pragma once
#ifndef PHYNETOUCH_BODY_LOAD_BC_HPP_
#define PHYNETOUCH_BODY_LOAD_BC_HPP_

#include "PHYNETOUCH/engine/conditions/load_curve_basic.hpp"
#include "PHYNETOUCH/engine/conditions/load_curve_factory.hpp"

#include <IMP/Vectors>
#include <vector>
#include <memory>
#include <limits>


namespace PNT {

/** \addtogroup Conditions \{ */

/**
 * \class BodyLoadBc
 * \author Konstantinos A. Mountris
 * \brief Body load boundary condition. If the boundary condition is applied on a scalar field its dimension is 1.
 *        Otherwise it is the same with the dimension of the geometry.
 * /tparam DIM the dimension of the body load boundary condition.
 */
template<short DIM>
class BodyLoadBc {

private:

    IMP::Vec<DIM,double> dir_;                      /**< The direction of the body load boundary condition application */

    double value_;                                  /**< The body load condition imposed value */

    std::shared_ptr<LoadCurveBasic> load_curve_;    /**< The load curve defining how the body load boundary condition is applied */

    LoadCurveType load_type_;

    double load_start_;

    double load_duration_;

public:

    /**
     * \brief The default constructor.
     */
    BodyLoadBc();


    /**
     * \brief The destructor.
     */
    virtual ~BodyLoadBc();


    /**
     * \brief Set the value of the External load condition.
     * \param [in] value The value of the External load condition.
     * \return [void]
     */
    inline void SetValue(double value);


    /**
     * \brief Set the application direction of the External load boundary condition.
     * \param [in] dir The application direction of the External load boundary condition application.
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
     * \brief Define the loading increments of the external load boundary condition.
     * \param [in] dt The time step of each loading increment.
     * \return [void]
     */
    inline void LoadingDefinition(double dt);


    /**
     * \brief Apply External load boundary condition on loads matrix.
     * \param [in] step The index of the integration's step.
     * \param [in] dt The time step to retrive the displament BC from the load curve.
     * \param [out] load The load matrix where the body load boundary condition will be applied.
     *  It is applied on the entries of the matrix that corresponds to the nodes indices of the boundary condition.
     * \return [void]
     */
    inline void Apply(int step, double dt, Eigen::RowVectorXd &load) const;


    /**
     * \brief Get the value of the External load condition.
     * \return [double] The value of the External load condition.
     */
    inline double Value() const { return this->value_; }


    /**
     * \brief Get the application direction of the External load boundary condition.
     * \return [const IMP::Vec<DIM,double>&] The application direction of the External load boundary condition.
     */
    inline auto & Direction() const { return this->dir_; }


    inline auto LoadType() const { return this->load_type_; }


    inline auto LoadStart() const { return this->load_start_; }


    inline auto LoadDuration() const { return this->load_duration_; }


};


/** \} End of Doxygen Groups */

} // End of namespace PNT

#endif //PHYNETOUCH_BODY_LOAD_BC_HPP_

#include "PHYNETOUCH/engine/conditions/body_load_bc.tpp"