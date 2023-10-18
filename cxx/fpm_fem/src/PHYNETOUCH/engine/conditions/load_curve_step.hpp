/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 */


/**
   \file step_load_curve.hpp
   \brief Header file for a step load curve describing time-varying loading conditions.
   \author Konstantinos A. Mountris
   \date 30/10/2020
*/

#ifndef PHYNETOUCH_CONDITIONS_LOAD_CURVE_STEP_HPP_
#define PHYNETOUCH_CONDITIONS_LOAD_CURVE_STEP_HPP_

#include "PHYNETOUCH/engine/conditions/load_curve_basic.hpp"
#include "PHYNETOUCH/engine/utilities/logger.hpp"

#include <IMP/IMP>

#include <cmath>
#include <limits>
#include <vector>
#include <iostream>
#include <exception>
#include <stdexcept>


namespace PNT {

/** \addtogroup Conditions \{ */

/**
 * \class StepLoadCurve
 * \author Konstantinos A. Mountris
 * \brief A step load curve to describe time-varying loading conditions.
 */
class LoadCurveStep : public LoadCurveBasic
{

private:

    std::vector<double> data_;      /**< The load curve data */


public:

    /**
     * \brief The default constructor of the LoadCurveStep.
     */
    LoadCurveStep();


    /**
     * \brief The destructor of the LoadCurveStep.
     */
    virtual ~LoadCurveStep();


    /**
     * \brief Define the load application at each timestep.
     * \param [in] load The maximum load of the loading curve.
     * \param [in] dt The time step of the loading curve discretization.
     * \return [void]
     */
    virtual void DefineLoad(double load, double dt);


    /**
     * \brief Get the data of the loading curve.
     * \return [const std::vector<double>&] The data of the loading curve.
     */
    inline const std::vector<double> & Data() const { return this->data_; }

};


/** \} End of Doxygen Groups */

} // End of namespace PNT

#endif //PHYNETOUCH_CONDITIONS_LOAD_CURVE_STEP_HPP_