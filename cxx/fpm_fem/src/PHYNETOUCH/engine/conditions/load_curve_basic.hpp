/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 */


/**
   \file load_curve.hpp
   \brief Header file for a load curve describing time-varying loading conditions.
   \author Konstantinos A. Mountris
   \date 29/07/2020
*/

#ifndef PHYNETOUCH_CONDITIONS_LOAD_CURVE_BASIC_HPP_
#define PHYNETOUCH_CONDITIONS_LOAD_CURVE_BASIC_HPP_

#include "PHYNETOUCH/engine/utilities/logger.hpp"

#include <IMP/IMP>

#include <string>
#include <limits>
#include <vector>
#include <iostream>
#include <exception>
#include <stdexcept>


namespace PNT {

/** \addtogroup Conditions \{ */


/**
 * \enum Type of load curve.
 * \author Konstantinos A. Mountris
 */
enum class LoadCurveType {
    unknown,   /**< Unknown type load curve */
    smooth,    /**< Smooth load curve */
    step       /**< Step load curve */
};


/**
 * \class LoadCurveBasic
 * \author Konstantinos A. Mountris
 * \brief A load curve to describe time-varying loading conditions.
 */
class LoadCurveBasic {

private:

    double start_;                  /**< The starting time of the load application */

    double duration_;               /**< The time duration of the load application */


public:

    /**
     * \brief The default constructor of the LoadCurve.
     */
    explicit LoadCurveBasic() : start_(0.), duration_(0.) {}


    /**
     * \brief The destructor of the LoadCurve.
     */
    virtual ~LoadCurveBasic() {}


    /**
     * \brief Set the starting time of the loading curve.
     * \param [in] time The starting time of the loading curve.
     * \return [void]
     */
    inline void SetStart(double time) { this->start_ = time; }


    /**
     * \brief Set the time duration of the loading curve.
     * \param [in] time The time duration of the loading curve.
     * \return [void]
     */
    inline void SetDuration(double time) { this->duration_ = time; }


    /**
     * \brief Define the load application at each timestep.
     * \param [in] load The maximum load of the loading curve.
     * \param [in] dt The time step of the loading curve discretization.
     * \return [void]
     */
    virtual void DefineLoad(double load, double dt) = 0;


    /**
     * \brief Get the starting time of the loading curve.
     * \return [double] The starting time of the loading curve.
     */
    inline double Start() const { return this->start_; }


    /**
     * \brief Get the time duration of the loading curve.
     * \return [double] The time duration of the loading curve.
     */
    inline double Duration() const { return this->duration_; }


    /**
     * \brief Get the data of the loading curve.
     * \return [const std::vector<double>&] The data of the loading curve.
     */
    virtual const std::vector<double> & Data() const = 0;

};


/** \} End of Doxygen Groups */

} // End of namespace PNT

#endif //PHYNETOUCH_CONDITIONS_LOAD_CURVE_BASIC_HPP_