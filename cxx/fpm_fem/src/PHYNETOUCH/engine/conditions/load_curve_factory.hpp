/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 */


/**
   \file load_curve_factory.hpp
   \brief Header file for a factory for load curve generation.
   \author Konstantinos A. Mountris
   \date 31/10/2021
*/

#ifndef PNT_CONDITIONS_LOAD_CURVE_FACTORY_HPP_
#define PNT_CONDITIONS_LOAD_CURVE_FACTORY_HPP_

#include "PHYNETOUCH/engine/conditions/load_curve_basic.hpp"
#include "PHYNETOUCH/engine/conditions/load_curve_smooth.hpp"
#include "PHYNETOUCH/engine/conditions/load_curve_step.hpp"

#include <memory>
#include <string>


namespace PNT {

/** \addtogroup Conditions \{ */

/**
 * \class LoadCurveFactory
 * \brief Class implementing a a factory for load curve generation.
 * \author Konstantinos A. Mountris
 */
class LoadCurveFactory
{
public:

    /**
     * \brief Creates a load curve according to the provided type.
     * \param [in] curve_type The type of the desired load curve.
     * \return [std::shared_ptr<LoadCurveBasic>] Unique pointer to the created load curve.
     */
    inline static std::shared_ptr<LoadCurveBasic> Create(LoadCurveType curve_type) {

        std::shared_ptr<LoadCurveBasic> curve_ptr;

        switch (curve_type) {
        case LoadCurveType::smooth :
            curve_ptr = std::make_unique<LoadCurveSmooth>();
            break;
        case LoadCurveType::step :
            curve_ptr = std::make_unique<LoadCurveStep>();
            break;
        default:
            std::string error_msg = Logger::Error("Could not create load curve. Not supported type.");
            throw std::invalid_argument(error_msg);
            break;
        }

        return curve_ptr;
    }
};


/** \} End of Doxygen Groups */

} // End of namespace PNT

#endif // PNT_CONDITIONS_LOAD_CURVE_FACTORY_HPP_