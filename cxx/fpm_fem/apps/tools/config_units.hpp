/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */

/**
   \file config_units.hpp
   \brief ConfigUnits class header file.
   \author Konstantinos A. Mountris
   \date 11/07/2021
*/

#pragma once
#ifndef PHYNETOUCH_APPS_TOOLS_CONFIG_UNITS_HPP_
#define PHYNETOUCH_APPS_TOOLS_CONFIG_UNITS_HPP_

#include "parser.hpp"

#include "PHYNETOUCH/engine/utilities/logger.hpp"
#include "PHYNETOUCH/engine/utilities/measure_units.hpp"
#include "PHYNETOUCH/engine/utilities/mpi_handler.hpp"

#include <termcolor/termcolor.hpp>

#include <string>
#include <iostream>

using namespace PNT;

namespace PNTSIM {

/** \addtogroup Application-Tools \{ */

/**
 * \class ConfigUnits
 * \brief Class to configure the measure units of a simulation.
 */
class ConfigUnits {

public:

    /**
     * \brief ConfigUnits object constructor.
     */
    ConfigUnits();


    /**
     * \brief ConfigUnits object destructor.
     */
    virtual ~ConfigUnits();


    /**
     * \brief Set up the reference scale of the units for proper unit conversion.
     * \param [in] parser The parser of the simulation file.
     * \param [in] mpi_handler The handler of the MPI variables.
     * \param [out] units The units after setting up the reference scale.
     * \param [out] stream The output logging stream.
     * \return [void]
     */
    void SetReferenceScale(const Parser &parser, const MpiHandler &mpi_handler,
        MeasureUnits &units, std::ostream &stream) const;

};

/** \} End of Doxygen Groups*/

} //end of namespace PNTSIM

#endif //PHYNETOUCH_APPS_TOOLS_CONFIG_UNITS_HPP_