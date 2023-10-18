/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */
/**
   \file config_sim.hpp
   \brief ConfigSim class header file.
   \author Konstantinos A. Mountris
   \date 11/07/2021
*/

#pragma once
#ifndef PHYNETOUCH_APPS_TOOLS_CONFIG_SIM_HPP_
#define PHYNETOUCH_APPS_TOOLS_CONFIG_SIM_HPP_

#include "PHYNETOUCH/PHYNETOUCH"

#include "config_approximation.hpp"
#include "config_geo.hpp"
#include "config_material.hpp"
#include "config_conditions.hpp"
#include "config_output.hpp"
#include "config_physics.hpp"
// #include "config_post_process.hpp"
#include "config_units.hpp"
#include "parser.hpp"

#include <CLOUDEA/CLOUDEA>
#include <IMP/IMP>

#include <termcolor/termcolor.hpp>

#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <memory>
#include <unordered_map>

using namespace PNT;

namespace PNTSIM {

/** \addtogroup Application-Tools \{ */

/**
 * \class ConfigSim
 * \brief Class to configure and execute a simulation wth the PhyNeTouchSim app.
 * \tparam DIM The dimensions of the geometry model of the simulation.
 * \tparam CELL_NODES The number of nodes of the geometry model's cells.
 */
template<short DIM, short CELL_NODES>
class ConfigSim {

public:

    /**
     * \brief ConfigSim object constructor.
     */
    inline ConfigSim();


    /**
     * \brief ConfigSim object destructor.
     */
    inline virtual ~ConfigSim();


    /**
     * \brief Check if a valid simulation file has been provided.
     * \param [in] parser The parser of the simulation file.
     * \return [void]
     */
    inline void CheckValid(const Parser &parser);


    /**
     * \brief Launch a simulation.
     * Both 2D and 3D tissues are supported.
     * \param [in] parser The parser of the simulation file.
     * \param [in] mpi_handler The handler of the MPI variables.
     * \param [out] stream The output logging stream.
     * \return [void].
     */
    inline void Launch(const Parser &parser, MpiHandler mpi_handler, std::ostream &stream);

};

/** \} End of Doxygen Groups*/

} //end of namespace PNTSIM

#endif // PHYNETOUCH_APPS_TOOLS_CONFIG_SIM_HPP_

#include "config_sim.tpp"