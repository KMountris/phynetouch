/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */

/**
   \file config_conditions.hpp
   \brief ConfigConditions class header file.
   \author Konstantinos A. Mountris
   \date 30/07/2021
*/

#pragma once
#ifndef PHYNETOUCH_APPS_TOOLS_CONFIG_CONDITIONS_HPP_
#define PHYNETOUCH_APPS_TOOLS_CONFIG_CONDITIONS_HPP_

#include "parser.hpp"

#include "PHYNETOUCH/engine/utilities/logger.hpp"
#include "PHYNETOUCH/engine/utilities/measure_units.hpp"
#include "PHYNETOUCH/engine/conditions/bound_conds.hpp"
#include "PHYNETOUCH/engine/conditions/dirichlet_bc.hpp"
#include "PHYNETOUCH/engine/conditions/initial_bc.hpp"

#include <IMP/IMP>
#include <termcolor/termcolor.hpp>

#include <string>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <filesystem>
#include <sstream>

using namespace PNT;

namespace PNTSIM {

/** \addtogroup Application-Tools \{ */

/**
 * \class ConfigConditions
 * \brief Class to configure the boundary conditions of a model.
 */
template<short DIM>
class ConfigConditions {


private:

    std::unordered_map<std::string,PNT::LoadCurveType> load_curve_types_;       /**< An unordered map for all possible loading types */


public:


    /**
     * \brief Constructor for Config Conditions.
     */
    ConfigConditions();


    /**
     * \brief Destructor for Config Conditions.
     */
    virtual ~ConfigConditions();


    /**
     * \brief Set up the boundary conditions of the model.
     * \param [in] parser The parser of the simulation input file.
     * \param [in] body_type The type of the model's body (tissue or catheter).
     * \param [in] nodesets The nodesets of the model to identify the nodes where boundary conditions will be applied.
     * \param [in] units The measure units.
     * \param [out] electrical_bc The boundary conditions for the electrical problem.
     * \param [out] bioheat_bc The boundary conditions for the bioheat problem.
     * \param [out] deform_bc The boundary conditions for the deformation problem.
     * \param [out] stream The message output stream.
     * \return [void]
     */
    void SetConditions(const Parser &parser, const MpiHandler &mpi_handler, const std::string &body_type,
        int body_nodes_num, const std::unordered_map<std::string, IMP::NodeSet> &nodesets, const MeasureUnits &units,
        BoundConds<1> &electrical_bc, BoundConds<1> &bioheat_bc, BoundConds<DIM> &deform_bc, std::ostream &stream) const;


protected:


    /**
     * \brief Set the initial boundary conditions of the model.
     * \param [in] parser The parser of the simulation input file.
     * \param [in] body_type The type of the model's body (tissue or catheter).
     * \param [in] units The measure units.
     * \param [out] electrical_bc The boundary conditions for the electrical problem where initial conditions will be set if any.
     * \param [out] bioheat_bc The boundary conditions for the bioheat problem where initial conditions will be set if any.
     * \param [out] deform_bc The boundary conditions for the deformation problem where initial conditions will be set if any.
     * \param [out] stream The message output stream.
     * \return [void]
     */
    void SetInitialBc(const Parser &parser, const MpiHandler &mpi_handler, const std::string &body_type,
        int body_nodesnum, const MeasureUnits &units, BoundConds<1> &electrical_bc, BoundConds<1> &bioheat_bc,
        BoundConds<DIM> &deform_bc, std::ostream &stream) const;


    /**
     * \brief Set the Dirichlet boundary conditions of the model.
     * \param [in] parser The parser of the simulation input file.
     * \param [in] body_type The type of the model's body (tissue or catheter).
     * \param [in] units The measure units.
     * \param [out] electrical_bc The boundary conditions for the electrical problem where Dirichlet conditions will be set if any.
     * \param [out] bioheat_bc The boundary conditions for the bioheat problem where Dirichlet conditions will be set if any.
     * \param [out] deform_bc The boundary conditions for the deformation problem where Dirichlet conditions will be set if any.
     * \param [out] stream The message output stream.
     * \return [void]
     */
    void SetDirichletBc(const Parser &parser, const MpiHandler &mpi_handler, const std::string &body_type,
        const std::unordered_map<std::string, IMP::NodeSet> &nodesets, const MeasureUnits &units,
        BoundConds<1> &electrical_bc, BoundConds<1> &bioheat_bc, BoundConds<DIM> &deform_bc, std::ostream &stream) const;


    /**
     * \brief Set the body load boundary conditions of the model.
     * \param [in] parser The parser of the simulation input file.
     * \param [in] body_type The type of the model's body (tissue or catheter).
     * \param [in] units The measure units.
     * \param [out] electrical_bc The boundary conditions for the electrical problem where body load conditions will be set if any.
     * \param [out] bioheat_bc The boundary conditions for the bioheat problem where body load conditions will be set if any.
     * \param [out] deform_bc The boundary conditions for the deformation problem where body load conditions will be set if any.
     * \param [out] stream The message output stream.
     * \return [void]
     */
    void SetBodyLoadBc(const Parser &parser, const MpiHandler &mpi_handler, const std::string &body_type,
        const MeasureUnits &units, BoundConds<1> &electrical_bc, BoundConds<1> &bioheat_bc,
        BoundConds<DIM> &deform_bc, std::ostream &stream) const;


};

/** \} End of Doxygen Groups*/

} //end of namespace PNTSIM

#endif //PHYNETOUCH_APPS_TOOLS_CONFIG_CONDITIONS_HPP_

#include "config_conditions.tpp"