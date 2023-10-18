/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */

/**
   \file config_geo.hpp
   \brief ConfigGeo class header file.
   \author Konstantinos A. Mountris
   \date 11/07/2021
*/

#pragma once
#ifndef PHYNETOUCH_APPS_TOOLS_CONFIG_GEO_HPP_
#define PHYNETOUCH_APPS_TOOLS_CONFIG_GEO_HPP_

#include "parser.hpp"

#include "PHYNETOUCH/engine/utilities/logger.hpp"

#include <IMP/IMP>
#include <termcolor/termcolor.hpp>

#include <string>
#include <filesystem>
#include <iostream>
#include <unordered_map>
#include <algorithm>

namespace PNTSIM {

/** \addtogroup Application-Tools \{ */

/**
 * \class ConfigGeo
 * \brief Class to configure the geometry model of a simulation.
 * \tparam DIM The dimensions of the geometry.
 * \tparam CELL_NODES The number of nodes of the geometry's cells.
 */
template<short DIM,short CELL_NODES>
class ConfigGeo {

public:

    /**
     * \brief ConfigGeo object constructor.
     */
    ConfigGeo();


    /**
     * \brief ConfigGeo object destructor.
     */
    virtual ~ConfigGeo();


    /**
     * \brief Set the geometry model.
     * \param [in] parser The parser of the simulation input script.
     * \param [in] body_type The type of the body that is represented by the geometry model (tissue or catheter).
     * \param [out] mesh The mesh topology of the geometry model.
     * \param [out] voro The voronoi tesselation of the geometry model.
     * \param [out] stream The logging stream.
     * \return [void]
     */
    void SetGeometryModel(const Parser &parser, const MpiHandler &mpi_handler, const std::string &body_type,
        IMP::Mesh<DIM, CELL_NODES> &mesh, IMP::Voronoi<DIM> &voro, std::ostream &stream) const;


    /**
     * \brief Extract nodesets from either the mesh topology or the voronoi tesselation of the geometry model
     * depending on the used numerical approximation.
     * \param [in] mesh The mesh topology of the geometry model.
     * \param [in] voro The voronoi tesselation of the geometry model.
     * \param [out] node_sets The nodesets that will be extracted.
     * \param [out] nodes_num The total number of nodes in the geometry model.
     * \return [void]
     */
    void ExtractNodeSets(const IMP::Mesh<DIM, CELL_NODES> &mesh, const IMP::Voronoi<DIM> &voro,
        std::unordered_map<std::string, IMP::NodeSet> &node_sets, int &nodes_num) const;
};

/** \} End of Doxygen Groups*/

} //end of namespace PNTSIM

#endif //PHYNETOUCH_APPS_TOOLS_CONFIG_GEO_HPP_

#include "config_geo.tpp"