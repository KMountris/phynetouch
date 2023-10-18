/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */

#ifndef PHYNETOUCH_APPS_TOOLS_CONFIG_GEO_TPP_
#define PHYNETOUCH_APPS_TOOLS_CONFIG_GEO_TPP_

#include "config_geo.hpp"

namespace PNTSIM
{

template<short DIM, short CELL_NODES>
ConfigGeo<DIM, CELL_NODES>::ConfigGeo()
{}


template<short DIM, short CELL_NODES>
ConfigGeo<DIM, CELL_NODES>::~ConfigGeo()
{}


template<short DIM, short CELL_NODES>
void ConfigGeo<DIM, CELL_NODES>::SetGeometryModel(const Parser &parser, const MpiHandler &mpi_handler,
    const std::string &body_type, IMP::Mesh<DIM, CELL_NODES> &mesh, IMP::Voronoi<DIM> &voro, std::ostream &stream) const
{
    // Check body type.
    if (body_type != "tissue" && body_type != "catheter") {
        std::string error_msg = "Could not set geometry model. Body type should be: tissue or catheter.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Resolve geometry file path.
    auto filename = parser.GetValue<std::string>(body_type+".geometry.mesh");
    parser.ResolveFilePath(filename);

    // Get type of numerical approximation method.
    auto method = parser.GetValue<std::string>("numerical approximation.method");
    std::transform(std::begin(method), std::end(method), std::begin(method), ::tolower);

    // Set geometry according to the numerical approximation method.
    if (method == "fem") {
        mesh.LoadFrom(filename);
        if (mpi_handler.rank_id == 0) {
            stream << Logger::Message("Loaded mesh: " + filename + "\n");
            stream << Logger::Message("Number of nodes: ") << mesh.NodesNum() << "\n";
            stream << Logger::Message("Number of cells: ") << mesh.CellsNum() << "\n";
        }
    } else if (method == "fpm") {
        // Temporary load of mesh for output.
        mesh.LoadFrom(filename);

        // Load the voronoi tesselation.
        auto voro_file = parser.GetValue<std::string>(body_type+".geometry.voronoi");
        parser.ResolveFilePath(voro_file);
        voro.LoadFrom(voro_file);
        voro.ComputeCellMeasures();

        if (mpi_handler.rank_id == 0) {
            stream << Logger::Message("Loaded mesh: " + filename + "\n");
            stream << Logger::Message("Loaded voronoi tesselation: " + voro_file + "\n");
            stream << Logger::Message("Number of nodes: ") << voro.NodesNum() << "\n";
            stream << Logger::Message("Number of voronoi points: ") << voro.PointsNum() << "\n";
            stream << Logger::Message("Number of voronoi facets: ") << voro.FacetsNum() << "\n";
            stream << Logger::Message("Number of voronoi cells: ") << voro.CellsNum() << "\n";
        }
    }
}


template<short DIM, short CELL_NODES>
void ConfigGeo<DIM, CELL_NODES>::ExtractNodeSets(const IMP::Mesh<DIM, CELL_NODES> &mesh, const IMP::Voronoi<DIM> &voro,
        std::unordered_map<std::string, IMP::NodeSet> &node_sets, int &nodes_num) const
{
    if (mesh.NodesNum() != 0) { // FEM case
        node_sets = mesh.NodeSets();
        nodes_num = mesh.NodesNum();
    } else if (voro.NodesNum() != 0) { // FPM case
        node_sets = voro.NodeSets();
        nodes_num = voro.NodesNum();
    } else {
        std::string error_msg = "Could not extract node sets. Simulation method is not supported. Supported: FEM | FPM";
        throw std::invalid_argument(Logger::Error(error_msg));
    }
}


} // end of namespace PNTSIM

#endif //PHYNETOUCH_APPS_TOOLS_CONFIG_UNITS_TPP_