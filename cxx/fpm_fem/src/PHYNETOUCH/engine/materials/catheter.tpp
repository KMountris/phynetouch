/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


#ifndef PHYNETOUCH_MATERIALS_CATHETER_TPP_
#define PHYNETOUCH_MATERIALS_CATHETER_TPP_

#include "PHYNETOUCH/engine/materials/catheter.hpp"

namespace PNT {

template <short DIM>
Catheter<DIM>::Catheter() : density_(), specific_heat_(), electrical_conductivity_(), thermal_conductivity_(),
        enthalpy_(), surface_node_ids_(), connector_node_ids_(), connected_node_ids_(), nodes_num_(0)
{}


template <short DIM>
Catheter<DIM>::~Catheter()
{}


template <short DIM>
void Catheter<DIM>::SetNodesNum(int nodes_num)
{
    // Check if given number of material nodes is positive.
    if (nodes_num < 0) {
        throw std::invalid_argument(Logger::Error("Could not set the number of nodes belonging to the cardiac tissue material. A negative number of material nodes was given."));
    }

    this->nodes_num_ = nodes_num;
}


template <short DIM>
void Catheter<DIM>::SetSurfaceNodeIds(const std::vector<int> &surface_node_ids)
{
    this->surface_node_ids_ = surface_node_ids;
    this->surface_node_ids_.clear();
    for (int i = 0; i != this->NodesNum(); ++i) {
        this->surface_node_ids_.emplace_back(i);
    }
    Logger::Warning("Replacing the surface nodeset with all the catheter nodes.");
}


template <short DIM>
void Catheter<DIM>::IdentifyTissueConnection(const std::vector<IMP::Vec<DIM,double>> &catheter_nodes,
    const std::vector<IMP::Vec<DIM,double>> &tissue_nodes, double influence_radius, int max_connected)
{
    if (this->surface_node_ids_.size() == 0) {
        auto error_msg = "Could not identify catheter connection with tissue. The surface node indices of the catheter have not been set.";
        std::runtime_error(Logger::Error(error_msg));
    }

    // Boost definitions.
    namespace bg = boost::geometry;
    namespace bgi = boost::geometry::index;
    typedef bg::model::point<double,DIM,bg::cs::cartesian> bg_point;
    typedef std::pair<bg_point,std::size_t> value;

    // Reset the connector and connected nodes containers.
    this->connector_node_ids_.clear(); this->connector_node_ids_.reserve(this->surface_node_ids_.size());
    this->connected_node_ids_.clear(); this->connected_node_ids_.reserve(this->surface_node_ids_.size());

    // Populate rtree with tissue nodes.
    auto rtree = bgi::rtree<value,bgi::quadratic<16>>{};
    auto query_point = bg_point{};
    for (std::size_t i = 0; i != tissue_nodes.size(); ++i) {
        if constexpr (DIM == 1) {
            bg::set<0>(query_point, tissue_nodes[i][0]);
        } else if constexpr (DIM == 2) {
            bg::set<0>(query_point, tissue_nodes[i][0]);
            bg::set<1>(query_point, tissue_nodes[i][1]);
        } else if constexpr (DIM == 3) {
            bg::set<0>(query_point, tissue_nodes[i][0]);
            bg::set<1>(query_point, tissue_nodes[i][1]);
            bg::set<2>(query_point, tissue_nodes[i][2]);
        }
        rtree.insert(std::make_pair(query_point, i));
    }

    // Search for the connected tissue nodes to each catheter surface node.
    auto connected_tissue_nodes = std::vector<value>{};
    auto surf_point = bg_point{};
    for (std::size_t i = 0; i != this->surface_node_ids_.size(); ++i) {
        // Get the index of the catheter surface node.
        auto sf_nid = this->surface_node_ids_[i];

        // Coordinates of the current surface point.
        if constexpr (DIM == 1) {
            bg::set<0>(surf_point, catheter_nodes[sf_nid][0]);
        } else if constexpr (DIM == 2) {
            bg::set<0>(surf_point, catheter_nodes[sf_nid][0]);
            bg::set<1>(surf_point, catheter_nodes[sf_nid][1]);
        } else if constexpr (DIM == 3) {
            bg::set<0>(surf_point, catheter_nodes[sf_nid][0]);
            bg::set<1>(surf_point, catheter_nodes[sf_nid][1]);
            bg::set<2>(surf_point, catheter_nodes[sf_nid][2]);
        }

        // Perform the search.
        rtree.query(bgi::satisfies([&](value const& v) { return bg::distance(v.first, surf_point) < influence_radius; }),
            std::back_inserter(connected_tissue_nodes));

        // Store the connector and connected nodes indices.
        if (connected_tissue_nodes.size() > 0) {
            this->connector_node_ids_.emplace_back(sf_nid);
            auto connected_ids = std::vector<int>{};
            for (const auto &cnode : connected_tissue_nodes) {
                connected_ids.emplace_back(cnode.second);
                if (connected_ids.size() >= static_cast<std::size_t>(max_connected))  { break; }
            }
            this->connected_node_ids_.emplace_back(connected_ids);
        }

        // Clear connected tissue nodes container for the next iteration.
        connected_tissue_nodes.clear();
    }
    this->connector_node_ids_.shrink_to_fit();
    this->connected_node_ids_.shrink_to_fit();

    if (this->connector_node_ids_.size() == 0) {
        std::string err_msg = "Could not set connector elements between tissue and catheter. Please, increase the connector radius.";
        throw std::invalid_argument(Logger::Error(err_msg));
    }
}


template <short DIM>
void Catheter<DIM>::ComputeEnthalpy()
{
    if (this->density_.size() == 0) {
        auto error_msg = "Could not compute catheter enthalpy. Set density first.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }
    if (this->specific_heat_.size() == 0) {
        auto error_msg = "Could not compute catheter enthalpy. Set specific heat first.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    this->enthalpy_.clear();
    this->enthalpy_.reserve(this->density_.size());
    std::transform(std::begin(this->density_), std::end(this->density_),
            std::begin(this->specific_heat_), std::begin(this->enthalpy_), std::multiplies<double>());
}


template <short DIM>
void Catheter<DIM>::SetDensity(const std::vector<double> &density)
{
    // Check if the nodes of the material have been set.
    if (this->nodes_num_ == 0) {
        throw std::runtime_error(Logger::Error("Could not set the density of the catheter. The material's nodes number has not been set."));
    }

    // Set density. Either different for each node or the same.
    if (static_cast<int>(density.size()) == this->nodes_num_) {
        this->density_ = density;
    } else if (density.size() == 1) {
        this->density_.assign(this->nodes_num_,density[0]);
        this->density_.shrink_to_fit();
    } else {
        auto error_msg = "Could not set density for catheter.\n"
                "The given density vector should contain either a single value or a value for each node of the material.";
        throw std::runtime_error(Logger::Error(error_msg));
    }
}


template <short DIM>
void Catheter<DIM>::SetSpecificHeat(const std::vector<double> &specific_heat)
{
    // Check if the nodes of the material have been set.
    if (this->nodes_num_ == 0) {
        throw std::runtime_error(Logger::Error("Could not set specific heat for catheter. The material's nodes number has not been set."));
    }

    // Set specific heat. Either different for each node or the same.
    if (static_cast<int>(specific_heat.size()) == this->nodes_num_) {
        this->specific_heat_ = specific_heat;
    } else if (specific_heat.size() == 1) {
        this->specific_heat_.assign(this->nodes_num_,specific_heat[0]);
        this->specific_heat_.shrink_to_fit();
    } else {
        auto error_msg = "Could not set specific heat for catheter.\n"
                "The given specific heat vector should contain either a single value or a value for each node of the material.";
        throw std::runtime_error(Logger::Error(error_msg));
    }
}


template <short DIM>
void Catheter<DIM>::SetElectricalConductivity(const std::vector<double> &conductivity)
{
    // Check if the nodes of the material have been set.
    if (this->nodes_num_ == 0) {
        throw std::runtime_error(Logger::Error("Could not set electrical conductivity for catheter. The material's nodes number has not been set."));
    }

    // Set electrical conductivity. Either different for each node or the same.
    if (static_cast<int>(conductivity.size()) == this->nodes_num_) {
        this->electrical_conductivity_ = conductivity;
    } else if (conductivity.size() == 1) {
        this->electrical_conductivity_.assign(this->nodes_num_,conductivity[0]);
        this->electrical_conductivity_.shrink_to_fit();
    } else {
        auto error_msg = "Could not set electrical conductivity for catheter.\n"
                "The given electrical conductivity vector should contain either a single value or a value for each node of the material.";
        throw std::runtime_error(Logger::Error(error_msg));
    }
}


template <short DIM>
void Catheter<DIM>::SetThermalConductivity(const std::vector<double> &conductivity)
{
    // Check if the nodes of the material have been set.
    if (this->nodes_num_ == 0) {
        throw std::runtime_error(Logger::Error("Could not set thermal conductivity for catheter. The material's nodes number has not been set."));
    }

    // Set thermal conductivity. Either different for each node or the same.
    if (static_cast<int>(conductivity.size()) == this->nodes_num_) {
        this->thermal_conductivity_ = conductivity;
    } else if (conductivity.size() == 1) {
        this->thermal_conductivity_.assign(this->nodes_num_,conductivity[0]);
        this->thermal_conductivity_.shrink_to_fit();
    } else {
        auto error_msg = "Could not set thermal conductivity for catheter material.\n"
                "The given thermal conductivity vector should contain either a single value or a value for each node of the material.";
        throw std::runtime_error(Logger::Error(error_msg));
    }
}


} // End of namespace PNT


#endif //PHYNETOUCH_MATERIALS_CATHETER_TPP_