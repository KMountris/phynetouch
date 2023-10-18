/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */


/**
   \file catheter.hpp
   \brief Catheter class header file.
   \author Konstantinos A. Mountris
   \date 23/07/2021
*/

#ifndef PHYNETOUCH_MATERIALS_CATHETER_HPP_
#define PHYNETOUCH_MATERIALS_CATHETER_HPP_


#include "PHYNETOUCH/engine/utilities/logger.hpp"

#include <Eigen/Dense>

#include <IMP/IMP>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <iostream>
#include <iterator>
#include <vector>
#include <numeric>
#include <string>
#include <stdexcept>
#include <exception>


namespace PNT {

/** \addtogroup Materials \{ */


/**
 * \class Catheter
 * \author Konstantinos A. Mountris
 * \brief Class implemmenting a catheter material.
 */
template <short DIM>
class Catheter
{

private:

    std::vector<double> density_;                           /**< The catheter's density */

    std::vector<double> specific_heat_;                     /**< The catheter's specific heat */

    std::vector<double> electrical_conductivity_;           /**< The catheter's electrical conductivity */

    std::vector<double> thermal_conductivity_;              /**< The catheter's thermal conductivity */

    std::vector<double> enthalpy_;                          /**< The catheter's enthalpy */

    std::vector<int> surface_node_ids_;                     /**< The node indices of the catheter surface */

    std::vector<int> connector_node_ids_;                   /**< The node indices of the catheter connected to another material */

    std::vector<std::vector<int>> connected_node_ids_;      /**< The connected node ids of another material that are connected to the connector node ids */

    int nodes_num_;                                         /**< The number of the material's nodes */


public:

    /**
     * \brief Catheter default constructor.
     */
    Catheter();


    /**
     * \brief Catheter default destructor.
     */
    virtual ~Catheter();


    /**
     * \brief Set the number of the catheter's nodes.
     * \param [in] nodes_num The number of the nodes that belong to the catheter.
     * \return [void]
     */
    inline void SetNodesNum(int nodes_num);


    /**
     * \brief Set the indices of the nodes at the surface of the catheter.
     * \param [in] surface_node_ids The indices of the catheter's surface nodes.
     * \return void
     */
    inline void SetSurfaceNodeIds(const std::vector<int> &surface_node_ids);


    /**
     * \brief Identify the connected nodes of the tissue to the nodes of the catheter surface.
     * \param [in] catheter_nodes The coordinates of the catheter material's geometry nodes.
     * \param [in] tissue_nodes The coordinates of the tissue material's geometry nodes.
     * \param [in] free_tissue_nodes_flag Flag to distinguish tissue nodes that are on the free boundary.
     * \param [in] influence_radius The radius in which the tissue nodes are considered connected to a catheter surface node.
     * \return [void]
     */
    inline void IdentifyTissueConnection(const std::vector<IMP::Vec<DIM,double>> &catheter_nodes,
        const std::vector<IMP::Vec<DIM,double>> &tissue_nodes, double influence_radius, int max_connected);


    /**
     * \brief Compute the enthalpy of the catheter.
     * \return [void]
     */
    inline void ComputeEnthalpy();


    /**
     * \brief Set the density of the catheter.
     * \param [in] density The density of the catheter.
     * \return [void]
     */
    inline void SetDensity(const std::vector<double> &density);


    /**
     * \brief Set the specific heat of the catheter.
     * \param [in] specific_heat The specific heat of the catheter.
     * \return [void]
     */
    inline void SetSpecificHeat(const std::vector<double> &specific_heat);


    /**
     * \brief Set the electrical conductivity of the catheter.
     * \param [in] conductivity The electrical conductivity of the catheter.
     * \return [void]
     */
    inline void SetElectricalConductivity(const std::vector<double> &conductivity);


    /**
     * \brief Set the thermal conductivity of the catheter.
     * \param [in] conductivity The thermal conductivity of the catheter.
     * \return [void]
     */
    inline void SetThermalConductivity(const std::vector<double> &conductivity);


    /**
     * \brief Get the number of nodes belonging to the electric material.
     * \return [int] The number of nodes belonging to the electric material.
     */
    inline auto NodesNum() const { return this->nodes_num_; }


    /**
     * \brief Get the density of the catheter.
     * \return [const std::vector<double>&] The density of the catheter.
     */
    inline auto & Density() const { return this->density_; }


    /**
     * \brief Get the density of a specific point of the catheter.
     * \param [in] id The index of the catheter's point.
     * \return [double] The density of a specific point of the catheter.
     */
    inline auto Density(std::size_t id) const { return this->density_[id]; }


    /**
     * \brief Get the specific heat of the catheter.
     * \return [const std::vector<double>&] The specific heat of the catheter.
     */
    inline auto & SpecificHeat() const { return this->specific_heat_; }


    /**
     * \brief Get the specific heat of a specific point of the catheter.
     * \param [in] id The index of the catheter's point.
     * \return [double] The specific heat of a specific point of the catheter.
     */
    inline auto SpecificHeat(std::size_t id) const { return this->specific_heat_[id]; }


    /**
     * \brief Get the electrical conductivity of the catheter.
     * \return [const std::vector<double>&] The electrical conductivity of the catheter.
     */
    inline auto & ElectricalConductivity() const { return this->electrical_conductivity_; }


    /**
     * \brief Get the electrical conductivity of a specific point of the catheter.
     * \param [in] id The index of the catheter's point.
     * \return [double] The electrical conductivity of a specific point of the catheter.
     */
    inline auto ElectricalConductivity(std::size_t id) const { return this->electrical_conductivity_[id]; }


    /**
     * \brief Get the thermal conductivity of the catheter.
     * \return [const std::vector<double>&] The thermal conductivity of the catheter.
     */
    inline auto & ThermalConductivity() const { return this->thermal_conductivity_; }


    /**
     * \brief Get the thermal conductivity of a specific point of the catheter.
     * \param [in] id The index of the catheter's point.
     * \return [double] The thermal conductivity of a specific point of the catheter.
     */
    inline auto ThermalConductivity(std::size_t id) const { return this->thermal_conductivity_[id]; }


    /**
     * \brief Get the enthalpy of the catheter.
     * \return [const std::vector<double>&] The enthalpy of the catheter.
     */
    inline auto & Enthalpy() const { return this->enthalpy_; }


    /**
     * \brief Get the enthalpy of a specific point of the catheter.
     * \param [in] id The index of the catheter's point.
     * \return [double] The enthalpy of a specific point of the catheter.
     */
    inline auto Enthalpy(std::size_t id) const { return this->enthalpy_[id]; }


    /**
     * \brief Get the indices of the connector nodes (nodes of the catheter that are connected to the tissue).
     * \return [const std::vector<int>&] The indices of the connector nodes (nodes of the catheter that are connected to the tissue).
     */
    inline auto & ConnectorNodeIds() const { return this->connector_node_ids_; }


    inline auto ConnectorNodeIds(std::size_t id) const { return this->connector_node_ids_[id]; }


    inline auto ConnectorNodesNum() const { return static_cast<int>(this->connector_node_ids_.size()); }


    /**
     * \brief Get the nodes of the tissue that are connected to the connector nodes of the catheter.
     * \return [const std::vector<std::vector<int>>&] The nodes of the tissue that are connected to the connector nodes of the catheter.
     */
    inline auto & ConnectedNodeIds() const { return this->connected_node_ids_; }


    inline auto & ConnectedNodeIds(std::size_t id) const { return this->connected_node_ids_[id]; }


};

/** \} End of Doxygen Groups */

} // End of namespace PNT

#endif //PHYNETOUCH_MATERIALS_CATHETER_HPP_

#include "PHYNETOUCH/engine/materials/catheter.tpp"