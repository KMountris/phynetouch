/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */

/**
   \file contact_handle.hpp
   \brief Header file for contact handling.
   \author Konstantinos A. Mountris
   \date 07/03/2022
*/

#ifndef PHYNETOUCH_ENGINE_PHYSICS_CONTACT_HANDLE_HPP_
#define PHYNETOUCH_ENGINE_PHYSICS_CONTACT_HANDLE_HPP_

#include "PHYNETOUCH/engine/utilities/kd_tree_impvec_vector_adaptor.hpp"
#include "PHYNETOUCH/engine/utilities/logger.hpp"

#include <IMP/IMP>

#include <nanoflann.hpp>

#include <vector>
#include <unordered_map>
#include <string>

namespace PNT {

/** \addtogroup Physics \{ */

/**
 * \class ContactHandle
 * \author Konstantinos A. Mountris
 * \brief Handling the contact between a soft and rigid body, such as the contact between the tissue and the catheter.
 * /tparam DIM the dimension of the Dirichlet boundary condition.
 */
template<short DIM, short CELL_NODES>
class ContactHandle {

private:

    std::vector<int> slave_node_ids_;                               /**< */

    std::vector<int> master_node_ids_;                              /**< */

    std::vector<IMP::Vec<DIM,double>> master_node_normals_;

    std::vector<std::vector<IMP::Vec<DIM,int>>> master_facets_conn_;        /**< */

    bool is_enabled_;


public:

    /**
     * \brief Construct a new Contact Handle object
     */
    ContactHandle();


    /**
     * \brief Destroy a Contact Handle object
     */
    virtual ~ContactHandle();


    /**
     * @brief
     */
    inline void Enable() { this->is_enabled_ = true; }


    /**
     * @brief
     */
    inline void Disable() { this->is_enabled_ = false; }


    /**
     * \brief Set the Slave object
     * \param [in] slave_mesh
     * \param [in] slave_nset_name
     */
    inline void SetSlave(const IMP::Mesh<DIM,CELL_NODES> &slave_mesh, const std::string slave_nset_name);


    /**
     * \brief Set the Slave object
     * \param [in] slave_voro
     * \param [in] slave_nset_name
     */
    inline void SetSlave(const IMP::Voronoi<DIM> &slave_voro, const std::string slave_nset_name);


    /**
     * \brief Set the Master object
     * \param [in] master_mesh
     * \param [in] master_nset_name
     */
    inline void SetMaster(const IMP::Mesh<DIM,CELL_NODES> &master_mesh, const std::string master_nset_name);


    /**
     * \brief Set the Master object
     * \param [in] master_voro
     * \param [in] master_nset_name
     */
    inline void SetMaster(const IMP::Voronoi<DIM> &master_voro, const std::string master_nset_name);


    /**
     * \brief
     * \param [in] master_nodes
     * \param [in] slave_nodes
     * \return [int]
     */
    inline void NearestMasterToSlave(const std::vector<IMP::Vec<DIM, double>> &master_nodes,
        const std::vector<IMP::Vec<DIM,double>> &slave_nodes, std::vector<int> &near_master_to_slave) const;


    /**
     * @brief
     * @return [true]
     * @return [false]
     */
    bool IsEnabled() const { return this->is_enabled_; }


    /**
     * @brief
     * @return int
     */
    const auto & MasterNodeIds() const { return this->master_node_ids_; }


    const auto & MasterNodeIds(std::size_t id) const { return this->master_node_ids_[id]; }


    const auto & MasterNodeIdsAt(std::size_t id) const { return this->master_node_ids_.at(id); }


    const auto & MasterNodeNormals() const { return this->master_node_normals_; }


    const auto & MasterNodeNormals(std::size_t id) const { return this->master_node_normals_[id]; }


    const auto & MasterNodeNormalsAt(std::size_t id) const { return this->master_node_normals_.at(id); }


    /**
     * @brief
     * @return int
     */
    auto MasterNodeIdsNum() const { return static_cast<int>(this->master_node_ids_.size()); }


    const auto & MasterFacetsConn() const { return this->master_facets_; }


    const auto & MasterFacetsConn(std::size_t id) const { return this->master_facets_conn_[id]; }


    const auto & MasterFacetsAt(std::size_t id) const { return this->master_facets_conn_.at(id); }


    /**
     * @brief
     * @return int
     */
    const auto & SlaveNodeIds() const { return this->slave_node_ids_; }


    const auto & SlaveNodeIds(std::size_t id) const { return this->slave_node_ids_[id]; }


    const auto & SlaveNodeIdsAt(std::size_t id) const { return this->slave_node_ids_.at(id); }


    /**
     * @brief
     * @return int
     */
    auto SlaveNodeIdsNum() const { return static_cast<int>(this->slave_node_ids_.size()); }


};


/** \} End of Doxygen Groups */

} // End of namespace PNT

#endif //PHYNETOUCH_ENGINE_PHYSICS_CONTACT_HANDLE_HPP_

#include "PHYNETOUCH/engine/physics/contact_handle.tpp"