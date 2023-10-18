/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */


#ifndef PHYNETOUCH_ENGINE_PHYSICS_CONTACT_HANDLE_TPP_
#define PHYNETOUCH_ENGINE_PHYSICS_CONTACT_HANDLE_TPP_

#include "PHYNETOUCH/engine/physics/contact_handle.hpp"

namespace PNT
{

template<short DIM, short CELL_NODES>
ContactHandle<DIM,CELL_NODES>::ContactHandle() : slave_node_ids_(), master_node_ids_(),
    master_node_normals_(), master_facets_conn_(), is_enabled_(false)
{}


template<short DIM, short CELL_NODES>
ContactHandle<DIM,CELL_NODES>::~ContactHandle()
{}


template<short DIM, short CELL_NODES>
void ContactHandle<DIM,CELL_NODES>::SetSlave(const IMP::Mesh<DIM,CELL_NODES> &slave_mesh, const std::string slave_nset_name)
{
    if (slave_mesh.NodeSets().find(slave_nset_name) != slave_mesh.NodeSets().end()) {
        this->slave_node_ids_ = slave_mesh.NodeSets().at(slave_nset_name).NodeIds();
    } else {
        auto err_msg = "Could not set slave in contact handle. The given set: "+slave_nset_name+" does not exist in the slave geometry.";
        throw std::runtime_error(Logger::Error(err_msg));
    }
}


template<short DIM, short CELL_NODES>
void ContactHandle<DIM,CELL_NODES>::SetSlave(const IMP::Voronoi<DIM> &slave_voro, const std::string slave_nset_name)
{
    if (slave_voro.NodeSets().find(slave_nset_name) != slave_voro.NodeSets().end()) {
        this->slave_node_ids_ = slave_voro.NodeSets().at(slave_nset_name).NodeIds();
    } else {
        auto err_msg = "Could not set slave in contact handle. The given set: "+slave_nset_name+" does not exist in the slave geometry.";
        throw std::runtime_error(Logger::Error(err_msg));
    }
}


template<short DIM, short CELL_NODES>
void ContactHandle<DIM,CELL_NODES>::SetMaster(const IMP::Mesh<DIM,CELL_NODES> &master_mesh, const std::string master_nset_name)
{
    if (master_mesh.NodeSets().find(master_nset_name) != master_mesh.NodeSets().end()) {
        this->master_node_ids_ = master_mesh.NodeSets().at(master_nset_name).NodeIds();
    } else {
        auto err_msg = "Could not set master in contact handle. The given node set: "+master_nset_name+" does not exist in the master geometry.";
        throw std::runtime_error(Logger::Error(err_msg));
    }

    // Set master nodes flag and position in the container.
    auto master_nodes_flag = std::vector<int>(master_mesh.NodesNum(), 0);
    auto master_nodes_it = std::unordered_map<int,int>{};
    auto it = int{0};
    for (const auto &nid : this->master_node_ids_) {
        master_nodes_flag[nid] = 1;
        master_nodes_it[nid] = it++;
    }

    // Get all the boundary facets of the master mesh.
    auto boundary_facets = master_mesh.FreeFacets();

    // Get the boundary facets that are attached on a master node.
    auto master_facets = std::vector<std::vector<IMP::HalfFacet>>(this->master_node_ids_.size());
    for (const auto &facet : boundary_facets) {

        // Check how many nodes of the facet are master nodes.
        auto master_node_counter = std::size_t{0};
        for (const auto &nid : master_mesh.Cells(facet.CellId()).FacetConnectivity(facet.FacetId())) {
            if (master_nodes_flag[nid] == 1) {
                master_node_counter++;
            }
        }

        // If all the nodes of the facet are master nodes, then we consider it as a master facet.
        if (master_node_counter == master_mesh.Cells(facet.CellId()).FacetConnectivity(facet.FacetId()).Dim()) {
            for (const auto &nid : master_mesh.Cells(facet.CellId()).FacetConnectivity(facet.FacetId())) {
                // Add the facet if the node is a master node.
                master_facets[master_nodes_it.at(nid)].emplace_back(facet);
            }
        }
    }

    // Store the master facets connectivity.
    this->master_facets_conn_ = std::vector<std::vector<IMP::Vec<DIM,int>>>(this->master_node_ids_.size());
    for (std::size_t i = 0; i != master_facets.size(); ++i) {
        for (const auto &facet : master_facets[i]) {
            // Get the facet connectivity.
            auto conn = master_mesh.Cells(facet.CellId()).FacetConnectivity(facet.FacetId());

            // Compute the centroid of the facet's cell.
            auto xc = IMP::Vec<DIM,double>{};
            for (const auto &n : master_mesh.Cells(facet.CellId()).Connectivity()) {
                for (short d = 0; d != DIM; ++d)  {
                    xc[d] += master_mesh.Nodes(n)[d];
                }
            }
            xc /= static_cast<double>(CELL_NODES);

            if constexpr (DIM == 2) {
                auto n0 = conn[0];
                auto n1 = conn[1];

                // Get position of the facet nodes.
                auto x0 = IMP::Vec<DIM,double>({master_mesh.Nodes(n0)[0], master_mesh.Nodes(n0)[1]});
                auto x1 = IMP::Vec<DIM,double>({master_mesh.Nodes(n1)[0], master_mesh.Nodes(n1)[1]});

                // Compute normal vector to the facet.
                auto en = IMP::Vec<DIM,double>({x1[1]-x0[1], -(x1[0]-x0[0])});

                // Ensure normal vector is outward.
                auto ec = IMP::Vec<DIM,double>({xc[0]-0.5*(x0[0]+x1[0]), xc[1]-0.5*(x0[1]+x1[1])});
                if (en.Dot(ec) > 0.) {
                    // Normal is found inward and we flip the facet connectivity.
                    auto n_temp = n0;
                    n0 = n1;
                    n1 = n_temp;
                }

                // Store the connectivity of the facet.
                this->master_facets_conn_[i].emplace_back(IMP::Vec<DIM,int>({n0,n1}));

            } else if constexpr (DIM == 3) {
                auto n0 = this->master_node_ids_[i];
                auto n1 = int{-1};
                auto n2 = int{-1};
                for (const auto &n : conn) {
                    if (n != n0 && n1 != -1) { n2 = n; }
                    if (n != n0 && n1 == -1) { n1 = n; }
                }

                // Get position of the facet nodes.
                auto x0 = IMP::Vec<DIM,double>({master_mesh.Nodes(n0)[0], master_mesh.Nodes(n0)[1], master_mesh.Nodes(n0)[2]});
                auto x1 = IMP::Vec<DIM,double>({master_mesh.Nodes(n1)[0], master_mesh.Nodes(n1)[1], master_mesh.Nodes(n1)[2]});
                auto x2 = IMP::Vec<DIM,double>({master_mesh.Nodes(n2)[0], master_mesh.Nodes(n2)[1], master_mesh.Nodes(n2)[2]});

                // Compute facet edges.
                auto e01 = x1-x0;
                auto e02 = x2-x0;

                // Compute normal vector to the facet.
                auto en = e01.Cross(e02);

                // Ensure normal vector is outward.
                auto ec = IMP::Vec<DIM,double>({xc[0]-(x0[0]+x1[0]+x2[0])/3, xc[1]-(x0[1]+x1[1]+x2[1])/3, xc[2]-(x0[2]+x1[2]+x2[2])/3});
                if (en.Dot(ec) > 0.) {
                    // Normal is found inward and we flip the facet connectivity.
                    auto n_temp = n1;
                    n1 = n2;
                    n2 = n_temp;
                }

                // Store the connectivity of the facet.
                this->master_facets_conn_[i].emplace_back(IMP::Vec<DIM,int>({n0,n1,n2}));
            }
        }
    }
}


template<short DIM, short CELL_NODES>
void ContactHandle<DIM,CELL_NODES>::SetMaster(const IMP::Voronoi<DIM> &master_voro, const std::string master_nset_name)
{
    if (master_voro.NodeSets().find(master_nset_name) != master_voro.NodeSets().end()) {
        this->master_node_ids_ = master_voro.NodeSets().at(master_nset_name).NodeIds();
        this->master_node_normals_.clear();
        this->master_node_normals_.resize(this->master_node_ids_.size());
    } else {
        auto err_msg = "Could not set master in contact handle. The node set: "+master_nset_name+" does not exist in the master geometry.";
        throw std::runtime_error(Logger::Error(err_msg));
    }

    // Get the indices of the boundary facets.
    auto bound_facets_ids = std::vector<std::vector<int>>(this->master_node_ids_.size());
    auto nid = int{0};
    for (auto i = std::size_t{0}; i != this->master_node_ids_.size(); ++i) {
        // Search boundary surfaces of the voronoi cell enclosing the master node.
        nid = this->master_node_ids_[i];
        for (const auto &fid : master_voro.Cells(nid).Connectivity()) {
            if (master_voro.Facets(fid).NeighCellId() == -1) {
                bound_facets_ids[i].emplace_back(fid);
            }
        }
    }

    // Extract connectivity of the facets.
    this->master_facets_conn_ = std::vector<std::vector<IMP::Vec<DIM,int>>>(this->master_node_ids_.size());
    if constexpr (DIM == 2) {

    } else if constexpr (DIM == 3) {
        // We store the connectivity of subfacets (triangles
        // constructed by the node and the segments of the facet)
        auto n0 = int{0};       // first vertex of triangle is the master node index.
        auto p1 = int{0};       // second vertex of triangle is the first point of the segment.
        auto p2 = int{0};       // third vertex of triangle is the second point of the segment.

        // Normal vector of the master node.
        auto master_normal = IMP::Vec<DIM,double>{};

        // Normal vector of the triangle.
        auto en = IMP::Vec<DIM,double>{};

        auto approx_zero = std::sqrt(std::numeric_limits<double>::epsilon());

        // Find master facets connectivity for each master node.
        for (auto i = std::size_t{0}; i != this->master_node_ids_.size(); ++i) {
            // First vertex of the triangle.
            n0 = this->master_node_ids_[i];

            // Initial master node normal.
            master_normal.SetZero();

            // Iterate over the facets of the master node.
            auto facet_normal = IMP::Vec<DIM,double>{};
            for (const auto &fid : bound_facets_ids[i]) {
                facet_normal.SetZero();
                for (auto j = std::size_t(0); j != master_voro.Facets(fid).Connectivity().size(); ++j) {
                    // Second and third vertices of the triangle.
                    p1 = master_voro.Facets(fid).C(j);
                    p2 = master_voro.Facets(fid).C( (j+1) % master_voro.Facets(fid).Connectivity().size() );

                    // Skip triangle if it is degenerated.
                    if (std::sqrt(master_voro.Nodes(n0).Distance2(master_voro.Points(p1))) < 100*approx_zero ||
                            std::sqrt(master_voro.Nodes(n0).Distance2(master_voro.Points(p2))) < 100*approx_zero) {
                        continue;
                    }

                    // Ensure that triangle normal is outward.
                    en = (master_voro.Points(p1) - master_voro.Nodes(n0)).Cross(master_voro.Points(p2) - master_voro.Nodes(n0));
                    if (en.Dot(master_voro.CellCentroid(n0) - master_voro.Nodes(n0)) > approx_zero) {
                        // en is inward and we switch connectivity.
                        auto p_temp = p1;
                        p1 = p2;
                        p2 = p_temp;

                        // We flip en.
                        en *= -1.;
                    }

                    // Add contribution to master node normal.
                    facet_normal += en;

                    // Store the triangle connectivity.
                    this->master_facets_conn_[i].emplace_back(IMP::Vec<DIM,int>({n0,p1,p2}));
                }

                facet_normal /= master_voro.Facets(fid).Connectivity().size();
                master_normal += facet_normal;
            }
            master_normal /= bound_facets_ids[i].size();
            this->master_node_normals_[i] = master_normal;
        }


    } else {
        auto err_msg = "Could not set master contact body. Expected geometry of either 2 or 3 dimensions.\n";
        throw std::runtime_error(Logger::Error(err_msg));
    }
}


template<short DIM, short CELL_NODES>
void ContactHandle<DIM,CELL_NODES>::NearestMasterToSlave(const std::vector<IMP::Vec<DIM, double>> &master_nodes,
    const std::vector<IMP::Vec<DIM,double>> &slave_nodes, std::vector<int> &near_master_to_slave) const
{
    near_master_to_slave.clear();
    near_master_to_slave.resize(slave_nodes.size());

	// construct a kd-tree index:
	// Dimensionality set at run-time (default: L2)
	// ------------------------------------------------------------
	typedef KDTreeIMPVecVectorAdaptor<std::vector<IMP::Vec<DIM, double>>, double> my_kd_tree_t;
	my_kd_tree_t   mat_index(DIM, master_nodes, 10 /* max leaf */ );
	mat_index.index->buildIndex();

	// do a knn search
	auto ret_indexes = std::vector<size_t>(1);
	auto out_dists_sqr = std::vector<double>(1);
	auto resultSet = nanoflann::KNNResultSet<double>(1);

    auto i = int{0};
    for (const auto &query : slave_nodes) {
        resultSet.init(&ret_indexes[0], &out_dists_sqr[0] );
        mat_index.index->findNeighbors(resultSet, &query[0], nanoflann::SearchParams(10));
        near_master_to_slave[i++] = ret_indexes[0];
    }
}


} // End of namespace PNT

#endif //PHYNETOUCH_ENGINE_PHYSICS_CONTACT_HANDLE_TPP_
