/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */


#include "PHYNETOUCH/engine/conditions/prescribed_bc.hpp"


namespace PNT {


PrescribedBc::PrescribedBc() :
    node_ids_(), values_()
{}


PrescribedBc::~PrescribedBc() {}


void PrescribedBc::ConnectBodies(const std::vector<int> &connector_node_ids,
    const std::vector<std::vector<int>> &connected_node_ids,
    const Eigen::VectorXd &connector_values, int connected_body_nodes_num)
{
    // Reset prescribed bc.
    this->node_ids_.clear();
    this->values_.clear();

    // Create tissue nodes values vector and occurrences counter vector.
    auto values_tis = std::vector<double>(connected_body_nodes_num, 0.0);
    auto tis_nodes_occur = std::vector<int>(connected_body_nodes_num, 0);

    // Prescribe values of catheter connector nodes to connected tissue nodes.
    auto val = double{0.0};
    auto id_cth = int{0};
    for (std::size_t i = 0; i != connector_node_ids.size(); ++i) {
        id_cth = connector_node_ids[i];
        val = connector_values.coeff(id_cth);

        // Assign catheter node value to connected tissue nodes.
        for (const auto &id_tis : connected_node_ids[i]) {
            values_tis[id_tis] += val;
            tis_nodes_occur[id_tis] += 1;
        }
    }

    // Set the prescribed bc.
    for (int i = 0; i != connected_body_nodes_num; ++i) {
        if (tis_nodes_occur[i] != 0) {
            this->node_ids_.emplace_back(i);
            this->values_.emplace_back(values_tis[i] / tis_nodes_occur[i]);
        }
    }
    this->node_ids_.shrink_to_fit();
    this->values_.shrink_to_fit();
}


} // End of namespace PNT
