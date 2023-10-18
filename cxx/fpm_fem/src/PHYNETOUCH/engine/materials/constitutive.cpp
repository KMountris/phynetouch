/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */


#include "PHYNETOUCH/engine/materials/constitutive.hpp"

namespace PNT {

Constitutive::Constitutive() : density_(), young_modulus_(), poisson_ratio_(),
    lame_lambda_(), lame_mu_(), bulk_modulus_(), wave_speed_(), nodes_num_(0), type_(ConstitutiveType::unknown)
{}


Constitutive::~Constitutive()
{}


void Constitutive::SetNodesNum(int nodes_num)
{
    // Check if given number of material nodes is positive.
    if (nodes_num < 0) {
        throw std::invalid_argument(Logger::Error("Could not set the number of nodes belonging to the material. A negative number of material nodes was given."));
    }

    this->nodes_num_ = nodes_num;
}


void Constitutive::SetDensity(const std::vector<double> &density)
{
    // Check if the nodes of the material have been set.
    if (this->nodes_num_ == 0) {
        throw std::runtime_error(Logger::Error("Could not set the density of the material. The material's nodes number has not been set."));
    }

    // Set density. Either different for each node or the same.
    if (static_cast<int>(density.size()) == this->nodes_num_) {
        this->density_ = density;
    } else if (density.size() == 1) {
        this->density_.assign(this->nodes_num_,density[0]);
        this->density_.shrink_to_fit();
    } else {
        std::string err_msg = "Could not set density for material.\n"
            "The given density vector should contain either a single value or a value for each node of the material.";
        throw std::runtime_error(Logger::Error(err_msg));
    }
}

} // End of namespace PNT