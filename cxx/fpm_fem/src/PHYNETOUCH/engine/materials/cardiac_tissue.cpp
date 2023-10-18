/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */


#include "PHYNETOUCH/engine/materials/cardiac_tissue.hpp"


namespace PNT {

CardiacTissue::CardiacTissue() : density_lp_(), density_gp_(), specific_heat_lp_not_(), specific_heat_lp_(), specific_heat_gp_(),
    electrical_conductivity_not_(), electrical_conductivity_(), thermal_conductivity_not_(), thermal_conductivity_(),
    enthalpy_(), latent_heat_(), water_content_(), nodes_num_(0), thread_looper_(), threads_num_(0)
{
    // Get the number of available threads.
    std::size_t available = std::thread::hardware_concurrency()-1;
    this->threads_num_ = std::max(available, 1ul);
}


CardiacTissue::~CardiacTissue()
{}


void CardiacTissue::SetNodesNum(int nodes_num)
{
    // Check if given number of material nodes is positive.
    if (nodes_num < 0) {
        throw std::invalid_argument(Logger::Error("Could not set the number of nodes belonging to the cardiac tissue material. A negative number of material nodes was given."));
    }

    this->nodes_num_ = nodes_num;
}


void CardiacTissue::ComputeTempDependent(const Eigen::VectorXd &temp)
{
    // Set conductivities-not to zero if either of them is missing.
    if (this->electrical_conductivity_not_.size() == 0) {
        this->electrical_conductivity_not_.assign(temp.size(), 0.);
    }
    if (this->thermal_conductivity_not_.size() == 0) {
        this->thermal_conductivity_not_.assign(temp.size(), 0.);
    }

    // Initialize temperature-dependent conductivities.
    this->electrical_conductivity_.clear();
    this->electrical_conductivity_.resize(temp.size());
    this->thermal_conductivity_.clear();
    this->thermal_conductivity_.resize(temp.size());
    this->specific_heat_lp_.clear();
    this->specific_heat_lp_.resize(temp.size());

    // Compute temperature-dependent conductivities.
    for (auto id = 0; id != temp.size(); ++id) {
        if (temp.coeff(id) < 0.) {
            // auto error_msg = "Could not compute conductivities for tissue material point: " + std::to_string(id) + ". Temperature is below 0. °C.";
            // throw std::invalid_argument(Logger::Error(error_msg));
            this->specific_heat_lp_[id] = this->specific_heat_lp_not_[id];
            this->electrical_conductivity_[id] = this->electrical_conductivity_not_[id];
            this->thermal_conductivity_[id] = this->thermal_conductivity_not_[id];
        } else {
            // Set linearly temperature-dependent parameters of the tissue.
            this->specific_heat_lp_[id] = this->specific_heat_lp_not_[id] * (1. - 0.0042*(temp.coeff(id) - 37.));
            this->electrical_conductivity_[id] = this->electrical_conductivity_not_[id] * (1. + 0.015*(temp.coeff(id) - 37.));
            this->thermal_conductivity_[id] = this->thermal_conductivity_not_[id] * (1. - 0.0005*(temp.coeff(id) - 37.));
        }
    }

}


void CardiacTissue::ComputeEnthalpy(const Eigen::VectorXd &temp)
{
    if (this->density_lp_.size() == 0 || this->density_gp_.size() == 0) {
        auto error_msg = "Could not compute tissue enthalpy. Set density first.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }
    if (this->specific_heat_lp_.size() == 0 || this->specific_heat_gp_.size() == 0) {
        auto error_msg = "Could not compute tissue enthalpy. Set specific heat first.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }
    if (this->latent_heat_.size() == 0) {
        auto error_msg = "Could not compute tissue enthalpy. Set latent heat.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }
    if (this->water_content_.size() == 0) {
        auto error_msg = "Could not compute tissue enthalpy. Set specific water content first.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Initialize enthalpy container.
    this->enthalpy_.clear();
    this->enthalpy_.resize(temp.size());

    for (auto id = 0; id != temp.size(); ++id) {
        this->enthalpy_[id] = this->DensityLp(id)*this->SpecificHeatLp(id);

        // if (temp.coeff(id) < 0.) {
        //     auto error_msg = "Could not compute enthalpy for tissue material point: " + std::to_string(id) + ". Temperature is below 0 °C.";
        //     throw std::invalid_argument(Logger::Error(error_msg));
        // }
        // if (temp.coeff(id) <= 99.) {
            // this->enthalpy_[id] = this->DensityLp(id)*this->SpecificHeatLp(id);
        // } else if (temp.coeff(id) <= 100.) {
            // this->enthalpy_[id] = this->LatentHeat(id)*this->WaterContent(id);
        // } else {
            // this->enthalpy_[id] = this->DensityGp(id)*this->SpecificHeatGp(id);
        // }
    }
}


void CardiacTissue::SetDensityLp(const std::vector<double> &density)
{
    // Check if the nodes of the material have been set.
    if (this->nodes_num_ == 0) {
        throw std::runtime_error(Logger::Error("Could not set the density of the cardiac tissue material at liquid phase. The material's nodes number has not been set."));
    }

    // Set density. Either different for each node or the same.
    if (static_cast<int>(density.size()) == this->nodes_num_) {
        this->density_lp_ = density;
    } else if (density.size() == 1) {
        this->density_lp_.assign(this->nodes_num_,density[0]);
        this->density_lp_.shrink_to_fit();
    } else {
        std::string error_msg = "Could not set density for cardiac tissue material at liquid phase.\n"
                                "The given density vector should contain either a single value or a value for each node of the material.";
        throw std::runtime_error(Logger::Error(error_msg));
    }
}


void CardiacTissue::SetDensityGp(const std::vector<double> &density)
{
    // Check if the nodes of the material have been set.
    if (this->nodes_num_ == 0) {
        throw std::runtime_error(Logger::Error("Could not set the density of the cardiac tissue material at gas phase. The material's nodes number has not been set."));
    }

    // Set density. Either different for each node or the same.
    if (static_cast<int>(density.size()) == this->nodes_num_) {
        this->density_gp_ = density;
    } else if (density.size() == 1) {
        this->density_gp_.assign(this->nodes_num_,density[0]);
        this->density_gp_.shrink_to_fit();
    } else {
        std::string error_msg = "Could not set density for cardiac tissue material at gas phase.\n"
                                "The given density vector should contain either a single value or a value for each node of the material.";
        throw std::runtime_error(Logger::Error(error_msg));
    }
}


void CardiacTissue::SetSpecificHeatLpNot(const std::vector<double> &specific_heat)
{
    // Check if the nodes of the material have been set.
    if (this->nodes_num_ == 0) {
        throw std::runtime_error(Logger::Error("Could not set specific heat not for cardiac tissue material at liquid phase. The material's nodes number has not been set."));
    }

    // Set specific heat. Either different for each node or the same.
    if (static_cast<int>(specific_heat.size()) == this->nodes_num_) {
        this->specific_heat_lp_not_ = specific_heat;
    } else if (specific_heat.size() == 1) {
        this->specific_heat_lp_not_.assign(this->nodes_num_,specific_heat[0]);
        this->specific_heat_lp_not_.shrink_to_fit();
    } else {
        auto error_msg = "Could not set specific heat not for cardiac tissue material at liquid phase.\n"
            "The given specific heat vector should contain either a single value or a value for each node of the material.";
        throw std::runtime_error(Logger::Error(error_msg));
    }
}


void CardiacTissue::SetSpecificHeatGp(const std::vector<double> &specific_heat)
{
    // Check if the nodes of the material have been set.
    if (this->nodes_num_ == 0) {
        throw std::runtime_error(Logger::Error("Could not set specific heat for cardiac tissue material at gas phase. The material's nodes number has not been set."));
    }

    // Set specific heat. Either different for each node or the same.
    if (static_cast<int>(specific_heat.size()) == this->nodes_num_) {
        this->specific_heat_gp_ = specific_heat;
    } else if (specific_heat.size() == 1) {
        this->specific_heat_gp_.assign(this->nodes_num_,specific_heat[0]);
        this->specific_heat_gp_.shrink_to_fit();
    } else {
        auto error_msg = "Could not set specific heat for cardiac tissue material at gas phase.\n"
            "The given specific heat vector should contain either a single value or a value for each node of the material.";
        throw std::runtime_error(Logger::Error(error_msg));
    }
}


void CardiacTissue::SetElectricalConductivityNot(const std::vector<double> &conductivity)
{
    // Check if the nodes of the material have been set.
    if (this->nodes_num_ == 0) {
        throw std::runtime_error(Logger::Error("Could not set electrical conductivity not for cardiac tissue material. The material's nodes number has not been set."));
    }

    // Set electrical conductivity not. Either different for each node or the same.
    if (static_cast<int>(conductivity.size()) == this->nodes_num_) {
        this->electrical_conductivity_not_ = conductivity;
    } else if (conductivity.size() == 1) {
        this->electrical_conductivity_not_.assign(this->nodes_num_,conductivity[0]);
        this->electrical_conductivity_not_.shrink_to_fit();
    } else {
        auto error_msg = "Could not set electrical conductivity not for cardiac tissue material.\n"
            "The given electrical conductivity not vector should contain either a single value or a value for each node of the material.";
        throw std::runtime_error(Logger::Error(error_msg));
    }
}


void CardiacTissue::SetThermalConductivityNot(const std::vector<double> &conductivity)
{
    // Check if the nodes of the material have been set.
    if (this->nodes_num_ == 0) {
        throw std::runtime_error(Logger::Error("Could not set thermal conductivity not for cardiac tissue material. The material's nodes number has not been set."));
    }

    // Set thermal conductivity. Either different for each node or the same.
    if (static_cast<int>(conductivity.size()) == this->nodes_num_) {
        this->thermal_conductivity_not_ = conductivity;
    } else if (conductivity.size() == 1) {
        this->thermal_conductivity_not_.assign(this->nodes_num_,conductivity[0]);
        this->thermal_conductivity_not_.shrink_to_fit();
    } else {
        auto error_msg = "Could not set thermal conductivity not for cardiac tissue material.\n"
            "The given thermal conductivity not vector should contain either a single value or a value for each node of the material.";
        throw std::runtime_error(Logger::Error(error_msg));
    }
}


void CardiacTissue::SetLatentHeat(const std::vector<double> &latent_heat)
{
    // Check if the nodes of the material have been set.
    if (this->nodes_num_ == 0) {
        throw std::runtime_error(Logger::Error("Could not set latent heat for cardiac tissue material. The material's nodes number has not been set."));
    }

    // Set latent heat. Either different for each node or the same.
    if (static_cast<int>(latent_heat.size()) == this->nodes_num_) {
        this->latent_heat_ = latent_heat;
    } else if (latent_heat.size() == 1) {
        this->latent_heat_.assign(this->nodes_num_,latent_heat[0]);
        this->latent_heat_.shrink_to_fit();
    } else {
        std::string error_msg = "Could not set latent heat for cardiac tissue material.\n"
            "The given latent heat vector should contain either a single value or a value for each node of the material.";
        throw std::runtime_error(Logger::Error(error_msg));
    }
}


void CardiacTissue::SetWaterContent(const std::vector<double> &water_content)
{
    // Check if the nodes of the material have been set.
    if (this->nodes_num_ == 0) {
        throw std::runtime_error(Logger::Error("Could not set water content for cardiac tissue material. The material's nodes number has not been set."));
    }

    // Set water content. Either different for each node or the same.
    if (static_cast<int>(water_content.size()) == this->nodes_num_) {
        this->water_content_ = water_content;
    } else if (water_content.size() == 1) {
        this->water_content_.assign(this->nodes_num_,water_content[0]);
        this->water_content_.shrink_to_fit();
    } else {
        std::string error_msg = "Could not set water content for cardiac tissue material.\n"
            "The given water content vector should contain either a single value or a value for each node of the material.";
        throw std::runtime_error(Logger::Error(error_msg));
    }
}


} // End of namespace PNT
