/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */


#include "PHYNETOUCH/engine/materials/neohookean.hpp"

namespace PNT {

Neohookean::Neohookean()
{
    this->EditType(ConstitutiveType::neohookean);
}


Neohookean::~Neohookean()
{}


void Neohookean::SetYoungModulus(const std::vector<double> &young_modulus)
{
    // Check if the nodes of the material have been set.
    if (this->NodesNum() == 0) {
        throw std::runtime_error(Logger::Error("Could not set Young modulus for neohookean material. The material's nodes number has not been set."));
    }

    // Set Young modulus. Either different for each node or the same.
    if (static_cast<int>(young_modulus.size()) == this->NodesNum()) {
        this->EditYoungModulus() = young_modulus;
    } else if (young_modulus.size() == 1) {
        this->EditYoungModulus().assign(this->NodesNum(),young_modulus[0]);
        this->EditYoungModulus().shrink_to_fit();
    } else {
        std::string err_msg = "Could not set Young modulus for neohookean material.\n"
            "The given Young modulus vector should contain either a single value or a value for each node of the material.";
        throw std::runtime_error(Logger::Error(err_msg));
    }
}


void Neohookean::SetPoissonRatio(const std::vector<double> &poisson_ratio)
{
    // Check if the nodes of the material have been set.
    if (this->NodesNum() == 0) {
        throw std::runtime_error(Logger::Error("Could not set Poisson's ratio for neohookean material. The material's nodes number has not been set."));
    }

    // Set Poisson's ratio. Either different for each node or the same.
    if (static_cast<int>(poisson_ratio.size()) == this->NodesNum()) {
        this->EditPoissonRatio() = poisson_ratio;
    } else if (poisson_ratio.size() == 1) {
        this->EditPoissonRatio().assign(this->NodesNum(),poisson_ratio[0]);
        this->EditPoissonRatio().shrink_to_fit();
    } else {
        std::string err_msg = "Could not set Poisson's ratio for neohookean material.\n"
            "The given Poisson's ratio vector should contain either a single value or a value for each node of the material.";
        throw std::runtime_error(Logger::Error(err_msg));
    }
}


void Neohookean::ComputeParameters()
{
    // Check if the nodes of the material and elastic parameters have been set.
    if (this->NodesNum() == 0)
        throw std::runtime_error(Logger::Error("Could not compute parameters. The material's nodes number has not been set."));
    if (static_cast<int>(this->YoungModulus().size()) != this->NodesNum())
        throw std::runtime_error(Logger::Error("Could not compute parameters. The material's Young modulus has not been set for each node."));
    if (static_cast<int>(this->PoissonRatio().size()) != this->NodesNum())
        throw std::runtime_error(Logger::Error("Could not compute parameters. The material's Poisson ratio has not been set for each node."));

    // Initialize the parameters.
    this->EditLameLambda().assign(this->NodesNum(), 0.);
    this->EditLameMu().assign(this->NodesNum(), 0.);
    this->EditBulkModulus().assign(this->NodesNum(), 0.);

    // Compute the parameters.
    for (int i = 0; i != this->NodesNum(); ++i) {
        this->EditLameLambda(i) = (this->YoungModulus(i)*this->PoissonRatio(i)) /
            ((1. + this->PoissonRatio(i))*(1. - 2.*this->PoissonRatio(i)));

        this->EditLameMu(i) = this->YoungModulus(i) / (2.*(1. + this->PoissonRatio(i)));

        this->EditBulkModulus(i) = this->YoungModulus(i) / (3.*(1. - 2.*this->PoissonRatio(i)));
    }
}


void Neohookean::ComputeWaveSpeed()
{
    // Check if the nodes of the material and elastic parameters have been set.
    if (this->NodesNum() == 0)
        throw std::runtime_error(Logger::Error("Could not compute wave speed. The material's nodes number has not been set."));
    if (static_cast<int>(this->LameLambda().size()) != this->NodesNum())
        throw std::runtime_error(Logger::Error("Could not compute parameters. The material's parameters have not been set for each node."));
    if (static_cast<int>(this->Density().size()) != this->NodesNum())
        throw std::runtime_error(Logger::Error("Could not compute parameters. The material's density has not been set for each node."));

    // Compute material wave speed
    this->EditWaveSpeed().assign(this->NodesNum(), 0.);
    for (int i = 0; i != this->NodesNum(); ++i) {
        this->EditWaveSpeed(i) = std::sqrt((this->LameLambda(i) + 2.*this->LameMu(i)) / this->Density(i));
    }
}


void Neohookean::ComputeStressSPK(const Eigen::MatrixXd &FT, const Eigen::MatrixXd &I,
    double mu, double kappa, Eigen::MatrixXd &stress) const
{
    double J = FT.determinant();
    Eigen::MatrixXd C = FT*FT.transpose();     // Right Cauchy Green
    Eigen::MatrixXd Cinv = C.inverse();

    stress = mu * std::pow(J, -(2./3.)) * (I - C.trace()/3.*Cinv) + kappa * J * (J - 1.) * Cinv;

}


} // End of namespace PNT