/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */


#include "PHYNETOUCH/engine/materials/rigid.hpp"

namespace PNT {

Rigid::Rigid()
{
    this->EditType(ConstitutiveType::rigid);
}


Rigid::~Rigid()
{}


void Rigid::SetYoungModulus(const std::vector<double> &young_modulus)
{
    // Check if the nodes of the material have been set.
    if (this->NodesNum() == 0) {
        throw std::runtime_error(Logger::Error("Could not set Young modulus for rigid material. The material's nodes number has not been set."));
    }

    // Set Young modulus to zero.
    this->EditYoungModulus().assign(this->NodesNum(),0.*young_modulus[0]);
}


void Rigid::SetPoissonRatio(const std::vector<double> &poisson_ratio)
{
    // Check if the nodes of the material have been set.
    if (this->NodesNum() == 0) {
        throw std::runtime_error(Logger::Error("Could not set Poisson's ratio for rigid material. The material's nodes number has not been set."));
    }

    // Set Poisson's ratio to zero.
    this->EditPoissonRatio().assign(this->NodesNum(),0.*poisson_ratio[0]);
}


void Rigid::ComputeParameters()
{
    // Check if the nodes of the material and elastic parameters have been set.
    if (this->NodesNum() == 0)
        throw std::runtime_error(Logger::Error("Could not compute parameters. The material's nodes number has not been set."));

    // Set the parameters to zero.
    this->EditLameLambda().assign(this->NodesNum(), 0.);
    this->EditLameMu().assign(this->NodesNum(), 0.);
    this->EditBulkModulus().assign(this->NodesNum(), 0.);
}


void Rigid::ComputeWaveSpeed()
{
    // Check if the nodes of the material and elastic parameters have been set.
    if (this->NodesNum() == 0)
        throw std::runtime_error(Logger::Error("Could not compute wave speed. The material's nodes number has not been set."));

    // Set material wave speed to zero.
    this->EditWaveSpeed().assign(this->NodesNum(), 0.);
}


void Rigid::ComputeStressSPK(const Eigen::MatrixXd &deform_grad, const Eigen::MatrixXd &identity,
    double mu, double kappa, Eigen::MatrixXd &stress) const
{
    // Set SPK stress to identity.
    stress = identity + 0.*mu*kappa*deform_grad;
    stress.setZero();
}


} // End of namespace PNT