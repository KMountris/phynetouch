/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */


#include "PHYNETOUCH/engine/materials/mooney_rivlin.hpp"

namespace PNT {

MooneyRivlin::MooneyRivlin()
{
    this->EditType(ConstitutiveType::MooneyRivlin);
}


MooneyRivlin::~MooneyRivlin()
{}


void MooneyRivlin::SetYoungModulus(const std::vector<double> &young_modulus)
{
    // Check if the nodes of the material have been set.
    if (this->NodesNum() == 0) {
        throw std::runtime_error(Logger::Error("Could not set Young modulus for MooneyRivlin material. The material's nodes number has not been set."));
    }

    // Set Young modulus. Either different for each node or the same.
    if (static_cast<int>(young_modulus.size()) == this->NodesNum()) {
        this->EditYoungModulus() = young_modulus;
    } else if (young_modulus.size() == 1) {
        this->EditYoungModulus().assign(this->NodesNum(),young_modulus[0]);
        this->EditYoungModulus().shrink_to_fit();
    } else {
        std::string err_msg = "Could not set Young modulus for MooneyRivlin material.\n"
            "The given Young modulus vector should contain either a single value or a value for each node of the material.";
        throw std::runtime_error(Logger::Error(err_msg));
    }
}


void MooneyRivlin::SetPoissonRatio(const std::vector<double> &poisson_ratio)
{
    // Check if the nodes of the material have been set.
    if (this->NodesNum() == 0) {
        throw std::runtime_error(Logger::Error("Could not set Poisson's ratio for MooneyRivlin material. The material's nodes number has not been set."));
    }

    // Set Poisson's ratio. Either different for each node or the same.
    if (static_cast<int>(poisson_ratio.size()) == this->NodesNum()) {
        this->EditPoissonRatio() = poisson_ratio;
    } else if (poisson_ratio.size() == 1) {
        this->EditPoissonRatio().assign(this->NodesNum(),poisson_ratio[0]);
        this->EditPoissonRatio().shrink_to_fit();
    } else {
        std::string err_msg = "Could not set Poisson's ratio for MooneyRivlin material.\n"
            "The given Poisson's ratio vector should contain either a single value or a value for each node of the material.";
        throw std::runtime_error(Logger::Error(err_msg));
    }
}


void MooneyRivlin::ComputeParameters()
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


void MooneyRivlin::ComputeWaveSpeed()
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


void MooneyRivlin::ComputeStressSPK(const Eigen::MatrixXd &deform_grad, const Eigen::MatrixXd &I,
    double mu, double kappa, Eigen::MatrixXd &stress) const
{
    if (deform_grad.rows() == 2) {
        auto C = Eigen::Matrix2d{};
        C.noalias() = deform_grad.transpose()*deform_grad;     // Right Cauchy Green
        auto Cinv = C.inverse();

        // Compute invariants
        auto I1 = C.trace();
        // I2 = 0.5*(trace(C)^2 - trace(C^2));
        auto I3 = C.determinant();

        // Compute reduced invariants
        //% J1 = I1*I3^(-1/3);
        //% J2 = I2*I3^(-2/3);
        auto J3 = std::sqrt(I3);

        // Calculate derivatives of invariants with respect to Lagrangian strain
        auto I1e = 2.*I;
        //% I2e = 2*(I1*eye(2)-C);
        auto I3e = 2.*I3*Cinv;

        // Calculate derivatives of reduced invariants
        auto J1e = std::pow(I3,(-1./3.))*I1e - (1./3.)*I1*std::pow(I3,(-4./3.))*I3e;
        // J2e = (I3^(-2/3))*I2e - ((2/3)*I2*I3^(-5/3))*I3e;
        auto J3e = 0.5*std::pow(I3,-0.5)*I3e;

        // Calculate Second Piola Kirchoff stress using Neo-Hookean material model
        stress = 0.5*mu*J1e + kappa*(J3-1.)*J3e;


    } else if (deform_grad.rows() == 3) {
        auto C = Eigen::Matrix3d{};
        C.noalias() = deform_grad.transpose()*deform_grad;     // Right Cauchy Green
        auto Cinv = C.inverse();

        // Compute invariants
        auto I1 = C.trace();
        // I2 = 0.5*(trace(C)^2 - trace(C^2));
        auto I3 = C.determinant();

        // Compute reduced invariants
        //% J1 = I1*I3^(-1/3);
        //% J2 = I2*I3^(-2/3);
        auto J3 = std::sqrt(I3);

        // Calculate derivatives of invariants with respect to Lagrangian strain
        auto I1e = 2.*I;
        //% I2e = 2*(I1*eye(2)-C);
        auto I3e = 2.*I3*Cinv;

        // Calculate derivatives of reduced invariants
        auto J1e = std::pow(I3,(-1./3.))*I1e - (1./3.)*I1*std::pow(I3,(-4./3.))*I3e;
        // J2e = (I3^(-2/3))*I2e - ((2/3)*I2*I3^(-5/3))*I3e;
        auto J3e = 0.5*std::pow(I3,-0.5)*I3e;

        // Calculate Second Piola Kirchoff stress using Neo-Hookean material model
        stress = 0.5*mu*J1e + kappa*(J3-1.)*J3e;

    } else {
        auto err_msg = "Could not compute StressSPK. Available only for either 2D or 3D.";
        throw std::runtime_error(Logger::Error(err_msg));
    }

}


} // End of namespace PNT