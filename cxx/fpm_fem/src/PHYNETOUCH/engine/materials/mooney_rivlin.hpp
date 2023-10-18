/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */


/**
   \file mooney_rivlin.hpp
   \brief MooneyRivlin class header file.
   \author Konstantinos A. Mountris
   \date 26/10/2021
*/

#ifndef PHYNETOUCH_MATERIALS_MOONEY_RIVLIN_HPP_
#define PHYNETOUCH_MATERIALS_MOONEY_RIVLIN_HPP_

#pragma once
#include "PHYNETOUCH/engine/materials/constitutive.hpp"
#include "PHYNETOUCH/engine/utilities/logger.hpp"

#include <Eigen/Dense>

#include <iostream>
#include <initializer_list>
#include <iterator>
#include <vector>
#include <numeric>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <exception>


namespace PNT {

/** \addtogroup Materials \{ */


/**
 * \class MooneyRivlin
 * \author Konstantinos A. Mountris
 * \brief Class implemmenting a Mooney Rivlin hyperelastic material.
 */
class MooneyRivlin : public Constitutive
{

public:

    /**
     * \brief The default constructor of the material.
     */
    MooneyRivlin();


    /**
     * \brief The default destructor of the material.
     */
    virtual ~MooneyRivlin();


    /**
     * \brief Set the Young modulus of the material.
     * \param [in] density The Young modulus of the material.
     * \return [void]
     */
    void SetYoungModulus(const std::vector<double> &young_modulus);


    /**
     * \brief Set the Poisson's ratio of the material.
     * \param [in] density The Poisson's ratio of the material.
     * \return [void]
     */
    void SetPoissonRatio(const std::vector<double> &poisson_ratio);


    /**
     * \brief Compute the material parameters - Lame parameters, bulk modulus.
     * \return [void]
     */
    void ComputeParameters();


    /**
     * \brief Compute the material wave speed.
     * \return [void]
     */
    void ComputeWaveSpeed();


    /**
     * \brief Compute the second Piola-Kirchhoff stress tensor.
     * \param [in] deform_grad The deformation gradient.
     * \param [in] identity The identity matrix.
     * \param [in] mu The Lame parameter mu.
     * \param [in] kappa The bulk modulus.
     * \param [out] stress The second Piola-Kirchhoff stress tensor.
     * \return [void]
     */
    void ComputeStressSPK(const Eigen::MatrixXd &deform_grad, const Eigen::MatrixXd &identity,
        double mu, double kappa, Eigen::MatrixXd &stress) const;


};

/** \} End of Doxygen Groups */

} // End of namespace PNT

#endif //PHYNETOUCH_MATERIALS_MOONEY_RIVLIN_HPP_