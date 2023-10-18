/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */


/**
   \file constitutive.hpp
   \brief Constitutive class header file.
   \author Konstantinos A. Mountris
   \date 26/10/2021
*/

#pragma once
#ifndef PHYNETOUCH_MATERIALS_CONSTITUTIVE_HPP_
#define PHYNETOUCH_MATERIALS_CONSTITUTIVE_HPP_


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
 * \enum ConstitutiveType
 * \author Konstantinos A. Mountris
 * \brief Enumeration to declare the type of the constitutive law.
 */
enum class ConstitutiveType {
    unknown,
    neohookean,        /**< NeoHookean hyperelastic body constitutive law */
    rigid              /**< Rigid body constitutive law */
};


/**
 * \class Constitutive
 * \author Konstantinos A. Mountris
 * \brief Abstract class implemmenting an interface to constitutive laws.
 */
class Constitutive
{

private:

    std::vector<double> density_;                       /**< The material's density */

    std::vector<double> young_modulus_;                 /**< The material's Young modulus */

    std::vector<double> poisson_ratio_;                 /**< The material's Poisson's ratio */

    std::vector<double> lame_lambda_;                   /**< The material's first Lame parameter */

    std::vector<double> lame_mu_;                       /**< The material's first Lame parameter */

    std::vector<double> bulk_modulus_;                  /**< The material's bulk modulus */

    std::vector<double> wave_speed_;                    /**< The material's wave speed */

    int nodes_num_;                                     /**< The number of the material's nodes */

    ConstitutiveType type_;                             /**< The type of the constitutive law */


public:

    /**
     * \brief The default constructor of the material.
     */
    Constitutive();


    /**
     * \brief The default destructor of the material.
     */
    virtual ~Constitutive();


    /**
     * \brief Set the number of the material's nodes.
     * \param [in] nodes_num The number of the nodes that belong to the material.
     * \return [void]
     */
    void SetNodesNum(int nodes_num);


    /**
     * \brief Set the density of the material.
     * \param [in] density The density of the material.
     * \return [void]
     */
    void SetDensity(const std::vector<double> &density);


    /**
     * \brief Set the Young modulus of the material.
     * \param [in] density The Young modulus of the material.
     * \return [void]
     */
    virtual void SetYoungModulus(const std::vector<double> &young_modulus) = 0;


    /**
     * \brief Set the Poisson's ratio of the material.
     * \param [in] density The Poisson's ratio of the material.
     * \return [void]
     */
    virtual void SetPoissonRatio(const std::vector<double> &poisson_ratio) = 0;


    /**
     * \brief Compute the material parameters - Lame parameters, bulk modulus.
     * \return [void]
     */
    virtual void ComputeParameters() = 0;


    /**
     * \brief Compute the material wave speed.
     * \return [void]
     */
    virtual void ComputeWaveSpeed() = 0;


    /**
     * \brief Compute the second Piola-Kirchhoff stress tensor.
     * \param [in] deform_grad The deformation gradient.
     * \param [in] identity The identity matrix.
     * \param [in] mu The Lame parameter mu.
     * \param [in] kappa The bulk modulus.
     * \param [out] stress The second Piola-Kirchhoff stress tensor.
     * \return [void]
     */
    virtual void ComputeStressSPK(const Eigen::MatrixXd &deform_grad, const Eigen::MatrixXd &identity,
        double mu, double kappa, Eigen::MatrixXd &stress) const = 0;


    /**
     * \brief Get the constitutive law type of the material.
     * \return [ConstitutiveType] The constitutive law type of the material. 
     */
    inline auto Type() const { return this->type_; }


    /**
     * \brief Get the number of nodes belonging to the material.
     * \return [int] The number of nodes belonging to the material.
     */
    inline int NodesNum() const { return this->nodes_num_; }


    /**
     * \brief Get the density of the material.
     * \return [const std::vector<double>&] The density of the material.
     */
    inline const std::vector<double> & Density() const { return this->density_; }


    /**
     * \brief Get the density of a specific point of the material.
     * \param [in] id The index of the material's point.
     * \return [double] The density of a specific point of the material.
     */
    inline double Density(std::size_t id) const { return this->density_[id]; }


    /**
     * \brief Get the Young modulus of the material.
     * \return [const std::vector<double>&] The Young modulus of the material.
     */
    inline const std::vector<double> & YoungModulus() const { return this->young_modulus_; }


    /**
     * \brief Get the Young modulus of a specific point of the material.
     * \param [in] id The index of the material's point.
     * \return [double] The Young modulus of a specific point of the material.
     */
    inline double YoungModulus(std::size_t id) const { return this->young_modulus_[id]; }


    /**
     * \brief Get the Poisson's ratio of the material.
     * \return [const std::vector<double>&] The Poisson's ratio of the material.
     */
    inline const std::vector<double> & PoissonRatio() const { return this->poisson_ratio_; }


    /**
     * \brief Get the Poisson's ratio of a specific point of the material.
     * \param [in] id The index of the material's point.
     * \return [double] The Poisson's ratio of a specific point of the material.
     */
    inline double PoissonRatio(std::size_t id) const { return this->poisson_ratio_[id]; }


    /**
     * \brief Get the Lame's first parameter of the material.
     * \return [const std::vector<double>&] The Lame's first parameter of the material.
     */
    inline const std::vector<double> & LameLambda() const { return this->lame_lambda_; }


    /**
     * \brief Get the Lame's first parameter of a specific point of the material.
     * \param [in] id The index of the material's point.
     * \return [double] The Lame's first parameter of a specific point of the material.
     */
    inline double LameLambda(std::size_t id) const { return this->lame_lambda_[id]; }


    /**
     * \brief Get the Lame's second parameter of the material.
     * \return [const std::vector<double>&] The Lame's second parameter of the material.
     */
    inline const std::vector<double> & LameMu() const { return this->lame_mu_; }


    /**
     * \brief Get the Lame's second parameter of a specific point of the material.
     * \param [in] id The index of the material's point.
     * \return [double] The Lame's second parameter of a specific point of the material.
     */
    inline double LameMu(std::size_t id) const { return this->lame_mu_[id]; }


    /**
     * \brief Get the bulk modulus of the material.
     * \return [const std::vector<double>&] The bulk modulus of the material.
     */
    inline const std::vector<double> & BulkModulus() const { return this->bulk_modulus_; }


    /**
     * \brief Get the bulk modulus of a specific point of the material.
     * \param [in] id The index of the material's point.
     * \return [double] The bulk modulus of a specific point of the material.
     */
    inline double BulkModulus(std::size_t id) const { return this->bulk_modulus_[id]; }


    /**
     * \brief Get the wave speed of the material.
     * \return [const std::vector<double>&] The wave speed of the material.
     */
    inline const std::vector<double> & WaveSpeed() const { return this->wave_speed_; }


    /**
     * \brief Get the wave speed of a specific point of the material.
     * \param [in] id The index of the material's point.
     * \return [double] The wave speed of a specific point of the material.
     */
    inline double WaveSpeed(std::size_t id) const { return this->wave_speed_[id]; }



protected:

    inline void EditType(ConstitutiveType type) { this->type_ = type; }
    

    /**
     * \brief Get the Young modulus of the material.
     * \return [const std::vector<double>&] The Young modulus of the material.
     */
    inline auto & EditYoungModulus() { return this->young_modulus_; }


    /**
     * \brief Get the Young modulus of a specific point of the material.
     * \param [in] id The index of the material's point.
     * \return [double] The Young modulus of a specific point of the material.
     */
    inline auto & EditYoungModulus(std::size_t id) { return this->young_modulus_[id]; }


    /**
     * \brief Get the Poisson's ratio of the material.
     * \return [const std::vector<double>&] The Poisson's ratio of the material.
     */
    inline auto & EditPoissonRatio() { return this->poisson_ratio_; }


    /**
     * \brief Get the Poisson's ratio of a specific point of the material.
     * \param [in] id The index of the material's point.
     * \return [double] The Poisson's ratio of a specific point of the material.
     */
    inline auto & EditPoissonRatio(std::size_t id) { return this->poisson_ratio_[id]; }


    /**
     * \brief Get the Lame's first parameter of the material.
     * \return [const std::vector<double>&] The Lame's first parameter of the material.
     */
    inline auto & EditLameLambda() { return this->lame_lambda_; }


    /**
     * \brief Get the Lame's first parameter of a specific point of the material.
     * \param [in] id The index of the material's point.
     * \return [double] The Lame's first parameter of a specific point of the material.
     */
    inline auto & EditLameLambda(std::size_t id) { return this->lame_lambda_[id]; }


    /**
     * \brief Get the Lame's second parameter of the material.
     * \return [const std::vector<double>&] The Lame's second parameter of the material.
     */
    inline auto & EditLameMu() { return this->lame_mu_; }


    /**
     * \brief Get the Lame's second parameter of a specific point of the material.
     * \param [in] id The index of the material's point.
     * \return [double] The Lame's second parameter of a specific point of the material.
     */
    inline auto & EditLameMu(std::size_t id) { return this->lame_mu_[id]; }


    /**
     * \brief Get the bulk modulus of the material.
     * \return [const std::vector<double>&] The bulk modulus of the material.
     */
    inline auto & EditBulkModulus() { return this->bulk_modulus_; }


    /**
     * \brief Get the bulk modulus of a specific point of the material.
     * \param [in] id The index of the material's point.
     * \return [double] The bulk modulus of a specific point of the material.
     */
    inline auto & EditBulkModulus(std::size_t id) { return this->bulk_modulus_[id]; }


    /**
     * \brief Get the wave speed of the material.
     * \return [const std::vector<double>&] The wave speed of the material.
     */
    inline auto & EditWaveSpeed() { return this->wave_speed_; }


    /**
     * \brief Get the wave speed of a specific point of the material.
     * \param [in] id The index of the material's point.
     * \return [double] The wave speed of a specific point of the material.
     */
    inline auto & EditWaveSpeed(std::size_t id) { return this->wave_speed_[id]; }


};

/** \} End of Doxygen Groups */

} // End of namespace PNT

#endif //PHYNETOUCH_MATERIALS_CONSTITUTIVE_HPP_