/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */


/**
   \file cardiac_tissue.hpp
   \brief CardiacTissue class header file.
   \author Konstantinos A. Mountris
   \date 13/07/2021
*/

#ifndef PHYNETOUCH_MATERIALS_CARDIAC_TISSUE_HPP_
#define PHYNETOUCH_MATERIALS_CARDIAC_TISSUE_HPP_


#include "PHYNETOUCH/engine/utilities/logger.hpp"

#include <CLOUDEA/CLOUDEA>
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
#include <thread>


namespace PNT {

/** \addtogroup Materials \{ */


/**
 * \class CardiacTissue
 * \author Konstantinos A. Mountris
 * \brief Class implemmenting a cardiac tissue material.
 */
class CardiacTissue
{

private:

    std::vector<double> density_lp_;                    /**< The cardiac tissue's density at liquid phase */

    std::vector<double> density_gp_;                    /**< The cardiac tissue's density at gas phase */

    std::vector<double> specific_heat_lp_not_;              /**< The cardiac tissue's specific heat not at liquid phase */

    std::vector<double> specific_heat_lp_;              /**< The cardiac tissue's specific heat at liquid phase */

    std::vector<double> specific_heat_gp_;              /**< The cardiac tissue's specific heat at gas phase */

    std::vector<double> electrical_conductivity_not_;   /**< The cardiac tissue's electrical conductivity not */

    std::vector<double> electrical_conductivity_;   /**< The cardiac tissue's electrical conductivity */

    std::vector<double> thermal_conductivity_not_;      /**< The cardiac tissue's thermal conductivity not */

    std::vector<double> thermal_conductivity_;      /**< The cardiac tissue's thermal conductivity */

    std::vector<double> enthalpy_;                      /**< The cardiac tissue's enthalpy */

    std::vector<double> latent_heat_;                   /**< The energy that is absorbed/released by the tissue during a change in its phase */

    std::vector<double> water_content_;                 /**< The content of the tissue in water */

    int nodes_num_;                                     /**< The number of the material's nodes */

    CLOUDEA::ThreadLoopManager thread_looper_;          /**< The managing object for the multithreaded loop execution */

    std::size_t threads_num_;                           /**< The number of threads for parallel execution of the monodomain model */


public:

    /**
     * \brief CardiacTissue default constructor.
     */
    CardiacTissue();


    /**
     * \brief CardiacTissue default destructor.
     */
    virtual ~CardiacTissue();


    /**
     * \brief Set the number of the electric material's nodes.
     * \param [in] nodes_num The number of the nodes that belong to the electric material.
     * \return [void]
     */
    void SetNodesNum(int nodes_num);


    /**
     * \brief Compute the temperature-dependent parameters of the tissue (specific heat at liquid phase, electrical & thermal conductivities).
     * \param [in] temp The temperature value of the tissue given in °C (degrees Celsius).
     * \return [void]
     */
    void ComputeTempDependent(const Eigen::VectorXd &temp);


    /**
     * \brief Compute the temperature-dependent enthalpy of the tissue.
     * \param [in] temp The temperature value of the tissue given in °C (degrees Celsius).
     * \return [void]
     */
    void ComputeEnthalpy(const Eigen::VectorXd &temp);


    /**
     * \brief Set the density of the cardiac tissue at liquid phase.
     * \param [in] density The density of the cardiac tissue at liquid phase.
     * \return [void]
     */
    void SetDensityLp(const std::vector<double> &density);


    /**
     * \brief Set the density of the cardiac tissue at gas phase.
     * \param [in] density The density of the cardiac tissue at gas phase.
     * \return [void]
     */
    void SetDensityGp(const std::vector<double> &density);


    /**
     * \brief Set the specific heat of the cardiac tissue at liquid phase.
     * \param [in] specific_heat The specific heat of the cardiac tissue at liquid phase.
     * \return [void]
     */
    void SetSpecificHeatLpNot(const std::vector<double> &specific_heat);


    /**
     * \brief Set the specific heat of the cardiac tissue at gas phase.
     * \param [in] specific_heat The specific heat of the cardiac tissue at gas phase.
     * \return [void]
     */
    void SetSpecificHeatGp(const std::vector<double> &specific_heat);


    /**
     * \brief Set the electrical conductivity not of the cardiac tissue.
     * \param [in] conductivity The electrical conductivity of the cardiac tissue.
     * \return [void]
     */
    void SetElectricalConductivityNot(const std::vector<double> &conductivity);


    /**
     * \brief Set the thermal conductivity not of the cardiac tissue.
     * \param [in] conductivity The thermal conductivity of the cardiac tissue.
     * \return [void]
     */
    void SetThermalConductivityNot(const std::vector<double> &conductivity);


    /**
     * \brief Set the latent heat of the cardiac tissue.
     * \param [in] latent_heat The latent heat of the cardiac tissue.
     * \return [void]
     */
    void SetLatentHeat(const std::vector<double> &latent_heat);


    /**
     * \brief Set the water content of the cardiac tissue.
     * \param [in] water_content The water content of the cardiac tissue.
     * \return [void]
     */
    void SetWaterContent(const std::vector<double> &water_content);


    /**
     * \brief Get the density at liquid phase of the material.
     * \return [const std::vector<double>&] The density at liquid phase of the material.
     */
    inline auto & DensityLp() const { return this->density_lp_; }


    /**
     * \brief Get the density at liquid phase of a specific point of the material.
     * \param [in] id The index of the material's point.
     * \return [double] The density at liquid phase of a specific point of the material.
     */
    inline auto DensityLp(std::size_t id) const { return this->density_lp_[id]; }


    /**
     * \brief Get the density at gas phase of the material.
     * \return [const std::vector<double>&] The density at gas phase of the material.
     */
    inline auto & DensityGp() const { return this->density_gp_; }


    /**
     * \brief Get the density at gas phase of a specific point of the material.
     * \param [in] id The index of the material's point.
     * \return [double] The density at gas phase of a specific point of the material.
     */
    inline auto DensityGp(std::size_t id) const { return this->density_gp_[id]; }


    /**
     * \brief Get the specific heat not at liquid phase of the material.
     * \return [const std::vector<double>&] The specific heat not at liquid phase value.
     */
    inline auto & SpecificHeatLpNot() const { return this->specific_heat_lp_not_; }


    /**
     * \brief Get the specific heat not at liquid phase of a specific point of the material.
     * \param [in] id The index of the material's point.
     * \return [double] The specific heat not at liquid phase of a specific point of the material.
     */
    inline auto SpecificHeatLpNot(std::size_t id) const { return this->specific_heat_lp_not_[id]; }


    /**
     * \brief Get the specific heat at liquid phase of the material.
     * \return [const std::vector<double>&] The specific heat at liquid phase value.
     */
    inline auto & SpecificHeatLp() const { return this->specific_heat_lp_; }


    /**
     * \brief Get the specific heat at liquid phase of a specific point of the material.
     * \param [in] id The index of the material's point.
     * \return [double] The specific heat at liquid phase of a specific point of the material.
     */
    inline auto SpecificHeatLp(std::size_t id) const { return this->specific_heat_lp_[id]; }


    /**
     * \brief Get the specific heat at gas phase of the material.
     * \return [const std::vector<double>&] The specific heat at gas phase value.
     */
    inline auto & SpecificHeatGp() const { return this->specific_heat_gp_; }


    /**
     * \brief Get the specific heat at gas phase of a specific point of the material.
     * \param [in] id The index of the material's point.
     * \return [double] The specific heat at gas phase of a specific point of the material.
     */
    inline auto SpecificHeatGp(std::size_t id) const { return this->specific_heat_gp_[id]; }


    /**
     * \brief Get the electrical conductivity not of the material.
     * \return [const std::vector<double>&] The electrical conductivity not value.
     */
    inline auto & ElectricalConductivityNot() const { return this->electrical_conductivity_not_; }


    /**
     * \brief Get the electrical conductivity not of a specific point of the material.
     * \param [in] id The index of the material's point.
     * \return [double] The electrical conductivity not of a specific point of the material.
     */
    inline auto ElectricalConductivityNot(std::size_t id) const { return this->electrical_conductivity_not_[id]; }


    /**
     * \brief Get the electrical conductivity of the material.
     * \return [const std::vector<double>&] The electrical conductivity value.
     */
    inline auto & ElectricalConductivity() const { return this->electrical_conductivity_; }


    /**
     * \brief Get the electrical conductivity of a specific point of the material.
     * \param [in] id The index of the material's point.
     * \return [double] The electrical conductivity of a specific point of the material.
     */
    inline auto ElectricalConductivity(std::size_t id) const { return this->electrical_conductivity_[id]; }


    /**
     * \brief Get the thermal conductivity not of the material.
     * \return [const std::vector<double>&] The thermal conductivity not value.
     */
    inline auto & ThermalConductivityNot() const { return this->thermal_conductivity_not_; }


    /**
     * \brief Get the thermal conductivity not of a specific point of the material.
     * \param [in] id The index of the material's point.
     * \return [double] The thermal conductivity not of a specific point of the material.
     */
    inline auto ThermalConductivityNot(std::size_t id) const { return this->thermal_conductivity_not_[id]; }


    /**
     * \brief Get the thermal conductivity not of the material.
     * \return [const std::vector<double>&] The thermal conductivity not value.
     */
    inline auto & ThermalConductivity() const { return this->thermal_conductivity_; }


    /**
     * \brief Get the thermal conductivity of a specific point of the material.
     * \param [in] id The index of the material's point.
     * \return [double] The thermal conductivity of a specific point of the material.
     */
    inline auto ThermalConductivity(std::size_t id) const { return this->thermal_conductivity_[id]; }


    /**
     * \brief Get the temperature-dependent enthalpy of the material.
     * \return [const std::vector<double>&] The temperature-dependent enthalpy value.
     */
    inline auto & Enthalpy() const { return this->enthalpy_; }


    /**
     * \brief Get the temperature-dependent enthalpy of a specific point of the material.
     * \param [in] id The index of the material's point.
     * \return [double] The temperature-dependent enthalpy of a specific point of the material.
     */
    inline auto Enthalpy(std::size_t id) const { return this->enthalpy_[id]; }


    /**
     * \brief Get the latent heat of the material.
     * \return [const std::vector<double>&] The latent heat value.
     */
    inline auto & LatentHeat() const { return this->latent_heat_; }


    /**
     * \brief Get the latent heat of a specific point of the material.
     * \param [in] id The index of the material's point.
     * \return [double] The latent heat of a specific point of the material.
     */
    inline auto LatentHeat(std::size_t id) const { return this->latent_heat_[id]; }


    /**
     * \brief Get the water content of the material.
     * \return [const std::vector<double>&] The water content value.
     */
    inline auto & WaterContent() const { return this->water_content_; }


    /**
     * \brief Get the latent heat of a specific point of the material.
     * \param [in] id The index of the material's point.
     * \return [double] The latent heat of a specific point of the material.
     */
    inline auto WaterContent(std::size_t id) const { return this->water_content_[id]; }


    /**
     * \brief Get the number of nodes belonging to the electric material.
     * \return [int] The number of nodes belonging to the electric material.
     */
    inline auto NodesNum() const { return this->nodes_num_; }

};

/** \} End of Doxygen Groups */

} // End of namespace PNT

#endif //PHYNETOUCH_MATERIALS_CARDIAC_TISSUE_HPP_