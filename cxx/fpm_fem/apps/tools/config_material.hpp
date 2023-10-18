/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */
/**
   \file config_material.hpp
   \brief ConfigMaterial class header file.
   \author Konstantinos A. Mountris
   \date 11/07/2021
*/

#pragma once
#ifndef PHYNETOUCH_APPS_TOOLS_CONFIG_MATERIAL_HPP_
#define PHYNETOUCH_APPS_TOOLS_CONFIG_MATERIAL_HPP_

#include "parser.hpp"

#include "PHYNETOUCH/engine/utilities/logger.hpp"
#include "PHYNETOUCH/engine/utilities/measure_units.hpp"
#include "PHYNETOUCH/engine/utilities/mpi_handler.hpp"
#include "PHYNETOUCH/engine/materials/cardiac_tissue.hpp"
#include "PHYNETOUCH/engine/materials/catheter.hpp"
#include "PHYNETOUCH/engine/materials/neohookean.hpp"

#include <IMP/IMP>
#include <termcolor/termcolor.hpp>

#include <Eigen/Dense>

#include <string>
#include <filesystem>
#include <utility>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <functional>
#include <memory>

namespace PNTSIM {

/** \addtogroup Application-Tools \{ */

/**
 * \class ConfigMaterial
 * \brief Class to configure the electric material properties of a model.
 * \tparam DIM The dimensions of the model.
 * \tparam CELL_NODES The number of nodes of the model's cells.
 */
template<short DIM, short CELL_NODES>
class ConfigMaterial {


private:
    std::unordered_map<std::string,PNT::ConstitutiveType> cs_type_map_;        /**< Map of constitutive law types */


public:

    /**
     * \brief ConfigMaterial object constructor.
     */
    ConfigMaterial();


    /**
     * \brief ConfigMaterial object destructor.
     */
    virtual ~ConfigMaterial();


    /**
     * \brief Set the material properties of the tissue material.
     * \param [in] parser The parser of the simulation input file.
     * \param [in] nodesets The nodesets of the tissue material geometry.
     * \param [in] units The measure units.
     * \param [out] material The tissue material.
     * \param [out] stream The message output stream.
     * \return [void]
     */
    void SetTissueMaterialProperties(const Parser &parser, const MpiHandler &mpi_handler,
        const std::unordered_map<std::string, IMP::NodeSet> &nodesets, const MeasureUnits &units,
        CardiacTissue &material, std::ostream &stream) const;


    /**
     * \brief Set the material properties of the catheter material.
     * \param [in] parser The parser of the simulation input file.
     * \param [in] nodesets The nodesets of the tiscathetersue material geometry.
     * \param [in] units The measure units.
     * \param [out] material The catheter material.
     * \param [out] stream The message output stream.
     * \return [void]
     */
    void SetCatheterMaterialProperties(const Parser &parser, const MpiHandler &mpi_handler,
        const std::unordered_map<std::string, IMP::NodeSet> &nodesets, const MeasureUnits &units,
        Catheter<DIM> &material, std::ostream &stream) const;


    /**
     * \brief Connect the catheter surface nodes with tissue nodes that are inside a radius.
     * \param [in] parser The parser of the simulation input file.
     * \param [in] tissue_nodes The nodes of the tissue material's geometry.
     * \param [in] free_tissue_nodes_flag Flag to distinguish tissue nodes that are on the free boundary.
     * \param [in] catheter_nodes The nodes of the catheter material's geometry.
     * \param [in] catheter_nodesets The nodesets of the catheter material's geometry.
     * \param [in] units The measure units.
     * \param [out] catheter The catheter material were the connecting node indices will be stored.
     * \param [out] stream The message output stream.
     * \return [void]
     */
    void ConnectCatheterTissue(const Parser &parser, const MpiHandler &mpi_handler,
        const std::vector<IMP::Vec<DIM,double>> &tissue_nodes, const std::vector<IMP::Vec<DIM,double>> &catheter_nodes,
        const std::unordered_map<std::string, IMP::NodeSet> &catheter_nodesets,
        const MeasureUnits &units, Catheter<DIM> &catheter, std::ostream &stream) const;


    /**
     * @brief Set the Hyperelastic Properties object
     * 
     * @param parser 
     * @param units 
     * @param stream 
     */
    void SetConstitutiveLaw(const Parser &parser, const MpiHandler &mpi_handler,const std::string &body_type,
        int nodes_num, const std::unordered_map<std::string, IMP::NodeSet> &nodesets, const MeasureUnits &units,
        std::shared_ptr<Constitutive> &hyperelastic, std::ostream &stream) const;


protected:

    /**
     * \brief Assign the density of a single-phase material.
     * \param [in] parser The parser of the simulation input file.
     * \param [in] body_type The type of the body represented by the material's geometry (tissue or catheter).
     * \param [in] nodesets The nodesets of the material's geometry.
     * \param [in] units The measure units.
     * \param [in] nodes_num The number of the nodes of the material.
     * \param [out] density The density values that will be assigned.
     * \return [void]
     */
    void AssignDensity(const Parser &parser, const std::string &body_type,
        const std::unordered_map<std::string, IMP::NodeSet> &nodesets, const MeasureUnits &units,
        int nodes_num, std::vector<double> &density) const;


    /**
     * \brief Assign the density of a dual-phase material.
     * \param [in] parser The parser of the simulation input file.
     * \param [in] body_type The type of the body represented by the material's geometry (tissue or catheter).
     * \param [in] nodesets The nodesets of the material's geometry.
     * \param [in] units The measure units.
     * \param [in] nodes_num The number of the nodes of the material.
     * \param [out] density_lp The liquid phase density values that will be assigned.
     * \param [out] density_gp The gas phase density values that will be assigned.
     * \return [void]
     */
    void AssignDualPhaseDensity(const Parser &parser, const std::string &body_type,
        const std::unordered_map<std::string, IMP::NodeSet> &nodesets, const MeasureUnits &units,
        int nodes_num, std::vector<double> &density_lp, std::vector<double> &density_gp) const;


    /**
     * \brief Assign the specific heat of a single-phase material.
     * \param [in] parser The parser of the simulation input file.
     * \param [in] body_type The type of the body represented by the material's geometry (tissue or catheter).
     * \param [in] nodesets The nodesets of the material's geometry.
     * \param [in] units The measure units.
     * \param [in] nodes_num The number of the nodes of the material.
     * \param [out] specific_heat The specific heat values that will be assigned.
     * \return [void]
     */
    void AssignSpecificHeat(const Parser &parser, const std::string &body_type,
        const std::unordered_map<std::string, IMP::NodeSet> &nodesets, const MeasureUnits &units,
        int nodes_num, std::vector<double> &specific_heat) const;


    /**
     * \brief Assign the specific heat of a dual-phase material.
     * \param [in] parser The parser of the simulation input file.
     * \param [in] body_type The type of the body represented by the material's geometry (tissue or catheter).
     * \param [in] nodesets The nodesets of the material's geometry.
     * \param [in] units The measure units.
     * \param [in] nodes_num The number of the nodes of the material.
     * \param [out] specific_heat_lp The liquid-phase specific heat values that will be assigned.
     * \param [out] specific_heat_gp The gas-phase specific heat values that will be assigned.
     * \return [void]
     */
    void AssignDualPhaseSpecificHeat(const Parser &parser, const std::string &body_type,
        const std::unordered_map<std::string, IMP::NodeSet> &nodesets, const MeasureUnits &units,
        int nodes_num, std::vector<double> &specific_heat_lp, std::vector<double> &specific_heat_gp) const;


    /**
     * \brief Assign the electrical conductivity of a material.
     * \param [in] parser The parser of the simulation input file.
     * \param [in] body_type The type of the body represented by the material's geometry (tissue or catheter).
     * \param [in] nodesets The nodesets of the material's geometry.
     * \param [in] units The measure units.
     * \param [in] nodes_num The number of the nodes of the material.
     * \param [out] electrical_conductivity The electrical conductivity values that will be assigned.
     * \return [void]
     */
    void AssignElectricalConductivity(const Parser &parser, const std::string &body_type,
        const std::unordered_map<std::string, IMP::NodeSet> &nodesets, const MeasureUnits &units,
        int nodes_num, std::vector<double> &electrical_conductivity) const;


    /**
     * \brief Assign the thermal conductivity of a material.
     * \param [in] parser The parser of the simulation input file.
     * \param [in] body_type The type of the body represented by the material's geometry (tissue or catheter).
     * \param [in] nodesets The nodesets of the material's geometry.
     * \param [in] units The measure units.
     * \param [in] nodes_num The number of the nodes of the material.
     * \param [out] thermal_conductivity The thermal conductivity values that will be assigned.
     * \return [void]
     */
    void AssignThermalConductivity(const Parser &parser, const std::string &body_type,
        const std::unordered_map<std::string, IMP::NodeSet> &nodesets, const MeasureUnits &units,
        int nodes_num, std::vector<double> &thermal_conductivity) const;


    /**
     * \brief Assign the latent heat of a material.
     * \param [in] parser The parser of the simulation input file.
     * \param [in] body_type The type of the body represented by the material's geometry (tissue or catheter).
     * \param [in] nodesets The nodesets of the material's geometry.
     * \param [in] units The measure units.
     * \param [in] nodes_num The number of the nodes of the material.
     * \param [out] latent_heat The latent heat values that will be assigned.
     * \return [void]
     */
    void AssignLatentHeat(const Parser &parser, const std::string &body_type,
        const std::unordered_map<std::string, IMP::NodeSet> &nodesets, const MeasureUnits &units,
        int nodes_num, std::vector<double> &latent_heat) const;


    /**
     * \brief Assign the water content of a material.
     * \param [in] parser The parser of the simulation input file.
     * \param [in] body_type The type of the body represented by the material's geometry (tissue or catheter).
     * \param [in] nodesets The nodesets of the material's geometry.
     * \param [in] units The measure units.
     * \param [in] nodes_num The number of the nodes of the material.
     * \param [out] water_content The water content values that will be assigned.
     * \return [void]
     */
    void AssignWaterContent(const Parser &parser, const std::string &body_type,
        const std::unordered_map<std::string, IMP::NodeSet> &nodesets, int nodes_num, std::vector<double> &water_content) const;


    /**
     * \brief Assign the young modulus of a material.
     * \param [in] parser The parser of the simulation input file.
     * \param [in] body_type The type of the body represented by the material's geometry (tissue or catheter).
     * \param [in] nodesets The nodesets of the material's geometry.
     * \param [in] units The measure units.
     * \param [in] nodes_num The number of the nodes of the material.
     * \param [out] young_modulus The young modulus values that will be assigned.
     * \return [void]
     */
    void AssignYoungModulus(const Parser &parser, const std::string &body_type,
        const std::unordered_map<std::string, IMP::NodeSet> &nodesets, const MeasureUnits &units,
        int nodes_num, std::vector<double> &young_modulus) const;


    /**
     * \brief Assign the poisson ratio of a material.
     * \param [in] parser The parser of the simulation input file.
     * \param [in] body_type The type of the body represented by the material's geometry (tissue or catheter).
     * \param [in] nodesets The nodesets of the material's geometry.
     * \param [in] nodes_num The number of the nodes of the material.
     * \param [out] poisson_ratio The poisson ratio values that will be assigned.
     * \return [void]
     */
    void AssignPoissonRatio(const Parser &parser, const std::string &body_type,
        const std::unordered_map<std::string, IMP::NodeSet> &nodesets, int nodes_num, std::vector<double> &poisson_ratio) const;


    /**
     * \brief Obtain the values of a field from an array of value-nodeset pairs.
     * \param [in] parser The parser of the simulation input file.
     * \param [in] attribute The attribute that for which we will obtain the values.
     * \param [in] nodesets The nodesets of the material's geometry.
     * \param [in] values_num The number of the values.
     * \param [out] values The values that will be assigned.
     * \return [void]
     */
    void ObtainValuesFromNodesets(const Parser &parser, const std::string &attribute,
        const std::unordered_map<std::string, IMP::NodeSet> &nodesets, int values_num, std::vector<double> &values) const;


};

/** \} End of Doxygen Groups*/

} //end of namespace PNTSIM

#endif //PHYNETOUCH_APPS_TOOLS_CONFIG_MATERIAL_HPP_

#include "config_material.tpp"