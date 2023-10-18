/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */

#ifndef PHYNETOUCH_APPS_TOOLS_CONFIG_MATERIAL_TPP_
#define PHYNETOUCH_APPS_TOOLS_CONFIG_MATERIAL_TPP_

#include "config_material.hpp"

namespace PNTSIM
{

///!  PUBLIC  ///

template<short DIM, short CELL_NODES>
ConfigMaterial<DIM, CELL_NODES>::ConfigMaterial() : cs_type_map_()
{
    this->cs_type_map_["neohookean"] = PNT::ConstitutiveType::neohookean;
    this->cs_type_map_["rigid"] = PNT::ConstitutiveType::rigid;
}


template<short DIM, short CELL_NODES>
ConfigMaterial<DIM, CELL_NODES>::~ConfigMaterial()
{}


template<short DIM, short CELL_NODES>
void ConfigMaterial<DIM, CELL_NODES>::SetTissueMaterialProperties(const Parser &parser, const MpiHandler &mpi_handler,
    const std::unordered_map<std::string,IMP::NodeSet> &nodesets, const MeasureUnits &units, CardiacTissue &material, std::ostream &stream) const
{

    // Assign density values to the tissue material.
    auto density_lp = std::vector<double>{};
    auto density_gp = std::vector<double>{};
    if (parser.HasAttribute("tissue.material.density.liquid phase value") &&
            parser.HasAttribute("tissue.material.density.gas phase value")) {
        if (mpi_handler.rank_id == 0)
            stream << Logger::Message("Assigned tissue material density...");
        this->AssignDualPhaseDensity(parser, "tissue", nodesets, units, material.NodesNum(), density_lp, density_gp);
        material.SetDensityLp(density_lp);
        material.SetDensityGp(density_gp);
        density_lp.clear(); density_lp.shrink_to_fit();
        density_gp.clear(); density_gp.shrink_to_fit();
        if (mpi_handler.rank_id == 0)
            stream << "OK\n";
    } else {
        if (mpi_handler.rank_id == 0)
            stream << Logger::Message("Assigned tissue material density...");
        this->AssignDensity(parser, "tissue", nodesets, units, material.NodesNum(), density_lp);
        // Set gas phase density equal to liquid phase since it was not provided.
        material.SetDensityLp(density_lp);
        material.SetDensityGp(density_lp);
        density_lp.clear(); density_lp.shrink_to_fit();
        density_gp.clear(); density_gp.shrink_to_fit();
        if (mpi_handler.rank_id == 0)
            stream << "OK\n";
    }

    // Assign specific heat values to the tissue material.
    auto specific_heat_lp = std::vector<double>{};
    auto specific_heat_gp = std::vector<double>{};
    if (parser.HasAttribute("tissue.material.specific heat")) {
        if (mpi_handler.rank_id == 0)
            stream << Logger::Message("Assigned tissue material specific heat...");
        this->AssignDualPhaseSpecificHeat(parser, "tissue", nodesets, units, material.NodesNum(), specific_heat_lp, specific_heat_gp);
        material.SetSpecificHeatLpNot(specific_heat_lp);
        material.SetSpecificHeatGp(specific_heat_gp);
        specific_heat_lp.clear(); specific_heat_lp.shrink_to_fit();
        specific_heat_gp.clear(); specific_heat_gp.shrink_to_fit();
        if (mpi_handler.rank_id == 0)
            stream << "OK\n";
    }

    // Assign electrical conductivity not to the tissue material.
    auto electrical_conductivity_not = std::vector<double>{};
    if (parser.HasAttribute("tissue.material.electrical conductivity")) {
        if (mpi_handler.rank_id == 0)
            stream << Logger::Message("Assigned tissue material electrical conductivity...");
        this->AssignElectricalConductivity(parser, "tissue", nodesets, units, material.NodesNum(), electrical_conductivity_not);
        material.SetElectricalConductivityNot(electrical_conductivity_not);
        electrical_conductivity_not.clear(); electrical_conductivity_not.shrink_to_fit();
        if (mpi_handler.rank_id == 0)
            stream << "OK\n";
    }

    // Assign thermal conductivity to the tissue material.
    auto thermal_conductivity_not = std::vector<double>{};
    if (parser.HasAttribute("tissue.material.thermal conductivity")) {
        if (mpi_handler.rank_id == 0)
            stream << Logger::Message("Assigned tissue material thermal conductivity...");
        this->AssignThermalConductivity(parser, "tissue", nodesets, units, material.NodesNum(), thermal_conductivity_not);
        material.SetThermalConductivityNot(thermal_conductivity_not);
        thermal_conductivity_not.clear(); thermal_conductivity_not.shrink_to_fit();
        if (mpi_handler.rank_id == 0)
            stream << "OK\n";
    }

    // Assign latent heat to the tissue material.
    auto latent_heat = std::vector<double>{};
    if (parser.HasAttribute("tissue.material.latent heat")) {
        if (mpi_handler.rank_id == 0)
            stream << Logger::Message("Assigned tissue material latent heat...");
        this->AssignLatentHeat(parser, "tissue", nodesets, units, material.NodesNum(), latent_heat);
        material.SetLatentHeat(latent_heat);
        latent_heat.clear(); latent_heat.shrink_to_fit();
        if (mpi_handler.rank_id == 0)
            stream << "OK\n";
    }

    // Assign water content to the tissue material.
    std::vector<double> water_content;
    if (parser.HasAttribute("tissue.material.water content")) {
        if (mpi_handler.rank_id == 0)
            stream << Logger::Message("Assigned tissue material water content...");
        this->AssignWaterContent(parser, "tissue", nodesets, material.NodesNum(), water_content);
        material.SetWaterContent(water_content);
        water_content.clear(); water_content.shrink_to_fit();
        if (mpi_handler.rank_id == 0)
            stream << "OK\n";
    }

    if (mpi_handler.rank_id == 0)
        stream << Logger::Message("Tissue material parameters assigned successfully\n");

}


template<short DIM, short CELL_NODES>
void ConfigMaterial<DIM, CELL_NODES>::SetCatheterMaterialProperties(const Parser &parser, const MpiHandler &mpi_handler,
    const std::unordered_map<std::string, IMP::NodeSet> &nodesets, const MeasureUnits &units, Catheter<DIM> &material, std::ostream &stream) const
{

    // Assign density values to the catheter material.
    auto density = std::vector<double>{};
    if (parser.HasAttribute("catheter.material.density")) {
        if (mpi_handler.rank_id == 0)
            stream << Logger::Message("Assigned catheter material density...");
        this->AssignDensity(parser, "catheter", nodesets, units, material.NodesNum(), density);
        material.SetDensity(density);
        density.clear(); density.shrink_to_fit();
        if (mpi_handler.rank_id == 0)
            stream << "OK\n";
    }

    // Assign specific heat values to the catheter material.
    auto specific_heat = std::vector<double>{};
    if (parser.HasAttribute("catheter.material.specific heat")) {
        if (mpi_handler.rank_id == 0)
            stream << Logger::Message("Assigned catheter material specific heat...");
        this->AssignSpecificHeat(parser, "catheter", nodesets, units, material.NodesNum(), specific_heat);
        material.SetSpecificHeat(specific_heat);
        specific_heat.clear(); specific_heat.shrink_to_fit();
        if (mpi_handler.rank_id == 0)
            stream << "OK\n";
    }

    // Assign electrical conductivity to the catheter material.
    auto electrical_conductivity = std::vector<double>{};
    if (parser.HasAttribute("catheter.material.electrical conductivity")) {
        if (mpi_handler.rank_id == 0)
            stream << Logger::Message("Assigned catheter material electrical conductivity...");
        this->AssignElectricalConductivity(parser, "catheter", nodesets, units, material.NodesNum(), electrical_conductivity);
        material.SetElectricalConductivity(electrical_conductivity);
        electrical_conductivity.clear(); electrical_conductivity.shrink_to_fit();
        if (mpi_handler.rank_id == 0)
            stream << "OK\n";
    }

    // Assign thermal conductivity to the catheter material.
    auto thermal_conductivity = std::vector<double>{};
    if (parser.HasAttribute("catheter.material.thermal conductivity")) {
        if (mpi_handler.rank_id == 0)
            stream << Logger::Message("Assigned catheter material thermal conductivity...");
        this->AssignThermalConductivity(parser, "catheter", nodesets, units, material.NodesNum(), thermal_conductivity);
        material.SetThermalConductivity(thermal_conductivity);
        thermal_conductivity.clear(); thermal_conductivity.shrink_to_fit();
        if (mpi_handler.rank_id == 0)
            stream << "OK\n";
    }

    if (mpi_handler.rank_id == 0)
        stream << Logger::Message("Catheter material parameters assigned successfully\n");

}


template<short DIM, short CELL_NODES>
void ConfigMaterial<DIM, CELL_NODES>::ConnectCatheterTissue(const Parser &parser, const MpiHandler &mpi_handler,
    const std::vector<IMP::Vec<DIM,double>> &tissue_nodes, const std::vector<IMP::Vec<DIM,double>> &catheter_nodes,
    const std::unordered_map<std::string, IMP::NodeSet> &catheter_nodesets,
    const MeasureUnits &units, Catheter<DIM> &catheter, std::ostream &stream) const
{
    // Get the connector radius.
    auto connect_radius = parser.GetValue<double>("catheter.connector radius.value");
    auto radius_unit = parser.GetValue<std::string>("catheter.connector radius.unit");

    // Get the maximum number of expected connectors.
    auto max_connected = parser.GetValue<int>("catheter.max connected");

    // Get the catheter surface node ids.
    auto surf_nset_name = parser.GetValue<std::string>("catheter.surface nodeset");
    if (catheter_nodesets.find(surf_nset_name) != catheter_nodesets.end()) {
        catheter.SetSurfaceNodeIds(catheter_nodesets.at(surf_nset_name).NodeIds());
    } else {
        auto error_msg = "Could not connect catheter and tissue geometries. The catheter surface nodeset [" + surf_nset_name + "] was not found.";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    if (mpi_handler.rank_id == 0)
        stream << Logger::Message("Identifying catheter-tissue connection...");
    catheter.IdentifyTissueConnection(catheter_nodes, tissue_nodes, connect_radius*units[radius_unit], max_connected);
    if (mpi_handler.rank_id == 0)  stream << "OK\n";

    auto con_num_min = std::numeric_limits<std::size_t>::max();
    auto con_num_max = std::numeric_limits<std::size_t>::min();
    for (const auto &connected : catheter.ConnectedNodeIds()) {
        if (connected.size() < con_num_min) con_num_min = connected.size();
        if (connected.size() > con_num_max) con_num_max = connected.size();
    }
    if (mpi_handler.rank_id == 0) {
        stream << Logger::Message("Catheter-tissue node connections number: ") << catheter.ConnectorNodeIds().size() << "\n";
        stream << Logger::Message("Tissue nodes per connection range: ") << con_num_min << " - " << con_num_max << "\n";
    }

    // stream << Logger::Message("DEBUG output in ConfigMaterial::ConnectCatheterTissue");
    // std::ofstream myfile;
    // myfile.open ("/home/mood/Desktop/phynetouch_sims/testing/connectorIds.txt");
    // for (const auto &cnId : catheter.ConnectorNodeIds()) {
    //     myfile << cnId << "\n";
    // }
    // myfile.close();
    // myfile.open ("/home/mood/Desktop/phynetouch_sims/testing/connectedIds.txt");
    // auto i = int{0};
    // for (const auto &connected : catheter.ConnectedNodeIds()) {
    //     for (const auto &cnId : connected) {
    //         myfile << i << " " << cnId << "\n";
    //     }
    //     i++;
    // }
    // myfile.close();

}


template<short DIM, short CELL_NODES>
void ConfigMaterial<DIM, CELL_NODES>::SetConstitutiveLaw(const Parser &parser, const MpiHandler &mpi_handler,
    const std::string &body_type, int nodes_num, const std::unordered_map<std::string, IMP::NodeSet> &nodesets,
    const MeasureUnits &units, std::shared_ptr<Constitutive> &hyperelastic, std::ostream &stream) const
{
    // Create the constitutive law.
    auto type = parser.GetValue<std::string>(body_type+".constitutive law.type");
    std::transform(std::begin(type), std::end(type), std::begin(type), ::tolower);
    hyperelastic = ConstitutiveFactory::Create(this->cs_type_map_.at(type));
    hyperelastic->SetNodesNum(nodes_num);

    // Assign density values to the hyperelastic material.
    if (mpi_handler.rank_id == 0)
        stream << Logger::Message("Assignment of hyperelastic material parameters...");
    auto density = std::vector<double>{};
    this->AssignDensity(parser, body_type, nodesets, units, hyperelastic->NodesNum(), density);
    hyperelastic->SetDensity(density);
    density.clear(); density.shrink_to_fit();

    // Assign Young modulus to the hyperelastic material.
    auto young_modulus = std::vector<double>{};
    if (type != "rigid") {
        this->AssignYoungModulus(parser, body_type+".constitutive law", nodesets, units, hyperelastic->NodesNum(), young_modulus);
    }
    hyperelastic->SetYoungModulus(young_modulus);
    young_modulus.clear(); young_modulus.shrink_to_fit();

    // Assign Poisson's ratio to the hyperelastic material.
    auto poisson_ratio = std::vector<double>{};
    if (type != "rigid") {
        this->AssignPoissonRatio(parser, body_type+".constitutive law", nodesets, hyperelastic->NodesNum(), poisson_ratio);
    }
    hyperelastic->SetPoissonRatio(poisson_ratio);
    poisson_ratio.clear(); poisson_ratio.shrink_to_fit();

    // Compute Lame parameters and wave speed.
    hyperelastic->ComputeParameters();
    hyperelastic->ComputeWaveSpeed();

    if (mpi_handler.rank_id == 0)
        stream << "OK\n";

}



///!  PROTECTED  ///


template<short DIM, short CELL_NODES>
void ConfigMaterial<DIM, CELL_NODES>::AssignDensity(const Parser &parser, const std::string &body_type,
    const std::unordered_map<std::string, IMP::NodeSet> &nodesets, const MeasureUnits &units,
    int nodes_num, std::vector<double> &density) const
{
    // Check body type.
    if (body_type != "tissue" && body_type != "catheter") {
        auto error_msg = "Could not assign density. Body type should be: tissue or catheter.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Get density of material.
    if (parser.IsSingleValue(body_type+".material.density.value")) {
        auto value = parser.GetValue<double>(body_type+".material.density.value");
        density.assign(nodes_num, value);
        density.shrink_to_fit();
    } else if (parser.IsArray(body_type+".material.density.value")) {
        auto object = parser.GetObject(body_type+".material.density.value");
        density = object.get<std::vector<double>>();

        // Check if the number of the parsed density values is the expected.
        if (static_cast<int>(density.size()) != nodes_num) {
            auto error_msg = "Could not assign density to material: " + body_type + ". The material has [" + std::to_string(nodes_num) +
                    "] and the parsed density values are [" + std::to_string(density.size()) + "]";
            throw std::runtime_error(Logger::Error(error_msg));
        }
    } else if (parser.IsMultiArray(body_type+".material.density.value")) {
        this->ObtainValuesFromNodesets(parser, body_type+".material.density.value", nodesets, nodes_num, density);

    } else {
        auto error_msg = "Could not assign density values to material: " + body_type + ". density attribute has unexpected format.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Get units of density.
    auto density_unit = parser.GetValue<std::string>(body_type+".material.density.unit");
    auto mass_unit = std::string{""};
    auto volume_unit = std::string{""};
    auto pos = density_unit.find('/');
    if (pos != std::string::npos) {
        mass_unit = density_unit.substr(0, pos);
        volume_unit = density_unit.substr(pos+1);
    } else {
        auto error_str = "Could not extract density unit for material: " + body_type + ". Expected unit format: mass / volume";
        throw std::invalid_argument(Logger::Error(error_str));
    }

    // Apply density unit corrections.
    double scale = units[mass_unit]/units[volume_unit];
    std::for_each(std::begin(density), std::end(density), [scale](double &val){ val *= scale; } );
}


template<short DIM, short CELL_NODES>
void ConfigMaterial<DIM, CELL_NODES>::AssignDualPhaseDensity(const Parser &parser, const std::string &body_type,
        const std::unordered_map<std::string, IMP::NodeSet> &nodesets, const MeasureUnits &units, int nodes_num,
        std::vector<double> &density_lp, std::vector<double> &density_gp) const
{
    // Check body type.
    if (body_type != "tissue" && body_type != "catheter") {
        auto error_msg = "Could not assign dual phase density. Body type should be: tissue or catheter.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Get density of material at liquid phase.
    if (parser.IsSingleValue(body_type+".material.density.liquid phase value")) {
        auto value = parser.GetValue<double>(body_type+".material.density.liquid phase value");
        density_lp.assign(nodes_num, value);
        density_lp.shrink_to_fit();
    } else if (parser.IsArray(body_type+".material.density.liquid phase value")) {
        auto object = parser.GetObject(body_type+".material.density.liquid phase value");
        density_lp = object.get<std::vector<double>>();

        // Check if the number of the parsed density values is the expected.
        if (static_cast<int>(density_lp.size()) != nodes_num) {
            std::string error_msg = "Could not assign liquid phase density to material: " + body_type + ". The material has [" + std::to_string(nodes_num) +
                    "] and the parsed density values are [" + std::to_string(density_lp.size()) + "]";
            throw std::runtime_error(Logger::Error(error_msg));
        }
    } else if (parser.IsMultiArray(body_type+".material.density.liquid phase value")) {
        this->ObtainValuesFromNodesets(parser, body_type+".material.density.liquid phase value", nodesets, nodes_num, density_lp);
    } else {
        std::string error_msg = "Could not assign liquid phase density values to material: " + body_type + ". density attribute has unexpected format.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }


    // Get density of material at gas phase.
    if (parser.IsSingleValue(body_type+".material.density.gas phase value")) {
        auto value = parser.GetValue<double>(body_type+".material.density.gas phase value");
        density_gp.assign(nodes_num, value);
        density_gp.shrink_to_fit();
    } else if (parser.IsArray(body_type+".material.density.gas phase value")) {
        auto object = parser.GetObject(body_type+".material.density.gas phase value");
        density_gp = object.get<std::vector<double>>();

        // Check if the number of the parsed density values is the expected.
        if (static_cast<int>(density_gp.size()) != nodes_num) {
            std::string error_msg = "Could not assign gas phase density to material: " + body_type + ". The material has [" + std::to_string(nodes_num) +
                    "] and the parsed density values are [" + std::to_string(density_gp.size()) + "]";
            throw std::runtime_error(Logger::Error(error_msg));
        }
    } else if (parser.IsMultiArray(body_type+".material.density.gas phase value")) {
        this->ObtainValuesFromNodesets(parser, body_type+".material.density.gas phase value", nodesets, nodes_num, density_gp);
    } else {
        std::string error_msg = "Could not assign gas phase density values to material: " + body_type + ". density attribute has unexpected format.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Get units of density.
    auto density_unit = parser.GetValue<std::string>(body_type+".material.density.unit");
    auto mass_unit = std::string{""};
    auto volume_unit = std::string{""};
    auto pos = density_unit.find('/');
    if (pos != std::string::npos) {
        mass_unit = density_unit.substr(0, pos);
        volume_unit = density_unit.substr(pos+1);
    } else {
        auto error_msg = "Could not extract density unit for material: " + body_type + ". Expected unit format: mass / volume";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Apply density unit corrections.
    double scale_lp = units[mass_unit]/units[volume_unit];
    std::for_each(std::begin(density_lp), std::end(density_lp), [scale_lp](double &val){ val *= scale_lp; } );

    double scale_gp = units[mass_unit]/units[volume_unit];
    std::for_each(std::begin(density_gp), std::end(density_gp), [scale_gp](double &val){ val *= scale_gp; } );

}


template<short DIM, short CELL_NODES>
void ConfigMaterial<DIM, CELL_NODES>::AssignSpecificHeat(const Parser &parser, const std::string &body_type,
        const std::unordered_map<std::string, IMP::NodeSet> &nodesets, const MeasureUnits &units, int nodes_num,
        std::vector<double> &specific_heat) const
{
    // Check body type.
    if (body_type != "tissue" && body_type != "catheter") {
        std::string error_msg = "Could not assign specific heat. Body type should be: tissue or catheter.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Get specific heat of material.
    if (parser.IsSingleValue(body_type+".material.specific heat.value")) {
        auto value = parser.GetValue<double>(body_type+".material.specific heat.value");
        specific_heat.assign(nodes_num, value);
        specific_heat.shrink_to_fit();
    } else if (parser.IsArray(body_type+".material.specific heat.value")) {
        auto object = parser.GetObject(body_type+".material.specific heat.value");
        specific_heat = object.get<std::vector<double>>();

        // Check if the number of the parsed specific heat values is the expected.
        if (static_cast<int>(specific_heat.size()) != nodes_num) {
            auto error_msg = "Could not assign specific heat to material: " + body_type + ". The material has [" + std::to_string(nodes_num) +
                    "] and the parsed specific heat values are [" + std::to_string(specific_heat.size()) + "]";
            throw std::runtime_error(Logger::Error(error_msg));
        }
    } else if (parser.IsMultiArray(body_type+".material.specific heat.value")) {
        this->ObtainValuesFromNodesets(parser, body_type+".material.specific heat.value", nodesets, nodes_num, specific_heat);
    } else {
        auto error_msg = "Could not assign specific heat values to material: " + body_type + ". specific heat attribute has unexpected format.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Get units of specific heat.
    auto specific_unit = parser.GetValue<std::string>(body_type+".material.specific heat.unit");
    auto energy_unit = std::string{""};
    auto mass_unit = std::string{""};
    auto temp_unit = std::string{""};
    auto pos1 = specific_unit.find('/');
    auto pos2 = specific_unit.find('g');
    if (pos1 != std::string::npos && pos2 != std::string::npos) {
        energy_unit = specific_unit.substr(0,pos1);
        mass_unit = specific_unit.substr(pos1+1,pos2-1);
        temp_unit = specific_unit.substr(pos2+1);
    } else {
        auto error_msg = "Could not extract density unit for material: " + body_type + ". Expected unit format: energy / (mass * temperature)";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Apply specific heat unit corrections.
    double scale = units[energy_unit]/(units[mass_unit]*units[temp_unit]);
    std::for_each(std::begin(specific_heat), std::end(specific_heat), [scale](double &val){ val *= scale; } );
}


template<short DIM, short CELL_NODES>
void ConfigMaterial<DIM, CELL_NODES>::AssignDualPhaseSpecificHeat(const Parser &parser, const std::string &body_type,
        const std::unordered_map<std::string, IMP::NodeSet> &nodesets, const MeasureUnits &units, int nodes_num,
        std::vector<double> &specific_heat_lp, std::vector<double> &specific_heat_gp) const
{
    // Check body type.
    if (body_type != "tissue" && body_type != "catheter") {
        std::string error_msg = "Could not assign dual phase specific heat. Body type should be: tissue or catheter.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Get specific heat of material at liquid phase.
    if (parser.IsSingleValue(body_type+".material.specific heat.liquid phase value")) {
        auto value = parser.GetValue<double>(body_type+".material.specific heat.liquid phase value");
        specific_heat_lp.assign(nodes_num, value);
        specific_heat_lp.shrink_to_fit();
    } else if (parser.IsArray(body_type+".material.specific heat.liquid phase value")) {
        auto object = parser.GetObject(body_type+".material.specific heat.liquid phase value");
        specific_heat_lp = object.get<std::vector<double>>();

        // Check if the number of the parsed specific heat values is the expected.
        if (static_cast<int>(specific_heat_lp.size()) != nodes_num) {
            auto error_msg = "Could not assign liquid phase specific heat to material: " + body_type + ". The material has [" + std::to_string(nodes_num) +
                    "] and the parsed specific heat values are [" + std::to_string(specific_heat_lp.size()) + "]";
            throw std::runtime_error(Logger::Error(error_msg));
        }
    } else if (parser.IsMultiArray(body_type+".material.specific heat.liquid phase value")) {
        this->ObtainValuesFromNodesets(parser, body_type+".material.specific heat.liquid phase value", nodesets, nodes_num, specific_heat_lp);
    } else {
        auto error_msg = "Could not assign liquid phase specific heat values to material: " + body_type + ". specific heat attribute has unexpected format.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }


    // Get specific heat of material at gas phase.
    if (parser.IsSingleValue(body_type+".material.specific heat.gas phase value")) {
        auto value = parser.GetValue<double>(body_type+".material.specific heat.gas phase value");
        specific_heat_gp.assign(nodes_num, value);
        specific_heat_gp.shrink_to_fit();
    } else if (parser.IsArray(body_type+".material.specific heat.gas phase value")) {
        auto object = parser.GetObject(body_type+".material.specific heat.gas phase value");
        specific_heat_gp = object.get<std::vector<double>>();

        // Check if the number of the parsed diffusivity values is the expected.
        if (static_cast<int>(specific_heat_gp.size()) != nodes_num) {
            auto error_msg = "Could not assign gas phase specific heat to material: " + body_type + ". The material has [" + std::to_string(nodes_num) +
                    "] and the parsed specific heat values are [" + std::to_string(specific_heat_gp.size()) + "]";
            throw std::runtime_error(Logger::Error(error_msg));
        }
    } else if (parser.IsMultiArray(body_type+".material.specific heat.gas phase value")) {
        this->ObtainValuesFromNodesets(parser, body_type+".material.specific heat.gas phase value", nodesets, nodes_num, specific_heat_gp);
    } else {
        auto error_msg = "Could not assign gas phase specific heat values to material: " + body_type + ". specific heat attribute has unexpected format.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Get units of specific heat.
    auto specific_unit = parser.GetValue<std::string>(body_type+".material.specific heat.unit");
    auto energy_unit = std::string{""};
    auto mass_unit = std::string{""};
    auto temp_unit = std::string{""};
    auto pos1 = specific_unit.find('/');
    auto pos2 = specific_unit.find('g');
    if (pos1 != std::string::npos && pos2 != std::string::npos) {
        energy_unit = specific_unit.substr(0,pos1);
        mass_unit = specific_unit.substr(pos1+1,pos2-1);
        temp_unit = specific_unit.substr(pos2+1);
    } else {
        auto error_msg = "Could not extract density unit for material: " + body_type + ". Expected unit format: energy / (mass * temperature)";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Apply specific heat unit corrections.
    double scale_lp = units[energy_unit]/(units[mass_unit]*units[temp_unit]);
    std::for_each(std::begin(specific_heat_lp), std::end(specific_heat_lp), [scale_lp](double &val){ val *= scale_lp; } );

    double scale_gp = units[energy_unit]/(units[mass_unit]*units[temp_unit]);
    std::for_each(std::begin(specific_heat_gp), std::end(specific_heat_gp), [scale_gp](double &val){ val *= scale_gp; } );

}


template<short DIM, short CELL_NODES>
void ConfigMaterial<DIM, CELL_NODES>::AssignElectricalConductivity(const Parser &parser, const std::string &body_type,
        const std::unordered_map<std::string, IMP::NodeSet> &nodesets, const MeasureUnits &units, int nodes_num,
        std::vector<double> &electrical_conductivity) const
{
    // Check body type.
    if (body_type != "tissue" && body_type != "catheter") {
        std::string error_msg = "Could not assign electrical conductivity. Body type should be: tissue or catheter.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Get electrical conductivity of tissue.
    if (parser.IsSingleValue(body_type+".material.electrical conductivity.value")) {
        auto value = parser.GetValue<double>(body_type+".material.electrical conductivity.value");
        electrical_conductivity.assign(nodes_num, value);
        electrical_conductivity.shrink_to_fit();
    } else if (parser.IsArray(body_type+".material.electrical conductivity.value")) {
        auto object = parser.GetObject(body_type+".material.electrical conductivity.value");
        electrical_conductivity = object.get<std::vector<double>>();

        // Check if the number of the parsed electrical conductivity values is the expected.
        if (static_cast<int>(electrical_conductivity.size()) != nodes_num) {
            auto error_msg = "Could not assign electrical conductivity to material: " + body_type + ". The material has [" + std::to_string(nodes_num) +
                    "] and the parsed electrical conductivity values are [" + std::to_string(electrical_conductivity.size()) + "]";
            throw std::runtime_error(Logger::Error(error_msg));
        }
    } else if (parser.IsMultiArray(body_type+".material.electrical conductivity.value")) {
        this->ObtainValuesFromNodesets(parser, body_type+".material.electrical conductivity.value", nodesets, nodes_num, electrical_conductivity);
    } else {
        auto error_msg = "Could not assign electrical conductivity values to material: " + body_type + ". electrical conductivity attribute has unexpected format.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Get units of electrical conductivity.
    auto conductivity_unit = parser.GetValue<std::string>(body_type+".material.electrical conductivity.unit");
    auto conductance_unit = std::string{""};
    auto length_unit = std::string{""};
    auto pos = conductivity_unit.find('/');
    if (pos != std::string::npos) {
        conductance_unit = conductivity_unit.substr(0, pos);
        length_unit = conductivity_unit.substr(pos+1);
    } else {
        auto error_msg = "Could not extract electrical conductivity unit for material: " + body_type + ". Expected unit format: conductance / length";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Apply electrical conductivity unit correction.
    double scale = units[conductance_unit]/units[length_unit];
    std::for_each(std::begin(electrical_conductivity), std::end(electrical_conductivity), [scale](double &val){ val *= scale; } );

}


template<short DIM, short CELL_NODES>
void ConfigMaterial<DIM, CELL_NODES>::AssignThermalConductivity(const Parser &parser, const std::string &body_type,
        const std::unordered_map<std::string,IMP::NodeSet> &nodesets, const MeasureUnits &units, int nodes_num,
        std::vector<double> &thermal_conductivity) const
{
    // Check body type.
    if (body_type != "tissue" && body_type != "catheter") {
        auto error_msg = "Could not assign thermal conductivity. Body type should be: tissue or catheter.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Get thermal conductivity of material.
    if (parser.IsSingleValue(body_type+".material.thermal conductivity.value")) {
        auto value = parser.GetValue<double>(body_type+".material.thermal conductivity.value");
        thermal_conductivity.assign(nodes_num, value);
        thermal_conductivity.shrink_to_fit();
    } else if (parser.IsArray(body_type+".material.thermal conductivity.value")) {
        auto object = parser.GetObject(body_type+".material.thermal conductivity.value");
        thermal_conductivity = object.get<std::vector<double>>();

        // Check if the number of the parsed thermal conductivity values is the expected.
        if (static_cast<int>(thermal_conductivity.size()) != nodes_num) {
            auto error_msg = "Could not assign thermal conductivity to material: " + body_type + ". The material has [" + std::to_string(nodes_num) +
                    "] and the parsed thermal conductivity values are [" + std::to_string(thermal_conductivity.size()) + "]";
            throw std::runtime_error(Logger::Error(error_msg));
        }
    } else if (parser.IsMultiArray(body_type+".material.thermal conductivity.value")) {
        this->ObtainValuesFromNodesets(parser, body_type+".material.thermal conductivity.value", nodesets, nodes_num, thermal_conductivity);
    } else {
        std::string error_msg = "Could not assign thermal conductivity values to material: " + body_type + ". thermal conductivity attribute has unexpected format.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Get units of thermal conductivity.
    auto conductivity_unit = parser.GetValue<std::string>(body_type+".material.thermal conductivity.unit");
    auto power_unit = std::string{""};
    auto temperature_unit = std::string{""};
    auto length_unit = std::string{""};
    auto pos1 = conductivity_unit.find('/');
    auto pos2 = conductivity_unit.find('K');
    if (pos1 != std::string::npos && pos2 != std::string::npos) {
        power_unit = conductivity_unit.substr(0,pos1);
        temperature_unit = conductivity_unit.substr(pos1+1,pos2-1);
        length_unit = conductivity_unit.substr(pos2+1);
    } else {
        auto error_msg = "Could not extract thermal conductivity unit for material: " + body_type + ". Expected unit format: power / (temperature * length)";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Apply thermal conductivity unit correction.
    double scale = units[power_unit]/(units[temperature_unit]*units[length_unit]);
    std::for_each(std::begin(thermal_conductivity), std::end(thermal_conductivity), [scale](double &val){ val *= scale; } );

}


template<short DIM, short CELL_NODES>
void ConfigMaterial<DIM, CELL_NODES>::AssignLatentHeat(const Parser &parser, const std::string &body_type,
        const std::unordered_map<std::string, IMP::NodeSet> &nodesets, const MeasureUnits &units, int nodes_num,
        std::vector<double> &latent_heat) const
{
    // Check body type.
    if (body_type != "tissue" && body_type != "catheter") {
        auto error_msg = "Could not assign latent heat. Body type should be: tissue or catheter.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Get latent heat of material.
    if (parser.IsSingleValue(body_type+".material.latent heat.value")) {
        auto value = parser.GetValue<double>(body_type+".material.latent heat.value");
        latent_heat.assign(nodes_num, value);
        latent_heat.shrink_to_fit();
    } else if (parser.IsArray(body_type+".material.latent heat.value")) {
        auto object = parser.GetObject(body_type+".material.latent heat.value");
        latent_heat = object.get<std::vector<double>>();

        // Check if the number of the parsed latent heat values is the expected.
        if (static_cast<int>(latent_heat.size()) != nodes_num) {
            auto error_msg = "Could not assign latent heat to material: " + body_type + ". The material has [" + std::to_string(nodes_num) +
                    "] and the parsed latent heat values are [" + std::to_string(latent_heat.size()) + "]";
            throw std::runtime_error(Logger::Error(error_msg));
        }
    } else if (parser.IsMultiArray(body_type+".material.latent heat.value")) {
        this->ObtainValuesFromNodesets(parser, body_type+".material.latent heat.value", nodesets, nodes_num, latent_heat);
    } else {
        std::string error_msg = "Could not assign latent heat values to material: " + body_type + ". latent heat attribute has unexpected format.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Get units of latent heat.
    auto lat_heat_unit = parser.GetValue<std::string>(body_type+".material.latent heat.unit");
    auto energy_unit = std::string{""};
    auto volume_unit = std::string{""};
    auto pos = lat_heat_unit.find('/');
    if (pos != std::string::npos) {
        energy_unit = lat_heat_unit.substr(0, pos);
        volume_unit = lat_heat_unit.substr(pos+1);
    } else {
        throw std::invalid_argument(Logger::Error("Could not extract electrical conductivity unit for material: " + body_type + ". Expected unit format: energy / volume"));
    }

    // Apply latent heat unit correction.
    double scale = units[energy_unit]/units[volume_unit];
    std::for_each(std::begin(latent_heat), std::end(latent_heat), [scale](double &val){ val *= scale; } );

}


template<short DIM, short CELL_NODES>
void ConfigMaterial<DIM, CELL_NODES>::AssignWaterContent(const Parser &parser, const std::string &body_type,
        const std::unordered_map<std::string, IMP::NodeSet> &nodesets, int nodes_num, std::vector<double> &water_content) const
{
    // Check body type.
    if (body_type != "tissue" && body_type != "catheter") {
        auto error_msg = "Could not assign water content. Body type should be: tissue or catheter.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Get water content of material.
    if (parser.IsSingleValue(body_type+".material.water content")) {
        auto value = parser.GetValue<double>(body_type+".material.water content");
        water_content.assign(nodes_num, value);
        water_content.shrink_to_fit();
    } else if (parser.IsArray(body_type+".material.water content")) {
        auto object = parser.GetObject(body_type+".material.water content");
        water_content = object.get<std::vector<double>>();

        // Check if the number of the parsed water content values is the expected.
        if (static_cast<int>(water_content.size()) != nodes_num) {
            auto error_msg = "Could not assign water content to material: " + body_type + ". The material has [" + std::to_string(nodes_num) +
                    "] and the parsed water content values are [" + std::to_string(water_content.size()) + "]";
            throw std::runtime_error(Logger::Error(error_msg));
        }
    } else if (parser.IsMultiArray(body_type+".material.water content")) {
        this->ObtainValuesFromNodesets(parser, body_type+".material.water content", nodesets, nodes_num, water_content);
    } else {
        auto error_msg = "Could not assign water content values to material: " + body_type + ". water content attribute has unexpected format.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }
}


template<short DIM, short CELL_NODES>
void ConfigMaterial<DIM, CELL_NODES>::AssignYoungModulus(const Parser &parser, const std::string &body_type,
    const std::unordered_map<std::string, IMP::NodeSet> &nodesets, const MeasureUnits &units,
    int nodes_num, std::vector<double> &young_modulus) const
{
    // Check body type.
    if (body_type != "tissue.constitutive law" && body_type != "catheter.constitutive law") {
        auto error_msg = "Could not assign Young modulus. Body type should be: tissue.constitutive law or catheter.constitutive law.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Get Young modulus of material.
    young_modulus.clear();
    if (parser.IsSingleValue(body_type + ".young modulus.value")) {
        auto value = parser.GetValue<double>(body_type + ".young modulus.value");
        young_modulus.assign(nodes_num, value);
        young_modulus.shrink_to_fit();
    } else if (parser.IsArray(body_type + ".young modulus.value")) {
        auto object = parser.GetObject(body_type + ".young modulus.value");
        young_modulus = object.get<std::vector<double>>();

        // Check if the number of the parsed Young modulus values is the expected.
        if (static_cast<int>(young_modulus.size()) != nodes_num) {
            auto error_msg = "Could not assign Young modulus to material: " + body_type + ". The material has [" + std::to_string(nodes_num) +
                "] and the parsed Young modulus values are [" + std::to_string(young_modulus.size()) + "]";
            throw std::runtime_error(Logger::Error(error_msg));
        }
    } else if (parser.IsMultiArray(body_type+".young modulus.value")) {
        this->ObtainValuesFromNodesets(parser, body_type + ".young modulus.value", nodesets, nodes_num, young_modulus);

    } else {
        auto error_msg = "Could not assign Young modulus values to material: " + body_type + ".young modulus attribute has unexpected format.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Get units of density.
    auto pressure_unit = parser.GetValue<std::string>(body_type + ".young modulus.unit");

    // Apply Young modulus unit corrections.
    double scale = units[pressure_unit];
    std::for_each(std::begin(young_modulus), std::end(young_modulus), [scale](double &val){ val *= scale; } );

}


template<short DIM, short CELL_NODES>
void ConfigMaterial<DIM, CELL_NODES>::AssignPoissonRatio(const Parser &parser, const std::string &body_type,
    const std::unordered_map<std::string, IMP::NodeSet> &nodesets, int nodes_num, std::vector<double> &poisson_ratio) const
{
    // Check body type.
    if (body_type != "tissue.constitutive law" && body_type != "catheter.constitutive law") {
        auto error_msg = "Could not assign Poisson ratio. Body type should be: tissue.constitutive law or catheter.constitutive law.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Get Poisson ratio of material.
    if (parser.IsSingleValue(body_type + ".poisson ratio")) {
        auto value = parser.GetValue<double>(body_type + ".poisson ratio");
        poisson_ratio.assign(nodes_num, value);
        poisson_ratio.shrink_to_fit();
    } else if (parser.IsArray(body_type + ".poisson ratio")) {
        auto object = parser.GetObject(body_type + ".poisson ratio");
        poisson_ratio = object.get<std::vector<double>>();

        // Check if the number of the parsed Poisson ratio values is the expected.
        if (static_cast<int>(poisson_ratio.size()) != nodes_num) {
            auto error_msg = "Could not assign Poisson ratio to material: " + body_type + ". The material has [" + std::to_string(nodes_num) +
                "] and the parsed Poisson ratio values are [" + std::to_string(poisson_ratio.size()) + "]";
            throw std::runtime_error(Logger::Error(error_msg));
        }
    } else if (parser.IsMultiArray(body_type + ".poisson ratio")) {
        this->ObtainValuesFromNodesets(parser, body_type + ".poisson ratio", nodesets, nodes_num, poisson_ratio);

    } else {
        auto error_msg = "Could not assign Poisson ratio values to material: " + body_type + ".poisson ratio attribute has unexpected format.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }
}


template<short DIM, short CELL_NODES>
void ConfigMaterial<DIM, CELL_NODES>::ObtainValuesFromNodesets(const Parser &parser, const std::string &attribute,
        const std::unordered_map<std::string, IMP::NodeSet> &nodesets, int values_num, std::vector<double> &values) const
{
    values.clear();
    values.resize(values_num, 0.);
    auto counter = int{0};
    auto value = double{0.};
    auto object = parser.GetObject(attribute);
    if (object[0][0].is_string() && object[0][1].is_string()) {
        // Case where values are given as strings together with nodeset names.
        auto entries = object.get<std::vector<std::pair<std::string,std::string>>>();

        // Assign values to the nodes of each nodeset.
        for (const auto &entry : entries) {
            value = std::stod(entry.first);
            for (const auto &id : nodesets.at(entry.second).NodeIds()) {
                values[id] = value;
                counter++;
            }
        }
    } else if (object[0][0].is_number() && object[0][1].is_string()) {
        // Case where values are given as numbers together with nodeset names.
        auto entries = object.get<std::vector<std::pair<double,std::string>>>();

        // Assign diffusivity values to the nodes of each nodeset.
        for (const auto &entry : entries) {
            value = entry.first;
            for (const auto &id : nodesets.at(entry.second).NodeIds()) {
                values[id] = value;
                counter++;
            }
        }
    }

    // Check that all nodes obtained a value.
    if (counter != values_num) {
        std::string error_str = "Could not obtain values from node sets correctly. The total number of nodes is [" + std::to_string(values_num) +
                "] while the number of the defined nodes in the node sets is [" + std::to_string(counter) + "]";
        throw std::runtime_error(Logger::Error(error_str));
    }

}


} // end of namespace PNTSIM

#endif //PHYNETOUCH_APPS_TOOLS_CONFIG_MATERIAL_TPP_