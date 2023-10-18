/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */

#ifndef PHYNETOUCH_APPS_TOOLS_CONFIG_OUTPUT_TPP_
#define PHYNETOUCH_APPS_TOOLS_CONFIG_OUTPUT_TPP_

#include "config_output.hpp"

namespace PNTSIM
{

////////////////////////////////////////////////
////!            P U B L I C              !////
//////////////////////////////////////////////


template<short DIM, short CELL_NODES>
ConfigOutput<DIM, CELL_NODES>::ConfigOutput() :
    states_types_(), field_names_() , field_types_(), field_steps_()
{}


template<short DIM, short CELL_NODES>
ConfigOutput<DIM, CELL_NODES>::~ConfigOutput()
{}


// template<short DIM, short CELL_NODES>
// void ConfigOutput<DIM, CELL_NODES>::OutputElectrical(const Parser &parser, const MpiHandler &mpi_handler,
//     const IMP::Mesh<DIM,CELL_NODES> &mesh_tis, const IMP::Mesh<DIM,CELL_NODES> &mesh_cath,
//     const Electrical<DIM, CELL_NODES> &electrical, std::ostream &stream) const
// {

//     // Output if ensight format was requested.
//     if (parser.HasAttribute("output.ensight.tissue")) {
//         this->OutputEnsightGeo(parser, mpi_handler, "tissue", mesh_tis, stream);
//         this->OutputEnsightScalarStates(parser, mpi_handler, "tissue", "voltage",
//             electrical.StoredVoltage(), stream);

//         auto states_types = std::vector<std::string>{"voltage"};
//         auto field_names = std::vector<std::string>{"Voltage"};
//         auto field_types = std::vector<std::string>{"scalar"};
//         auto field_steps = std::vector<std::size_t>{electrical.StoredVoltage().size()};
//         this->OutputEnsightAnimation(parser, mpi_handler, "tissue", states_types, field_names, field_types, field_steps, 1., stream);
//     }

//     if (parser.HasAttribute("output.ensight.catheter")) {
//         this->OutputEnsightGeo(parser, mpi_handler, "catheter", mesh_cath, stream);
//         this->OutputEnsightScalarStates(parser, mpi_handler, "catheter", "voltage",
//             electrical.StoredVoltage(), stream);

//         auto states_types = std::vector<std::string>{"voltage"};
//         auto field_names = std::vector<std::string>{"Voltage"};
//         auto field_types = std::vector<std::string>{"scalar"};
//         auto field_steps = std::vector<std::size_t>{electrical.StoredVoltage().size()};
//         this->OutputEnsightAnimation(parser, mpi_handler, "catheter", states_types, field_names, field_types, field_steps, 1., stream);
//     }

// //     // Output if ascii format was requested.
// //     if (parser.HasAttribute("output.ascii")) {
// //         this->OutputToAscii(parser, post_process, stream);
// //     }

// }

// template<short DIM, short CELL_NODES>
// void ConfigOutput<DIM, CELL_NODES>::OutputBioheat(const Parser &parser, const MpiHandler &mpi_handler,
//     const IMP::Mesh<DIM,CELL_NODES> &mesh_tis, const IMP::Mesh<DIM,CELL_NODES> &mesh_cath,
//     const Bioheat<DIM, CELL_NODES> &bioheat, std::ostream &stream) const
// {

//     // Output if ensight format was requested.
//     if (parser.HasAttribute("output.ensight.tissue")) {
//         this->OutputEnsightGeo(parser, mpi_handler, "tissue", mesh_tis, stream);
//         this->OutputEnsightScalarStates(parser, mpi_handler, "tissue", "temperature",
//             bioheat.StoredTemperature(), stream);

//         auto states_types = std::vector<std::string>{"temperature"};
//         auto field_names = std::vector<std::string>{"Temperature"};
//         auto field_types = std::vector<std::string>{"scalar"};
//         auto field_steps = std::vector<std::size_t>{bioheat.StoredTemperature().size()};
//         this->OutputEnsightAnimation(parser, mpi_handler, "tissue", states_types, field_names, field_types, field_steps, bioheat.Dt(), stream);
//     }

//     if (parser.HasAttribute("output.ensight.catheter")) {
//         this->OutputEnsightGeo(parser, mpi_handler, "catheter", mesh_cath, stream);
//         this->OutputEnsightScalarStates(parser, mpi_handler, "catheter", "temperature",
//             bioheat.StoredTemperature(), stream);

//         auto states_types = std::vector<std::string>{"temperature"};
//         auto field_names = std::vector<std::string>{"Temperature"};
//         auto field_types = std::vector<std::string>{"scalar"};
//         auto field_steps = std::vector<std::size_t>{bioheat.StoredTemperature().size()};
//         this->OutputEnsightAnimation(parser, mpi_handler, "catheter", states_types, field_names, field_types, field_steps, bioheat.Dt(), stream);
//     }
// }


template<short DIM, short CELL_NODES>
void ConfigOutput<DIM, CELL_NODES>::OutputDeformation(const Parser &parser, const MpiHandler &mpi_handler,
    const IMP::Mesh<DIM,CELL_NODES> &mesh_tis, const IMP::Mesh<DIM,CELL_NODES> &mesh_cath,
    const Deformation<DIM, CELL_NODES> &deformation, std::ostream &stream)
{
    this->Reset();

    this->field_steps_.emplace_back(deformation.DispHistory().size());
    this->states_types_.emplace_back("displacement");
    this->field_names_.emplace_back("Displacement");
    this->field_types_.emplace_back("vector");

    // Output if ensight format was requested.
    if (parser.HasAttribute("output.ensight.tissue")) {
        this->OutputEnsightGeo(parser, mpi_handler, "tissue", mesh_tis, stream);
        this->OutputEnsightVectorStatesDeprecated(parser, mpi_handler, "tissue", "displacement",
            deformation.DispHistory(), mesh_tis.NodesNum(), stream);

        auto time_interv = parser.GetValue<double>("physics.deformation.output interval.value");
        this->OutputEnsightAnimation(parser, mpi_handler, "tissue", time_interv, stream);
    }
    if (parser.HasAttribute("output.ensight.catheter")) {
        this->OutputEnsightGeo(parser, mpi_handler, "catheter", mesh_cath, stream);
        this->OutputEnsightVectorStatesDeprecated(parser, mpi_handler, "catheter", "displacement",
            deformation.DispHistory(), mesh_cath.NodesNum(), stream);

        auto time_interv = parser.GetValue<double>("physics.deformation.output interval.value");
        this->OutputEnsightAnimation(parser, mpi_handler, "catheter", time_interv, stream);
    }
}


template<short DIM, short CELL_NODES>
void ConfigOutput<DIM, CELL_NODES>::OutputMultiphysics(const Parser &parser, const MpiHandler &mpi_handler,
    const IMP::Mesh<DIM,CELL_NODES> &mesh_tis, const IMP::Mesh<DIM,CELL_NODES> &mesh_cath,
    const IMP::Voronoi<DIM> &voro_tis, const IMP::Voronoi<DIM> &voro_cath,
    const CardiacTissue &mat_tis, const Catheter<DIM> &mat_cath, const CLOUDEA::Fpm<DIM> &fpm_tis, const CLOUDEA::Fpm<DIM> &fpm_cath,
    const Electrical<DIM, CELL_NODES> &electrical_tis, const Electrical<DIM, CELL_NODES> &electrical_cath,
    const Bioheat<DIM, CELL_NODES> &bioheat_tis, const Bioheat<DIM, CELL_NODES> &bioheat_cath, std::ostream &stream)
{
    // Ensight output for tissue.
    if (parser.HasAttribute("tissue") && parser.HasAttribute("output.ensight.tissue")) {
        this->OutputMultiphysicsToEnsight(parser, mpi_handler, "tissue", mesh_tis, voro_tis,
            fpm_tis, mat_tis.ElectricalConductivity(), electrical_tis, bioheat_tis, stream);
    }

    // Ensight output for catheter.
    if (parser.HasAttribute("catheter") && parser.HasAttribute("output.ensight.catheter")) {
        this->OutputMultiphysicsToEnsight(parser, mpi_handler, "catheter", mesh_cath, voro_cath,
            fpm_cath, mat_cath.ElectricalConductivity(), electrical_cath, bioheat_cath, stream);
    }
}


////////////////////////////////////////////////
////!          P R O T E C T E D          !////
//////////////////////////////////////////////


template<short DIM, short CELL_NODES>
void ConfigOutput<DIM, CELL_NODES>::Reset()
{
    this->states_types_.clear();
    this->field_names_.clear();
    this->field_types_.clear();
    this->field_steps_.clear();
}


template<short DIM, short CELL_NODES>
void ConfigOutput<DIM, CELL_NODES>::OutputMultiphysicsToEnsight(const Parser &parser, const MpiHandler &mpi_handler,
    const std::string &body_type, const IMP::Mesh<DIM,CELL_NODES> &mesh, const IMP::Voronoi<DIM> &voro,
    const CLOUDEA::Fpm<DIM> &fpm, const std::vector<double> &electrical_conductivity,
    const Electrical<DIM, CELL_NODES> &electrical, const Bioheat<DIM, CELL_NODES> &bioheat, std::ostream &stream)
{
    this->Reset();

    if (parser.HasAttribute("output.ensight."+body_type+".geometry"))
        this->OutputEnsightGeo(parser, mpi_handler, body_type, mesh, stream);

    if (parser.HasAttribute("output.ensight."+body_type+".temperature"))
        this->OutputEnsightTemperature(parser, mpi_handler, body_type, bioheat.StoredTemperature(), stream);

    if (parser.HasAttribute("output.ensight."+body_type+".voltage"))
        this->OutputEnsightVoltage(parser, mpi_handler, body_type, electrical.StoredVoltage(), stream);

    if (parser.HasAttribute("output.ensight."+body_type+".heat source"))
        this->OutputEnsightHeatSource(parser, mpi_handler, body_type, mesh, voro,
            fpm, electrical_conductivity, electrical.StoredVoltage(), stream);

    if (parser.HasAttribute("output.ensight."+body_type+".electric field"))
        this->OutputEnsightElectricField(parser, mpi_handler, body_type, mesh, voro,
            fpm, electrical.StoredVoltage(), stream);

    if (parser.HasAttribute("output.ensight."+body_type+".animation")) {
        auto time_interv = bioheat.Dt() * parser.GetValue<double>("physics.bioheat.output interval");
        this->OutputEnsightAnimation(parser, mpi_handler, body_type, time_interv, stream);
    }

}


template<short DIM, short CELL_NODES>
void ConfigOutput<DIM, CELL_NODES>::OutputEnsightGeo(const Parser &parser, const MpiHandler &mpi_handler,
    const std::string &body_type, const IMP::Mesh<DIM,CELL_NODES> &mesh, std::ostream &stream)
{
    // Check body type.
    if (body_type != "tissue" && body_type != "catheter") {
        std::string error_msg = "Could not output Ensight geometry file. Body type should be: tissue or catheter.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Initialize ensight exporter.
    EnsightExporter<DIM,CELL_NODES> ens_export;

    // Save geometry to file.
    std::string geo_file = parser.GetValue<std::string>("output.ensight."+body_type+".geometry");
    std::string geo_ext = std::filesystem::path(geo_file).extension();
    if (geo_ext != ".geo")  geo_file += ".geo";
    ens_export.SaveGeo(mesh.Nodes(), mesh.Cells(), geo_file);

    if (mpi_handler.rank_id == 0)
        stream << Logger::Message("Saved Ensight geometry: " + geo_file + "\n");
}


template<short DIM, short CELL_NODES>
void ConfigOutput<DIM, CELL_NODES>::OutputEnsightTemperature(const Parser &parser, const MpiHandler &mpi_handler,
    const std::string &body_type, const std::vector<Eigen::VectorXd> &temperatures, std::ostream &stream)
{
    // Output body temperature.
    this->OutputEnsightScalarStates(parser, mpi_handler, body_type, "temperature", temperatures, stream);

    this->field_steps_.emplace_back(temperatures.size());
    this->states_types_.emplace_back("temperature");
    this->field_names_.emplace_back("Temperature");
    this->field_types_.emplace_back("scalar");
}


template<short DIM, short CELL_NODES>
void ConfigOutput<DIM, CELL_NODES>::OutputEnsightVoltage(const Parser &parser, const MpiHandler &mpi_handler,
    const std::string &body_type, const std::vector<Eigen::VectorXd> &voltages, std::ostream &stream)
{
    // Output body voltage.
    this->OutputEnsightScalarStates(parser, mpi_handler, body_type, "voltage", voltages, stream);

    this->field_steps_.emplace_back(voltages.size());
    this->states_types_.emplace_back("voltage");
    this->field_names_.emplace_back("Voltage");
    this->field_types_.emplace_back("scalar");
}


template<short DIM, short CELL_NODES>
void ConfigOutput<DIM, CELL_NODES>::OutputEnsightHeatSource(const Parser &parser, const MpiHandler &mpi_handler,
    const std::string &body_type, const IMP::Mesh<DIM,CELL_NODES> &mesh, const IMP::Voronoi<DIM> &voro,
    const CLOUDEA::Fpm<DIM> &fpm, const std::vector<double> &electrical_conductivity,
    const std::vector<Eigen::VectorXd> &voltages, std::ostream &stream)
{
    auto method = parser.GetValue<std::string>("numerical approximation.method");
    std::transform(std::begin(method), std::end(method), std::begin(method), ::tolower);

    // Output body heat source.
    auto heat_sources = std::vector<Eigen::VectorXd>(voltages.size());
    Eigen::VectorXd heat_source;

    Electrical<DIM, CELL_NODES> electrical;
    for (std::size_t id = 0; id != voltages.size(); ++id) {
        electrical.SetVoltage(voltages[id]);
        if (method == "fem") {
            electrical.ComputeElectricField(mesh);
        } else {
            electrical.ComputeElectricField(voro, fpm);
        }

        // Update the heat source.
        heat_source = Eigen::Map<const Eigen::VectorXd>(electrical_conductivity.data(), electrical_conductivity.size());
        heat_source = heat_source.cwiseProduct(electrical.ElectricField().rowwise().squaredNorm());
        heat_sources[id] = heat_source;
    }

    this->OutputEnsightScalarStates(parser, mpi_handler, body_type, "heat source", heat_sources, stream);

    this->field_steps_.emplace_back(heat_sources.size());
    this->states_types_.emplace_back("heat source");
    this->field_names_.emplace_back("HeatSource");
    this->field_types_.emplace_back("scalar");
}


template<short DIM, short CELL_NODES>
void ConfigOutput<DIM, CELL_NODES>::OutputEnsightElectricField(const Parser &parser, const MpiHandler &mpi_handler,
    const std::string &body_type, const IMP::Mesh<DIM,CELL_NODES> &mesh, const IMP::Voronoi<DIM> &voro,
    const CLOUDEA::Fpm<DIM> &fpm, const std::vector<Eigen::VectorXd> &voltages, std::ostream &stream)
{
    auto method = parser.GetValue<std::string>("numerical approximation.method");
    std::transform(std::begin(method), std::end(method), std::begin(method), ::tolower);

    // Output body electric field.
    auto electric_fields = std::vector<Eigen::MatrixXd>(voltages.size());

    Electrical<DIM, CELL_NODES> electrical;
    for (std::size_t id = 0; id != voltages.size(); ++id) {
        electrical.SetVoltage(voltages[id]);
        if (method == "fem") {
            electrical.ComputeElectricField(mesh);
        } else {
            electrical.ComputeElectricField(voro, fpm);
        }

        // Update the electric field.
        electric_fields[id] = electrical.ElectricField();
    }

    this->OutputEnsightVectorStates(parser, mpi_handler, body_type, "electric field", electric_fields, stream);

    this->field_steps_.emplace_back(electric_fields.size());
    this->states_types_.emplace_back("electric field");
    this->field_names_.emplace_back("ElectricField");
    this->field_types_.emplace_back("vector");
}


template<short DIM, short CELL_NODES>
void ConfigOutput<DIM, CELL_NODES>::OutputEnsightScalarStates(const Parser &parser, const MpiHandler &mpi_handler,
    const std::string &body_type, const std::string &states_type, const std::vector<Eigen::VectorXd> &scalars, std::ostream &stream)
{
    // Check body type.
    if (body_type != "tissue" && body_type != "catheter") {
        std::string error_msg = "Could not output scalar field states file. Body type should be: tissue or catheter.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Initialize ensight exporter.
    EnsightExporter<DIM,CELL_NODES> ens_export;

    // Get states file path and strip down the extension.
    std::string state_file = parser.GetValue<std::string>("output.ensight."+body_type+"."+states_type);
    std::filesystem::path states_path(state_file);
    if (states_path.has_extension()) {
        std::filesystem::path e = "";
        states_path.replace_extension(e);
        state_file = states_path.string();
    }

    // Number of solution states.
    std::size_t states_num = scalars.size();

    // Create wildcard string.
    std::size_t wildcard_size = std::to_string(states_num).size();
    std::string wildcard(wildcard_size, '0');

    // Export vector states.
    std::string state_id = "";
    std::size_t replace_start = 0;
    for (std::size_t i = 0; i != states_num; ++i) {
        state_id = std::to_string(i);
        replace_start = wildcard_size - state_id.size();
        wildcard.replace(replace_start,std::string::npos,state_id);
        ens_export.SaveScalarField(scalars[i], state_file+wildcard+".ens");
    }

    if (mpi_handler.rank_id == 0)
        stream << Logger::Message("Saved Ensight "+states_type+": "+state_file+".ens\n");

}


template<short DIM, short CELL_NODES>
void ConfigOutput<DIM, CELL_NODES>::OutputEnsightVectorStates(const Parser &parser, const MpiHandler &mpi_handler,
    const std::string &body_type, const std::string &states_type, const std::vector<Eigen::MatrixXd> &vectors,
    std::ostream &stream)
{
    // Check body type.
    if (body_type != "tissue" && body_type != "catheter") {
        std::string error_msg = "Could not output vector field states file. Body type should be: tissue or catheter.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Initialize ensight exporter.
    EnsightExporter<DIM,CELL_NODES> ens_export;

    // Get states file path and strip down the extension.
    std::string state_file = parser.GetValue<std::string>("output.ensight."+body_type+"."+states_type);
    std::filesystem::path states_path(state_file);
    if (states_path.has_extension()) {
        std::filesystem::path e = "";
        states_path.replace_extension(e);
        state_file = states_path.string();
    }

    // Number of solution states.
    std::size_t states_num = vectors.size();

    // Create wildcard string.
    std::size_t wildcard_size = std::to_string(states_num).size();
    std::string wildcard(wildcard_size, '0');

    // Export vector states.
    std::string state_id = "";
    std::size_t replace_start = 0;
    for (std::size_t i = 0; i != states_num; ++i) {
        state_id = std::to_string(i);
        replace_start = wildcard_size - state_id.size();
        wildcard.replace(replace_start,std::string::npos,state_id);
        ens_export.SaveVectorField(vectors[i], state_file+wildcard+".ens");
    }

    if (mpi_handler.rank_id == 0)
        stream << Logger::Message("Saved Ensight "+states_type+": "+state_file+".ens\n");

}


template<short DIM, short CELL_NODES>
void ConfigOutput<DIM, CELL_NODES>::OutputEnsightVectorStatesDeprecated(const Parser &parser, const MpiHandler &mpi_handler,
    const std::string &body_type, const std::string &states_type, const std::vector<Eigen::MatrixXd> &vectors,
    int nodes_num, std::ostream &stream)
{
    // Check body type.
    if (body_type != "tissue" && body_type != "catheter") {
        std::string error_msg = "Could not output vector field states file. Body type should be: tissue or catheter.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Initialize ensight exporter.
    EnsightExporter<DIM,CELL_NODES> ens_export;

    // Get states file path and strip down the extension.
    std::string state_file = parser.GetValue<std::string>("output.ensight."+body_type+"."+states_type);
    std::filesystem::path states_path(state_file);
    if (states_path.has_extension()) {
        std::filesystem::path e = "";
        states_path.replace_extension(e);
        state_file = states_path.string();
    }

    // Number of solution states.
    std::size_t states_num = vectors.size();

    // Create wildcard string.
    std::size_t wildcard_size = std::to_string(states_num).size();
    std::string wildcard(wildcard_size, '0');

    // Export vector states.
    std::string state_id = "";
    std::size_t replace_start = 0;
    for (std::size_t i = 0; i != states_num; ++i) {
        state_id = std::to_string(i);
        replace_start = wildcard_size - state_id.size();
        wildcard.replace(replace_start,std::string::npos,state_id);
        if (body_type == "tissue") {
            ens_export.SaveVectorField(vectors[i].topRows(nodes_num), state_file+wildcard+".ens");
        } else {
            ens_export.SaveVectorField(vectors[i].bottomRows(nodes_num), state_file+wildcard+".ens");
        }
    }

    if (mpi_handler.rank_id == 0)
        stream << Logger::Message("Saved Ensight "+states_type+": "+state_file+".ens\n");

}


template<short DIM, short CELL_NODES>
void ConfigOutput<DIM, CELL_NODES>::OutputEnsightAnimation(const Parser &parser, const MpiHandler &mpi_handler,
    const std::string &body_type, double time_inc, std::ostream &stream)
{
    std::string anim_file = parser.GetValue<std::string>("output.ensight."+body_type+".animation");
    std::string anim_path = std::filesystem::path(anim_file).parent_path();
    if (std::filesystem::path(anim_file).extension() != ".case") { anim_file += ".case"; }

    std::string geo_file = parser.GetValue<std::string>("output.ensight."+body_type+".geometry");
    std::string geo_path = std::filesystem::path(geo_file).parent_path();
    if (geo_path != anim_path) {
        std::string error_msg = "Could not output ensight animation. Geometry *.geo file must have the same path with animation *.case file.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }
    if (std::filesystem::path(geo_file).extension() != ".geo")  geo_file += ".geo";
    geo_file = std::filesystem::path(geo_file).filename();

    std::vector<std::string> state_files(this->states_types_.size());
    std::string state_file = "";
    int it = 0;
    for (const auto &states_type : this->states_types_) {
        // Get the path of the states type without extension.
        state_file = parser.GetValue<std::string>("output.ensight."+body_type+"."+states_type);
        if (std::filesystem::path(state_file).parent_path() != anim_path) {
            std::string error_msg = "Could not output ensight animation. State *.ens file must have the same path with animation *.case file.";
            throw std::invalid_argument(Logger::Error(error_msg));
        }
        state_file = std::filesystem::path(state_file).stem();

        // Store the states type file path with wildcard and extension.
        std::string digits = std::to_string(this->field_steps_[it]);
        std::string wildcard(digits.size(), '*');
        state_files[it++] = state_file+wildcard+".ens";
    }

    // Save the animation file.
    EnsightExporter<DIM,CELL_NODES> ens_export;
    ens_export.SaveAnimation(anim_file, geo_file, state_files,
        this->field_names_, this->field_types_, this->field_steps_, time_inc);
    if (mpi_handler.rank_id == 0)
        stream << Logger::Message("Saved Ensight animation: " + anim_file + "\n");
}




// template<short DIM, short CELL_NODES>
// void ConfigOutput<DIM, CELL_NODES>::OutputToAscii(const Parser &parser, const PostProcess &post_process, std::ostream &stream) const
// {
//     AsciiExporter ascii_export;

//     if (parser.HasAttribute("output.ascii.action potential")) {
//         std::string ap_file = parser.GetValue<std::string>("output.ascii.action potential");
//         std::string ap_ext = std::filesystem::path(ap_file).extension();
//         if (ap_ext != ".txt") { ap_file += ".txt"; }

//         ascii_export.WriteActionPotentials(post_process.ActionPotentials(), ap_file);
//         stream << Logger::Message("Saved action potential: "+ap_file+"\n");
//     }

// }


} // end of namespace PNTSIM

#endif //PHYNETOUCH_APPS_TOOLS_CONFIG_OUTPUT_TPP_