/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */

#ifndef PHYNETOUCH_APPS_TOOLS_CONFIG_CONDITIONS_TPP_
#define PHYNETOUCH_APPS_TOOLS_CONFIG_CONDITIONS_TPP_

#include "config_conditions.hpp"

namespace PNTSIM
{

///!  PUBLIC  ///

template<short DIM>
ConfigConditions<DIM>::ConfigConditions() : load_curve_types_()
{
    // Create an unordered map for all possible loading types.
    this->load_curve_types_["smooth"] = PNT::LoadCurveType::smooth;
    this->load_curve_types_["step"] = PNT::LoadCurveType::step;
}


template<short DIM>
ConfigConditions<DIM>::~ConfigConditions()
{}


template<short DIM>
void ConfigConditions<DIM>::SetConditions(const Parser &parser, const MpiHandler &mpi_handler, const std::string &body_type,
    int body_nodes_num, const std::unordered_map<std::string, IMP::NodeSet> &nodesets, const MeasureUnits &units,
    BoundConds<1> &electrical_bc, BoundConds<1> &bioheat_bc, BoundConds<DIM> &deform_bc, std::ostream &stream) const
{
    // Check body type.
    if (body_type != "tissue" && body_type != "catheter") {
        auto error_msg = "Could not assign density. Body type should be tissue or catheter.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Set the Initial boundary conditions.
    if (parser.HasAttribute(body_type+".boundary conditions.initial")) {
        this->SetInitialBc(parser, mpi_handler, body_type, body_nodes_num, units, electrical_bc, bioheat_bc, deform_bc, stream);
    }

    // Set the Dirichlet boundary conditions.
    if (parser.HasAttribute(body_type+".boundary conditions.dirichlet")) {
        this->SetDirichletBc(parser, mpi_handler, body_type, nodesets, units, electrical_bc, bioheat_bc, deform_bc, stream);
    }

    // Set the body load boundary conditions.
    if (parser.HasAttribute(body_type+".boundary conditions.body loads")) {
        this->SetBodyLoadBc(parser, mpi_handler, body_type, units, electrical_bc, bioheat_bc, deform_bc, stream);
    }

    // Set the prescribed values boundary condition.
    electrical_bc.SetPrescribed(PrescribedBc());
}



///!  PROTECTED  ///

template<short DIM>
void ConfigConditions<DIM>::SetInitialBc(const Parser &parser, const MpiHandler &mpi_handler,
    const std::string &body_type, int body_nodesnum, const MeasureUnits &units, BoundConds<1> &electrical_bc,
    BoundConds<1> &bioheat_bc, BoundConds<DIM> &deform_bc, std::ostream &stream) const
{
    auto bc_num = parser.GetValue<unsigned short>(body_type+".boundary conditions.initial.conditions number");

    // Set initial boundary conditions.
    auto condition = InitialBc{};
    Eigen::VectorXd init_values;
    for (unsigned short i = 1; i <= bc_num; ++i) {
        auto condition_path = body_type+".boundary conditions.initial.initial-" + std::to_string(i);
        auto physics_type = parser.GetValue<std::string>(condition_path+".physics");
        auto condition_unit = parser.GetValue<std::string>(condition_path+".unit");

        bool cond_value = false;
        // Set same condition value to all nodes.
        if (parser.HasAttribute(condition_path+".value")) {
            auto condition_value = parser.GetValue<double>(condition_path+".value");
                
            if (physics_type == "bioheat") {
                init_values = IMP::ALGORITHMS::InCelsius(condition_value, condition_unit) * Eigen::VectorXd::Ones(body_nodesnum);
            } else {
                init_values = (condition_value * units[condition_unit]) * Eigen::VectorXd::Ones(body_nodesnum);
            }

            condition.SetValues(init_values);
            cond_value = true;
        }
        
        bool cond_preinit = false;
        // Set preinitialized condition value for each node.
        if (parser.HasAttribute(condition_path+".preinit")) {

            auto filename = parser.GetValue<std::string>(condition_path+".preinit");

            // Check if file is in Ensight format (.ens).
            if (fs::path(filename).extension().string() != ".ens") {
                auto err_msg = "Expected file in Ensight format (*.ens)";
                throw std::invalid_argument(Logger::Error(err_msg));
            }

            // Open file.
            auto file = std::ifstream(filename, std::ios::in);
            if (!file.is_open()) {
                auto err_msg = "Could not open the initial BC preinit file.";
                throw std::runtime_error(Logger::Error(err_msg));
            }

            // Read preinit file line by line.
            auto line = std::string{""};
            auto ss = std::stringstream{};
            auto bc_val = double{0.};
            auto bc_preinit = std::vector<double>{};

            // Skip 4-line header.
            for (auto i = 0; i < 4; ++i) {
                std::getline(file, line);
            }
            // Read values.
            while (std::getline(file, line)) {
                ss.str(line);
                if (!(ss >> bc_val)) { break; }
                ss.str(std::string{});
                ss.clear();

                if (physics_type == "bioheat") {
                    bc_preinit.emplace_back(IMP::ALGORITHMS::InCelsius(bc_val, condition_unit));
                } else {
                    bc_preinit.emplace_back(bc_val*units[condition_unit]);
                }
            }
            // Close the preinit file.
            file.close();

            init_values = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(bc_preinit.data(), bc_preinit.size());
            condition.SetValues(init_values);

            cond_preinit = true;
        } 
        
        if (!cond_value && !cond_preinit) {
            std::string error_msg = "Could not set initial BC. Requires a single value or a preinitialization file.";
            throw std::invalid_argument(Logger::Error(error_msg));
        }

        // Assign initial condition to proper container according to physics type.
        if (physics_type == "bioheat") {
            bioheat_bc.SetInitial(condition);
        } else if (physics_type == "electrical") {
            electrical_bc.SetInitial(condition);
        } else if (physics_type == "deformation") {
            deform_bc.SetInitial(condition);
        } else {
            auto error_msg = "Could not set initial-" + std::to_string(i) + ". Supported physics: [bioheat | electrical | deformation]";
            throw std::invalid_argument(Logger::Error(error_msg));
        }

        // Print information about the current initial.
        if (mpi_handler.rank_id == 0) {
            stream << Logger::Message("Set up initial-") << i << "\n";
            stream << Logger::Message("       - physics: ") << physics_type << "\n";
            stream << Logger::Message("       - max value:     ") << init_values.maxCoeff() << " " << condition_unit << "\n";
        }
    }

}


template<short DIM>
void ConfigConditions<DIM>::SetDirichletBc(const Parser &parser, const MpiHandler &mpi_handler, const std::string &body_type,
    const std::unordered_map<std::string, IMP::NodeSet> &nodesets, const MeasureUnits &units,
    BoundConds<1> &electrical_bc, BoundConds<1> &bioheat_bc, BoundConds<DIM> &deform_bc, std::ostream &stream) const
{
    auto bc_num = parser.GetValue<unsigned short>(body_type + ".boundary conditions.dirichlet.conditions number");

    // Set dirichlet boundary conditions.
    auto scalar_dirichlet = DirichletBc<1>{};
    auto vector_dirichlet = DirichletBc<DIM>{};
    auto electrical_dirichlet = std::vector<DirichletBc<1>>{};
    auto bioheat_dirichlet = std::vector<DirichletBc<1>>{};
    auto deform_dirichlet = std::vector<DirichletBc<DIM>>{};
    for (unsigned short i = 1; i <= bc_num; ++i) {
        auto condition_path = body_type + ".boundary conditions.dirichlet.dirichlet-" + std::to_string(i);
        auto physics_type = parser.GetValue<std::string>(condition_path + ".physics");

        // Set condition value.
        auto value = parser.GetValue<double>(condition_path + ".value");
        auto unit = parser.GetValue<std::string>(condition_path + ".unit");
        if (physics_type == "deformation") {
            vector_dirichlet.SetValue(value * units[unit]);
        } else if (physics_type == "bioheat") {
            // scalar_dirichlet.SetValue(IMP::ALGORITHMS::InKelvin(value, unit));
            scalar_dirichlet.SetValue(IMP::ALGORITHMS::InCelsius(value, unit));
        } else if ((physics_type == "electrical")) {
            scalar_dirichlet.SetValue(value * units[unit]);
        }

        // Set loading curve.
        auto load_type = std::string{"step"}, load_start_unit = std::string{""}, load_duration_unit = std::string{""};
        auto load_start = double{0.}, load_duration = double{0.};
        if (parser.HasAttribute(condition_path + ".loading type")) {
            // Type of the loading curve.
            load_type = parser.GetValue<std::string>(condition_path + ".loading type");
            std::transform(std::begin(load_type), std::end(load_type), std::begin(load_type), ::tolower);

            // Starting time of the loading curve.
            load_start_unit = parser.GetValue<std::string>(condition_path + ".loading start.unit");
            load_start = parser.GetValue<double>(condition_path + ".loading start.value") * units[load_start_unit];

            // Duration of the loading curve.
            load_duration_unit = parser.GetValue<std::string>(condition_path + ".loading duration.unit");
            load_duration = parser.GetValue<double>(condition_path + ".loading duration.value") * units[load_duration_unit];
        }

        // Set loading curve for the boundary condition.
        if (physics_type == "deformation") {
            vector_dirichlet.SetLoadingCurve(this->load_curve_types_.at(load_type), load_start, load_duration);
        } else {
            scalar_dirichlet.SetLoadingCurve(this->load_curve_types_.at(load_type), load_start, load_duration);
        }

        // Set the condition node ids.
        auto nset_name = parser.GetValue<std::string>(condition_path + ".nodeset");
        if (nodesets.find(nset_name) != nodesets.end()) {
            if (physics_type == "deformation") {
                vector_dirichlet.SetNodeIds(nodesets.at(nset_name).NodeIds());
            } else {
                scalar_dirichlet.SetNodeIds(nodesets.at(nset_name).NodeIds());
            }
        } else {
            auto err_msg = "Could not set dirichlet-" + std::to_string(i) + ". The nodeset [" + nset_name + "] does not exist.";
            throw std::runtime_error(Logger::Error(err_msg));
        }

        // Set the condition direction.
        if (parser.HasAttribute(condition_path + ".direction")) {
            // Get the direction entry.
            auto dir = parser.GetObject(condition_path + ".direction").get<std::vector<double>>();
            if (dir.size() != static_cast<std::size_t>(DIM)) {
                auto err_msg = "Could not set direction of dirichlet-"+std::to_string(i)+
                    ". The direction should have the same number of components as the model's dimensions.";
                throw std::invalid_argument(Logger::Error(err_msg));
            }
            vector_dirichlet.SetDirection(IMP::Vec<DIM,double>(dir));
        }

        // Assign dirichlet condition to proper container according to physics type.
        if (physics_type == "deformation") {
            deform_dirichlet.emplace_back(vector_dirichlet);
        } else if (physics_type == "bioheat") {
            bioheat_dirichlet.emplace_back(scalar_dirichlet);
        } else if (physics_type == "electrical") {
            electrical_dirichlet.emplace_back(scalar_dirichlet);
        } else {
            auto error_msg = "Could not set dirichlet-" + std::to_string(i) + ". Supported physics: [bioheat | electrical | deformation]";
            throw std::invalid_argument(Logger::Error(error_msg));
        }

        // Print information about the current dirichlet.
        if (mpi_handler.rank_id == 0) {
            stream << Logger::Message("Set up dirichlet-") << i << "\n";
            stream << Logger::Message("       - physics:              ") << physics_type << "\n";
            stream << Logger::Message("       - nodeset:              ") << nset_name << "\n";
            stream << Logger::Message("       - value:                ") << value << " " << unit << "\n";
            if (parser.HasAttribute(condition_path+".loading type")) {
                stream << Logger::Message("       - loading type:     ") << load_type << "\n";
                stream << Logger::Message("       - loading start:    ") << load_start << " " << load_start_unit << "\n";
                stream << Logger::Message("       - loading duration: ") << load_duration << " " << load_duration_unit << "\n";
            }
            if (physics_type == "deformation") {
                stream << Logger::Message("       - direction:  [");
                for (short d = 0; d != DIM-1; ++d) {
                    stream <<  vector_dirichlet.Direction()[d] << ", ";
                }
                stream << vector_dirichlet.Direction()[DIM-1] << "]\n";
            }
        }

    } // End of Set dirichlet boundary conditions.

    // Set the boundary conditions containers.
    if (deform_dirichlet.size() != 0)  deform_bc.SetDirichlet(deform_dirichlet);
    if (bioheat_dirichlet.size() != 0)  bioheat_bc.SetDirichlet(bioheat_dirichlet);
    if (electrical_dirichlet.size() != 0)  electrical_bc.SetDirichlet(electrical_dirichlet);
}


template<short DIM>
void ConfigConditions<DIM>::SetBodyLoadBc(const Parser &parser, const MpiHandler &mpi_handler, const std::string &body_type,
    const MeasureUnits &units, BoundConds<1> &electrical_bc, BoundConds<1> &bioheat_bc,
    BoundConds<DIM> &deform_bc, std::ostream &stream) const
{
    auto bc_num = parser.GetValue<unsigned short>(body_type + ".boundary conditions.body loads.conditions number");

    // Set external load boundary conditions.
    auto scalar_load = BodyLoadBc<1>{};
    auto vector_load = BodyLoadBc<DIM>{};
    for (unsigned short i = 1; i <= bc_num; ++i) {
        auto condition_path = body_type + ".boundary conditions.body loads.load-" + std::to_string(i);
        auto physics_type = parser.GetValue<std::string>(condition_path + ".physics");

        // Set condition value.
        auto value = parser.GetValue<double>(condition_path + ".value");
        auto unit = parser.GetValue<std::string>(condition_path + ".unit");
        if (physics_type == "deformation") {
            if (unit == "g" || unit == "G") {
                // vector_load.SetValue(IMP::ALGORITHMS::InNewton(value, "g") * units["N"]);
                vector_load.SetValue(IMP::ALGORITHMS::InPascal(value, "g") * units["kPa"]);

            } else {
                vector_load.SetValue(value * units[unit]);
            }
        } else {
            scalar_load.SetValue(value * units[unit]);
        }

        // Set loading curve.
        auto load_type = std::string{"step"}, load_start_unit = std::string{""}, load_duration_unit = std::string{""};
        auto load_start = double{0.}, load_duration = double{0.};
        if (parser.HasAttribute(condition_path + ".loading type")) {
            // Type of the loading curve.
            load_type = parser.GetValue<std::string>(condition_path + ".loading type");
            std::transform(std::begin(load_type), std::end(load_type), std::begin(load_type), ::tolower);

            // Starting time of the loading curve.
            load_start_unit = parser.GetValue<std::string>(condition_path + ".loading start.unit");
            load_start = parser.GetValue<double>(condition_path + ".loading start.value") * units[load_start_unit];

            // Duration of the loading curve.
            load_duration_unit = parser.GetValue<std::string>(condition_path + ".loading duration.unit");
            load_duration = parser.GetValue<double>(condition_path + ".loading duration.value") * units[load_duration_unit];
        }

        // Set loading curve for the boundary condition.
        if (physics_type == "deformation") {
            vector_load.SetLoadingCurve(this->load_curve_types_.at(load_type), load_start, load_duration);
        } else {
            scalar_load.SetLoadingCurve(this->load_curve_types_.at(load_type), load_start, load_duration);
        }

        // Set the condition direction.
        if (parser.HasAttribute(condition_path + ".direction")) {
            // Get the direction entry.
            auto dir = parser.GetObject(condition_path + ".direction").get<std::vector<double>>();
            if (dir.size() != static_cast<std::size_t>(DIM)) {
                auto err_msg = "Could not set direction of load-"+std::to_string(i)+
                    ". The direction should have the same number of components as the model's dimensions.";
                throw std::invalid_argument(Logger::Error(err_msg));
            }
            vector_load.SetDirection(IMP::Vec<DIM,double>(dir));
        }

        // Set the boundary conditions containers.
        if (physics_type == "deformation") {
            deform_bc.SetBodyLoad(vector_load);
        } else if (physics_type == "bioheat") {
            bioheat_bc.SetBodyLoad(scalar_load);
        } else if (physics_type == "electrical") {
            electrical_bc.SetBodyLoad(scalar_load);
        } else {
            auto error_msg = "Could not set load-" + std::to_string(i) + ". Supported physics: [bioheat | electrical | deformation]";
            throw std::invalid_argument(Logger::Error(error_msg));
        }

        // Print information about the current external load.
        if (mpi_handler.rank_id == 0) {
            stream << Logger::Message("Body loads") << i << "\n";
            stream << Logger::Message("Set up load-") << i << "\n";
            stream << Logger::Message("       - physics:              ") << physics_type << "\n";
            stream << Logger::Message("       - value:                ") << value << " " << unit << "\n";
            if (parser.HasAttribute(condition_path+".loading type")) {
                stream << Logger::Message("       - loading type:     ") << load_type << "\n";
                stream << Logger::Message("       - loading start:    ") << load_start << " " << load_start_unit << "\n";
                stream << Logger::Message("       - loading duration: ") << load_duration << " " << load_duration_unit << "\n";
            }
            if (physics_type == "deformation") {
                stream << Logger::Message("       - direction:  [");
                for (short d = 0; d != DIM-1; ++d) {
                    stream <<  vector_load.Direction()[d] << ", ";
                }
                stream << vector_load.Direction()[DIM-1] << "]\n";
            }
        }

    } // End of Set external loads boundary conditions.

}

} // end of namespace PNTSIM

#endif //PHYNETOUCH_APPS_TOOLS_CONFIG_CONDITIONS_TPP_