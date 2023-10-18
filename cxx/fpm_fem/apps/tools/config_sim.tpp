/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */

#ifndef PHYNETOUCH_APPS_TOOLS_CONFIG_SIM_TPP_
#define PHYNETOUCH_APPS_TOOLS_CONFIG_SIM_TPP_

#include "config_sim.hpp"

namespace PNTSIM {


template<short DIM, short CELL_NODES>
ConfigSim<DIM, CELL_NODES>::ConfigSim()
{}


template<short DIM, short CELL_NODES>
ConfigSim<DIM, CELL_NODES>::~ConfigSim()
{}


template<short DIM, short CELL_NODES>
void ConfigSim<DIM, CELL_NODES>::CheckValid(const Parser &parser)
{
    auto software = parser.GetValue<std::string>("application");
    auto author = parser.GetValue<std::string>("author");
    auto email = parser.GetValue<std::string>("email");
    auto licence = parser.GetValue<std::string>("license");

    if (software.find("PhyNeTouchSim") == std::string::npos ||
        author != "Konstantinos A. Mountris" ||
        email != "konstantinos.mountris@gmail.com" ||
        licence != "all rights reserved") {
        throw std::invalid_argument(Logger::Error("Header info in configuration file is not consistent with PhyNeTouchSim."));
    }

    std::string version = software.substr(software.find("v")+1);
    if (version != PHYNETOUCH_VERSION) {
        std::string error_str = "Could not start simulation. The configuration file version ["
                                + version + "] does not match the PhyNeTouchSim version [" + PHYNETOUCH_VERSION + "]";
        throw std::runtime_error(Logger::Error(error_str));
    }
}


template<short DIM, short CELL_NODES>
void ConfigSim<DIM, CELL_NODES>::Launch(const Parser &parser, MpiHandler mpi_handler, std::ostream &stream)
{
    // Check validity of the provided configuration file.
    this->CheckValid(parser);

    // Initialize time recorder.
    auto timer = Timer{};

    // Print relevant information.
    if (mpi_handler.rank_id == 0) {
        stream << termcolor::green << termcolor::bold << "[*** SIMULATION SET UP ***]\n" << termcolor::reset;
        stream << Logger::Message("Name:  " + parser.GetValue<std::string>("simulation name") + "\n");
        stream << Logger::Message("Numerical method: " + parser.GetValue<std::string>("numerical approximation.method") + "\n");
    }

    // Set up measure units.
    auto units = MeasureUnits{};
    auto config_units = ConfigUnits{};
    config_units.SetReferenceScale(parser, mpi_handler, units, stream);

    if (mpi_handler.rank_id == 0)
        std::cout << Logger::Message("Elapsed time: ") << termcolor::cyan << termcolor::bold << timer.PrintElapsedTime() << termcolor::reset << "\n\n";


    // Set up model geometry.
    timer.Reset();
    if (mpi_handler.rank_id == 0)
        stream << termcolor::green << termcolor::bold << "[*** GEOMETRY SET UP ***]\n" << termcolor::reset;
    auto config_geo = ConfigGeo<DIM, CELL_NODES>{};

    // Set tissue geometry if required.
    auto mesh_tis = IMP::Mesh<DIM, CELL_NODES>{};
    auto voro_tis = IMP::Voronoi<DIM>{};
    auto tis_nsets = std::unordered_map<std::string,IMP::NodeSet>{};
    auto tis_nodes_num = int{0};
    if (parser.HasAttribute("tissue")) {
        if (mpi_handler.rank_id == 0)
            stream << termcolor::magenta << "** Tissue geometry:\n" << termcolor::reset;
        config_geo.SetGeometryModel(parser, mpi_handler, "tissue", mesh_tis, voro_tis, stream);
        config_geo.ExtractNodeSets(mesh_tis, voro_tis, tis_nsets, tis_nodes_num);
    }

    // Set catheter geometry if required.
    auto mesh_cath = IMP::Mesh<DIM, CELL_NODES>{};
    auto voro_cath = IMP::Voronoi<DIM>{};
    auto cath_nsets = std::unordered_map<std::string,IMP::NodeSet>{};
    auto cath_nodes_num = int{0};
    if (parser.HasAttribute("catheter")) {
        if (mpi_handler.rank_id == 0)
            stream << termcolor::magenta << "\n** Catheter geometry:\n" << termcolor::reset;
        config_geo.SetGeometryModel(parser, mpi_handler, "catheter", mesh_cath, voro_cath, stream);
        config_geo.ExtractNodeSets(mesh_cath, voro_cath, cath_nsets, cath_nodes_num);
    }
    if (mpi_handler.rank_id == 0)
        stream << Logger::Message("Elapsed time: ") << termcolor::cyan << termcolor::bold << timer.PrintElapsedTime() << termcolor::reset << "\n\n";


    // Set up materials.
    timer.Reset();
    if (mpi_handler.rank_id == 0)
        stream << termcolor::green << termcolor::bold << "[*** MATERIAL SET UP ***]\n" << termcolor::reset;
    auto config_material = ConfigMaterial<DIM, CELL_NODES>{};

    // Set tissue material.
    auto tis_mat = CardiacTissue{};
    if (parser.HasAttribute("tissue.material")) {
        if (mpi_handler.rank_id == 0)
            stream << termcolor::magenta << "** Tissue material:\n" << termcolor::reset;
        tis_mat.SetNodesNum(tis_nodes_num);
        config_material.SetTissueMaterialProperties(parser, mpi_handler, tis_nsets, units, tis_mat, stream);
    }

    // Set catheter material.
    auto cath_mat = Catheter<DIM>{};
    if (parser.HasAttribute("catheter.material")) {
        if (mpi_handler.rank_id == 0)
            stream << termcolor::magenta << "\n** Catheter material:\n" << termcolor::reset;

        cath_mat.SetNodesNum(cath_nodes_num);
        config_material.SetCatheterMaterialProperties(parser, mpi_handler, cath_nsets, units, cath_mat, stream);

        if (parser.HasAttribute("physics.electrical")) {
            config_material.ConnectCatheterTissue(parser, mpi_handler, mesh_tis.Nodes(), mesh_cath.Nodes(),
                cath_nsets, units, cath_mat, stream);
        }
    }

    // Set tissue and catheter hyperelastic materials.
    auto tis_elastic = std::shared_ptr<Constitutive>{};
    auto cath_elastic = std::shared_ptr<Constitutive>{};
    if (parser.HasAttribute("tissue.constitutive law")) {
        if (mpi_handler.rank_id == 0)
            stream << termcolor::magenta << "\n** Tissue hyperelastic material:\n" << termcolor::reset;
        config_material.SetConstitutiveLaw(parser, mpi_handler, "tissue", tis_nodes_num, tis_nsets,
            units, tis_elastic, stream);
    }
    if (parser.HasAttribute("catheter.constitutive law")) {
        if (mpi_handler.rank_id == 0)
            stream << termcolor::magenta << "\n** Catheter hyperelastic material:\n" << termcolor::reset;
        config_material.SetConstitutiveLaw(parser, mpi_handler, "catheter", cath_nodes_num, cath_nsets,
            units, cath_elastic, stream);
    }
    if (mpi_handler.rank_id == 0)
        stream << Logger::Message("Elapsed time: ") << termcolor::cyan << termcolor::bold << timer.PrintElapsedTime() << termcolor::reset << "\n\n";


    // Set up numerical approximation.
    timer.Reset();
    if (mpi_handler.rank_id == 0)
        stream << termcolor::green << termcolor::bold << "[*** NUMERICAL APPROXIMATION SET UP ***]\n" << termcolor::reset;
    auto config_approx = ConfigApproximation<DIM,CELL_NODES>{};
    auto fpm_tis = CLOUDEA::Fpm<DIM>{};
    auto fpm_cath = CLOUDEA::Fpm<DIM>{};

    // Matrices to collect FEM shape functions, derivatives,
    // and jacobians over all mesh quadrature points.
    auto fem_tis  = CLOUDEA::FemMats<DIM,CELL_NODES>{};
    auto fem_cath = CLOUDEA::FemMats<DIM,CELL_NODES>{};

    // Set numerical approximation.
    auto approx_type = parser.GetValue<std::string>("numerical approximation.method");
    std::transform(std::begin(approx_type), std::end(approx_type), std::begin(approx_type), ::tolower);
    if (approx_type == "fem") {
        if (parser.HasAttribute("tissue")) {
            config_approx.SetFemApproximation(mpi_handler, mesh_tis, fem_tis, stream);
        }
        if (parser.HasAttribute("catheter")) {
            config_approx.SetFemApproximation(mpi_handler, mesh_cath, fem_cath, stream);
        }
    } else if (approx_type == "fpm") {
        if (parser.HasAttribute("tissue")) {
            config_approx.SetFpmApproximation(parser, mpi_handler, voro_tis, fpm_tis, stream);
        }
        if (parser.HasAttribute("catheter")) {
            config_approx.SetFpmApproximation(parser, mpi_handler, voro_cath, fpm_cath, stream);
        }
    } else {
        throw std::invalid_argument(Logger::Error("Invalid simulation method in configuration file. Expected: FEM | FPM"));
    }
    if (mpi_handler.rank_id == 0)
        stream << Logger::Message("Elapsed time: ") << termcolor::cyan << termcolor::bold << timer.PrintElapsedTime() << termcolor::reset << "\n\n";


    // Set up boundary conditions.
    if (mpi_handler.rank_id == 0)
        stream << termcolor::green << termcolor::bold << "[*** BOUNDARY CONDITIONS SET UP ***]\n" << termcolor::reset;
    auto config_conds = ConfigConditions<DIM>{};

    // Set tissue boundary conditions.
    auto tis_electrical_bc = BoundConds<1>{};
    auto tis_bioheat_bc = BoundConds<1>{};
    auto tis_deform_bc = BoundConds<DIM>{};
    if (parser.HasAttribute("tissue.boundary conditions")) {
        if (mpi_handler.rank_id == 0)
            stream << termcolor::magenta << "** Tissue boundary conditions:\n" << termcolor::reset;
        config_conds.SetConditions(parser, mpi_handler, "tissue", tis_nodes_num, tis_nsets, units, tis_electrical_bc, tis_bioheat_bc, tis_deform_bc, stream);
    }

    // Set catheter boundary conditions.
    auto cath_electrical_bc = BoundConds<1>{};
    auto cath_bioheat_bc = BoundConds<1>{};
    auto cath_deform_bc = BoundConds<DIM>{};
    if (parser.HasAttribute("catheter.boundary conditions")) {
        if (mpi_handler.rank_id == 0)
            stream << termcolor::magenta << "\n** Catheter boundary conditions:\n" << termcolor::reset;
        config_conds.SetConditions(parser, mpi_handler, "catheter", cath_nodes_num, cath_nsets, units, cath_electrical_bc, cath_bioheat_bc, cath_deform_bc, stream);
    }
    if (mpi_handler.rank_id == 0)
        stream << Logger::Message("Elapsed time: ") << termcolor::cyan << termcolor::bold << timer.PrintElapsedTime() << termcolor::reset << "\n\n";

    // Set up physics.
    timer.Reset();
    if (mpi_handler.rank_id == 0)
        stream << termcolor::green << termcolor::bold << "[*** PHYSICS SOLUTION ***]\n" << termcolor::reset;
    auto config_physics = ConfigPhysics<DIM,CELL_NODES>{};
    auto electrical_tis = Electrical<DIM,CELL_NODES>{};
    auto electrical_cath = Electrical<DIM,CELL_NODES>{};
    auto bioheat_tis = Bioheat<DIM,CELL_NODES>{};
    auto bioheat_cath = Bioheat<DIM,CELL_NODES>{};
    auto deformation = Deformation<DIM,CELL_NODES>{};
    if (parser.HasAttribute("physics.electrical")) {
        auto electrical_status = parser.GetValue<std::string>("physics.electrical.status");
        std::transform(std::begin(electrical_status), std::end(electrical_status), std::begin(electrical_status), ::tolower);

        auto bioheat_status = std::string{""};
        if (parser.HasAttribute("physics.bioheat")) {
            bioheat_status = parser.GetValue<std::string>("physics.bioheat.status");
            std::transform(std::begin(bioheat_status), std::end(bioheat_status), std::begin(bioheat_status), ::tolower);
        }

        if (bioheat_status == "on" && electrical_status == "on") {
            config_physics.SolveMultiphysics(parser, mpi_handler, mesh_tis, voro_tis, mesh_cath, voro_cath,
                tis_mat, cath_mat, fpm_tis, fpm_cath, tis_electrical_bc, cath_electrical_bc,
                tis_bioheat_bc, cath_bioheat_bc, units, electrical_tis, electrical_cath, bioheat_tis, bioheat_cath, stream);
        } else if (bioheat_status != "on" && electrical_status == "on") {
            // config_physics.SolveElectrical(parser, mpi_handler, mesh_tis, voro_tis, mesh_cath, voro_cath,
            //     tis_mat, cath_mat, fpm_tis, fpm_cath, tis_electrical_bc, cath_electrical_bc,
            //     electrical_tis, electrical_cath, bioheat_tis, bioheat_cath, stream);
        } else if (bioheat_status == "on" && electrical_status != "on") {
            config_physics.SolveBioheat(parser, mpi_handler, mesh_tis, voro_tis, mesh_cath, voro_cath,
                tis_mat, cath_mat, fpm_tis, fpm_cath, tis_bioheat_bc, cath_bioheat_bc, units,
                bioheat_tis, bioheat_cath, stream);
        }
    }
    if (parser.HasAttribute("physics.deformation")) {
        auto deform_status = parser.GetValue<std::string>("physics.deformation.status");
        std::transform(std::begin(deform_status), std::end(deform_status), std::begin(deform_status), ::tolower);
        if (deform_status == "on") {
            config_physics.SolveDeformation(parser, mpi_handler, mesh_tis, voro_tis, mesh_cath, voro_cath,
                tis_elastic, cath_elastic, fem_tis, fem_cath, fpm_tis, fpm_cath, units,
                tis_deform_bc, cath_deform_bc, deformation, stream);
        }
    }
    if (mpi_handler.rank_id == 0)
        stream << Logger::Message("Elapsed time: ") << termcolor::cyan << termcolor::bold << timer.PrintElapsedTime() << termcolor::reset << "\n\n";

    // Set up output.
    timer.Reset();
    if (mpi_handler.rank_id == 0) {
        stream << termcolor::green << termcolor::bold << "[*** OUTPUT ***]\n" << termcolor::reset;
        auto config_output = ConfigOutput<DIM,CELL_NODES>{};
        if (parser.HasAttribute("physics.electrical") && parser.HasAttribute("physics.bioheat")) {
            config_output.OutputMultiphysics(parser,  mpi_handler, mesh_tis, mesh_cath,
                voro_tis, voro_cath, tis_mat, cath_mat, fpm_tis, fpm_cath,
                electrical_tis, electrical_cath, bioheat_tis, bioheat_cath, stream);
        } /*else if (parser.HasAttribute("physics.electrical")) {
            // config_output.OutputElectrical(parser, mpi_handler, mesh_tis, mesh_cath, electrical_tis, electrical_cath, stream);
        } else if (parser.HasAttribute("physics.bioheat")) {
            // config_output.OutputBioheat(parser, mpi_handler, mesh_tis, mesh_cath, bioheat_tis, bioheat_cath, stream);
        }*/
        if (parser.HasAttribute("physics.deformation")) {
            config_output.OutputDeformation(parser, mpi_handler, mesh_tis, mesh_cath, deformation, stream);
        }
        stream << Logger::Message("Elapsed time: ") << termcolor::cyan << termcolor::bold << timer.PrintElapsedTime() << termcolor::reset << "\n\n";
    }
}


} // End of namespace PNTSIM

#endif //PHYNETOUCH_APPS_TOOLS_CONFIG_SIM_TPP_