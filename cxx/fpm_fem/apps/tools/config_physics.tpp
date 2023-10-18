/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */

#ifndef PHYNETOUCH_APPS_TOOLS_CONFIG_PHYSICS_TPP_
#define PHYNETOUCH_APPS_TOOLS_CONFIG_PHYSICS_TPP_

#include "config_physics.hpp"

namespace PNTSIM
{

///!  PUBLIC  ///


template<short DIM, short CELL_NODES>
ConfigPhysics<DIM,CELL_NODES>::ConfigPhysics()
{}


template<short DIM, short CELL_NODES>
ConfigPhysics<DIM,CELL_NODES>::~ConfigPhysics()
{}


template<short DIM, short CELL_NODES>
void ConfigPhysics<DIM,CELL_NODES>::SetContactHandle(const Parser &parser, const MpiHandler &mpi_handler,
    const IMP::Mesh<DIM,CELL_NODES> &tis_mesh, const IMP::Voronoi<DIM> &tis_voro, const IMP::Mesh<DIM,CELL_NODES> &cath_mesh,
    const IMP::Voronoi<DIM> &cath_voro, ContactHandle<DIM,CELL_NODES> &contact_handle, std::ostream &stream) const
{
    if (mpi_handler.rank_id == 0)
        stream << Logger::Message("Contact deformation...ENABLED\n");
    contact_handle.Enable();

    auto method = parser.GetValue<std::string>("numerical approximation.method");
    if (!IMP::ALGORITHMS::ExistWord(method, "fem") && !IMP::ALGORITHMS::ExistWord(method, "fpm")) {
        auto err_msg = "Can not set contact handle for unknown numerical approximation method. Supported: [FEM | FPM]";
        throw std::invalid_argument(Logger::Error(err_msg));
    }

    // Get the master and slave bodies.
    auto master_body = parser.GetValue<std::string>("physics.deformation.contact.master body");
    auto slave_body = parser.GetValue<std::string>("physics.deformation.contact.slave body");

    // Get the master and slave nodesets.
    auto master_nset_name = parser.GetValue<std::string>("physics.deformation.contact.master nodeset");
    auto slave_nset_name = parser.GetValue<std::string>("physics.deformation.contact.slave nodeset");

    // Set master and slave for FEM.
    if (IMP::ALGORITHMS::ExistWord(method, "fem")) {
        // Master is the tissue and slave the catheter.
        if (IMP::ALGORITHMS::ExistWord(master_body, "tissue")) {
            contact_handle.SetMaster(tis_mesh, master_nset_name);
            contact_handle.SetSlave(cath_mesh, slave_nset_name);

            if (mpi_handler.rank_id == 0) {
                stream << Logger::Message("Contact master body: tissue\n");
                stream << Logger::Message("Contact master nodeset: "+master_nset_name+"\n");
                stream << Logger::Message("Contact slave body: catheter\n");
                stream << Logger::Message("Contact slave nodeset: "+slave_nset_name+"\n");
            }
        } else {
            // Master is the catheter and slave the tissue.
            contact_handle.SetMaster(cath_mesh, master_nset_name);
            contact_handle.SetSlave(tis_mesh, slave_nset_name);

            if (mpi_handler.rank_id == 0) {
                stream << Logger::Message("Contact master body: catheter\n");
                stream << Logger::Message("Contact master nodeset: "+master_nset_name+"\n");
                stream << Logger::Message("Contact slave body: tissue\n");
                stream << Logger::Message("Contact slave nodeset: "+slave_nset_name+"\n");
            }
        }
    } else {  // Set master and slave for FPM.
        // Master is the tissue and slave the catheter.
        if (IMP::ALGORITHMS::ExistWord(master_body, "tissue")) {
            contact_handle.SetMaster(tis_voro, master_nset_name);
            contact_handle.SetSlave(cath_voro, slave_nset_name);

            if (mpi_handler.rank_id == 0) {
                stream << Logger::Message("Contact master body: tissue\n");
                stream << Logger::Message("Contact master nodeset: "+master_nset_name+"\n");
                stream << Logger::Message("Contact slave body: catheter\n");
                stream << Logger::Message("Contact slave nodeset: "+slave_nset_name+"\n");
            }
        } else {
            // Master is the catheter and slave the tissue.
            contact_handle.SetMaster(cath_voro, master_nset_name);
            contact_handle.SetSlave(tis_voro, slave_nset_name);

            if (mpi_handler.rank_id == 0) {
                stream << Logger::Message("Contact master body: catheter\n");
                stream << Logger::Message("Contact master nodeset: "+master_nset_name+"\n");
                stream << Logger::Message("Contact slave body: tissue\n");
                stream << Logger::Message("Contact slave nodeset: "+slave_nset_name+"\n");
            }
        }
    }

}


template<short DIM, short CELL_NODES>
void ConfigPhysics<DIM,CELL_NODES>::SolveDeformation(const Parser &parser, const MpiHandler &mpi_handler,
    const IMP::Mesh<DIM,CELL_NODES> &tis_mesh, const IMP::Voronoi<DIM> &tis_voro, const IMP::Mesh<DIM,CELL_NODES> &cath_mesh,
    const IMP::Voronoi<DIM> &cath_voro, const std::shared_ptr<Constitutive> &tis_elastic,
    const std::shared_ptr<Constitutive> &cath_elastic, const CLOUDEA::FemMats<DIM,CELL_NODES> &tis_fem,
    const CLOUDEA::FemMats<DIM,CELL_NODES> &cath_fem, const CLOUDEA::Fpm<DIM> &tis_fpm, const CLOUDEA::Fpm<DIM> &cath_fpm,
    const MeasureUnits &units, BoundConds<DIM> &tis_deform_bc, BoundConds<DIM> &cath_deform_bc,
    Deformation<DIM,CELL_NODES> &deformation, std::ostream &stream) const
{
    if (mpi_handler.rank_id == 0)
        stream << Logger::Message("Physics model: Deformation\n");

    auto method = parser.GetValue<std::string>("numerical approximation.method");
    if (!IMP::ALGORITHMS::ExistWord(method,"fem") && !IMP::ALGORITHMS::ExistWord(method,"fpm")) {
        auto err_msg = "Can not solve deformation for unknown numerical approximation method. Supported: [FEM | FPM]";
        throw std::invalid_argument(Logger::Error(err_msg));
    }

    // Apply mass scaling if requested.
    bool mass_scaling_status = false;
    if (parser.HasAttribute("physics.deformation.mass scaling")) {
        if (IMP::ALGORITHMS::ExistWord(parser.GetValue<std::string>("physics.deformation.mass scaling"), "yes")) {
            mass_scaling_status = true;
        }
    }

    // Establish contact handling if requested.
    ContactHandle<DIM,CELL_NODES> contact_handle;
    if (parser.HasAttribute("physics.deformation.contact")) {
        if (IMP::ALGORITHMS::ExistWord(parser.GetValue<std::string>("physics.deformation.contact.enable"), "yes")) {
            this->SetContactHandle(parser, mpi_handler, tis_mesh, tis_voro, cath_mesh, cath_voro, contact_handle, stream);
        }
    }

    // Set simulation time.
    auto ref_time_unit = parser.GetValue<std::string>("reference units.time");
    auto sim_time_unit = parser.GetValue<std::string>("physics.deformation.simulation time.unit");
    auto sim_time = parser.GetValue<double>("physics.deformation.simulation time.value");
    deformation.SetSimulationTime(sim_time * units[sim_time_unit]);
    if (mpi_handler.rank_id == 0)
        stream << Logger::Message("Simulation time: "+std::to_string(deformation.SimulationTime())+" "+ref_time_unit+"\n");

    // Set simulation time step.
    auto sim_dt_unit = parser.GetValue<std::string>("physics.deformation.dt.unit");
    auto sim_dt = parser.GetValue<double>("physics.deformation.dt.value");

    // Assemble lumped mass vector.
    if (mpi_handler.rank_id == 0)
        stream << Logger::Message("Form lumped mass...");
    if (IMP::ALGORITHMS::ExistWord(method,"fem")) {
        deformation.FormLumpedMass(tis_mesh, cath_mesh, tis_fem, cath_fem,
            tis_elastic, cath_elastic, mass_scaling_status);
    } else {
        deformation.FormLumpedMass(tis_voro, cath_voro, tis_fpm, cath_fpm,
            tis_elastic, cath_elastic, mass_scaling_status);
    }
    if (mpi_handler.rank_id == 0)
        stream << "OK\n";

    deformation.SetDt(sim_dt * units[sim_dt_unit]);
    if (mpi_handler.rank_id == 0) {
        stream << Logger::Message("Critical time step: "+std::to_string(deformation.DtCritical())+" "+ref_time_unit+"\n");
        stream << Logger::Message("Simulation time step: "+std::to_string(deformation.Dt())+" "+ref_time_unit+"\n");
    }

    // Define the load curves of the dirichlet boundary conditions.
    for (auto &dirichlet : tis_deform_bc.Dirichlet()) {
        dirichlet.LoadingDefinition(deformation.Dt());
    }
    for (auto &dirichlet : cath_deform_bc.Dirichlet()) {
        dirichlet.LoadingDefinition(deformation.Dt());
    }

    // Define the load curves of the body load boundary conditions.
    if (tis_deform_bc.HasBodyLoad()) {
        tis_deform_bc.BodyLoad().LoadingDefinition(deformation.Dt());
    }
    if (cath_deform_bc.HasBodyLoad()) {
        cath_deform_bc.BodyLoad().LoadingDefinition(deformation.Dt());
    }

    // Solve for deformation using TLED.
    if (mpi_handler.rank_id == 0)
        stream << Logger::Message("Solving deformation problem...\n");
    auto steps_num = static_cast<int>(std::round(deformation.SimulationTime()/deformation.Dt()));
    if (steps_num < 1) {
        auto err_msg = "Not sufficient steps for deformation problem solution. Check given simulation time and time step.";
        throw std::invalid_argument(Logger::Error(err_msg));
    }

    // Set ouput interval.
    auto output_steps_num = steps_num;
    if (parser.HasAttribute("physics.deformation.output interval")) {
        auto output_time = parser.GetValue<double>("physics.deformation.output interval.value");
        auto output_unit = parser.GetValue<std::string>("physics.deformation.output interval.unit");
        output_steps_num = static_cast<int>(std::round(deformation.SimulationTime()/(output_time*units[output_unit])));
    }

    // Initialize displacements progress matrices.
    auto nnum = tis_mesh.NodesNum()+cath_mesh.NodesNum();
    Eigen::MatrixXd u = Eigen::MatrixXd::Zero(nnum,DIM);
    Eigen::MatrixXd u_old = u;
    Eigen::MatrixXd u_new = u;

    auto save_interval = steps_num / output_steps_num;
    if (save_interval*output_steps_num == steps_num) {
        deformation.InitHistory(output_steps_num, u);
    } else {
        deformation.InitHistory(output_steps_num+1, u);
    }

    // Initialize forces.
    Eigen::MatrixXd int_forces = Eigen::MatrixXd(nnum,DIM);
    int_forces.setZero();

    Eigen::MatrixXd ext_forces = Eigen::MatrixXd(nnum,DIM);
    ext_forces.setZero();

    // Compute nodal volumes for FEM non-locking elements.
    if (IMP::ALGORITHMS::ExistWord(method,"fem")) {
        deformation.ComputeCellVolumes(tis_mesh, cath_mesh, tis_fem, cath_fem);
        deformation.ComputeNodalVolumes(tis_mesh, cath_mesh);
    }

    // Perform time integration with TLED.
    auto store_cnt = int{0};
    bool disp_error_check = false;
    for (int step = 0; step != steps_num; ++step) {
        // Guard previous displacements.
        u_old = u;
        u = u_new;

        // Compute new displacements.
        if (IMP::ALGORITHMS::ExistWord(method,"fem")) {
            deformation.UpdateInternalForces(tis_mesh, cath_mesh, tis_fem, cath_fem, tis_elastic, cath_elastic, u, int_forces);
            deformation.UpdateExternalForces(tis_mesh, cath_mesh, tis_fem, cath_fem, tis_deform_bc, cath_deform_bc, step, ext_forces);
            deformation.UpdateDisplacements(tis_deform_bc, cath_deform_bc, step, tis_mesh.NodesNum(), ext_forces, int_forces, u_old, u, u_new, disp_error_check);
            if (contact_handle.IsEnabled()) {
                deformation.ApplyContact(contact_handle, cath_mesh, tis_mesh, tis_mesh.NodesNum(), 0, u_new);
            }
        } else {
            // FPM solution.
            deformation.UpdateInternalForces(tis_voro, cath_voro, tis_fpm, cath_fpm, tis_elastic, cath_elastic, u, int_forces);
            deformation.UpdateExternalForces(tis_voro, cath_voro, tis_fpm, cath_fpm, tis_deform_bc, cath_deform_bc, step, ext_forces);
            deformation.UpdateDisplacements(tis_deform_bc, cath_deform_bc, step, tis_voro.NodesNum(), ext_forces, int_forces, u_old, u, u_new, disp_error_check);
            if (contact_handle.IsEnabled()) {
                deformation.ApplyContact(contact_handle, cath_voro, tis_voro, tis_voro.NodesNum(), 0, u_new);
            }
        }

        if (output_steps_num < steps_num) {
            if ((step+1) % save_interval == 0 || (step+1) == steps_num) {
                deformation.AddInHistory(store_cnt++,u_new);
            }
        } else {
            deformation.AddInHistory(store_cnt++,u_new);
        }

        if (disp_error_check) {
            deformation.AddInHistory(store_cnt++,u_new);
            break;
        }

        if (mpi_handler.rank_id == 0)
            stream << Logger::Message("Deformation solution step: ") << step+1 << "/" << steps_num << "\r" << std::flush;
    }
    if (mpi_handler.rank_id == 0)
        stream << "\n";
}


template<short DIM, short CELL_NODES>
void ConfigPhysics<DIM,CELL_NODES>::SolveElectrical(const Parser &parser, const MpiHandler &mpi_handler,
    const IMP::Mesh<DIM,CELL_NODES> &mesh_tis, const IMP::Voronoi<DIM> &voro_tis,
    const IMP::Mesh<DIM,CELL_NODES> &mesh_cath, const IMP::Voronoi<DIM> &voro_cath,
    CardiacTissue &mat_tis, const Catheter<DIM> &mat_cath,
    const CLOUDEA::Fpm<DIM> &fpm_tis, const CLOUDEA::Fpm<DIM> &fpm_cath,
    const BoundConds<1> &tissue_bc, const BoundConds<1> &catheter_bc,
    Electrical<DIM,CELL_NODES> &electrical, std::ostream &stream) const
{
    if (mpi_handler.rank_id == 0)
        stream << Logger::Message("Physics model: Electrical\n");

    auto method = parser.GetValue<std::string>("numerical approximation.method");
    std::transform(std::begin(method), std::end(method), std::begin(method), ::tolower);

    // Number of total nodes.
    auto nnum_tis = mat_tis.NodesNum();
    auto nnum_cath = mat_cath.NodesNum();
    auto nnum = nnum_tis + nnum_cath;

    // Compute tissue conductivities. Interested for the electrical one.
    auto temp = Eigen::VectorXd{Eigen::VectorXd::Ones(nnum_tis)};
    mat_tis.ComputeTempDependent(temp);

    // Assign initial voltage value.
    auto init_volt = Eigen::VectorXd{Eigen::VectorXd::Zero(nnum)};
    init_volt.array().head(nnum_tis) += tissue_bc.Initial().Value();
    init_volt.array().tail(nnum_cath) += catheter_bc.Initial().Value();
    electrical.SetVoltage(init_volt);

    // Initialize voltage storage container.
    electrical.InitStoredVoltage(2, Eigen::VectorXd::Zero(nnum));
    electrical.StoreCurrentVoltageAt(0);

    // Assemble electrical matrix.
    if (mpi_handler.rank_id == 0) { stream << Logger::Message("Assemble electrical algebraic system..."); }
    if (method == "fem") {
        electrical.AssembleSystem(mesh_tis, mesh_cath, mat_tis, mat_cath);
    } else if (method == "fpm") {
        electrical.AssembleSystem(voro_tis, voro_cath, mat_tis, mat_cath, fpm_tis, fpm_cath);
    } else {
        auto error_msg = "Could not assemble algebraic system for electrical problem. Supported method: [FEM | FPM]";
        throw std::invalid_argument(Logger::Error(error_msg));
    }
    if (mpi_handler.rank_id == 0) { stream << "OK\n"; }

    // Compute electrical voltage.
    if (mpi_handler.rank_id == 0) { stream << Logger::Message("Computing solution..."); }
    electrical.ComputeVoltage(tissue_bc.Dirichlet(), catheter_bc.Dirichlet(), nnum_tis);
    electrical.StoreCurrentVoltageAt(1);
    if (mpi_handler.rank_id == 0) { stream << "OK\n"; }
}


template<short DIM, short CELL_NODES>
void ConfigPhysics<DIM,CELL_NODES>::SolveBioheat(const Parser &parser, const MpiHandler &mpi_handler,
    const IMP::Mesh<DIM,CELL_NODES> &mesh_tis, const IMP::Voronoi<DIM> &voro_tis,
    const IMP::Mesh<DIM,CELL_NODES> &mesh_cath, const IMP::Voronoi<DIM> &voro_cath,
    CardiacTissue &mat_tis, Catheter<DIM> &mat_cath,
    const CLOUDEA::Fpm<DIM> &fpm_tis, const CLOUDEA::Fpm<DIM> &fpm_cath,
    BoundConds<1> &bioheat_bc_tis, BoundConds<1> &bioheat_bc_cath, const MeasureUnits &units,
    Bioheat<DIM,CELL_NODES> &bioheat_tis, Bioheat<DIM,CELL_NODES> &bioheat_cath, std::ostream &stream) const
{
    if (mpi_handler.rank_id == 0)
        stream << Logger::Message("Physics model: Bioheat\n");

    // Set bioheat simulation parameters.
    auto sim_time_unit = parser.GetValue<std::string>("physics.bioheat.simulation time.unit");
    auto sim_time = parser.GetValue<double>("physics.bioheat.simulation time.value")*units[sim_time_unit];

    auto time_step_unit = parser.GetValue<std::string>("physics.bioheat.time step.unit");
    auto time_step = parser.GetValue<double>("physics.bioheat.time step.value")*units[time_step_unit];
    bioheat_tis.SetDt(time_step);
    bioheat_cath.SetDt(time_step);

    auto steps_num = static_cast<int>(std::ceil(sim_time / time_step));
    auto output_interval = parser.GetValue<int>("physics.bioheat.output interval");

    auto method = parser.GetValue<std::string>("numerical approximation.method");
    std::transform(std::begin(method), std::end(method), std::begin(method), ::tolower);

    auto temp_depend_val = parser.GetValue<std::string>("tissue.material.temperature-dependent");
    std::transform(std::begin(temp_depend_val), std::end(temp_depend_val), std::begin(temp_depend_val), ::tolower);
    bool temp_depend_status = false;
    if (temp_depend_val == "yes" || temp_depend_val == "on") { temp_depend_status = true; }

    if (mpi_handler.rank_id == 0) {
        stream << Logger::Message("Bioheat simulation time: ") << sim_time << " " << sim_time_unit << "\n";
        stream << Logger::Message("Bioheat total time steps: ") << steps_num << " steps\n";
        stream << Logger::Message("Bioheat output interval: ") << output_interval << " steps\n";
        stream << Logger::Message("Bioheat time step: ") << time_step << " " << time_step_unit << "\n";
    }

    // Maximum tissue temperature value recordings.
    auto max_tis_temp_old = double{37.0};
    auto max_tis_temp_new = double{37.0};

    // Initial tissue.
    auto nnum_tis  = int{0};
    if (parser.HasAttribute("tissue")) {

        nnum_tis = mesh_tis.NodesNum();

        // Assign initial temperature value.
        bioheat_tis.SetTemperature(bioheat_bc_tis.Initial().Values());

        // Initialize temperature storage container.
        bioheat_tis.InitStoredTemperature(1+(steps_num/output_interval), Eigen::VectorXd::Zero(nnum_tis));
        bioheat_tis.StoreCurrentTemperatureAt(0);

        // Compute properties and enthalpy of the tissue.
        mat_tis.ComputeTempDependent(bioheat_tis.Temperature());
        mat_tis.ComputeEnthalpy(bioheat_tis.Temperature());

        // Define the load curves of the dirichlet boundary conditions.
        for (auto &dirichlet : bioheat_bc_tis.Dirichlet()) {
            dirichlet.LoadingDefinition(bioheat_tis.Dt());
        }
    }

    // Initial catheter.
    auto nnum_cath = int{0};
    if (parser.HasAttribute("catheter")) {

        nnum_cath = mesh_cath.NodesNum();

        // Assign initial temperature value.
        bioheat_cath.SetTemperature(bioheat_bc_cath.Initial().Values());

        // Initialize temperature storage container.
        bioheat_cath.InitStoredTemperature(1+(steps_num/output_interval), Eigen::VectorXd::Zero(nnum_cath));
        bioheat_cath.StoreCurrentTemperatureAt(0);

        // Compute the enthalpy of the catheter material.
        mat_cath.ComputeEnthalpy();

        // Define the load curves of the dirichlet boundary conditions.
        for (auto &dirichlet : bioheat_bc_cath.Dirichlet()) {
            dirichlet.LoadingDefinition(bioheat_cath.Dt());
        }
    }

    // Initialize the internal heat source of the tissue and catheter.
    // Considered always zero for catheter for now.
    Eigen::VectorXd rf_heat_tis  = Eigen::VectorXd::Zero(nnum_tis);
    Eigen::VectorXd rf_heat_cath = Eigen::VectorXd::Zero(nnum_cath);

    // Assemble electrical and bioheat system matrices.
    IMP::Timer timer;
    timer.Reset();
    if (method == "fem") {
        if (parser.HasAttribute("catheter")) {
            bioheat_cath.AssembleSystem(mesh_cath, mat_cath.ThermalConductivity(), mat_cath.Enthalpy(), rf_heat_cath);
        }
        if (parser.HasAttribute("tissue")) {
            bioheat_tis.AssembleSystem(mesh_tis, mat_tis.ThermalConductivity(), mat_tis.Enthalpy(), rf_heat_tis);
        }
    } else if (method == "fpm") {
        if (parser.HasAttribute("catheter")) {
            bioheat_cath.AssembleSystem(voro_cath, fpm_cath, mat_cath.ThermalConductivity(), mat_cath.Enthalpy(), rf_heat_cath);
        }
        if (parser.HasAttribute("tissue")) {
            bioheat_tis.AssembleSystem(voro_tis, fpm_tis, mat_tis.ThermalConductivity(), mat_tis.Enthalpy(), rf_heat_tis);
        }
    }
    auto assemble_time = timer.PrintElapsedTime();

    // Iterate over time.
    auto steps_counter = int{0}, pos = int{1};
    auto temp_increment = double{0.};
    const auto update_thres = parser.GetValue<double>("physics.bioheat.update threshold");
    for (int step = 0; step != steps_num; ++step) {
        // Increase the steps counter.
        steps_counter++;

        // Update time dependent BC values.
        if (parser.HasAttribute("catheter")) {
            bioheat_bc_cath.UpdateDirichlet(step, bioheat_cath.Dt());
        }
        if (parser.HasAttribute("tissue")) {
            bioheat_bc_tis.UpdateDirichlet(step, bioheat_tis.Dt());
        }

        // Solve for catheter.
        if (parser.HasAttribute("catheter")) {
            bioheat_cath.ComputeTemperature(bioheat_bc_cath.Dirichlet());
        }

        // Check temperature increment for temperature dependent variables recalculation.
        if (temp_depend_status && step != 0) {
            double temp_val = 0.;
            int temp_counter = 0;
            for (int i = 0; i != nnum_tis; ++i) {
                if (bioheat_tis.Temperature().coeff(i) > 38.) {
                    temp_val += bioheat_tis.Temperature().coeff(i);
                    temp_counter++;
                }
            }
            if (temp_counter != 0) {
                max_tis_temp_new = temp_val / static_cast<double>(temp_counter);
            }

            if  (max_tis_temp_new > max_tis_temp_old) {
                temp_increment = max_tis_temp_new - max_tis_temp_old;
            } else {
                temp_increment = max_tis_temp_old - max_tis_temp_new;
            }

            if (mpi_handler.rank_id == 0) {
                stream << Logger::Message("temp increment: ") << temp_increment << "\n";
            }

            if (temp_increment > update_thres) {
                if (mpi_handler.rank_id == 0) {
                    stream << Logger::Message("Will update systems\n");
                    stream << Logger::Message("Old temp: ") << max_tis_temp_old << " new temp: " << max_tis_temp_new << "\n";
                }
                if (max_tis_temp_new > max_tis_temp_old) {
                    max_tis_temp_old = max_tis_temp_new;
                }
            }
        }

        // Update tissue material time dependent properties.
        if (temp_increment > update_thres) {
            mat_tis.ComputeTempDependent(bioheat_tis.Temperature());
            mat_tis.ComputeEnthalpy(bioheat_tis.Temperature());
        }

        // Compute temperature. Assembling in each step since fields are temperature-dependent.
        if (method == "fem") {
            if (temp_increment > update_thres) {
                if (parser.HasAttribute("tissue")) {
                    bioheat_tis.UpdateSystem(mesh_tis, mat_tis.ThermalConductivity(), mat_tis.Enthalpy(), rf_heat_tis);
                }
            } else {
                if (parser.HasAttribute("tissue")) {
                    bioheat_tis.UpdateOnlyHeatSource(mesh_tis, rf_heat_tis);
                }
            }
        } else if (method == "fpm") {
            if (temp_increment > update_thres) {
                if (parser.HasAttribute("tissue")) {
                    bioheat_tis.UpdateSystem(voro_tis, fpm_tis,
                        mat_tis.ThermalConductivity(), mat_tis.Enthalpy(), rf_heat_tis);
                }
            } else {
                if (parser.HasAttribute("tissue")) {
                    bioheat_tis.UpdateOnlyHeatSource(voro_tis, fpm_tis, rf_heat_tis);
                }
            }
        }

        // Compute temperature.
        if (parser.HasAttribute("tissue")) {
            bioheat_tis.ComputeTemperature(bioheat_bc_tis.Dirichlet());
        }

        if ((step+1) % output_interval == 0) {
            // Store voltage and temperature.
            if (parser.HasAttribute("catheter")) {
                bioheat_cath.StoreCurrentTemperatureAt(pos);
            }
            if (parser.HasAttribute("tissue")) {
                bioheat_tis.StoreCurrentTemperatureAt(pos);
            }
            pos++;
        }

        if (mpi_handler.rank_id == 0)
            stream << Logger::Message("Completed bioheat simulation steps: " + std::to_string(step+1) + "/" + std::to_string(steps_num) + "\n");

    } // End of Iterate over time.
    if (mpi_handler.rank_id == 0) { stream << std::endl; }

    // Store final state of current potential if the last step was not saved in the loop.
    if (output_interval != 0) {
        if (steps_counter % output_interval != 0) {
            if (parser.HasAttribute("catheter")) {
                bioheat_cath.StoreCurrentTemperatureAt(pos+1);
            }
            if (parser.HasAttribute("tissue")) {
                bioheat_tis.StoreCurrentTemperatureAt(pos+1);
            }
        }
    }
    else {  // Store final state if output_interval is zero.
        if (parser.HasAttribute("catheter")) {
            bioheat_cath.StoreCurrentTemperatureAt(1);
        }
        if (parser.HasAttribute("tissue")) {
            bioheat_tis.StoreCurrentTemperatureAt(1);
        }
    }

    if (mpi_handler.rank_id == 0)
        stream << Logger::Message("Assemble time: "+assemble_time+"\n");
}


template<short DIM, short CELL_NODES>
void ConfigPhysics<DIM,CELL_NODES>::SolveMultiphysics(const Parser &parser, const MpiHandler &mpi_handler,
    const IMP::Mesh<DIM,CELL_NODES> &mesh_tis, const IMP::Voronoi<DIM> &voro_tis,
    const IMP::Mesh<DIM,CELL_NODES> &mesh_cath, const IMP::Voronoi<DIM> &voro_cath,
    CardiacTissue &mat_tis, Catheter<DIM> &mat_cath,
    const CLOUDEA::Fpm<DIM> &fpm_tis, const CLOUDEA::Fpm<DIM> &fpm_cath,
    BoundConds<1> &electrical_bc_tis, BoundConds<1> &electrical_bc_cath, BoundConds<1> &bioheat_bc_tis,
    BoundConds<1> &bioheat_bc_cath, const MeasureUnits &units,
    Electrical<DIM,CELL_NODES> &electrical_tis, Electrical<DIM,CELL_NODES> &electrical_cath,
    Bioheat<DIM,CELL_NODES> &bioheat_tis, Bioheat<DIM,CELL_NODES> &bioheat_cath, std::ostream &stream) const
{
    // Check physics problems to be solved in multiphysics solution.
    if (mpi_handler.rank_id == 0)
        stream << Logger::Message("Physics models:\n");
    auto has_electrical = false, has_bioheat = false;
    if (parser.HasAttribute("physics.electrical.status")) {
        auto status = parser.GetValue<std::string>("physics.electrical.status");
        std::transform(std::begin(status), std::end(status), std::begin(status), ::tolower);
        if (status == "on") {
            has_electrical = true;
            if (mpi_handler.rank_id == 0)
                stream << Logger::Message("      - electrical\n");
        }
    }
    if (parser.HasAttribute("physics.bioheat.status")) {
        std::string status = parser.GetValue<std::string>("physics.bioheat.status");
        std::transform(std::begin(status), std::end(status), std::begin(status), ::tolower);
        if (status == "on") {
            has_bioheat = true;
            if (mpi_handler.rank_id == 0)
                stream << Logger::Message("      - bioheat\n");
        }
    }

    // Solve coupled bioheat and electrical problems.
    if (has_bioheat && has_electrical) {
        this->SolveCoupledBioheat(parser, mpi_handler, mesh_tis, voro_tis, mesh_cath, voro_cath, mat_tis, mat_cath,
            fpm_tis, fpm_cath, electrical_bc_tis, electrical_bc_cath,
            bioheat_bc_tis, bioheat_bc_cath, units, electrical_tis, electrical_cath, bioheat_tis, bioheat_cath, stream);
    } else {
        std::string error_msg = "Multiphysics problem requires at least the bioheat and electrical problems.";
        throw std::invalid_argument(error_msg);
    }

}



///!  PROTECTED  ///

template<short DIM, short CELL_NODES>
void ConfigPhysics<DIM,CELL_NODES>::SolveCoupledBioheat(const Parser &parser, const MpiHandler &mpi_handler,
    const IMP::Mesh<DIM,CELL_NODES> &mesh_tis, const IMP::Voronoi<DIM> &voro_tis,
    const IMP::Mesh<DIM,CELL_NODES> &mesh_cath, const IMP::Voronoi<DIM> &voro_cath,
    CardiacTissue &mat_tis, Catheter<DIM> &mat_cath,
    const CLOUDEA::Fpm<DIM> &fpm_tis, const CLOUDEA::Fpm<DIM> &fpm_cath,
    BoundConds<1> &electrical_bc_tis, BoundConds<1> &electrical_bc_cath, BoundConds<1> &bioheat_bc_tis,
    BoundConds<1> &bioheat_bc_cath, const MeasureUnits &units,
    Electrical<DIM,CELL_NODES> &electrical_tis, Electrical<DIM,CELL_NODES> &electrical_cath,
    Bioheat<DIM,CELL_NODES> &bioheat_tis, Bioheat<DIM,CELL_NODES> &bioheat_cath, std::ostream &stream) const
{
    // Set bioheat simulation parameters.
    auto sim_time_unit = parser.GetValue<std::string>("physics.bioheat.simulation time.unit");
    auto sim_time = parser.GetValue<double>("physics.bioheat.simulation time.value")*units[sim_time_unit];

    auto time_step_unit = parser.GetValue<std::string>("physics.bioheat.time step.unit");
    auto time_step = parser.GetValue<double>("physics.bioheat.time step.value")*units[time_step_unit];
    bioheat_tis.SetDt(time_step);
    bioheat_cath.SetDt(time_step);

    auto steps_num = static_cast<int>(std::ceil(sim_time / time_step));
    auto output_interval = parser.GetValue<int>("physics.bioheat.output interval");

    auto method = parser.GetValue<std::string>("numerical approximation.method");
    std::transform(std::begin(method), std::end(method), std::begin(method), ::tolower);

    auto temp_depend_val = parser.GetValue<std::string>("tissue.material.temperature-dependent");
    std::transform(std::begin(temp_depend_val), std::end(temp_depend_val), std::begin(temp_depend_val), ::tolower);
    bool temp_depend_status = false;
    if (temp_depend_val == "yes" || temp_depend_val == "on") { temp_depend_status = true; }

    if (mpi_handler.rank_id == 0) {
        stream << Logger::Message("Bioheat simulation time: ") << sim_time << " " << sim_time_unit << "\n";
        stream << Logger::Message("Bioheat total time steps: ") << steps_num << " steps\n";
        stream << Logger::Message("Bioheat output interval: ") << output_interval << " steps\n";
        stream << Logger::Message("Bioheat time step: ") << time_step << " " << time_step_unit << "\n";
    }

    // Maximum tissue temperature value recordings.
    auto max_tis_temp_old = double{37.0};
    auto max_tis_temp_new = double{37.0};

    // Initial tissue.
    auto nnum_tis  = int{0};
    if (parser.HasAttribute("tissue")) {

        nnum_tis = mesh_tis.NodesNum();

        electrical_tis.SetVoltage(electrical_bc_tis.Initial().Values());

        // Initialize voltage storage container.
        electrical_tis.InitStoredVoltage(1+(steps_num/output_interval), Eigen::VectorXd::Zero(nnum_tis));
        electrical_tis.StoreCurrentVoltageAt(0);

        // Assign initial temperature value.
        bioheat_tis.SetTemperature(bioheat_bc_tis.Initial().Values());

        // Initialize temperature storage container.
        bioheat_tis.InitStoredTemperature(1+(steps_num/output_interval), Eigen::VectorXd::Zero(nnum_tis));
        bioheat_tis.StoreCurrentTemperatureAt(0);

        // Compute properties and enthalpy of the tissue.
        mat_tis.ComputeTempDependent(bioheat_tis.Temperature());
        mat_tis.ComputeEnthalpy(bioheat_tis.Temperature());

        // Define the load curves of the dirichlet boundary conditions.
        for (auto &dirichlet : electrical_bc_tis.Dirichlet()) {
            dirichlet.LoadingDefinition(bioheat_tis.Dt());
        }
        for (auto &dirichlet : bioheat_bc_tis.Dirichlet()) {
            dirichlet.LoadingDefinition(bioheat_tis.Dt());
        }
    }

    // Initial catheter.
    auto nnum_cath = int{0};
    if (parser.HasAttribute("catheter")) {

        nnum_cath = mesh_cath.NodesNum();

        electrical_cath.SetVoltage(electrical_bc_cath.Initial().Values());

        // Initialize voltage storage container.
        electrical_cath.InitStoredVoltage(1+(steps_num/output_interval), Eigen::VectorXd::Zero(nnum_cath));
        electrical_cath.StoreCurrentVoltageAt(0);

        // Assign initial temperature value.
        bioheat_cath.SetTemperature(bioheat_bc_cath.Initial().Values());

        // Initialize temperature storage container.
        bioheat_cath.InitStoredTemperature(1+(steps_num/output_interval), Eigen::VectorXd::Zero(nnum_cath));
        bioheat_cath.StoreCurrentTemperatureAt(0);

        // Compute the enthalpy of the catheter material.
        mat_cath.ComputeEnthalpy();

                // Define the load curves of the dirichlet boundary conditions.
        for (auto &dirichlet : electrical_bc_cath.Dirichlet()) {
            dirichlet.LoadingDefinition(bioheat_cath.Dt());
        }
        for (auto &dirichlet : bioheat_bc_cath.Dirichlet()) {
            dirichlet.LoadingDefinition(bioheat_cath.Dt());
        }
    }

    // Initialize the internal heat source of the tissue and catheter.
    // Considered always zero for catheter for now.
    Eigen::VectorXd rf_heat_tis  = Eigen::VectorXd::Zero(nnum_tis);
    Eigen::VectorXd rf_heat_cath = Eigen::VectorXd::Zero(nnum_cath);

    // Assemble electrical and bioheat system matrices.
    IMP::Timer timer;
    timer.Reset();
    if (method == "fem") {
        if (parser.HasAttribute("catheter")) {
            electrical_cath.AssembleSystem(mesh_cath, mat_cath.ElectricalConductivity());
            bioheat_cath.AssembleSystem(mesh_cath, mat_cath.ThermalConductivity(), mat_cath.Enthalpy(), rf_heat_cath);
        }
        if (parser.HasAttribute("tissue")) {
            electrical_tis.AssembleSystem(mesh_tis, mat_tis.ElectricalConductivity());
            bioheat_tis.AssembleSystem(mesh_tis, mat_tis.ThermalConductivity(), mat_tis.Enthalpy(), rf_heat_tis);
        }
    } else if (method == "fpm") {
        if (parser.HasAttribute("catheter")) {
            electrical_cath.AssembleSystem(voro_cath, mat_cath.ElectricalConductivity(), fpm_cath);
            bioheat_cath.AssembleSystem(voro_cath, fpm_cath, mat_cath.ThermalConductivity(), mat_cath.Enthalpy(), rf_heat_cath);
        }
        if (parser.HasAttribute("tissue")) {
            electrical_tis.AssembleSystem(voro_tis, mat_tis.ElectricalConductivity(), fpm_tis);
            bioheat_tis.AssembleSystem(voro_tis, fpm_tis, mat_tis.ThermalConductivity(), mat_tis.Enthalpy(), rf_heat_tis);
        }
    }
    auto assemble_time = timer.PrintElapsedTime();

    // Iterate over time.
    auto steps_counter = int{0}, pos = int{1};
    auto temp_increment = double{0.};
    const auto update_thres = parser.GetValue<double>("physics.bioheat.update threshold");
    for (int step = 0; step != steps_num; ++step) {
        // Increase the steps counter.
        steps_counter++;

        // Update time dependent BC values.
        if (parser.HasAttribute("catheter")) {
            electrical_bc_cath.UpdateDirichlet(step, bioheat_cath.Dt());
            bioheat_bc_cath.UpdateDirichlet(step, bioheat_cath.Dt());
        }
        if (parser.HasAttribute("tissue")) {
            bioheat_bc_tis.UpdateDirichlet(step, bioheat_tis.Dt());
            electrical_bc_tis.UpdateDirichlet(step, bioheat_tis.Dt());
        }

        // Solve for catheter.
        if (parser.HasAttribute("catheter")) {
            electrical_cath.ComputeVoltage(electrical_bc_cath);
            electrical_cath.ComputeElectricField(mesh_cath);
            bioheat_cath.ComputeTemperature(bioheat_bc_cath.Dirichlet());

            // If solving also for tissue make electrical connection between catheter and tissue.
            if (parser.HasAttribute("tissue")) {
                electrical_bc_tis.Prescribed().ConnectBodies(mat_cath.ConnectorNodeIds(),
                    mat_cath.ConnectedNodeIds(), electrical_cath.Voltage(), nnum_tis);
            }
        }

        // Check temperature increment for temperature dependent variables recalculation.
        if (temp_depend_status && step != 0) {
            double temp_val = 0.;
            int temp_counter = 0;
            for (int i = 0; i != nnum_tis; ++i) {
                if (bioheat_tis.Temperature().coeff(i) > 38.) {
                    temp_val += bioheat_tis.Temperature().coeff(i);
                    temp_counter++;
                }
            }
            if (temp_counter != 0) {
                max_tis_temp_new = temp_val / static_cast<double>(temp_counter);
            }

            if  (max_tis_temp_new > max_tis_temp_old) {
                temp_increment = max_tis_temp_new - max_tis_temp_old;
            } else {
                temp_increment = max_tis_temp_old - max_tis_temp_new;
            }

            if (mpi_handler.rank_id == 0) {
                stream << Logger::Message("temp increment: ") << temp_increment << "\n";
            }

            if (temp_increment > update_thres) {
                if (mpi_handler.rank_id == 0) {
                    stream << Logger::Message("Will update systems\n");
                    stream << Logger::Message("Old temp: ") << max_tis_temp_old << " new temp: " << max_tis_temp_new << "\n";
                }
                if (max_tis_temp_new > max_tis_temp_old) {
                    max_tis_temp_old = max_tis_temp_new;
                }
            }
        }

        // Update tissue material time dependent properties.
        if (temp_increment > update_thres) {
            mat_tis.ComputeTempDependent(bioheat_tis.Temperature());
            mat_tis.ComputeEnthalpy(bioheat_tis.Temperature());
        }

        // Compute potential and electric field. Assembling in each step since fields are temperature-dependent.
        if (method == "fem") {
            if (temp_increment > update_thres) {
                if (parser.HasAttribute("tissue")) {
                    electrical_tis.UpdateSystem(mesh_tis, mat_tis.ElectricalConductivity());
                }
            }
            
            if (parser.HasAttribute("tissue")) {
                electrical_tis.ComputeVoltage(electrical_bc_tis);
                electrical_tis.ComputeElectricField(mesh_tis);
            }
            
        } else if (method == "fpm") {
            if (temp_increment > update_thres) {
                if (parser.HasAttribute("tissue")) {
                    electrical_tis.UpdateSystem(voro_tis, mat_tis.ElectricalConductivity(), fpm_tis);
                }
            }
            
            if (parser.HasAttribute("tissue")) {
                electrical_tis.ComputeVoltage(electrical_bc_tis);
                electrical_tis.ComputeElectricField(voro_tis, fpm_tis);
            }
            
        }

        // Update the tissue heat source.
        if (parser.HasAttribute("tissue")) {
            rf_heat_tis = Eigen::Map<const Eigen::VectorXd>(mat_tis.ElectricalConductivity().data(),
                mat_tis.ElectricalConductivity().size());
            rf_heat_tis = rf_heat_tis.cwiseProduct(electrical_tis.ElectricField().rowwise().squaredNorm());

            if (mpi_handler.rank_id == 0) {
                stream << Logger::Message("sigma: ") << mat_tis.ElectricalConductivity().data()[0] << std::endl;
                stream << Logger::Message("rho: ") << mat_tis.DensityLp().data()[0] << std::endl;
                stream << Logger::Message("k: ") << mat_tis.ThermalConductivity().data()[0] << std::endl;
                stream << Logger::Message("c: ") << mat_tis.SpecificHeatLp().data()[0] << std::endl;
                if (method == "fem") {
                    stream << Logger::Message("Total Power: ") << electrical_tis.TotalPower(mesh_tis, mat_tis.ElectricalConductivity()) << std::endl;
                } else {
                    stream << Logger::Message("Total Power: ") << electrical_tis.TotalPower(voro_tis, mat_tis.ElectricalConductivity()) << std::endl;
                }
            }
        }

        // Compute temperature. Assembling in each step since fields are temperature-dependent.
        if (method == "fem") {
            if (temp_increment > update_thres) {
                if (parser.HasAttribute("tissue")) {
                    bioheat_tis.UpdateSystem(mesh_tis, mat_tis.ThermalConductivity(), mat_tis.Enthalpy(), rf_heat_tis);
                }
            } else {
                if (parser.HasAttribute("tissue")) {
                    bioheat_tis.UpdateOnlyHeatSource(mesh_tis, rf_heat_tis);
                }
            }
        } else if (method == "fpm") {
            if (temp_increment > update_thres) {
                if (parser.HasAttribute("tissue")) {
                    bioheat_tis.UpdateSystem(voro_tis, fpm_tis,
                        mat_tis.ThermalConductivity(), mat_tis.Enthalpy(), rf_heat_tis);
                }
            } else {
                if (parser.HasAttribute("tissue")) {
                    bioheat_tis.UpdateOnlyHeatSource(voro_tis, fpm_tis, rf_heat_tis);
                }
            }
        }

        // Compute temperature.
        if (parser.HasAttribute("tissue")) {
            bioheat_tis.ComputeTemperature(bioheat_bc_tis.Dirichlet());
        }

        if ((step+1) % output_interval == 0) {
            // Store voltage and temperature.
            if (parser.HasAttribute("catheter")) {
                electrical_cath.StoreCurrentVoltageAt(pos);
                bioheat_cath.StoreCurrentTemperatureAt(pos);
            }
            if (parser.HasAttribute("tissue")) {
                electrical_tis.StoreCurrentVoltageAt(pos);
                bioheat_tis.StoreCurrentTemperatureAt(pos);
            }
            pos++;
        }

        if (mpi_handler.rank_id == 0)
            stream << Logger::Message("Completed bioheat simulation steps: " + std::to_string(step+1) + "/" + std::to_string(steps_num) + "\n");

    } // End of Iterate over time.
    if (mpi_handler.rank_id == 0) { stream << std::endl; }

    // Store final state of current potential if the last step was not saved in the loop.
    if (output_interval != 0) {
        if (steps_counter % output_interval != 0) {
            if (parser.HasAttribute("catheter")) {
                electrical_cath.StoreCurrentVoltageAt(pos+1);
                bioheat_cath.StoreCurrentTemperatureAt(pos+1);
            }
            if (parser.HasAttribute("tissue")) {
                electrical_tis.StoreCurrentVoltageAt(pos+1);
                bioheat_tis.StoreCurrentTemperatureAt(pos+1);
            }
        }
    }
    else {  // Store final state if output_interval is zero.
        if (parser.HasAttribute("catheter")) {
            electrical_cath.StoreCurrentVoltageAt(1);
            bioheat_cath.StoreCurrentTemperatureAt(1);
        }
        if (parser.HasAttribute("tissue")) {
            electrical_tis.StoreCurrentVoltageAt(1);
            bioheat_tis.StoreCurrentTemperatureAt(1);
        }
    }

    if (mpi_handler.rank_id == 0)
        stream << Logger::Message("Assemble time: "+assemble_time+"\n");
}

} // end of namespace PNTSIM

#endif //PHYNETOUCH_APPS_TOOLS_CONFIG_PHYSICS_TPP_