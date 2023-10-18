/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */

#ifndef PHYNETOUCH_PHYSICS_DEFORMATION_TPP_
#define PHYNETOUCH_PHYSICS_DEFORMATION_TPP_


#include "PHYNETOUCH/engine/physics/deformation.hpp"


namespace PNT {


template<short DIM, short CELL_NODES>
Deformation<DIM, CELL_NODES>::Deformation() : disp_history_(), lumped_mass_(),
    cell_volumes_(), nodal_volumes_(), attached_cells_to_nodes_(), sim_time_(0.), dt_(-1.), dt_crit_(-1.),
    thread_looper_(), threads_num_(), tmutex_()
{
    // Get the number of available threads.
    std::size_t available = std::thread::hardware_concurrency()-1;
    this->threads_num_ = std::max(available, 1ul);

    // Solve laplacian for free nodes with BICGSTAB solver.
    Eigen::initParallel();
    Eigen::setNbThreads(this->threads_num_);
}


template<short DIM, short CELL_NODES>
Deformation<DIM, CELL_NODES>::~Deformation()
{}


template<short DIM, short CELL_NODES>
void Deformation<DIM, CELL_NODES>::FormLumpedMass(const IMP::Mesh<DIM,CELL_NODES> &tis_mesh, const IMP::Mesh<DIM,CELL_NODES> &cath_mesh,
    const CLOUDEA::FemMats<DIM,CELL_NODES> &tis_fem, const CLOUDEA::FemMats<DIM,CELL_NODES> &cath_fem,
    const std::shared_ptr<Constitutive> &tis_elastic, const std::shared_ptr<Constitutive> &cath_elastic, bool mass_scaling)
{

    // Get quadrature points for tissue and catheter.
    auto tis_qnum = int{0};
    auto tis_quad_per_cell = int{0};
    if (tis_mesh.CellsNum() != 0) {
        tis_qnum = static_cast<int>(tis_fem.GaussWeights().size());
        tis_quad_per_cell = tis_fem.GaussPerCell();
    }
    auto cath_qnum = int{0};
    auto cath_quad_per_cell = int{0};
    if (cath_mesh.CellsNum() != 0) {
        cath_qnum = static_cast<int>(cath_fem.GaussWeights().size());
        cath_quad_per_cell = cath_fem.GaussPerCell();
    }
    auto qnum = tis_qnum+cath_qnum;

    auto tis_nnum = tis_mesh.NodesNum();
    auto cath_nnum = cath_mesh.NodesNum();
    auto nnum = tis_nnum + cath_nnum;

    // Compute critical time step per quadrature point.
    auto qcrit_dts = std::vector<double>(qnum,0.);
    auto min_tis_qcrit_dt = std::numeric_limits<double>::max();
    // Tissue quadrature points.
    for (int qid = 0; qid != tis_qnum; ++qid) {
        auto cid = qid/tis_quad_per_cell;

        // Compute the sum of the wave speed over the cell nodes.
        auto wave_speed = double{0.};
        for (const auto &nid : tis_mesh.Cells(cid).Connectivity()) {
            wave_speed += tis_elastic->WaveSpeed(nid);
        }
        wave_speed /= static_cast<double>(CELL_NODES);

        // Compute integration point critical time step.
        auto lmax = wave_speed*std::sqrt(static_cast<double>(CELL_NODES) * (tis_fem.Grads(qid).array()*tis_fem.Grads(qid).array()).sum());
        qcrit_dts[qid] = 2./lmax;
        if (qcrit_dts[qid] < min_tis_qcrit_dt) {
            min_tis_qcrit_dt = qcrit_dts[qid];
        }
    }
    // Catheter quadrature points.
    if (cath_qnum > 0) {
        if (cath_elastic->Type() == ConstitutiveType::rigid) {
            for (int qid = 0; qid != cath_qnum; ++qid) {
                qcrit_dts[tis_qnum+qid] = min_tis_qcrit_dt;
            }
        } else {
            for (int qid = 0; qid != cath_qnum; ++qid) {
                auto cid = qid/cath_quad_per_cell;

                // Compute the sum of the wave speed over the cell nodes.
                auto wave_speed = double{0.};
                for (const auto &nid : cath_mesh.Cells(cid).Connectivity()) {
                    wave_speed += cath_elastic->WaveSpeed(nid);
                }
                wave_speed /= static_cast<double>(CELL_NODES);

                // Compute integration point critical time step.
                auto lmax = wave_speed*std::sqrt(static_cast<double>(CELL_NODES) * (cath_fem.Grads(qid).array()*cath_fem.Grads(qid).array()).sum());
                qcrit_dts[tis_qnum+qid] = 2./lmax;

            }
        }
    }

    // Compute the critical time step and the mass scaling coefficients.
    auto scaling_coeffs = std::vector<double>(qnum,1.);
    if (mass_scaling) {
        this->dt_crit_ = *std::max_element(std::begin(qcrit_dts), std::end(qcrit_dts));

        for (std::size_t i = 0; i != scaling_coeffs.size(); ++i) {
            scaling_coeffs[i] = std::pow(this->dt_crit_/qcrit_dts[i], 2.);
        }
    } else {
        this->dt_crit_ = *std::min_element(std::begin(qcrit_dts), std::end(qcrit_dts));
    }
    // Apply a safety coefficient to the critical time step.
    this->dt_crit_ *= 0.9;

    // Compute lumped mass vector
    auto cell_mass = Eigen::MatrixXd{};
    auto mass_vec = Eigen::VectorXd(nnum);
    mass_vec.setZero();
    // Tissue mesh quadrature points.
    for (int qid = 0; qid != tis_qnum; ++qid) {
        // Mean cell density.
        auto cid = qid/tis_quad_per_cell;
        auto cell_density = double{0.};
        for (const auto &nid : tis_mesh.Cells(cid).Connectivity()) {
            cell_density += tis_elastic->Density(nid);
        }
        cell_density /= static_cast<double>(CELL_NODES);

        // Compute cell mass matrix.
        cell_mass = scaling_coeffs[qid] * cell_density * tis_fem.GaussWeights(qid) *
            tis_fem.JacobianDets(qid)*(tis_fem.ShapeFuncs(qid)*tis_fem.ShapeFuncs(qid).transpose());

        // Distribute cell mass to nodes.
        for (short i = 0; i != CELL_NODES; ++i) {
            mass_vec.coeffRef(tis_mesh.Cells(cid).N(i)) += cell_mass.row(i).sum();
        }
    }
    // Catheter mesh quadrature points.
    for (int qid = 0; qid != cath_qnum; ++qid) {
        // Mean cell density.
        auto cid = qid/cath_quad_per_cell;
        auto cell_density = double{0.};
        for (const auto &nid : cath_mesh.Cells(cid).Connectivity()) {
            cell_density += cath_elastic->Density(nid);
        }
        cell_density /= static_cast<double>(CELL_NODES);

        // Compute cell mass matrix.
        cell_mass = scaling_coeffs[tis_qnum+qid] * cell_density * cath_fem.GaussWeights(qid) *
            cath_fem.JacobianDets(qid)*(cath_fem.ShapeFuncs(qid)*cath_fem.ShapeFuncs(qid).transpose());

        // Distribute cell mass to nodes.
        for (short i = 0; i != CELL_NODES; ++i) {
            mass_vec.coeffRef(tis_nnum+cath_mesh.Cells(cid).N(i)) += cell_mass.row(i).sum();
        }
    }

    // Create a lumped mass matrix.
    this->lumped_mass_ = Eigen::MatrixXd::Zero(nnum,DIM);
    for (short d = 0; d != DIM; ++d) { this->lumped_mass_.col(d) = mass_vec; }
}


template<short DIM, short CELL_NODES>
void Deformation<DIM, CELL_NODES>::FormLumpedMass(const IMP::Voronoi<DIM> &tis_voro, const IMP::Voronoi<DIM> &cath_voro,
    const CLOUDEA::Fpm<DIM> &tis_fpm, const CLOUDEA::Fpm<DIM> &cath_fpm,
    const std::shared_ptr<Constitutive> &tis_elastic, const std::shared_ptr<Constitutive> &cath_elastic, bool mass_scaling)
{
    // Get nodes number for tissue and catheter.
    auto tis_nnum = tis_voro.NodesNum();
    auto cath_nnum = cath_voro.NodesNum();
    auto nnum = tis_nnum + cath_nnum;

    // Compute critical time step per node.
    auto ncrit_dts = std::vector<double>(nnum,0.);
    auto neighs_num = std::size_t{0};
    auto min_tis_ncrit_dt = std::numeric_limits<double>::max();
    auto lmax = double{0.};
    // Tissue nodes.
    for (int nid = 0; nid != tis_nnum; ++nid) {
        neighs_num = tis_fpm.Support().InfluenceNodeIds(nid).size();
        lmax = tis_elastic->WaveSpeed(nid)*std::sqrt(neighs_num * (tis_fpm.PhiGrad(nid).array()*tis_fpm.PhiGrad(nid).array()).sum());
        ncrit_dts[nid] = 2./lmax;
        if (ncrit_dts[nid] < min_tis_ncrit_dt) {
            min_tis_ncrit_dt = ncrit_dts[nid];
        }
    }
    // Catheter nodes.
    if (cath_nnum > 0) {
        if (cath_elastic->Type() == ConstitutiveType::rigid) {
            for (int nid = 0; nid != cath_nnum; ++nid) {
                ncrit_dts[tis_nnum+nid] = min_tis_ncrit_dt;
            }
        } else {
            for (int nid = 0; nid != cath_nnum; ++nid) {
                neighs_num = cath_fpm.Support().InfluenceNodeIds(nid).size();
                lmax = cath_elastic->WaveSpeed(nid)*std::sqrt(neighs_num * (cath_fpm.PhiGrad(nid).array()*cath_fpm.PhiGrad(nid).array()).sum());
                ncrit_dts[tis_nnum+nid] = 2./lmax;
            }
        }
    }

    auto scaling_coeffs = std::vector<double>(nnum,1.);
    if (mass_scaling) {
        // Set final critical time step to max since there is mass scaling.
        this->dt_crit_ = *std::max_element(std::begin(ncrit_dts), std::end(ncrit_dts));

        // Compute the mass scaling coefficients.
        for (int nid = 0; nid != nnum; ++nid) {
            scaling_coeffs[nid] = std::pow(this->dt_crit_/ncrit_dts[nid], 2.);
        }
    } else {
        this->dt_crit_ = *std::min_element(std::begin(ncrit_dts), std::end(ncrit_dts));
    }
    // Apply a safety coefficient to the critical time step.
    this->dt_crit_ *= 0.9 / std::sqrt(tis_fpm.Penalty());


    // Compute lumped mass vector
    auto mass_vec = Eigen::VectorXd(nnum);
    mass_vec.setZero();
    auto N1 = Eigen::VectorXd{};
    auto N1_t = Eigen::RowVectorXd{};
    auto cell_mass = Eigen::MatrixXd{};
    auto cell_centroid = IMP::Vec<DIM,double>{};
    // Tissue nodes.
    for (int nid = 0; nid != tis_nnum; ++nid) {

        // Compute shape functions for parent and neighbor cells.
        cell_centroid = tis_voro.CellCentroid(nid);
        N1 = tis_fpm.PhiGrad(nid) * (cell_centroid - tis_voro.Nodes(nid)).CopyToEigen();
        N1.coeffRef(0) = N1.coeff(0) + 1.;
        N1_t = N1.transpose();

        // Compute cell mass matrix.
        cell_mass = scaling_coeffs[nid] * tis_elastic->Density(nid) * tis_voro.CellMeasures(nid) * (N1*N1_t);

        // Distribute cell mass to nodes.
        auto i = int{0};
        for (const auto &neigh_id : tis_fpm.Support().InfluenceNodeIds(nid)) {
            mass_vec.coeffRef(neigh_id) += cell_mass.row(i++).sum();
        }
    }
    // Catheter nodes.
    for (int nid = 0; nid != cath_nnum; ++nid) {

        // Compute shape functions for parent and neighbor cells.
        cell_centroid = cath_voro.CellCentroid(nid);
        N1 = cath_fpm.PhiGrad(nid) * (cell_centroid - cath_voro.Nodes(nid)).CopyToEigen();
        N1.coeffRef(0) = N1.coeff(0) + 1.;
        N1_t = N1.transpose();

        // Compute cell mass matrix.
        cell_mass = scaling_coeffs[tis_nnum+nid] * cath_elastic->Density(nid) * cath_voro.CellMeasures(nid) * (N1*N1_t);

        // Distribute cell mass to nodes.
        auto i = int{0};
        for (const auto &neigh_id : cath_fpm.Support().InfluenceNodeIds(nid)) {
            mass_vec.coeffRef(tis_nnum+neigh_id) += cell_mass.row(i++).sum();
        }
    }

    // Create a lumped mass matrix.
    this->lumped_mass_ = Eigen::MatrixXd::Zero(nnum,DIM);
    for (short d = 0; d != DIM; ++d) { this->lumped_mass_.col(d) = mass_vec; }

}


template<short DIM, short CELL_NODES>
void Deformation<DIM, CELL_NODES>::SetSimulationTime(double sim_time)
{
    this->sim_time_ = sim_time;
}


template<short DIM, short CELL_NODES>
void Deformation<DIM, CELL_NODES>::SetDt(double dt)
{
    if (this->dt_crit_ < 0.) {
        auto err_msg = "Critical time step is not available. Form lump mass before setting the time step to ensure that critical time step is available.";
        throw std::runtime_error(Logger::Warning(err_msg));
    }

    if (dt > this->dt_crit_) {
        this->dt_ = this->dt_crit_;
        auto wrn_msg = "Selected deformation time step is larger than the critical time step. The critical time step will be used instead.\n";
        std::cout << Logger::Warning(wrn_msg);
    } else {
        this->dt_ = dt;
    }
}


template<short DIM, short CELL_NODES>
void Deformation<DIM, CELL_NODES>::ComputeCellVolumes(const IMP::Mesh<DIM,CELL_NODES> &tis_mesh,
    const IMP::Mesh<DIM,CELL_NODES> &cath_mesh, const CLOUDEA::FemMats<DIM,CELL_NODES> &tis_fem,
    const CLOUDEA::FemMats<DIM,CELL_NODES> &cath_fem)
{
    auto cnum_tis  = tis_mesh.CellsNum();
    auto cnum_cath = cath_mesh.CellsNum();

    this->cell_volumes_.clear();
    this->cell_volumes_.resize(cnum_tis+cnum_cath, 0.);

    auto cell_vol = double{0.};

    // Compute volume for each cell of the tissue mesh.
    for (int cid = 0; cid != cnum_tis; ++cid) {
        cell_vol = 0.;
        for (auto q = cid*tis_fem.GaussPerCell(); q != cid*tis_fem.GaussPerCell() + tis_fem.GaussPerCell(); ++q) {
            cell_vol += tis_fem.GaussWeights(q) * tis_fem.JacobianDets(q);
        }
        this->cell_volumes_[cid] = cell_vol;
    }

    // Compute volume for each cell of the catheter mesh.
    for (int cid = 0; cid != cnum_cath; ++cid) {
        cell_vol = 0.;
        for (auto q = cid*cath_fem.GaussPerCell(); q != cid*cath_fem.GaussPerCell() + cath_fem.GaussPerCell(); ++q) {
            cell_vol += cath_fem.GaussWeights(q) * cath_fem.JacobianDets(q);
        }
        this->cell_volumes_[cid+cnum_tis] = cell_vol;
    }

}


template<short DIM, short CELL_NODES>
void Deformation<DIM, CELL_NODES>::ComputeNodalVolumes(const IMP::Mesh<DIM,CELL_NODES> &tis_mesh,
    const IMP::Mesh<DIM,CELL_NODES> &cath_mesh)
{
    auto cnum_tis = tis_mesh.CellsNum();
    auto cnum_cath = cath_mesh.CellsNum();
    auto cnum  = static_cast<std::size_t>(cnum_tis + cnum_cath);
    if (this->cell_volumes_.size() != cnum) {
        auto err_msg = "Could not compute nodal volumes. Compute cell volumes first.";
        throw std::runtime_error(Logger::Error(err_msg));
    }

    auto nnum_tis = tis_mesh.NodesNum();
    auto nnum_cath = cath_mesh.NodesNum();

    this->nodal_volumes_.clear();
    this->nodal_volumes_.resize(nnum_tis+nnum_cath, 0.);

    this->attached_cells_to_nodes_.clear();
    this->attached_cells_to_nodes_.resize(nnum_tis+nnum_cath);

    auto attached_cell_ids = std::vector<int>{};

    // Compute nodal volumes for tissue mesh.
    for (int nid = 0; nid != nnum_tis; ++nid) {
        attached_cell_ids = tis_mesh.AttachedCellIdsToNode(nid);
        for (const auto &cid : attached_cell_ids) {
            this->nodal_volumes_[nid] += (this->cell_volumes_[cid] / static_cast<double>(CELL_NODES));
        }
        this->attached_cells_to_nodes_[nid] = attached_cell_ids;
    }

    // Nodal volumes for catheter mesh.
    for (int nid = 0; nid != nnum_cath; ++nid) {
        attached_cell_ids = cath_mesh.AttachedCellIdsToNode(nid);
        for (const auto &cid : attached_cell_ids) {
            this->nodal_volumes_[nid+nnum_tis] += (this->cell_volumes_[cid+cnum_tis] / static_cast<double>(CELL_NODES));
        }
        this->attached_cells_to_nodes_[nid+nnum_tis] = attached_cell_ids;
    }
}


template<short DIM, short CELL_NODES>
void Deformation<DIM, CELL_NODES>::UpdateInternalForces(const IMP::Mesh<DIM,CELL_NODES> &tis_mesh, const IMP::Mesh<DIM,CELL_NODES> &cath_mesh,
    const CLOUDEA::FemMats<DIM,CELL_NODES> &tis_fem, const CLOUDEA::FemMats<DIM,CELL_NODES> &cath_fem,
    const std::shared_ptr<Constitutive> &tis_elastic, const std::shared_ptr<Constitutive> &cath_elastic,
    const Eigen::MatrixXd &disp, Eigen::MatrixXd &forces)
{

    auto CellPressureCallback = [](std::size_t thread_id, std::size_t gauss_per_cell,
        const CLOUDEA::ThreadLoopManager &thread_looper, const IMP::Mesh<DIM,CELL_NODES> &mesh,
        const std::shared_ptr<Constitutive> &hyperelastic, const CLOUDEA::FemMats<DIM,CELL_NODES> &fem,
        const Eigen::MatrixXd &disps, Eigen::VectorXd &tpressures) {

        // Initialize pressures in the thread.
        tpressures = Eigen::VectorXd::Zero(mesh.CellsNum()*gauss_per_cell);

        // Initialize matrices to compute cell contribution to forces.
        Eigen::MatrixXd u_cell = Eigen::MatrixXd::Zero(CELL_NODES,DIM);     // Cell displacements matrix.
        Eigen::MatrixXd FT = Eigen::MatrixXd::Zero(DIM,DIM);                // Deformation gradient.
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(DIM,DIM);             // Identity matrix.

        // Compute pressure of each cell in the thread.
        for (auto cid = thread_looper.LoopStartId(thread_id);
             cid != thread_looper.LoopEndId(thread_id); ++cid) {

            // Compute average bulk modulus of the cell.
            auto kappa_cell = double{0.};
            for (short i = 0; i != CELL_NODES; ++i) {
                u_cell.row(i) = disps.row(mesh.Cells(cid).N(i));
                kappa_cell += hyperelastic->BulkModulus(mesh.Cells(cid).N(i));
            }
            kappa_cell /= static_cast<double>(CELL_NODES);

            // Compute cell pressure on gauss points.
            for (auto q = cid*gauss_per_cell; q != cid*gauss_per_cell + gauss_per_cell; ++q) {
                // Compute deformation gradient of the cell gauss point.
                FT = fem.Grads(q).transpose()*u_cell + I;
                auto J = FT.determinant();

                // Compute pressure of the cell gauss point.
                tpressures.coeffRef(q) = kappa_cell * (J - 1.);
            }
        }

    };


    // Lambda callback for multithreaded force computation per mesh structure.
    auto ForceCallback = [](std::size_t thread_id, std::size_t gauss_per_cell, const CLOUDEA::ThreadLoopManager &thread_looper,
        const IMP::Mesh<DIM,CELL_NODES> &mesh, const std::shared_ptr<Constitutive> &hyperelastic,
        const CLOUDEA::FemMats<DIM,CELL_NODES> &fem, const std::vector<double> &nodal_vols,
        const std::vector<double> &cell_vols, const Eigen::VectorXd &cell_pres,
        const std::vector<std::vector<int>> &attached_cells_to_nodes,
        const Eigen::MatrixXd &u, Eigen::MatrixXd &tforces) {

        // Initialize forces in the thread.
        tforces = Eigen::MatrixXd::Zero(hyperelastic->NodesNum(),DIM);

        // Initialize matrices to compute cell contribution to forces.
        Eigen::MatrixXd u_cell = Eigen::MatrixXd::Zero(CELL_NODES,DIM);     // Cell displacements matrix.
        Eigen::MatrixXd f_cell = Eigen::MatrixXd::Zero(CELL_NODES,DIM);     // Cell forces matrix.
        Eigen::MatrixXd FT = Eigen::MatrixXd::Zero(DIM,DIM);                // Deformation gradient.
        Eigen::MatrixXd S = Eigen::MatrixXd::Zero(DIM,DIM);                 // Second Piola-Kirchhoff stress.
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(DIM,DIM);             // Identity matrix.

        // Compute force contribution of each cell of the mesh.
        for (auto cid = thread_looper.LoopStartId(thread_id);
             cid != thread_looper.LoopEndId(thread_id); ++cid) {

            // Get displacements at cell nodes
            // and mean hyperelastic parameters.
            auto mu_cell = double{0.};
            auto kappa_cell = double{0.};
            for (short i = 0; i != CELL_NODES; ++i) {
                u_cell.row(i) = u.row(mesh.Cells(cid).N(i));
                mu_cell += hyperelastic->LameMu(mesh.Cells(cid).N(i));
                kappa_cell += hyperelastic->BulkModulus(mesh.Cells(cid).N(i));
            }
            mu_cell /= static_cast<double>(CELL_NODES);
            kappa_cell /= static_cast<double>(CELL_NODES);

            // Iterate over quadrature point weights of the cell.
            f_cell.setZero();
            if (hyperelastic->Type() != ConstitutiveType::rigid) {
                auto q_it = int{0};
                for (auto q = cid*gauss_per_cell; q != cid*gauss_per_cell + gauss_per_cell; ++q) {

                    // Compute deformation gradient of the cell at gauss point.
                    FT = fem.Grads(q).transpose()*u_cell + I;
                    auto J = FT.determinant();

                    // Compute ANP cell pressure at gauss point.
                    auto pe_bar = double{0.};
                    for (const auto &nid : mesh.Cells(cid).Connectivity()) {

                        // Compute nodal pressure using the attached cells to the node.
                        auto pa = double{0.};
                        for (const auto &att_id : attached_cells_to_nodes[nid]) {
                            pa += cell_pres.coeff(att_id*gauss_per_cell+q_it) * cell_vols[att_id];
                        }
                        pa /= nodal_vols[nid] * static_cast<double>(CELL_NODES);
                        pe_bar += pa;
                    }
                    pe_bar /= static_cast<double>(CELL_NODES);
                    auto J_bar = (pe_bar / kappa_cell) + 1.;

                    // Modify deformation gradient.
                    if (J > std::sqrt(std::numeric_limits<double>::epsilon())) {
                        FT = std::pow(J_bar,1./3.)*std::pow(J,-(1./3.))*FT;
                    }

                    // Compute second Piola-Kirchhoff stress.
                    hyperelastic->ComputeStressSPK(FT, I, mu_cell, kappa_cell, S);

                    // Compute force contribution of the quadrature point.
                    f_cell += fem.GaussWeights(q) * fem.JacobianDets(q) * fem.Grads(q) * S*FT;
                    q_it++;
                }

                // Add contribution from the cell to the forces in the thread.
                for (short i = 0; i != CELL_NODES; ++i) {
                    tforces.row(mesh.Cells(cid).N(i)) += f_cell.row(i);
                }
            }
        } // End of Compute force contribution of each cell of the mesh.
    };

    auto nnum_tis = tis_mesh.NodesNum();
    auto cnum_tis = tis_mesh.CellsNum();

    auto nnum_cath = cath_mesh.NodesNum();
    auto cnum_cath = cath_mesh.CellsNum();

    if (this->nodal_volumes_.size() != static_cast<std::size_t>(nnum_tis+nnum_cath)) {
        auto err_msg = "Could not compute internal forces. Compute nodal volumes first.";
        throw std::runtime_error(Logger::Error(err_msg));
    }

    std::vector<std::vector<int>>  tis_attached_cells_to_nodes =
        {this->attached_cells_to_nodes_.begin(), this->attached_cells_to_nodes_.begin() + nnum_tis};

    std::vector<std::vector<int>> cath_attached_cells_to_nodes =
        {this->attached_cells_to_nodes_.begin() + nnum_tis, this->attached_cells_to_nodes_.end()};

    std::vector<double>  tis_nodal_volumes = {this->nodal_volumes_.begin(), this->nodal_volumes_.begin() + nnum_tis};
    std::vector<double> cath_nodal_volumes = {this->nodal_volumes_.begin() + nnum_tis, this->nodal_volumes_.end()};

    std::vector<double>  tis_cell_volumes = {this->cell_volumes_.begin(), this->cell_volumes_.begin() + cnum_tis};
    std::vector<double> cath_cell_volumes = {this->cell_volumes_.begin() + cnum_tis, this->cell_volumes_.end()};

    Eigen::VectorXd  tis_cell_pressures = Eigen::VectorXd::Zero(cnum_tis*tis_fem.GaussPerCell());
    Eigen::VectorXd cath_cell_pressures = Eigen::VectorXd::Zero(cnum_cath*cath_fem.GaussPerCell());

    // Compute tissue mesh cells pressure.
    if (tis_mesh.NodesNum() != 0) {
        if (tis_elastic->Type() != ConstitutiveType::rigid) {
            auto workers = std::vector<std::thread>{};
            workers.reserve(this->threads_num_);

            auto th_pressures = std::vector<Eigen::VectorXd>{};
            th_pressures.resize(this->threads_num_);

            this->thread_looper_.SetLoopRanges(static_cast<std::size_t>(tis_mesh.CellsNum()), this->threads_num_);

            // Get displacements for tissue nodes and volumes for tissue cells.
            Eigen::MatrixXd tis_disp = disp.topRows(tis_mesh.NodesNum());

            // Multithreaded pressure computation.
            for (std::size_t t = 0; t != this->threads_num_; ++t) {
                workers.emplace_back(std::thread(CellPressureCallback, t, tis_fem.GaussPerCell(),
                std::cref(this->thread_looper_), std::cref(tis_mesh), std::cref(tis_elastic),
                std::cref(tis_fem), std::cref(tis_disp), std::ref(th_pressures[t])));
            }
            std::for_each(workers.begin(), workers.end(), std::mem_fn(&std::thread::join));

            // Add tissue cell pressures contribution.
            tis_cell_pressures = th_pressures[0];
            for (std::size_t t = 1; t != this->threads_num_; ++t) {
                tis_cell_pressures += th_pressures[t];
            }
        }
    }

    // Compute catheter mesh cells pressure.
    if (cath_mesh.NodesNum() != 0) {
        if (cath_elastic->Type() != ConstitutiveType::rigid) {
            auto workers = std::vector<std::thread>{};
            workers.reserve(this->threads_num_);

            auto th_pressures = std::vector<Eigen::VectorXd>{};
            th_pressures.resize(this->threads_num_);

            this->thread_looper_.SetLoopRanges(static_cast<std::size_t>(cath_mesh.CellsNum()), this->threads_num_);

            // Get displacements for tissue nodes and volumes for tissue cells.
            Eigen::MatrixXd cath_disp = disp.bottomRows(cath_mesh.NodesNum());

            // Multithreaded pressure computation.
            for (std::size_t t = 0; t != this->threads_num_; ++t) {
                workers.emplace_back(std::thread(CellPressureCallback, t, cath_fem.GaussPerCell(),
                std::cref(this->thread_looper_), std::cref(cath_mesh), std::cref(cath_elastic),
                std::cref(cath_fem), std::cref(cath_disp), std::ref(th_pressures[t])));
            }
            std::for_each(workers.begin(), workers.end(), std::mem_fn(&std::thread::join));

            // Add tissue cell pressures contribution.
            cath_cell_pressures = th_pressures[0];
            for (std::size_t t = 1; t != this->threads_num_; ++t) {
                cath_cell_pressures += th_pressures[t];
            }
        }
    }

    // Initialize forces.
    forces = Eigen::MatrixXd::Zero(nnum_tis+nnum_cath,DIM);

    // Compute tissue forces.
    if (tis_mesh.NodesNum() != 0) {
        if (tis_elastic->Type() != ConstitutiveType::rigid) {
            auto workers = std::vector<std::thread>{};
            workers.reserve(this->threads_num_);

            auto th_forces = std::vector<Eigen::MatrixXd>{};
            th_forces.resize(this->threads_num_);

            this->thread_looper_.SetLoopRanges(static_cast<std::size_t>(tis_mesh.CellsNum()), this->threads_num_);

            // Multithreaded force computation.
            Eigen::MatrixXd tis_disp = disp.topRows(tis_mesh.NodesNum());
            for (std::size_t t = 0; t != this->threads_num_; ++t) {
                workers.emplace_back(std::thread(ForceCallback, t, tis_fem.GaussPerCell(), std::cref(this->thread_looper_),
                std::cref(tis_mesh), std::cref(tis_elastic), std::cref(tis_fem), std::cref(tis_nodal_volumes),
                std::cref(tis_cell_volumes), std::cref(tis_cell_pressures),
                std::cref(tis_attached_cells_to_nodes), std::cref(tis_disp), std::ref(th_forces[t])));
            }
            std::for_each(workers.begin(), workers.end(), std::mem_fn(&std::thread::join));

            // Add tissue forces contribution.
            forces.topRows(tis_mesh.NodesNum()) = th_forces[0];
            for (std::size_t t = 1; t != this->threads_num_; ++t) {
                forces.topRows(tis_mesh.NodesNum()) += th_forces[t];
            }
        }
    }

    // Compute catheter forces.
    if (cath_mesh.NodesNum() != 0) {
        if (cath_elastic->Type() != ConstitutiveType::rigid) {
            auto workers = std::vector<std::thread>{};
            workers.reserve(this->threads_num_);

            auto th_forces = std::vector<Eigen::MatrixXd>{};
            th_forces.resize(this->threads_num_);

            this->thread_looper_.SetLoopRanges(static_cast<std::size_t>(cath_mesh.CellsNum()), this->threads_num_);

            // Multithreaded force computation.
            Eigen::MatrixXd cath_disp = disp.bottomRows(cath_mesh.NodesNum());
            for (std::size_t t = 0; t != this->threads_num_; ++t) {
                workers.emplace_back(std::thread(ForceCallback, t, cath_fem.GaussPerCell(), std::cref(this->thread_looper_),
                std::cref(cath_mesh), std::cref(cath_elastic), std::cref(cath_fem), std::cref(cath_nodal_volumes),
                std::cref(cath_cell_volumes), std::cref(cath_cell_pressures),
                std::cref(cath_attached_cells_to_nodes), std::cref(cath_disp), std::ref(th_forces[t])));
            }
            std::for_each(workers.begin(), workers.end(), std::mem_fn(&std::thread::join));

            // Add catheter forces contribution.
            forces.bottomRows(cath_mesh.NodesNum()) = th_forces[0];
            for (std::size_t t = 1; t != this->threads_num_; ++t) {
                forces.bottomRows(cath_mesh.NodesNum()) += th_forces[t];
            }
        }
    }

}


template<short DIM, short CELL_NODES>
void Deformation<DIM, CELL_NODES>::UpdateInternalForces(const IMP::Voronoi<DIM> &tis_voro, const IMP::Voronoi<DIM> &cath_voro,
    const CLOUDEA::Fpm<DIM> &tis_fpm, const CLOUDEA::Fpm<DIM> &cath_fpm,
    const std::shared_ptr<Constitutive> &tis_elastic, const std::shared_ptr<Constitutive> &cath_elastic,
    const Eigen::MatrixXd &disp, Eigen::MatrixXd &forces)
{
    // Lambda callback for multithreaded force computation.
    auto ForceCallback = [](std::size_t thread_id, const CLOUDEA::ThreadLoopManager &thread_looper,
        const IMP::Voronoi<DIM> &voro, const std::shared_ptr<Constitutive> &hyperelastic, const CLOUDEA::Fpm<DIM> &fpm,
        const Eigen::MatrixXd &u, double &tpenalty, double &tdomain_measure, Eigen::MatrixXd &tforces) {

        // Initialize forces in the thread.
        tforces = Eigen::MatrixXd::Zero(hyperelastic->NodesNum(),DIM);

        // Initialize matrices to compute cell contribution to forces.
        auto u_neighs = Eigen::MatrixXd{};                                  // Neighbor nodes displacements matrix.
        auto f_neighs = Eigen::MatrixXd{};                                  // Neighbor nodes forces matrix.
        Eigen::MatrixXd FT = Eigen::MatrixXd::Zero(DIM,DIM);                 // Deformation gradient.
        Eigen::MatrixXd S = Eigen::MatrixXd::Zero(DIM,DIM);                 // Second Piola-Kirchhoff stress.
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(DIM,DIM);             // Identity matrix.

        // Compute force contribution of each node of the voronoi tesselation.
        tpenalty = 0.;
        tdomain_measure = 0.;
        auto neighs_num = std::size_t{0};
        for (auto nid = thread_looper.LoopStartId(thread_id);
             nid != thread_looper.LoopEndId(thread_id); ++nid) {

            // Get displacements at neighbor nodes.
            neighs_num = fpm.Support().InfluenceNodeIds(nid).size();
            u_neighs = Eigen::MatrixXd::Zero(neighs_num,DIM);
            for (std::size_t i = 0; i != neighs_num; ++i) {
                u_neighs.row(i) = u.row(fpm.Support().InfluenceNodeIds(nid)[i]);
            }

            // Compute deformation gradient of the cell.
            FT = fpm.PhiGrad(nid).transpose()*u_neighs + I;

            // Compute second Piola-Kirchhoff stress.
            hyperelastic->ComputeStressSPK(FT, I, hyperelastic->LameMu(nid), hyperelastic->BulkModulus(nid), S);

            // Compute force contribution on the neighbor nodes.
            f_neighs = voro.CellMeasures(nid) * fpm.PhiGrad(nid) * (S*FT);

            // Add force contribution to the forces in the thread.
            for (std::size_t i = 0; i != neighs_num; ++i) {
                tforces.row(fpm.Support().InfluenceNodeIds(nid)[i]) += f_neighs.row(i);
            }

            // Add contribution to penalty and domain measure.
            tpenalty += voro.CellMeasures(nid) * hyperelastic->YoungModulus(nid);
            tdomain_measure += voro.CellMeasures(nid);
        }
    };

    // Initialize forces.
    forces = Eigen::MatrixXd::Zero(tis_voro.NodesNum()+cath_voro.NodesNum(),DIM);

    // Compute tissue forces.
    if (tis_voro.NodesNum() != 0) {
        if (tis_elastic->Type() != ConstitutiveType::rigid) {
            auto workers = std::vector<std::thread>{};
            workers.reserve(this->threads_num_);

            this->thread_looper_.SetLoopRanges(static_cast<std::size_t>(tis_voro.NodesNum()), this->threads_num_);

            // Initialize thread forces, penalty, and domain measure.
            auto th_forces = std::vector<Eigen::MatrixXd>{};
            th_forces.resize(this->threads_num_);

            auto th_penal = std::vector<double>{};
            th_penal.resize(this->threads_num_);

            auto th_measure = std::vector<double>{};
            th_measure.resize(this->threads_num_);

            // Multithreaded force computation.
            Eigen::MatrixXd tis_disp = disp.topRows(tis_voro.NodesNum());
            for (std::size_t t = 0; t != this->threads_num_; ++t) {
                workers.emplace_back(std::thread(ForceCallback, t, std::cref(this->thread_looper_),
                    std::cref(tis_voro), std::cref(tis_elastic), std::cref(tis_fpm), std::cref(tis_disp),
                    std::ref(th_penal[t]), std::ref(th_measure[t]), std::ref(th_forces[t])));
            }
            std::for_each(workers.begin(), workers.end(), std::mem_fn(&std::thread::join));

            // Gather forces and penalty from the threads.
            auto tis_forces = th_forces[0];
            auto penal = th_penal[0];
            auto measure = th_measure[0];
            for (std::size_t t = 1; t != this->threads_num_; ++t) {
                tis_forces += th_forces[t];
                penal += th_penal[t];
                measure += th_measure[t];
            }
            penal = (penal/measure) * tis_fpm.Penalty();

            // Apply forces correction.
            auto forces_correction = Eigen::MatrixXd{};
            this->ComputeFpmForceCorrection(tis_voro, tis_fpm, tis_elastic, tis_disp, penal, forces_correction);
            tis_forces = tis_forces + forces_correction;

            // Add tissue forces contribution.
            forces.topRows(tis_voro.NodesNum()) = tis_forces;
        }
    }

    // Compute catheter forces.
    if (cath_voro.NodesNum() != 0) {
        if (cath_elastic->Type() != ConstitutiveType::rigid) {
            auto workers = std::vector<std::thread>{};
            workers.reserve(this->threads_num_);

            this->thread_looper_.SetLoopRanges(static_cast<std::size_t>(cath_voro.CellsNum()), this->threads_num_);

            // Initialize thread forces, penalty, and domain measure.
            auto th_forces = std::vector<Eigen::MatrixXd>{};
            th_forces.resize(this->threads_num_);

            auto th_penal = std::vector<double>{};
            th_penal.resize(this->threads_num_);

            auto th_measure = std::vector<double>{};
            th_measure.resize(this->threads_num_);

            // Multithreaded force computation.
            Eigen::MatrixXd cath_disp = disp.bottomRows(cath_voro.NodesNum());
            for (std::size_t t = 0; t != this->threads_num_; ++t) {
                workers.emplace_back(std::thread(ForceCallback, t, std::cref(this->thread_looper_),
                    std::cref(cath_voro), std::cref(cath_elastic), std::cref(cath_fpm), std::cref(cath_disp),
                    std::ref(th_penal[t]), std::ref(th_measure[t]), std::ref(th_forces[t])));
            }
            std::for_each(workers.begin(), workers.end(), std::mem_fn(&std::thread::join));

            // Gather forces and penalty from the threads.
            auto cath_forces = th_forces[0];
            auto penal = th_penal[0];
            auto measure = th_measure[0];
            for (std::size_t t = 1; t != this->threads_num_; ++t) {
                cath_forces += th_forces[t];
                penal += th_penal[t];
                measure += th_measure[t];
            }
            penal = (penal/measure) * cath_fpm.Penalty();

            // Apply forces correction.
            auto forces_correction = Eigen::MatrixXd{};
            this->ComputeFpmForceCorrection(cath_voro, cath_fpm, cath_elastic, cath_disp, penal, forces_correction);
            cath_forces += forces_correction;

            // Add tissue forces contribution.
            forces.bottomRows(cath_voro.NodesNum()) = cath_forces;
        }
    }
}


template<short DIM, short CELL_NODES>
void Deformation<DIM, CELL_NODES>::UpdateExternalForces(const IMP::Mesh<DIM,CELL_NODES> &tis_mesh, const IMP::Mesh<DIM,CELL_NODES> &cath_mesh,
    const CLOUDEA::FemMats<DIM,CELL_NODES> &tis_fem, const CLOUDEA::FemMats<DIM,CELL_NODES> &cath_fem,
    const BoundConds<DIM> &tis_deform_bc, const BoundConds<DIM> &cath_deform_bc, int step, Eigen::MatrixXd &forces)
{
    // Reset external forces.
    forces.setZero();

    // Apply body loads on tissue.
    if (tis_deform_bc.HasBodyLoad()) {

        // Update body loads conditions.
        auto body_load = Eigen::RowVectorXd(DIM);
        body_load.setZero();
        tis_deform_bc.BodyLoad().Apply(step, this->dt_, body_load);

        auto force_cell = Eigen::MatrixXd(CELL_NODES,DIM);
        for (auto cid = int{0}; cid != tis_mesh.CellsNum(); ++cid) {

            // Compute the load of the cell.
            force_cell.setZero();
            for (auto qid = cid*tis_fem.GaussPerCell(); qid != cid*tis_fem.GaussPerCell() + tis_fem.GaussPerCell(); ++qid) {
                force_cell += tis_fem.GaussWeights(qid)*tis_fem.JacobianDets(qid) * (tis_fem.ShapeFuncs(qid) * body_load);
            }

            // Assemble load forces.
            auto gi = int{0};
            for (short li = 0; li != CELL_NODES; ++li) {
                gi = tis_mesh.Cells(cid).N(li);
                forces.row(gi) += force_cell.row(li);
            }
        }
    }

    // Apply body loads on catheter.
    if (cath_deform_bc.HasBodyLoad()) {

        // Update body loads conditions.
        auto body_load = Eigen::RowVectorXd(DIM);
        body_load.setZero();
        cath_deform_bc.BodyLoad().Apply(step, this->dt_, body_load);

        // std::cout << body_load << std::endl;

        auto force_cell = Eigen::MatrixXd(CELL_NODES,DIM);
        for (auto cid = int{0}; cid != cath_mesh.CellsNum(); ++cid) {

            // Compute the load of the cell.
            force_cell.setZero();
            for (auto qid = cid*cath_fem.GaussPerCell(); qid != cid*cath_fem.GaussPerCell() + cath_fem.GaussPerCell(); ++qid) {
                force_cell += cath_fem.GaussWeights(qid)* cath_fem.JacobianDets(qid) * (cath_fem.ShapeFuncs(qid) * body_load);
                // std::cout << cath_fem.ShapeFuncs(qid) * body_load << std::endl << std::endl;
            }

            // Assemble load forces.
            auto gi = int{0};
            for (short li = 0; li != CELL_NODES; ++li) {
                gi = tis_mesh.NodesNum() + cath_mesh.Cells(cid).N(li);
                forces.row(gi) += force_cell.row(li);
            }
        }
    }
}


template<short DIM, short CELL_NODES>
void Deformation<DIM, CELL_NODES>::UpdateExternalForces(const IMP::Voronoi<DIM> &tis_voro, const IMP::Voronoi<DIM> &cath_voro,
    const CLOUDEA::Fpm<DIM> &tis_fpm, const CLOUDEA::Fpm<DIM> &cath_fpm,
    const BoundConds<DIM> &tis_deform_bc, const BoundConds<DIM> &cath_deform_bc, int step, Eigen::MatrixXd &forces)
{
    // Reset external forces.
    forces.setZero();

    // Apply body loads on tissue.
    auto force_cell = Eigen::MatrixXd{};
    if (tis_deform_bc.HasBodyLoad()) {

        // Update body loads conditions.
        auto body_load = Eigen::RowVectorXd(DIM);
        body_load.setZero();
        tis_deform_bc.BodyLoad().Apply(step, this->dt_, body_load);

        for (auto nid = int{0}; nid != tis_voro.NodesNum(); ++nid) {

            // Compute the load of the node's voronoi.
            auto cell_centroid = tis_voro.CellCentroid(nid);
            Eigen::VectorXd N1 = tis_fpm.PhiGrad(nid) * (cell_centroid - tis_voro.Nodes(nid)).CopyToEigen();
            N1.coeffRef(0) = N1.coeff(0) + 1.;

            force_cell = tis_voro.CellMeasures(nid) * (N1 * body_load);

            // Assemble load forces.
            auto gi = int{0};
            for (Eigen::Index li = 0; li != force_cell.rows(); ++li) {
                gi = tis_fpm.Support().InfluenceNodeIds(nid)[li];
                forces.row(gi) += force_cell.row(li);
            }
        }
    }

    // Apply body loads on catheter.
    if (cath_deform_bc.HasBodyLoad()) {

        // Update body loads conditions.
        auto body_load = Eigen::RowVectorXd(DIM);
        body_load.setZero();
        cath_deform_bc.BodyLoad().Apply(step, this->dt_, body_load);


        for (auto nid = int{0}; nid != cath_voro.NodesNum(); ++nid) {

            // Compute the load of the cell.
            auto cell_centroid = cath_voro.CellCentroid(nid);
            Eigen::VectorXd N1 = cath_fpm.PhiGrad(nid) * (cell_centroid - cath_voro.Nodes(nid)).CopyToEigen();
            N1.coeffRef(0) = N1.coeff(0) + 1.;
            force_cell = cath_voro.CellMeasures(nid) * (N1 * body_load);

            // Assemble load forces.
            auto gi = int{0};
            for (Eigen::Index li = 0; li != force_cell.rows(); ++li) {
                gi = tis_voro.NodesNum() + cath_fpm.Support().InfluenceNodeIds(nid)[li];
                forces.row(gi) += force_cell.row(li);
            }
        }
    }
}


template<short DIM, short CELL_NODES>
void Deformation<DIM, CELL_NODES>::UpdateDisplacements(const BoundConds<DIM> &tis_deform_bc, const BoundConds<DIM> &cath_deform_bc,
    int step, int tis_nodes_num, const Eigen::MatrixXd &ext_forces, const Eigen::MatrixXd &int_forces, const Eigen::MatrixXd &old_disp,
    const Eigen::MatrixXd &disp, Eigen::MatrixXd &new_disp, bool &disp_error_check)
{
    // Update displacements with central finite differences.
    new_disp = (this->dt_*this->dt_) * ((ext_forces-int_forces).array() / this->lumped_mass_.array()).matrix() + 2.*disp - old_disp;

    // Apply tissue Dirichlet boundary conditions.
    for (const auto &dirichlet : tis_deform_bc.Dirichlet()) {
        dirichlet.Apply(0, step, this->dt_, new_disp);

        // Ensure that solution is not diverging.
        auto abs_bc_val = std::fabs(dirichlet.Value());
        if (abs_bc_val > 0. && new_disp.cwiseAbs().maxCoeff() > 1.5*abs_bc_val) {
            auto wrn_msg = "Solution has become unbounded at step: "+std::to_string(step)+". Please reduce the time step.";
            disp_error_check = true;
            std::cout << termcolor::bold << termcolor::yellow << Logger::Warning(wrn_msg) << termcolor::reset << std::endl;
        }
    }

    // Apply catheter Dirichlet boundary conditions.
    for (const auto &dirichlet : cath_deform_bc.Dirichlet()) {
        dirichlet.Apply(tis_nodes_num, step, this->dt_, new_disp);

        // Ensure that solution is not diverging.
        auto abs_bc_val = std::fabs(dirichlet.Value());
        if (abs_bc_val > 0. && new_disp.cwiseAbs().maxCoeff() > 1.5*abs_bc_val) {
            auto wrn_msg = "Solution has become unbounded at step: "+std::to_string(step)+". Please reduce the time step.";
            disp_error_check = true;
            std::cout << termcolor::bold << termcolor::yellow << Logger::Warning(wrn_msg) << termcolor::reset << std::endl;
        }
    }

}


template<short DIM, short CELL_NODES>
void Deformation<DIM, CELL_NODES>::ApplyContact(const ContactHandle<DIM,CELL_NODES> &contact_handle, const IMP::Mesh<DIM,CELL_NODES> &master_mesh,
    const IMP::Mesh<DIM,CELL_NODES> &slave_mesh, int master_pad, int slave_pad, Eigen::MatrixXd &disp)
{
    // Get the master nodes initial coordinates and the displaced positions.
    auto master_pos = std::vector<IMP::Vec<DIM, double>>(contact_handle.MasterNodeIdsNum());
    for (int i = 0; i != contact_handle.MasterNodeIdsNum(); ++i) {
        for (short d = 0; d != DIM; ++d) {
            master_pos[i][d] = master_mesh.Nodes(contact_handle.MasterNodeIds(i))[d] + disp.coeff(contact_handle.MasterNodeIds(i)+master_pad,d);
        }
    }

    // Get the slave nodes initial coordinates and the displaced positions.
    auto slave_coords = std::vector<IMP::Vec<DIM, double>>(contact_handle.SlaveNodeIdsNum());
    auto slave_pos = slave_coords;
    for (int i = 0; i != contact_handle.SlaveNodeIdsNum(); ++i) {
        for (short d = 0; d != DIM; ++d) {
            slave_coords[i][d] = slave_mesh.Nodes(contact_handle.SlaveNodeIds(i))[d];
            slave_pos[i][d] = slave_coords[i][d] + disp.coeff(contact_handle.SlaveNodeIds(i)+slave_pad,d);
        }
    }

    // Update slave nodes position for contact.
    if constexpr (DIM == 2) {
        this->CorrectSlavePos2D(contact_handle, master_mesh, disp, master_pad, master_pos, slave_pos);
    } else if constexpr (DIM == 3) {
        this->CorrectSlavePos3D(contact_handle, master_mesh, disp, master_pad, master_pos, slave_pos);
    } else {
        auto err_msg = "Could not apply contact. Contact correction is available only for 2D and 3D problems.";
        throw std::runtime_error(Logger::Error(err_msg));
    }

    // Apply displacement correction to the slave nodes.
    for (int i = 0; i != contact_handle.SlaveNodeIdsNum(); ++i) {
        for (short d = 0; d != DIM; ++d) {
            disp.coeffRef(contact_handle.SlaveNodeIds(i)+slave_pad, d) = slave_pos[i][d] - slave_coords[i][d];
        }
    }

}


template<short DIM, short CELL_NODES>
void Deformation<DIM, CELL_NODES>::ApplyContact(ContactHandle<DIM,CELL_NODES> &contact_handle, const IMP::Voronoi<DIM> &master_voro,
        const IMP::Voronoi<DIM> &slave_voro, int master_pad, int slave_pad, Eigen::MatrixXd &disp)
{
    // Get the master nodes initial coordinates and the displaced positions.
    auto master_pos = std::vector<IMP::Vec<DIM, double>>(contact_handle.MasterNodeIdsNum());
    for (int i = 0; i != contact_handle.MasterNodeIdsNum(); ++i) {
        for (short d = 0; d != DIM; ++d) {
            master_pos[i][d] = master_voro.Nodes(contact_handle.MasterNodeIds(i))[d] + disp.coeff(contact_handle.MasterNodeIds(i)+master_pad,d);
        }
    }

    // Get the slave nodes initial coordinates and the displaced positions.
    auto slave_coords = std::vector<IMP::Vec<DIM, double>>(contact_handle.SlaveNodeIdsNum());
    auto slave_pos = slave_coords;
    for (int i = 0; i != contact_handle.SlaveNodeIdsNum(); ++i) {
        for (short d = 0; d != DIM; ++d) {
            slave_coords[i][d] = slave_voro.Nodes(contact_handle.SlaveNodeIds(i))[d];
            slave_pos[i][d] = slave_coords[i][d] + disp.coeff(contact_handle.SlaveNodeIds(i)+slave_pad,d);
        }
    }

    // Update slave nodes position for contact.
    if constexpr (DIM == 2) {
        // this->CorrectSlavePos2D(contact_handle, master_mesh, disp, master_pad, master_pos, slave_pos);
    } else if constexpr (DIM == 3) {
        this->CorrectSlavePos3D_v2(contact_handle, master_voro, disp, master_pad, master_pos, slave_pos);
    } else {
        auto err_msg = "Could not apply contact. Contact correction is available only for 2D and 3D problems.";
        throw std::runtime_error(Logger::Error(err_msg));
    }

    // Apply displacement correction to the slave nodes.
    for (int i = 0; i != contact_handle.SlaveNodeIdsNum(); ++i) {
        for (short d = 0; d != DIM; ++d) {
            disp.coeffRef(contact_handle.SlaveNodeIds(i)+slave_pad, d) = slave_pos[i][d] - slave_coords[i][d];
        }
    }

}


template<short DIM, short CELL_NODES>
void Deformation<DIM, CELL_NODES>::InitHistory(int instants_num, const Eigen::MatrixXd &disp)
{
    this->disp_history_.clear();
    this->disp_history_.resize(instants_num, disp);
}


template<short DIM, short CELL_NODES>
void Deformation<DIM, CELL_NODES>::AddInHistory(std::size_t pos, const Eigen::MatrixXd &disp)
{
    this->disp_history_[pos] = disp;
}



///!  PROTECTED  ///


template<short DIM, short CELL_NODES>
void Deformation<DIM, CELL_NODES>::CorrectSlavePos2D(const ContactHandle<DIM,CELL_NODES> &contact_handle, const IMP::Mesh<DIM,CELL_NODES> &master_mesh,
    const Eigen::MatrixXd &disp, int master_pad, const std::vector<IMP::Vec<DIM, double>> &master_pos, std::vector<IMP::Vec<DIM, double>> &slave_pos)
{
    auto x0 = IMP::Vec<DIM, double>{};
    auto x1 = IMP::Vec<DIM, double>{};
    auto tangent = IMP::Vec<DIM, double>{};
    auto et = IMP::Vec<DIM, double>{};
    auto en = IMP::Vec<DIM, double>{};
    auto vnode = IMP::Vec<DIM, double>{};

    auto near_master_to_slave = std::vector<int>{};
    contact_handle.NearestMasterToSlave(master_pos, slave_pos, near_master_to_slave);

    for (std::size_t i = 0; i != slave_pos.size(); ++i) {

        // Process attached facets of the master node.
        auto slave_corrected = false;
        auto penetrations = short{0};
        for (const auto &conn : contact_handle.MasterFacetsConn(near_master_to_slave[i])) {
            // Stop loop facets loop if the slave node got corrected.
            if (slave_corrected) {
                break;
            }

            // Get node indices of the facets connectivity.
            auto n0 = conn[0];
            auto n1 = conn[1];

            // Compute current position of the facet nodes.
            x0.Set({master_mesh.Nodes(n0)[0]+disp.coeff(n0+master_pad,0), master_mesh.Nodes(n0)[1]+disp.coeff(n0+master_pad,1)});
            x1.Set({master_mesh.Nodes(n1)[0]+disp.coeff(n1+master_pad,0), master_mesh.Nodes(n1)[1]+disp.coeff(n1+master_pad,1)});

            // Compute tangential and normal vectors to the facet.
            tangent = x1 - x0;
            et = tangent / tangent.Norm();
            en.Set({et[1], -et[0]});

            // Check penetration in facet.
            vnode = slave_pos[i] - x0;
            if (en.Dot(vnode) <= 0.) {
                // Compute natural coordinate.
                auto ksi = et.Dot(vnode) / tangent.Norm();
                if (ksi >= 0. && ksi <= 1.) {
                    slave_pos[i] = x0 + ksi*tangent;
                    slave_corrected = true;
                } else {
                    penetrations++;
                }
            }
        }

        if (penetrations == 2) {
            // Target slave position is the master node position.
            slave_pos[i] = master_pos[near_master_to_slave[i]];
        }
    }
}


template<short DIM, short CELL_NODES>
void Deformation<DIM, CELL_NODES>::CorrectSlavePos3D(const ContactHandle<DIM,CELL_NODES> &contact_handle, const IMP::Mesh<DIM,CELL_NODES> &master_mesh,
    const Eigen::MatrixXd &disp, int master_pad, const std::vector<IMP::Vec<DIM, double>> &master_pos, std::vector<IMP::Vec<DIM, double>> &slave_pos)
{
    auto x0 = IMP::Vec<DIM, double>{};
    auto x1 = IMP::Vec<DIM, double>{};
    auto x2 = IMP::Vec<DIM, double>{};

    auto r = IMP::Vec<DIM, double>{};
    auto w = IMP::Vec<DIM, double>{};

    auto e01 = IMP::Vec<DIM, double>{};
    auto e02 = IMP::Vec<DIM, double>{};

    auto en = IMP::Vec<DIM, double>{};

    auto approx_zero = std::sqrt(std::numeric_limits<double>::epsilon());

    auto near_master_to_slave = std::vector<int>{};
    contact_handle.NearestMasterToSlave(master_pos, slave_pos, near_master_to_slave);

    for (std::size_t i = 0; i != slave_pos.size(); ++i) {

        // Process attached facets of the master node.
        for (const auto &conn : contact_handle.MasterFacetsConn(near_master_to_slave[i])) {

            // Get node indices of the facet.
            auto n0 = conn[0];
            auto n1 = conn[1];
            auto n2 = conn[2];

            // Compute current position of the facet nodes.
            x0.Set({master_mesh.Nodes(n0)[0]+disp.coeff(n0+master_pad,0),
                master_mesh.Nodes(n0)[1]+disp.coeff(n0+master_pad,1), master_mesh.Nodes(n0)[2]+disp.coeff(n0+master_pad,2)});

            x1.Set({master_mesh.Nodes(n1)[0]+disp.coeff(n1+master_pad,0),
                master_mesh.Nodes(n1)[1]+disp.coeff(n1+master_pad,1), master_mesh.Nodes(n1)[2]+disp.coeff(n1+master_pad,2)});

            x2.Set({master_mesh.Nodes(n2)[0]+disp.coeff(n2+master_pad,0),
                master_mesh.Nodes(n2)[1]+disp.coeff(n2+master_pad,1), master_mesh.Nodes(n2)[2]+disp.coeff(n2+master_pad,2)});

            // Compute facet edges.
            e01 = x1-x0;
            e02 = x2-x0;

            // Compute normal vector to the facet.
            en = e01.Cross(e02);

            // // Compute the barycentric coordinates of the slave point in the facet.
            w = slave_pos[i] - x0;
            auto one_over_4area_squared = 1. / en.Dot(en);
            auto gamma = e01.Cross(w).Dot(en) * one_over_4area_squared;
            auto beta = w.Cross(e02).Dot(en) * one_over_4area_squared;
            auto alpha = 1. - beta - gamma;

            // Check if slave point is inside the facet.
            if (alpha > approx_zero && alpha < 1. && beta > approx_zero && beta < 1. && gamma > approx_zero && gamma < 1.) {
                // The projection of the slave node.
                r = alpha*x0 + beta*x1 + gamma*x2;
                if ((r-slave_pos[i]).Dot(en) > approx_zero) {
                    slave_pos[i] += (r-slave_pos[i]);
                    break;
                }
            }
        }
    }

}


template<short DIM, short CELL_NODES>
void Deformation<DIM, CELL_NODES>::CorrectSlavePos3D(ContactHandle<DIM,CELL_NODES> &contact_handle, const IMP::Voronoi<DIM> &master_voro,
    const Eigen::MatrixXd &disp, int master_pad, const std::vector<IMP::Vec<DIM, double>> &master_pos, std::vector<IMP::Vec<DIM, double>> &slave_pos)
{
    auto x0 = IMP::Vec<DIM, double>{};
    auto x1 = IMP::Vec<DIM, double>{};
    auto x2 = IMP::Vec<DIM, double>{};

    auto r = IMP::Vec<DIM, double>{};
    auto w = IMP::Vec<DIM, double>{};

    auto e01 = IMP::Vec<DIM, double>{};
    auto e02 = IMP::Vec<DIM, double>{};

    auto en = IMP::Vec<DIM, double>{};

    auto TL3 = double{1.01};
    auto approx_zero = 10.*(std::numeric_limits<double>::epsilon());

    // Find the iterator index of the nearest master nodes to the slave nodes.
    auto near_master_to_slave = std::vector<int>{};
    contact_handle.NearestMasterToSlave(master_pos, slave_pos, near_master_to_slave);

    int near_master_id = -1;
    for (std::size_t i = 0; i != slave_pos.size(); ++i) {

        near_master_id = near_master_to_slave[i];

        // Process attached facets of the master node.
        for (const auto &conn : contact_handle.MasterFacetsConn(near_master_id)) {

            // Get node indices of the facet.
            auto n0 = conn[0];
            auto p1 = conn[1];
            auto p2 = conn[2];

            // Compute current position of the facet vertices.
            x0.Set({master_voro.Nodes(n0)[0]+disp.coeff(n0+master_pad,0),
                master_voro.Nodes(n0)[1]+disp.coeff(n0+master_pad,1), master_voro.Nodes(n0)[2]+disp.coeff(n0+master_pad,2)});

            // We consider that the points of the facet have the same displacement as the node since the body is assumed rigid.
            x1.Set({master_voro.Points(p1)[0]+disp.coeff(n0+master_pad,0),
                master_voro.Points(p1)[1]+disp.coeff(n0+master_pad,1), master_voro.Points(p1)[2]+disp.coeff(n0+master_pad,2)});

            x2.Set({master_voro.Points(p2)[0]+disp.coeff(n0+master_pad,0),
                master_voro.Points(p2)[1]+disp.coeff(n0+master_pad,1), master_voro.Points(p2)[2]+disp.coeff(n0+master_pad,2)});

            // Compute facet edges.
            e01 = x1-x0;
            e02 = x2-x0;

            // Compute normal vector to the facet.
            en = e01.Cross(e02);

            // // Compute the barycentric coordinates of the slave point in the facet.
            w = slave_pos[i] - x0;
            auto one_over_4area_squared = 1. / en.Dot(en);
            auto gamma = e01.Cross(w).Dot(en) * one_over_4area_squared;
            auto beta = w.Cross(e02).Dot(en) * one_over_4area_squared;
            auto alpha = 1. - beta - gamma;

            // if((beta < -TL3) || (beta > TL3) || (gamma < -TL3) || (gamma > TL3))  continue;

            // r = alpha*x0 + beta*x1 + gamma*x2;
            // if ((r-slave_pos[i]).Dot(en) > approx_zero) {
            //     slave_pos[i] = r;
            //     break;
            // }

            // Check if slave point is inside the facet.
            if (alpha > approx_zero && alpha < 1. && beta > approx_zero && beta < 1. && gamma > approx_zero && gamma < 1.) {
                // The projection of the slave node.
                r = alpha*x0 + beta*x1 + gamma*x2;
                if ((r-slave_pos[i]).Dot(en) > approx_zero) {
                    slave_pos[i] = r;
                    break;
                }
            }
        }
    }
}


template<short DIM, short CELL_NODES>
void Deformation<DIM, CELL_NODES>::CorrectSlavePos3D_v2(ContactHandle<DIM,CELL_NODES> &contact_handle, const IMP::Voronoi<DIM> &master_voro,
    const Eigen::MatrixXd &disp, int master_pad, const std::vector<IMP::Vec<DIM, double>> &master_pos, std::vector<IMP::Vec<DIM, double>> &slave_pos)
{
    const auto EPS = double{1.0e-6};
    const auto TL1 = double{1.01};
    const auto TL2 = double{2.01};
    const auto TL3 = double{0.1};
    auto GN = double{0.0};

    auto x0 = IMP::Vec<DIM,double>{};
    auto x1 = IMP::Vec<DIM,double>{};
    auto x2 = IMP::Vec<DIM,double>{};

    auto XS = IMP::Vec<DIM,double>{};
    auto XC = IMP::Vec<DIM,double>{};
    auto XN = IMP::Vec<DIM,double>{};

    auto T1 = IMP::Vec<DIM,double>{};
    auto T2 = IMP::Vec<DIM,double>{};

    Eigen::Matrix2d A;
    Eigen::Vector2d B, DXI;

    // Lambda for mapping to isoparametric space.
    auto Cutl = [](const IMP::Vec<DIM,double> &slave_pos, const IMP::Vec<DIM,double> &x0,
        const IMP::Vec<DIM,double> &x1, const IMP::Vec<DIM,double> &x2, double ksi1, double ksi2,
        IMP::Vec<DIM,double> &T1, IMP::Vec<DIM,double> &T2, IMP::Vec<DIM,double> &XS, IMP::Vec<DIM,double> &XC) {

        XC = ksi1*x0 + ksi2*x1 + (1.0-ksi1-ksi2)*x2;        // projection of slave to triangle
        XS = slave_pos - XC;                                // vector from projection to actual position of slave
        T1 = x0 - x2;                                       // tangent vector no1 of triangle
        T2 = x1 - x2;                                       // tangent vector no2 of triangle
    };

    // Find the iterator index of the nearest master nodes to the slave nodes.
    auto near_master_to_slave = std::vector<int>{};
    contact_handle.NearestMasterToSlave(master_pos, slave_pos, near_master_to_slave);

    int near_master_id = -1;
    for (std::size_t i = 0; i != slave_pos.size(); ++i) {

        near_master_id = near_master_to_slave[i];

        // Process attached facets of the master node.
        for (const auto &conn : contact_handle.MasterFacetsConn(near_master_id)) {

            // Get node indices of the facet.
            auto n0 = conn[0];
            auto p1 = conn[1];
            auto p2 = conn[2];

            // Compute current position of the facet vertices.
            x0.Set({master_voro.Nodes(n0)[0]+disp.coeff(n0+master_pad,0),
                master_voro.Nodes(n0)[1]+disp.coeff(n0+master_pad,1), master_voro.Nodes(n0)[2]+disp.coeff(n0+master_pad,2)});

            // We consider that the points of the facet have the same displacement as the node since the body is assumed rigid.
            x1.Set({master_voro.Points(p1)[0]+disp.coeff(n0+master_pad,0),
                master_voro.Points(p1)[1]+disp.coeff(n0+master_pad,1), master_voro.Points(p1)[2]+disp.coeff(n0+master_pad,2)});

            x2.Set({master_voro.Points(p2)[0]+disp.coeff(n0+master_pad,0),
                master_voro.Points(p2)[1]+disp.coeff(n0+master_pad,1), master_voro.Points(p2)[2]+disp.coeff(n0+master_pad,2)});

            // Compute two tangent vectors and approximated contact point at the center.
            double ksi1 = 0.;
            double ksi2 = 0.;
            Cutl(slave_pos[i], x0, x1, x2, ksi1, ksi2, T1, T2, XS, XC);

            // Compute normal vector to the facet.
            XN = T1.Cross(T2);
            XN /= XN.Norm();
            GN = XN.Dot(XS);

            // Update isoparametric coords.
            ksi1 = XS.Dot(T1) / (2.*std::pow(T1.Norm(),2.));       //check this
            ksi2 = XS.Dot(T2) / (2.*std::pow(T2.Norm(),2.));       //check this

            // If natural coord is out of bound, no contact is assumed.
            if ((ksi1 < -TL1) || (ksi1 > TL2) || (ksi2 < -TL1) || (ksi2 > TL2) || GN > TL3) continue;

            // Finde exact contact point through Newton-Raphson method.
            int icount = 0;
            while (icount < 20) {
                Cutl(slave_pos[i], x0, x1, x2, ksi1, ksi2, T1, T2, XS, XC);
                A << -T1.Dot(T1),  -T2.Dot(T1),
                     -T2.Dot(T1),  -T2.Dot(T2);

                B << -XS.Dot(T1), -XS.Dot(T2);

                DXI = A.colPivHouseholderQr().solve(B);
                ksi1 += DXI.coeff(0);
                ksi2 += DXI.coeff(1);

                if(DXI.norm() < EPS)  break;
                icount++;
            }

            // Check for reference coordinate within range
            if ((ksi1 < 0.01) || (ksi1 > TL1) || (ksi2 < 0.01) || (ksi2 > TL1)) continue;

            // Normal gap function
            XN = T1.Cross(T2);
            XN /= XN.Norm();
            GN = XN.Dot(XS);
            if (GN > 0.00001)  continue;

            // std::cout << slave_pos[i] << std::endl;
            slave_pos[i] = XC*(1.+10.*std::sqrt(std::numeric_limits<double>::epsilon()));
            // std::cout << slave_pos[i] << std::endl << std::endl;
            break;
        }
    }
}


template<short DIM, short CELL_NODES>
void Deformation<DIM, CELL_NODES>::ComputeFpmForceCorrection(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm,
    const std::shared_ptr<Constitutive> &hyperelastic, const Eigen::MatrixXd &disp, double penalty, Eigen::MatrixXd &forces)
{
    // Lambda callback for multithreaded force correction computation.
    auto ForceCorrectionCallback = [](std::size_t thread_id, const CLOUDEA::ThreadLoopManager &thread_looper,
        const IMP::Voronoi<DIM> &voro, const std::shared_ptr<Constitutive> &hyperelastic,const CLOUDEA::Fpm<DIM> &fpm,
        const Eigen::MatrixXd &u, double penalty, Eigen::MatrixXd &tforces) {

        // Initialize matrices to compute internal surface force corrections.
        auto u_neighs1 = Eigen::MatrixXd{};                                 // Displacements of support domain neighs for parent node -> 1
        auto u_neighs2 = Eigen::MatrixXd{};                                 // and neighbor node -> 2
        Eigen::MatrixXd FT1 = Eigen::MatrixXd::Zero(DIM,DIM);                // Deformation gradient.
        auto FT2 = FT1;
        auto S1 = FT1;                                                       // Second Piola-Kirchhoff stress tensor.
        auto S2 = FT1;
        const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(DIM,DIM);       // Identity matrix.
        auto N1 = Eigen::VectorXd{};                                        // Shape function of the interface quadrature point.
        auto N2 = N1;
        auto traction = Eigen::VectorXd{};                                  // Traction at the subdomains interface.
        auto alpha = Eigen::MatrixXd{};                                     // Penalty matrix to impose continuity.
        auto force_correct1 = Eigen::MatrixXd{};                            // Force correction contribution of the support domain.
        auto force_correct2 = force_correct1;

        // Initialize the correction forces.
        tforces = Eigen::MatrixXd::Zero(voro.NodesNum(),DIM);

        // Correct forces for each internal facet.
        auto fid = int{0}, n1 = int{0}, n2 = int{0};
        for (auto i = thread_looper.LoopStartId(thread_id);
             i != thread_looper.LoopEndId(thread_id); ++i) {
            fid = fpm.FluxCorrector().FacetIds(i);

            // Indices of cells sharing the facet.
            n1 = voro.Facets(fid).ParentCellId();
            n2 = voro.Facets(fid).NeighCellId();
            if (n2 == -1) {
                continue;
            }

            // Get neighbors displacements in the two support domain.
            u_neighs1 = Eigen::MatrixXd::Zero(fpm.Support().InfluenceNodeIds(n1).size(), DIM);
            for (std::size_t i = 0; i != fpm.Support().InfluenceNodeIds(n1).size(); ++i) {
                u_neighs1.row(i) = u.row(fpm.Support().InfluenceNodeIds(n1)[i]);
            }

            u_neighs2 = Eigen::MatrixXd::Zero(fpm.Support().InfluenceNodeIds(n2).size(), DIM);
            for (std::size_t i = 0; i != fpm.Support().InfluenceNodeIds(n2).size(); ++i) {
                u_neighs2.row(i) = u.row(fpm.Support().InfluenceNodeIds(n2)[i]);
            }

            // Shape functions for each support domain.
            N1 = fpm.PhiGrad(n1) * (fpm.FluxCorrector().FacetCentroids(i) - voro.Nodes(n1).CopyToEigen());
            N1.coeffRef(0) = N1.coeff(0) + 1.;
            N2 = fpm.PhiGrad(n2) * (fpm.FluxCorrector().FacetCentroids(i) - voro.Nodes(n2).CopyToEigen());
            N2.coeffRef(0) = N2.coeff(0) + 1.;

            // Compute second Piola-Kirchhoff stress of the n1 cell.
            FT1 = fpm.PhiGrad(n1).transpose()*u_neighs1 + I;
            hyperelastic->ComputeStressSPK(FT1, I, hyperelastic->LameMu(n1), hyperelastic->BulkModulus(n1), S1);

            // Compute second Piola-Kirchhoff stress of the n2 cell.
            FT2 = fpm.PhiGrad(n2).transpose()*u_neighs2 + I;
            hyperelastic->ComputeStressSPK(FT2, I, hyperelastic->LameMu(n2), hyperelastic->BulkModulus(n2), S2);

            // Compute interface traction to impose continuity.
            if (DIM > 1) {
                alpha = (penalty / std::pow(fpm.FluxCorrector().FacetMeasures(i), 1./(DIM-1))) * I;
            } else {
                alpha = (penalty / fpm.FluxCorrector().FacetMeasures(i)) * I;
            }

            traction = 0.5 * (S1*FT1*fpm.FluxCorrector().FacetNormals(i) + S2*FT2*fpm.FluxCorrector().FacetNormals(i));
            traction = traction - (alpha * (u_neighs1.transpose()*N1 - u_neighs2.transpose()*N2));

            // Compute force corrections.
            force_correct1 = fpm.FluxCorrector().FacetMeasures(i) * (N1 * traction.transpose());
            force_correct2 = fpm.FluxCorrector().FacetMeasures(i) * (N2 * traction.transpose());

            // Apply interface traction
            for (std::size_t j = 0; j != fpm.Support().InfluenceNodeIds(n1).size(); ++j) {
                tforces.row(fpm.Support().InfluenceNodeIds(n1)[j]) =
                    tforces.row(fpm.Support().InfluenceNodeIds(n1)[j]) - force_correct1.row(j);
            }
            for (std::size_t j = 0; j != fpm.Support().InfluenceNodeIds(n2).size(); ++j) {
                tforces.row(fpm.Support().InfluenceNodeIds(n2)[j]) =
                    tforces.row(fpm.Support().InfluenceNodeIds(n2)[j]) + force_correct2.row(j);
            }
        }
    };

    // Initialize thread looper over the mesh cells.
    this->thread_looper_.SetLoopRanges(static_cast<std::size_t>(fpm.FluxCorrector().FacetsNum()), this->threads_num_);

    // Initialize thread forces, penalty, and domain measure.
    auto th_forces = std::vector<Eigen::MatrixXd>{};
    th_forces.resize(this->threads_num_);

    // Compute force correction in worker threads.
    auto workers = std::vector<std::thread>{};
    workers.reserve(this->threads_num_);
    for (std::size_t t = 0; t != this->threads_num_; ++t) {
        workers.emplace_back(std::thread(ForceCorrectionCallback, t, std::cref(this->thread_looper_),
        std::cref(voro), std::cref(hyperelastic), std::cref(fpm), std::cref(disp),
        penalty, std::ref(th_forces[t])));
    }
    std::for_each(workers.begin(), workers.end(), std::mem_fn(&std::thread::join));

    // Gather forces correction from the threads.
    forces = th_forces[0];
    for (std::size_t t = 1; t != this->threads_num_; ++t) {
        forces += th_forces[t];
    }
}


} // End of namespace PNT


#endif //PHYNETOUCH_PHYSICS_ELECTRICAL_TPP_