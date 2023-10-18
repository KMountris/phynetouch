/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */

#ifndef PHYNETOUCH_PHYSICS_ELECTRICAL_TPP_
#define PHYNETOUCH_PHYSICS_ELECTRICAL_TPP_


#include "PHYNETOUCH/engine/physics/electrical.hpp"


namespace PNT {


////////////////////////////////////////////////
////!            P U B L I C              !////
//////////////////////////////////////////////


template<short DIM, short CELL_NODES>
Electrical<DIM, CELL_NODES>::Electrical() :
    conduct_mat_(PETSC_NULL), electric_field_(), stored_voltage_(), voltage_()
{}


template<short DIM, short CELL_NODES>
Electrical<DIM, CELL_NODES>::~Electrical()
{
    if (this->conduct_mat_ != PETSC_NULL)
        MatDestroy(&this->conduct_mat_);
}


template<short DIM, short CELL_NODES>
void Electrical<DIM, CELL_NODES>::SetVoltage(const Eigen::VectorXd &voltage)
{
    this->voltage_ = voltage;
}


template<short DIM, short CELL_NODES>
void Electrical<DIM, CELL_NODES>::InitStoredVoltage(int instants_num, const Eigen::VectorXd &voltage)
{
    this->stored_voltage_.clear();
    this->stored_voltage_.resize(instants_num, voltage);
}


template<short DIM, short CELL_NODES>
void Electrical<DIM, CELL_NODES>::StoreCurrentVoltageAt(std::size_t pos)
{
    this->stored_voltage_[pos] = this->voltage_;
}


template<short DIM, short CELL_NODES>
PetscErrorCode Electrical<DIM,CELL_NODES>::AssembleSystem(const IMP::Mesh<DIM,CELL_NODES> &mesh,
    const std::vector<double> &electrical_conductivity)
{
    // Number of nodes in the geometries.
    PetscInt nnum = mesh.NodesNum();

    // Create conduction matrix.
    PetscCall(MatCreate(PETSC_COMM_WORLD, &this->conduct_mat_));
    PetscCall(MatSetType(this->conduct_mat_, MATAIJ));
    PetscCall(MatSetSizes(this->conduct_mat_, PETSC_DECIDE, PETSC_DECIDE, nnum, nnum));
    PetscCall(MatSetFromOptions(this->conduct_mat_));

    // Preallocate memory for the matrix.
    PetscCall(MatSeqAIJSetPreallocation(this->conduct_mat_, 500, nullptr));
    PetscCall(MatMPIAIJSetPreallocation(this->conduct_mat_, 500, nullptr, 500, nullptr));

    // Set up the matrix.
    PetscCall(MatSetOption(this->conduct_mat_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
    PetscCall(MatSetUp(this->conduct_mat_));

    // Assemble the matrix of the system.
    this->Assembly(mesh, electrical_conductivity);

    PetscFunctionReturn(0);
}


template<short DIM, short CELL_NODES>
PetscErrorCode Electrical<DIM,CELL_NODES>::AssembleSystem(const IMP::Voronoi<DIM> &voro,
    const std::vector<double> &electrical_conductivity, const CLOUDEA::Fpm<DIM> &fpm)
{
    // Number of nodes in the geometries.
    PetscInt nnum = voro.NodesNum();

    // Create conduction matrix.
    PetscCall(MatCreate(PETSC_COMM_WORLD, &this->conduct_mat_));
    PetscCall(MatSetType(this->conduct_mat_, MATAIJ));
    PetscCall(MatSetSizes(this->conduct_mat_, PETSC_DECIDE, PETSC_DECIDE, nnum, nnum));
    PetscCall(MatSetFromOptions(this->conduct_mat_));

    // Preallocate memory for the matrix.
    PetscCall(MatSeqAIJSetPreallocation(this->conduct_mat_, 500, nullptr));
    PetscCall(MatMPIAIJSetPreallocation(this->conduct_mat_, 500, nullptr, 500, nullptr));

    // Set up the matrix.
    PetscCall(MatSetOption(this->conduct_mat_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
    PetscCall(MatSetUp(this->conduct_mat_));

    // Assemble the matrix of the system.
    this->Assembly(voro, electrical_conductivity, fpm);

    PetscFunctionReturn(0);
}


template<short DIM, short CELL_NODES>
PetscErrorCode Electrical<DIM,CELL_NODES>::UpdateSystem(const IMP::Mesh<DIM,CELL_NODES> &mesh,
    const std::vector<double> &electrical_conductivity)
{
    // Reset the matrix of the system.
    PetscCall(MatZeroEntries(this->conduct_mat_));

    // Assemble the matrix of the system.
    this->Assembly(mesh, electrical_conductivity);

    PetscFunctionReturn(0);
}


template<short DIM, short CELL_NODES>
PetscErrorCode Electrical<DIM,CELL_NODES>::UpdateSystem(const IMP::Voronoi<DIM> &voro,
    const std::vector<double> &electrical_conductivity, const CLOUDEA::Fpm<DIM> &fpm)
{
    // Reset the matrix of the system.
    PetscCall(MatZeroEntries(this->conduct_mat_));

    // Assemble the matrix of the system.
    this->Assembly(voro, electrical_conductivity, fpm);

    PetscFunctionReturn(0);
}


template<short DIM, short CELL_NODES>
PetscErrorCode Electrical<DIM, CELL_NODES>::ApplyConstraints(const BoundConds<1> &bc, Mat *A, Vec *b)
{
    // Get total number of nodes in the geometries.
    PetscInt nnum = static_cast<PetscInt>(this->voltage_.size());

    // Get range of processor.
    PetscInt i_start, i_end;
    PetscCall(MatGetOwnershipRange(*A, &i_start, &i_end));

    // Initialize boundary condition flag and values for nodes in the processor's range.
    auto bc_flag = std::vector<short>(nnum, 0);
    auto bc_vals = std::vector<double>(nnum, 0.);
    auto bc_global_ids = std::vector<PetscInt>{};
    bc_global_ids.reserve(nnum);

    // Get Dirichlet boundary condition information.
    for (const auto &dirichlet : bc.Dirichlet()) {
        for (const auto &nid : dirichlet.NodeIds()) {
            bc_flag[nid] = 1;
            bc_vals[nid] = dirichlet.Value();
            bc_global_ids.emplace_back(nid);
        }
    }
    // Get Prescribed boundary condition information.
    for (std::size_t i = 0; i != bc.Prescribed().NodeIds().size(); ++i) {
        auto nid = bc.Prescribed().NodeIds()[i];
        bc_flag[nid] = 1;
        bc_vals[nid] = bc.Prescribed().Values()[i];
        bc_global_ids.emplace_back(nid);
    }
    bc_global_ids.shrink_to_fit();
    PetscInt nbcs = static_cast<PetscInt>(bc_global_ids.size());

    // Substract scaled BC cols (rows due to symmetry) from RHS vector.
    PetscInt ncols;
    const PetscInt *cols;
    const PetscScalar *vals;
    for (PetscInt i = i_start; i != i_end; ++i) {
        if (bc_flag[i] == 1) {
            // Get the row of the BC node.
            PetscCall(MatGetRow(*A, i, &ncols, &cols, &vals));

            // Modify the row values and add to RHS;
            PetscScalar *vals_upd = new PetscScalar[ncols];
            for (PetscInt j = 0; j != ncols; ++j) {
                vals_upd[j] = -bc_vals[i]*vals[j];
            }
            PetscCall(VecSetValues(*b, ncols, cols, vals_upd, ADD_VALUES));
            PetscCall(MatRestoreRow(*A, i, &ncols, &cols, &vals));
            delete [] vals_upd;
        }
    }
    PetscCall(MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY));
    PetscCall(VecAssemblyBegin(*b));
    PetscCall(VecAssemblyEnd(*b));

    // Set the BC value to the RHS vector.
    for (PetscInt i = i_start; i != i_end; ++i) {
        if (bc_flag[i] == 1) {
            PetscCall(VecSetValues(*b, 1, &i, &bc_vals[i], INSERT_VALUES));
        }
    }
    PetscCall(VecAssemblyBegin(*b));
    PetscCall(VecAssemblyEnd(*b));

    // Zero out the row and column values of the matrix.
    PetscCall(MatZeroRowsColumns(*A, nbcs, bc_global_ids.data(), 1., 0, 0));
    PetscCall(MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY));

    PetscFunctionReturn(0);
}


template<short DIM, short CELL_NODES>
PetscErrorCode Electrical<DIM, CELL_NODES>::ComputeVoltage(const BoundConds<1> &bc)
{
    // Get total number of nodes in the geometries.
    PetscInt nnum = static_cast<PetscInt>(this->voltage_.size());

    Mat A;
    PetscCall(MatDuplicate(this->conduct_mat_, MAT_COPY_VALUES, &A));
    PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
    PetscCall(MatSetOption(A, MAT_SYMMETRIC, PETSC_TRUE));

    // Create RHS & solution vectors.
    Vec b, u;
    PetscCall(MatCreateVecs(A, &b, &u));
    PetscCall(VecSet(b, 0.0));

    // Assemble solution & RHS vectors.
    PetscCall(VecAssemblyBegin(b));
    PetscCall(VecAssemblyEnd(b));

    // Apply boundary conditions.
    PetscCall(this->ApplyConstraints(bc, &A, &b));

    // Set up Krylov space solver.
    KSP ksp;
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCall(KSPSetOperators(ksp, A, A));
    PetscCall(KSPSetType(ksp, KSPCG));
    PetscCall(KSPSetTolerances(ksp, 1.e-8, PETSC_DEFAULT, PETSC_DEFAULT, 200));
    PetscCall(KSPSetFromOptions(ksp));

    // Set up the preconditioner.
    PC pc;
    PetscCall(KSPGetPC(ksp,&pc));
    PetscCall(PCSetType(pc, PCGAMG));
    PetscCall(PCSetFromOptions(pc));

    PetscCall(KSPSolve(ksp, b, u));

    PetscCall(VecAssemblyBegin(u));
    PetscCall(VecAssemblyEnd(u));

    // Check convergence.
    PetscReal rnorm;
    PetscInt its;
    Vec r;
    PetscCall(VecDuplicate(b, &r));
    PetscCall(MatMult(A, u, r));
    PetscCall(VecAYPX(r, -1., b));
    PetscCall(VecNorm(r, NORM_2, &rnorm));

    PetscCall(KSPGetIterationNumber(ksp, &its));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Electrical -> Norm of error %g Iterations %d\n", (double)(rnorm), its));

    // Scatter the solution's part of each rank to the others.
    VecScatter ctx;
    Vec u_seq;
    PetscCall(VecScatterCreateToAll(u, &ctx, &u_seq));
    PetscCall(VecScatterBegin(ctx, u, u_seq, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterEnd(ctx, u, u_seq, INSERT_VALUES, SCATTER_FORWARD));

    // Copy the scattered solution to the voltage Eigen vector.
    PetscScalar *u_vals;
    PetscCall(VecGetArray(u_seq, &u_vals));
    for (PetscInt i = 0; i < nnum; i++) {
      this->voltage_.coeffRef(i) = u_vals[i];
    }
    PetscCall(VecRestoreArray(u_seq, &u_vals));
    PetscCall(VecScatterDestroy(&ctx));
    PetscCall(VecDestroy(&u_seq));

    // Clean up.
    PetscCall(KSPDestroy(&ksp));
    PetscCall(MatDestroy(&A));
    PetscCall(VecDestroy(&u));
    PetscCall(VecDestroy(&b));

    PetscFunctionReturn(0);
}


template<short DIM, short CELL_NODES>
double Electrical<DIM,CELL_NODES>::TotalPower(const IMP::Mesh<DIM,CELL_NODES> &mesh,
    const std::vector<double> &electrical_conductivity)
{

    // Create a finite element of the given type.
    auto elem = CLOUDEA::FemFactory<DIM>::Create(mesh.CellsShape());
    switch (mesh.CellsShape()) {
        case IMP::CellShape::tri :
            elem->SetQuadrature({3});  break;
        case IMP::CellShape::quad :
            elem->SetQuadrature({2,2});  break;
        case IMP::CellShape::tet :
            elem->SetQuadrature({4});  break;
        case IMP::CellShape::hex :
            elem->SetQuadrature({2,2,2});  break;
        default:
            auto error_msg = "Could not compute local conduction matrix. It is required a mesh with cells of type: [tri | quad | tet | hex].";
            throw std::invalid_argument(Logger::Error(error_msg));
        break;
    }

    // Compute natural derivatives and shape functions of the element.
    elem->ComputeDerivsNatural();
    elem->ComputeShapeFunctions();

    double total_power = 0.;

    // Iterate over mesh cells (tissue and catheter).
    auto cell_node_coords = std::vector<IMP::Vec<DIM, double>>{CELL_NODES};
    Eigen::RowVectorXd sigma = Eigen::RowVectorXd::Zero(CELL_NODES);
    Eigen::VectorXd grad_phi_norm_square = Eigen::VectorXd::Zero(CELL_NODES);
    for (const auto& cell : mesh.Cells()) {
        auto i = int{0};
        for (const auto &nid : cell.Connectivity()) {
            for (short d = 0; d != DIM; ++d) { cell_node_coords[i][d] = mesh.Nodes(nid)[d]; }
            sigma.coeffRef(i) = electrical_conductivity[nid];
            grad_phi_norm_square.coeffRef(i) = this->ElectricField().row(nid).squaredNorm();
            i++;
        }

        // Compute the jacobian and physical derivatives for the current cell.
        elem->ComputeJacobians(cell_node_coords);

        // Construct the local conduction cell matrix over the element's quadrature.
        auto qid = int{0};
        for (const auto &qweight : elem->Quadrature().Weights()) {
            auto q_power = (qweight*elem->DetJacobians(qid)/static_cast<double>(CELL_NODES)) * (sigma*grad_phi_norm_square);
            total_power += q_power.coeff(0);
            qid++;
        }
    }

    return total_power;

}


template<short DIM, short CELL_NODES>
double Electrical<DIM,CELL_NODES>::TotalPower(const IMP::Voronoi<DIM> &voro,
    const std::vector<double> &electrical_conductivity)
{
    auto total_power = double{0.};
    for (int nid = 0; nid != voro.NodesNum(); nid++) {
        total_power += electrical_conductivity[nid] *
            this->ElectricField().row(nid).squaredNorm() * voro.CellMeasures(nid);
    }
    return total_power;
}

////////////////////////////////////////////////
////!          P R O T E C T E D          !////
//////////////////////////////////////////////


template<short DIM, short CELL_NODES>
PetscErrorCode Electrical<DIM,CELL_NODES>::Assembly(const IMP::Mesh<DIM,CELL_NODES> &mesh,
    const std::vector<double> &electrical_conductivity)
{
    // Number of cells in the geometries.
    PetscInt cnum = mesh.CellsNum();

    // Obtain number of CPUs and current CPU rank.
    PetscMPIInt cpu_num, cpu_rank;
    PetscCallMPI(MPI_Comm_size(MPI_COMM_WORLD, &cpu_num));
    PetscCallMPI(MPI_Comm_rank(MPI_COMM_WORLD, &cpu_rank));

    // Determine cells range in cpu rank.
    PetscIdxRange local_ids;
    local_ids.Compute(cnum, cpu_rank, cpu_num);

    // Assemble local cell contributions.
    EigRowMatrixXd Ke;
    PetscInt idx[CELL_NODES];
    PetscInt it;
    for (PetscInt cid = local_ids.IdStart(); cid != local_ids.IdEnd(); ++cid) {
        // Compute conductivity matrix associated with cell i.
        it = 0;
        PetscCall(this->ComputeLocalConductMat(mesh, electrical_conductivity, cid, Ke));
        for (const auto &nid : mesh.Cells(cid).Connectivity()) {
            idx[it++] = nid;
        }

        // Flatten in row major order the local matrix.
        PetscCall(MatSetValues(this->conduct_mat_, CELL_NODES, idx, CELL_NODES, idx, Ke.data(), ADD_VALUES));
    }

    // Assemble the final matrix.
    PetscCall(MatAssemblyBegin(this->conduct_mat_, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(this->conduct_mat_, MAT_FINAL_ASSEMBLY));

    PetscFunctionReturn(0);
}


template<short DIM, short CELL_NODES>
PetscErrorCode Electrical<DIM,CELL_NODES>::Assembly(const IMP::Voronoi<DIM> &voro,
    const std::vector<double> &electrical_conductivity, const CLOUDEA::Fpm<DIM> &fpm)
{
    // Number of nodes in the geometries.
    PetscInt nnum = voro.NodesNum();

    // Obtain number of CPUs and current CPU rank.
    PetscMPIInt cpu_num, cpu_rank;
    PetscCallMPI(MPI_Comm_size(MPI_COMM_WORLD, &cpu_num));
    PetscCallMPI(MPI_Comm_rank(MPI_COMM_WORLD, &cpu_rank));

    // Compute FPM correction penalty.
    auto penal = double{0.};
    auto vol = double{0.};
    for (auto i = 0; i != voro.NodesNum(); ++i) {
        penal += electrical_conductivity[i]*voro.CellMeasures(i);
        vol += voro.CellMeasures(i);
    }
    penal = fpm.Penalty() * (penal/vol);

    // Range of indices of either nodes or internal facets, local to the cpu rank.
    PetscIdxRange local_ids;
    local_ids.Compute(nnum, cpu_rank, cpu_num);

    // Assemble nodal contributions.
    EigRowMatrixXd Ke, K11, K12, K21, K22;
    PetscInt n2, neighs1_num, neighs2_num;
    auto neighs1 = std::vector<int>{};
    auto neighs2 = std::vector<int>{};
    // PetscInt *idx_n1, *idx_n2;
    for (PetscInt n1 = local_ids.IdStart(); n1 != local_ids.IdEnd(); ++n1) {
        // Compute conductivity matrix.
        PetscCall(this->ComputeLocalConductMat(voro, fpm, electrical_conductivity, n1, Ke));

        neighs1 = fpm.Support().InfluenceNodeIds(n1);
        neighs1_num = static_cast<PetscInt>(neighs1.size());
        // idx_n1 = new PetscInt[neighs1_num];
        // for (PetscInt i = 0; i != neighs1_num; ++i) {
        //     idx_n1[i] = neighs1[i];
        // }

        // Ke
        PetscCall(MatSetValues(this->conduct_mat_, neighs1_num, neighs1.data(),
            neighs1_num, neighs1.data(), Ke.data(), ADD_VALUES));

        // Compute facets correction.
        for (const auto& fid : voro.FacetsOnNodes(n1)) {

            // Correction if facets has neighbor cell (it is internal).
            if (voro.Facets(fid).NeighCellId() != -1) {
                n2 = voro.Facets(fid).NeighCellId();

                PetscCall(this->ComputeConductMatCorrection(voro, fpm, electrical_conductivity, penal, fid, K11, K12, K21, K22));
                neighs2 = fpm.Support().InfluenceNodeIds(n2);
                neighs2_num = static_cast<PetscInt>(neighs2.size());
                // idx_n2 = new PetscInt[neighs2_num];
                // for (PetscInt j = 0; j != neighs2_num; ++j) {
                //     idx_n2[j] = neighs2[j];
                // }

                // K11
                PetscCall(MatSetValues(this->conduct_mat_, neighs1_num, neighs1.data(),
                    neighs1_num, neighs1.data(), K11.data(), ADD_VALUES));
                // K12
                PetscCall(MatSetValues(this->conduct_mat_, neighs1_num, neighs1.data(),
                    neighs2_num, neighs2.data(), K12.data(), ADD_VALUES));
                // K21
                PetscCall(MatSetValues(this->conduct_mat_, neighs2_num, neighs2.data(),
                    neighs1_num, neighs1.data(), K21.data(), ADD_VALUES));
                // K22
                PetscCall(MatSetValues(this->conduct_mat_, neighs2_num, neighs2.data(),
                    neighs2_num, neighs2.data(), K22.data(), ADD_VALUES));

                // delete [] idx_n2;
            }
        }

        // delete [] idx_n1;
    }

    // Assemble the final matrix.
    PetscCall(MatAssemblyBegin(this->conduct_mat_, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(this->conduct_mat_, MAT_FINAL_ASSEMBLY));

    PetscFunctionReturn(0);
}


template<short DIM, short CELL_NODES>
PetscErrorCode Electrical<DIM,CELL_NODES>::ComputeLocalConductMat(const IMP::Mesh<DIM,CELL_NODES> &mesh,
    const std::vector<double> &conductivity, int cell_id, EigRowMatrixXd &conduct_mat) const
{
    // Create a finite element of the given type.
    auto elem = CLOUDEA::FemFactory<DIM>::Create(mesh.CellsShape());
    switch (mesh.CellsShape()) {
        case IMP::CellShape::tri :
            elem->SetQuadrature({3});  break;
        case IMP::CellShape::quad :
            elem->SetQuadrature({2,2});  break;
        case IMP::CellShape::tet :
            elem->SetQuadrature({4});  break;
        case IMP::CellShape::hex :
            elem->SetQuadrature({2,2,2});  break;
        default:
            auto error_msg = "Could not compute local conduction matrix. It is required a mesh with cells of type: [tri | quad | tet | hex].";
            throw std::invalid_argument(Logger::Error(error_msg));
        break;
    }

    // Compute natural derivatives and shape functions of the element.
    elem->ComputeDerivsNatural();
    elem->ComputeShapeFunctions();

    // Iterate over mesh cells (tissue and catheter).
    auto cell_node_coords = std::vector<IMP::Vec<DIM, double>>{CELL_NODES};
    auto mean_conductivity = double{0.};
    auto i = int{0};
    for (const auto &nid : mesh.Cells(cell_id).Connectivity()) {
        for (short d = 0; d != DIM; ++d) { cell_node_coords[i][d] = mesh.Nodes(nid)[d]; }
        mean_conductivity += conductivity[nid];
        i++;
    }
    mean_conductivity /= static_cast<double>(CELL_NODES);

    // Compute the jacobian and physical derivatives for the current cell.
    elem->ComputeJacobians(cell_node_coords);
    elem->ComputeDerivs();

    // Construct the local conduction cell matrix over the element's quadrature.
    auto qid = int{0};
    conduct_mat = EigRowMatrixXd::Zero(CELL_NODES, CELL_NODES);
    for (const auto &qweight : elem->Quadrature().Weights()) {
        conduct_mat += qweight*elem->DetJacobians(qid) * (elem->Derivs(qid)*mean_conductivity*elem->Derivs(qid).transpose());
        qid++;
    }
    PetscFunctionReturn(0);
}


template<short DIM, short CELL_NODES>
PetscErrorCode Electrical<DIM, CELL_NODES>::ComputeLocalConductMat(const IMP::Voronoi<DIM> &voro,
    const CLOUDEA::Fpm<DIM> &fpm, const std::vector<double> &conductivity, int node_id, EigRowMatrixXd &conduct_mat) const
{
    // Compute local conductivity matrix for node_id.
    conduct_mat = fpm.PhiGrad(node_id) * conductivity[node_id] * fpm.PhiGrad(node_id).transpose() * voro.CellMeasures(node_id);
    PetscFunctionReturn(0);
}


template<short DIM, short CELL_NODES>
PetscErrorCode Electrical<DIM, CELL_NODES>::ComputeConductMatCorrection(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm,
    const std::vector<double> &conductivity, double penalty, int fid, EigRowMatrixXd &corr_mat11,
    EigRowMatrixXd &corr_mat12, EigRowMatrixXd &corr_mat21, EigRowMatrixXd &corr_mat22) const
{
    // Indices of cells sharing the facet.
    if (fid != fpm.FluxCorrector().FacetIds(fid)) {
        std::string err_msg = "Could not compute bioheat mat correction. Facet id mismatch.";
        throw std::runtime_error(Logger::Error(err_msg));
    }

    auto n1 = voro.Facets(fid).ParentCellId();
    auto n2 = voro.Facets(fid).NeighCellId();

    // Boundary dependent parameter with unit of length
    auto he = std::sqrt(voro.Nodes(n2).Distance2(voro.Nodes(n1)));

    // Set facet normals pointing towards the parent and neighbor cell.
    Eigen::VectorXd ne1 = fpm.FluxCorrector().FacetNormals(fid);
    Eigen::VectorXd ne2 = -ne1;
    Eigen::RowVectorXd ne1_t = ne1.transpose();
    Eigen::RowVectorXd ne2_t = ne2.transpose();


    // Get fpm gradients for the cells.
    Eigen::MatrixXd B1_t = fpm.PhiGrad(n1).transpose();
    Eigen::MatrixXd B2_t = fpm.PhiGrad(n2).transpose();

    // Compute shape functions for parent and neighbor cells.
    Eigen::VectorXd aux1 = fpm.FluxCorrector().FacetCentroids(fid);
    Eigen::VectorXd aux2 = fpm.FluxCorrector().FacetCentroids(fid);
    for (short d = 0; d != DIM; ++d) {
        aux1.coeffRef(d) = aux1.coeff(d) - voro.Nodes(n1)[d];
        aux2.coeffRef(d) = aux2.coeff(d) - voro.Nodes(n2)[d];
    }
    Eigen::VectorXd N1 = fpm.PhiGrad(n1)*aux1;
    Eigen::VectorXd N2 = fpm.PhiGrad(n2)*aux2;

    N1.coeffRef(0) = N1.coeff(0) + 1.;
    N2.coeffRef(0) = N2.coeff(0) + 1.;

    Eigen::RowVectorXd N1_t = N1.transpose();
    Eigen::RowVectorXd N2_t = N2.transpose();

    // Compute facet correction matrices.
    corr_mat11 = fpm.FluxCorrector().FacetMeasures(fid) *
        (-0.5*(N1*ne1_t*conductivity[n1]*B1_t + fpm.PhiGrad(n1)*conductivity[n1]*ne1*N1_t) + (penalty/he)*(N1*N1_t));

    corr_mat12 = fpm.FluxCorrector().FacetMeasures(fid) *
        (-0.5*(N1*ne1_t*conductivity[n2]*B2_t + fpm.PhiGrad(n1)*conductivity[n1]*ne2*N2_t) - (penalty/he)*(N1*N2_t));

    corr_mat21 = fpm.FluxCorrector().FacetMeasures(fid) *
        (-0.5*(N2*ne2_t*conductivity[n1]*B1_t + fpm.PhiGrad(n2)*conductivity[n2]*ne1*N1_t) - (penalty/he)*(N2*N1_t));

    corr_mat22 = fpm.FluxCorrector().FacetMeasures(fid) *
        (-0.5*(N2*ne2_t*conductivity[n2]*B2_t + fpm.PhiGrad(n2)*conductivity[n2]*ne2*N2_t) + (penalty/he)*(N2*N2_t));

    PetscFunctionReturn(0);
}


template<short DIM, short CELL_NODES>
void Electrical<DIM, CELL_NODES>::ComputeElectricField(const IMP::Mesh<DIM,CELL_NODES> &mesh)
{
    // Compute voltage gradient.
    Eigen::MatrixXd voltage_grad;
    this->ComputeGradient(mesh, this->voltage_, voltage_grad);

    // Compute the electric field.
    this->electric_field_ = -voltage_grad;
}


template<short DIM, short CELL_NODES>
void Electrical<DIM, CELL_NODES>::ComputeElectricField(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm)
{
    // Compute voltage gradient.
    Eigen::MatrixXd voltage_grad;
    this->ComputeGradient(voro, fpm, this->voltage_, voltage_grad);

    // Compute the total electric field.
    this->electric_field_ = -voltage_grad;
}


template<short DIM, short CELL_NODES>
void Electrical<DIM, CELL_NODES>::ComputeGradient(const IMP::Mesh<DIM,CELL_NODES> &mesh, const Eigen::VectorXd &nodal_scalars, Eigen::MatrixXd &nodal_gradients)
{
    // Set cell faces connectivity.
    auto faces = Eigen::MatrixXd{};
    if (DIM == 3 && CELL_NODES == 4) {
        faces.resize(4,3);
        faces << 0, 2, 1,
                 0, 3, 2,
                 1, 2, 3,
                 0, 1, 3;
    } else {
        auto error_msg = "Could not compute gradient. Currently supports only tetrahedral meshes.";
        throw std::runtime_error(Logger::Error(error_msg));
    }
    // } else if (CELL_NODES == 8) {
    //     faces << 0, 3, 2, 1,
    //              4, 5, 6, 7,
    //              0, 1, 5, 4,
    //              2, 3, 7, 6,
    //              1, 2, 6, 5,
    //              0, 4, 7, 3;
    // }

    // Tetrahedron volume lambda.
    auto TetVolume = [](IMP::Vec<DIM,double> v0, IMP::Vec<DIM,double> v1,
            IMP::Vec<DIM,double> v2, IMP::Vec<DIM,double> v3)
    {
        auto a = v0 - v3;
        auto b = v1 - v3;
        auto c = v2 - v3;
        auto det = a[0] * (b[1]*c[2] - c[1]*b[2]) -
                   b[0] * (a[1]*c[2] - c[1]*a[2]) +
                   c[0] * (a[1]*b[2] - b[1]*a[2]);
        return std::abs(det)/6.;
    };

    // // Hexahedron volume lambda.
    // double HexVolume = [](IMP::Vec<DIM,double> v0, IMP::Vec<DIM,double> v1,
    //                       IMP::Vec<DIM,double> v2, IMP::Vec<DIM,double> v3,
    //                       IMP::Vec<DIM,double> v4, IMP::Vec<DIM,double> v5,
    //                       IMP::Vec<DIM,double> v6, IMP::Vec<DIM,double> v7)
    // { };

    // Compute gradient at cell centroids.
    auto cell_gradients = Eigen::MatrixXd{Eigen::MatrixXd::Zero(mesh.CellsNum(),DIM)};
    auto cell_grad = Eigen::VectorXd{};
    auto cell_centroids = std::vector<IMP::Vec<DIM,double>>(mesh.CellsNum());
    auto attached_cells_to_nodes = std::vector<std::vector<int>>(mesh.NodesNum());
    auto face_centroid = IMP::Vec<DIM,double>{};
    auto face_normal = IMP::Vec<DIM,double>{};
    auto cell_vol = double{0.};
    auto face_scalar = double{0.};
    auto cid = int{0};
    for (const auto &cell : mesh.Cells()) {

        // Compute the centroid of the cell.
        cell_centroids[cid].SetZero();
        for (const auto &nid : cell.Connectivity()) {
            // Add the cell to the container of attached cell to each of its nodes.
            attached_cells_to_nodes[nid].emplace_back(cid);
            cell_centroids[cid] += mesh.Nodes(nid);
        }
        cell_centroids[cid] /= static_cast<double>(CELL_NODES);

        // Compute the volume of the cell.
        cell_vol = TetVolume(mesh.Nodes(cell.N(0)), mesh.Nodes(cell.N(1)),
                mesh.Nodes(cell.N(2)), mesh.Nodes(cell.N(3)));

        // Compute the gradient at the cell's centroid.
        cell_grad = Eigen::VectorXd::Zero(DIM);
        for (int f=0; f!=faces.rows(); ++f) {

            // Compute face normal vector with Newell's method.
            face_normal.SetZero();
            face_centroid.SetZero();
            face_scalar = 0.;
            for (int i=0; i!=faces.cols(); ++i) {
                auto v = mesh.Nodes( cell.N(faces(f,i)) );
                auto w = mesh.Nodes( cell.N(faces(f,(i+1)%faces.cols())) );

                // Compute face normal by cross product of sequential nodes
                face_normal[0] += v[1]*w[2] - v[2]*w[1];
                face_normal[1] += v[2]*w[0] - v[0]*w[2];
                face_normal[2] += v[0]*w[1] - v[1]*w[0];

                // Compute face centroid.
                face_centroid += mesh.Nodes( cell.N(faces(f,i)) );

                // Compute scalar at face centroid.
                face_scalar += nodal_scalars.coeffRef( cell.N(faces(f,i)) );
            }

            // Normalize face centroid and face scalar.
            face_centroid /= static_cast<double>(faces.cols());
            face_scalar /= static_cast<double>(faces.cols());

            // Check face normal direction.
            if (face_normal.CwiseMul(face_centroid-cell_centroids[cid]).Sum() < 0.)  face_normal = -face_normal;

            // Add the gradient contribution from the face centroid.
            cell_grad += face_scalar*0.5*face_normal.CopyToEigen();
        }
        cell_gradients.row(cid++) = cell_grad / cell_vol;

    } // End of Iterate over the mesh cells.

    // Compute nodal gradients.
    auto nid = int{0};
    auto weight = double{0.};
    auto weights_sum = double{0.};
    auto nodal_grad = Eigen::RowVectorXd{Eigen::RowVectorXd::Zero(DIM)};
    nodal_gradients = Eigen::MatrixXd::Zero(nodal_scalars.rows(),DIM);
    for (const auto &attached_cells : attached_cells_to_nodes) {
        weights_sum = 0.;
        nodal_grad = Eigen::RowVectorXd::Zero(DIM);
        for (const auto &cid : attached_cells) {
            weight = 1. / std::sqrt(mesh.Nodes(nid).Distance2(cell_centroids[cid]));
            weights_sum += weight;
            nodal_grad += weight*cell_gradients.row(cid);
        }
        nodal_gradients.row(nid++) = nodal_grad/weights_sum;
    }
}


template<short DIM, short CELL_NODES>
void Electrical<DIM, CELL_NODES>::ComputeGradient(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm, const Eigen::VectorXd &nodal_scalars, Eigen::MatrixXd &nodal_gradients)
{
    // Compute nodal gradients.
    auto neighs = std::vector<int>{};
    auto neighs_scalar = Eigen::RowVectorXd{};
    nodal_gradients = Eigen::MatrixXd::Zero(nodal_scalars.rows(),DIM);
    for (int nid = 0; nid != voro.NodesNum(); ++nid) {

        // Get neighbor nodes.
        neighs = fpm.Support().InfluenceNodeIds(nid);

        // Collect the nodal scalar values of the neighbor nodes.
        neighs_scalar = Eigen::RowVectorXd::Zero(neighs.size());
        for (std::size_t i = 0; i != neighs.size(); ++i) {
            neighs_scalar.coeffRef(i) = nodal_scalars.coeff(neighs[i]);
        }

        // Compute nodal gradient.
        nodal_gradients.row(nid) = neighs_scalar*fpm.PhiGrad(nid);
    }
}


template<short DIM, short CELL_NODES>
void Electrical<DIM, CELL_NODES>::PrintMatrixAscii(Mat *A, const std::string &filename)
{
    PetscViewer viewer;
    PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
    PetscViewerSetType(viewer, PETSCVIEWERASCII);
    PetscViewerFileSetName(viewer, filename.c_str());
    PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);

    PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
    MatView(*A, viewer);
    PetscViewerPopFormat(viewer);

    PetscViewerDestroy(&viewer);
}


} // End of namespace PNT


#endif //PHYNETOUCH_PHYSICS_ELECTRICAL_TPP_