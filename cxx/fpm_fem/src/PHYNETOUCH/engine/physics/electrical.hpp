/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */

/**
   \file electrical.hpp
   \brief Electrical class header file.
   \author Konstantinos A. Mountris
   \date 28/07/2021
*/

#ifndef PHYNETOUCH_PHYSICS_ELECTRICAL_HPP_
#define PHYNETOUCH_PHYSICS_ELECTRICAL_HPP_

#include "PHYNETOUCH/engine/conditions/bound_conds.hpp"
#include "PHYNETOUCH/engine/conditions/dirichlet_bc.hpp"
#include "PHYNETOUCH/engine/conditions/initial_bc.hpp"
#include "PHYNETOUCH/engine/materials/cardiac_tissue.hpp"
#include "PHYNETOUCH/engine/materials/catheter.hpp"
#include "PHYNETOUCH/engine/utilities/logger.hpp"
#include "PHYNETOUCH/engine/utilities/measure_units.hpp"
#include "PHYNETOUCH/engine/utilities/mpi_handler.hpp"
#include "PHYNETOUCH/engine/utilities/petsc_idx_range.hpp"

#include <CLOUDEA/CLOUDEA>
#include <IMP/IMP>

#include <Eigen/Eigen>

#include <petscksp.h>

#include <vector>
#include <string>
#include <memory>
#include <numeric>
#include <utility>
#include <iterator>
#include <stdexcept>
#include <exception>
#include <thread>
#include <mutex>
#include <cmath>
#include <algorithm>


namespace PNT {

/** \addtogroup Physics \{ */

/**
 * \class Electrical
 * \brief Class implemmenting the Electrical equation for the solution of the electrical problem. *
 * \tparam DIM The dimensions of the monodomain model's domain. Supported: [1 | 2 | 3].
 * \tparam CELL_NODES The number of nodes composing the cells of the monodomain model's domain.
 */
template<short DIM, short CELL_NODES=1>
class Electrical
{

typedef Eigen::Matrix<double,-1,-1,Eigen::RowMajor> EigRowMatrixXd;

private:


    Mat conduct_mat_;                                           /**< The PetSC conduction matrix */

    Eigen::MatrixXd electric_field_;                                  /**< The electric vector field */

    std::vector<Eigen::VectorXd> stored_voltage_;                     /**< The stored values of the voltage scalar field */

    Eigen::VectorXd voltage_;                                         /**< The voltage scalar field */


public:

    /**
     * \brief Default constructor of Electrical.
     */
    Electrical();


    /**
     * \brief Default destructor of Electrical.
     *
     */
    virtual ~Electrical();


        /**
     * \brief Set the voltage to specified values.
     * \param [in] voltage The voltage values to set.
     * \return [void]
     */
    inline void SetVoltage(const Eigen::VectorXd &voltage);


    /**
     * \brief Initialize the container for storing several instants of the voltage scalar field.
     * \param [in] instants_num The number of instants to store.
     * \param [in] voltage The initializing values for the container.
     * \return [void]
     */
    inline void InitStoredVoltage(int instants_num, const Eigen::VectorXd &voltage);


    /**
     * \brief Store the voltage in the specified position of the storing container.
     * \param [in] pos The position of the container where the voltage is stored.
     * \return [void]
     */
    inline void StoreCurrentVoltageAt(std::size_t pos);


    /**
     * \brief Assemble the matrix systems for the solution of the electrical problem using Finite Element Method (FEM).
     * \param [in] mesh_tis The mesh topology of the tissue.
     * \param [in] mesh_cath The mesh topology of the catheter.
     * \param [in] mat_tis The material of the tissue.
     * \param [in] mat_cath The material of the catheter.
     * \return [void]
     */
    inline PetscErrorCode AssembleSystem(const IMP::Mesh<DIM,CELL_NODES> &mesh,
        const std::vector<double> &electrical_conductivity);


    /**
     * \brief Assemble the matrix systems for the solution of the electrical problem using Fragile Points Method (FPM).
     * \param [in] voro_tis The voronoi tesselation of the tissue.
     * \param [in] voro_cath The voronoi tesselation of the catheter.
     * \param [in] mat_tis The material of the tissue.
     * \param [in] mat_cath The material of the catheter.
     * \param [in] fpm_tis The Fragile Points Method approximant of the tissue.
     * \param [in] fpm_cath The Fragile Points Method approximant of the catheter.
     * \return [void]
     */
    inline PetscErrorCode AssembleSystem(const IMP::Voronoi<DIM> &voro,
        const std::vector<double> &electrical_conductivity, const CLOUDEA::Fpm<DIM> &fpm);


    inline PetscErrorCode UpdateSystem(const IMP::Mesh<DIM,CELL_NODES> &mesh,
        const std::vector<double> &electrical_conductivity);


    inline PetscErrorCode UpdateSystem(const IMP::Voronoi<DIM> &voro,
        const std::vector<double> &electrical_conductivity, const CLOUDEA::Fpm<DIM> &fpm);


    /**
     * \brief Apply constraints related to boundary conditions on the system matrix.
     * \param [in] dirichlet_tis The tissue Dirichlet boundary conditions.
     * \param [in] dirichlet_cath The catheter Dirichlet boundary conditions.
     * \param [in] nnum_tis The number of nodes in the tissue model.
     * \param [out] A The sparse system matrix.
     * \param [out] b The system vector.
     */
    inline PetscErrorCode ApplyConstraints(const BoundConds<1> &bc, Mat *A, Vec *b);


    /**
     * \brief
     * \param nsets
     */
    inline PetscErrorCode ComputeVoltage(const BoundConds<1> &bc);


    /**
     * \brief Compute the electric field that corresponds to the voltage scalar field.
     * \param [in] sigma
     * \return [void]
     */
    inline void ComputeElectricField(const IMP::Mesh<DIM,CELL_NODES> &mesh);


    /**
     * \brief Compute the electric field that corresponds to the voltage scalar field.
     * \param [in] sigma
     * \return [void]
     */
    inline void ComputeElectricField(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm);


    /**
     * \brief Get the voltage vector containing the voltage value at the tissue and catheter nodes.
     * \return [const Eigen::VectorXd&] The voltage vector containing the voltage value at the tissue and catheter nodes.
     */
    inline auto & Voltage() const { return this->voltage_; }


    inline auto & StoredVoltage() const { return this->stored_voltage_; }


    /**
     * \brief Get the electric field.
     * \return [const Eigen::MatrixXd&] The electric field.
     */
    inline auto & ElectricField() const { return this->electric_field_; }


    inline double TotalPower(const IMP::Mesh<DIM,CELL_NODES> &mesh,
        const std::vector<double> &electrical_conductivity);


    inline double TotalPower(const IMP::Voronoi<DIM> &voro,
        const std::vector<double> &electrical_conductivity);


protected:

    inline PetscErrorCode Assembly(const IMP::Mesh<DIM,CELL_NODES> &mesh,
        const std::vector<double> &electrical_conductivity);


    inline PetscErrorCode Assembly(const IMP::Voronoi<DIM> &voro,
        const std::vector<double> &electrical_conductivity, const CLOUDEA::Fpm<DIM> &fpm);


    inline PetscErrorCode ComputeLocalConductMat(const IMP::Mesh<DIM,CELL_NODES> &mesh,
        const std::vector<double> &conductivity, int cell_id, EigRowMatrixXd &conduct_mat) const;


    inline PetscErrorCode ComputeLocalConductMat(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm,
        const std::vector<double> &conductivity, int node_id, EigRowMatrixXd &conduct_mat) const;


    inline PetscErrorCode ComputeConductMatCorrection(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm,
        const std::vector<double> &conductivity, double penalty, int facet_pos, EigRowMatrixXd &corr_mat11,
        EigRowMatrixXd &corr_mat12, EigRowMatrixXd &corr_mat21, EigRowMatrixXd &corr_mat22) const;


    /**
     * @brief Compute the gradient of a scalar field in a given mesh.
     * The gradient is computed using the node-based Green-Gauss gradient scheme at the centroid of each cell
     * and mapped to the nodes by computing the weighted sum of the gradient in the attached cells to each node.
     * @param [in] mesh The mesh on which the gradient is computed.
     * @param [in] scalars The scalar field for which we compute the gradient.
     * @param [in] gradient The matrix of the scalar field's gradient.
     */
    inline void ComputeGradient(const IMP::Mesh<DIM,CELL_NODES> &mesh, const Eigen::VectorXd &scalars,
        Eigen::MatrixXd &gradient);


    /**
     * @brief Compute the gradient of a scalar field in a given mesh.
     * The gradient is computed using the node-based Green-Gauss gradient scheme at the centroid of each cell
     * and mapped to the nodes by computing the weighted sum of the gradient in the attached cells to each node.
     * @param [in] mesh The mesh on which the gradient is computed.
     * @param [in] scalars The scalar field for which we compute the gradient.
     * @param [in] gradient The matrix of the scalar field's gradient.
     */
    inline void ComputeGradient(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm,
        const Eigen::VectorXd &scalars, Eigen::MatrixXd &gradient);

    
    inline void PrintMatrixAscii(Mat *A, const std::string &filename);

};

/** \} End of Doxygen Groups */

} // End of namespace PNT

#endif //PHYNETOUCH_PHYSICS_ELECTRICAL_HPP_

#include "PHYNETOUCH/engine/physics/electrical.tpp"