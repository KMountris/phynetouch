/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */

/**
   \file bioheat.hpp
   \brief Bioheat class header file.
   \author Konstantinos A. Mountris
   \date 28/07/2021
*/

#ifndef PHYNETOUCH_PHYSICS_BIOHEAT_HPP_
#define PHYNETOUCH_PHYSICS_BIOHEAT_HPP_

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
#include <IMP/Tesselations>
#include <IMP/Topology>
#include <IMP/Utilities>

#include <Eigen/Eigen>

#include <petscksp.h>

#include <vector>
#include <string>
#include <memory>
#include <numeric>
#include <limits>
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
 * \class Bioheat
 * \brief Class implemmenting the solution to the bioheat problem using FEM and meshfree methods.
 * \tparam DIM The dimensions of the model's. Supported: [1 | 2 | 3].
 * \tparam CELL_NODES The number of nodes composing the cells of the model's mesh.
 */
template<short DIM, short CELL_NODES=1>
class Bioheat
{

typedef Eigen::Matrix<double,-1,-1,Eigen::RowMajor> EigRowMatrixXd;

private:

    Mat conduct_mat_;

    Mat enthalpy_mat_;

    Eigen::MatrixXd thermal_field_;                                    /**< The thermal vector field */

    std::vector<Eigen::VectorXd> stored_temperature_;                  /**< Stored values of the temperature scalar field */

    Eigen::VectorXd temperature_;                                      /**< The temperature scalar field */

    Vec body_heat_;

    double dt_;


public:

    /**
     * \brief Default constructor of Laplacian.
     */
    Bioheat();


    /**
     * \brief Default destructor of Laplacian.
     *
     */
    virtual ~Bioheat();


    /**
     * \brief Set the time step of the integration.
     * \param [in] dt The time step of the integration.
     * \return [void]
     */
    inline void SetDt(double dt);


    /**
     * \brief Assemble the conductivity matrix using FEM.
     * \param [in] mesh The mesh discretization of the cardiac geometry.
     */
    inline PetscErrorCode AssembleSystem(const IMP::Mesh<DIM,CELL_NODES> &mesh,
        const std::vector<double> &thermal_conductivity, const std::vector<double> &enthalpy,
        const Eigen::VectorXd &rf_heat);


    inline PetscErrorCode AssembleSystem(const IMP::Voronoi<DIM> &voro,
        const CLOUDEA::Fpm<DIM> &fpm, const std::vector<double> &thermal_conductivity,
        const std::vector<double> &enthalpy, const Eigen::VectorXd &rf_heat);


    inline PetscErrorCode UpdateSystem(const IMP::Mesh<DIM,CELL_NODES> &mesh,
        const std::vector<double> &thermal_conductivity, const std::vector<double> &enthalpy,
        const Eigen::VectorXd &rf_heat);


    inline PetscErrorCode UpdateSystem(const IMP::Voronoi<DIM> &voro,
        const CLOUDEA::Fpm<DIM> &fpm, const std::vector<double> &thermal_conductivity,
        const std::vector<double> &enthalpy, const Eigen::VectorXd &rf_heat);


    inline PetscErrorCode UpdateOnlyHeatSource(const IMP::Mesh<DIM,CELL_NODES> &mesh,
        const Eigen::VectorXd &rf_heat);


    inline PetscErrorCode UpdateOnlyHeatSource(const IMP::Voronoi<DIM> &voro,
        const CLOUDEA::Fpm<DIM> &fpm, const Eigen::VectorXd &rf_heat);


    /**
     * \brief Apply constraints related to boundary conditions on the system matrix.
     * \param [in] dirichlet_tis The tissue Dirichlet boundary conditions.
     * \param [in] dirichlet_cath The catheter Dirichlet boundary conditions.
     * \param [in] nnum_tis The number of nodes in the tissue model.
     * \param [out] A The sparse system matrix.
     * \param [out] b The system vector.
     */
    inline PetscErrorCode ApplyConstraints(const std::vector<DirichletBc<1>> &dirichlet, Mat *A, Vec *b);


    /**
     * \brief
     * \param nsets
     */
    inline PetscErrorCode ComputeTemperature(const std::vector<DirichletBc<1>> &dirichlet);


    /**
     * \brief Set the temperature to specified values.
     * \param [in] temperature The temperature values to set.
     * \return [void]
     */
    inline void SetTemperature(const Eigen::VectorXd &temperature);


    /**
     * \brief Initialize the container for storing several instants of the temperature scalar field.
     * \param [in] instants_num The number of instants to store.
     * \param [in] temperature The initializing values for the container.
     * \return [void]
     */
    inline void InitStoredTemperature(int instants_num, const Eigen::VectorXd &temperature);


    /**
     * \brief Store the temperature in the specified position of the storing container.
     * \param [in] pos The position of the container where the temperature is stored.
     * \return [void]
     */
    inline void StoreCurrentTemperatureAt(int pos);


    /**
     * @brief
     * @return [double]
     */
    inline auto Dt() const { return this->dt_; }


    /**
     * \brief Get the temperature scalar field at the tissue and catheter nodes.
     * \return [const Eigen::VectorXd&] The temperature scalar field at the tissue and catheter nodes.
     */
    inline auto & Temperature() const { return this->temperature_; }


    /**
     * @brief
     *
     * @return const std::vector<Eigen::VectorXd>&
     */
    inline auto & StoredTemperature() const { return this->stored_temperature_; }


    /**
     * \brief Get the temperature stored in the container at a given position.
     * \param [in] id The index of the position to retrieve the temperature.
     * \return [const Eigen::VectorXd&] The temperature stored in the container at a given position.
     */
    inline auto & StoredTemperature(std::size_t id) const { return this->stored_temperature_[id]; }


protected:


    inline PetscErrorCode Assembly(const IMP::Mesh<DIM,CELL_NODES> &mesh,
        const std::vector<double> &thermal_conductivity, const std::vector<double> &enthalpy,
        const Eigen::VectorXd &rf_heat);


    inline PetscErrorCode Assembly(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm,
        const std::vector<double> &thermal_conductivity, const std::vector<double> &enthalpy,
        const Eigen::VectorXd &rf_heat);


    inline PetscErrorCode AssemblyHeatSource(const IMP::Mesh<DIM,CELL_NODES> &mesh, const Eigen::VectorXd &rf_heat);


    inline PetscErrorCode AssemblyHeatSource(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm,
        const Eigen::VectorXd &rf_heat);


    inline PetscErrorCode ComputeLocalMats(const IMP::Mesh<DIM,CELL_NODES> &mesh,
        const std::vector<double> &thermal_conductivity, const std::vector<double> &enthalpy,
        const Eigen::VectorXd rf_heat, int cell_id, EigRowMatrixXd &conduct_mat,
        EigRowMatrixXd &enthalpy_mat, Eigen::VectorXd &body_heat_vec) const;


    inline PetscErrorCode ComputeLocalMats(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm,
        const std::vector<double> &thermal_conductivity, const std::vector<double> &enthalpy,
        const Eigen::VectorXd rf_heat, int node_id, EigRowMatrixXd &conduct_mat,
        EigRowMatrixXd &enthalpy_mat, Eigen::VectorXd &body_heat_vec) const;


    inline PetscErrorCode ComputeLocalHeatVec(const IMP::Mesh<DIM,CELL_NODES> &mesh,
        const Eigen::VectorXd rf_heat, int cell_id, Eigen::VectorXd &body_heat_vec) const;


    inline PetscErrorCode ComputeLocalHeatVec(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm,
        const Eigen::VectorXd rf_heat, int node_id, Eigen::VectorXd &body_heat_vec) const;


    inline PetscErrorCode ComputeConductMatCorrection(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm,
        const std::vector<double> &thermal_conductivity, double penalty, int facet_id, EigRowMatrixXd &corr_mat11,
        EigRowMatrixXd &corr_mat12, EigRowMatrixXd &corr_mat21, EigRowMatrixXd &corr_mat22) const;


    /**
     * @brief Compute the gradient of a scalar field in a given mesh.
     *
     * The gradient is computed using the node-based Green-Gauss gradient scheme at the centroid of each cell
     * and mapped to the nodes by computing the weighted sum of the gradient in the attached cells to each node.
     *
     * @param [in] mesh The mesh on which the gradient is computed.
     * @param [in] scalars The scalar field for which we compute the gradient.
     * @param [in] gradient The matrix of the scalar field's gradient.
     */
    inline void ComputeGradient(const IMP::Mesh<DIM,CELL_NODES> &mesh, const Eigen::VectorXd &scalars, Eigen::MatrixXd &gradient);


    /**
     * @brief Compute the gradient of a scalar field in a given mesh.
     *
     * The gradient is computed using the node-based Green-Gauss gradient scheme at the centroid of each cell
     * and mapped to the nodes by computing the weighted sum of the gradient in the attached cells to each node.
     *
     * @param [in] mesh The mesh on which the gradient is computed.
     * @param [in] scalars The scalar field for which we compute the gradient.
     * @param [in] gradient The matrix of the scalar field's gradient.
     */
    inline void ComputeGradient(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm, const Eigen::VectorXd &scalars, Eigen::MatrixXd &gradient);


    inline void PrintMatrixAscii(Mat *A, const std::string &filename);
};

/** \} End of Doxygen Groups */

} // End of namespace PNT

#endif //PHYNETOUCH_PHYSICS_BIOHEAT_HPP_

#include "PHYNETOUCH/engine/physics/bioheat.tpp"