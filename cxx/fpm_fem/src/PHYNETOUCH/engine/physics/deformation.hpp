/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */

/**
   \file deformation.hpp
   \brief Deformation class header file.
   \author Konstantinos A. Mountris
   \date 01/11/2021
*/

#ifndef PHYNETOUCH_PHYSICS_DEFORMATION_HPP_
#define PHYNETOUCH_PHYSICS_DEFORMATION_HPP_

#include "PHYNETOUCH/engine/physics/contact_handle.hpp"
#include "PHYNETOUCH/engine/conditions/bound_conds.hpp"
#include "PHYNETOUCH/engine/conditions/dirichlet_bc.hpp"
#include "PHYNETOUCH/engine/conditions/initial_bc.hpp"
#include "PHYNETOUCH/engine/materials/constitutive.hpp"
#include "PHYNETOUCH/engine/utilities/logger.hpp"
#include "PHYNETOUCH/engine/utilities/measure_units.hpp"

#include <CLOUDEA/CLOUDEA>
#include <IMP/IMP>

#include <Eigen/Dense>
#include <termcolor/termcolor.hpp>

#include <vector>
#include <string>
#include <memory>
#include <numeric>
#include <limits>
#include <utility>
#include <iterator>
#include <thread>
#include <mutex>
#include <stdexcept>
#include <exception>
#include <cmath>


namespace PNT {

/** \addtogroup Physics \{ */

/**
 * \class Deformation
 * \brief Class implemmenting the solution to the deformation problem using FEM and meshfree methods.
 * \tparam DIM The dimensions of the model's. Supported: [1 | 2 | 3].
 * \tparam CELL_NODES The number of nodes composing the cells of the model's mesh.
 */
template<short DIM, short CELL_NODES=1>
class Deformation
{

private:

    std::vector<Eigen::MatrixXd> disp_history_;     /**< Container of the displacements at each time step */

    Eigen::MatrixXd lumped_mass_;                   /**< The lumped mass vector */

    std::vector<double> cell_volumes_;              /**< The volume of each cell of the model's geometry */

    std::vector<double> nodal_volumes_;             /**< The volume associated with each node of the model's geometry */

    std::vector<std::vector<int>> attached_cells_to_nodes_; 

    double sim_time_;                               /**< The simulation time */

    double dt_;                                     /**< The simulation time step */

    double dt_crit_;                                /**< The critical time step for explicit integration */

    CLOUDEA::ThreadLoopManager thread_looper_;      /**< The managing object for the multithreaded loop execution */

    std::size_t threads_num_;                       /**< The number of threads for parallel execution of the monodomain model */

    std::mutex tmutex_;                             /**< The mutex for the multithreaded execution */


public:

    /**
     * \brief Default constructor of Deformation.
     */
    Deformation();


    /**
     * \brief Default destructor of Deformation.
     *
     */
    virtual ~Deformation();


    /**
     * \brief Assemble the conductivity matrix using FEM.
     * \param [in] mesh The mesh discretization of the cardiac geometry.
     */
    inline void FormLumpedMass(const IMP::Mesh<DIM,CELL_NODES> &tis_mesh, const IMP::Mesh<DIM,CELL_NODES> &cath_mesh,
        const CLOUDEA::FemMats<DIM,CELL_NODES> &tis_fem, const CLOUDEA::FemMats<DIM,CELL_NODES> &cath_fem,
        const std::shared_ptr<Constitutive> &tis_elastic, const std::shared_ptr<Constitutive> &cath_elastic, bool mass_scaling);


    inline void FormLumpedMass(const IMP::Voronoi<DIM> &tis_voro, const IMP::Voronoi<DIM> &cath_voro,
        const CLOUDEA::Fpm<DIM> &tis_fpm, const CLOUDEA::Fpm<DIM> &cath_fpm,
        const std::shared_ptr<Constitutive> &tis_elastic, const std::shared_ptr<Constitutive> &cath_elastic, bool mass_scaling);


    /**
     * \brief Set the simulation time for the deformation problem.
     * \param [in] sim_time The simulation time.
     * \return [void]
     */
    inline void SetSimulationTime(double sim_time);


    /**
     * \brief Set the simulation time step for the deformation problem.
     * If explicit integration is chosen and the provided time step is larger than the critical one,
     * then the critical time step is used instead.
     * \param [in] dt The simulation time step.
     * \return [void]
     */
    inline void SetDt(double dt);


    inline void ComputeCellVolumes(const IMP::Mesh<DIM,CELL_NODES> &tis_mesh, const IMP::Mesh<DIM,CELL_NODES> &cath_mesh,
        const CLOUDEA::FemMats<DIM,CELL_NODES> &tis_fem, const CLOUDEA::FemMats<DIM,CELL_NODES> &cath_fem);


    inline void ComputeNodalVolumes(const IMP::Mesh<DIM,CELL_NODES> &tis_mesh, const IMP::Mesh<DIM,CELL_NODES> &cath_mesh);


    inline void UpdateInternalForces(const IMP::Mesh<DIM,CELL_NODES> &tis_mesh, const IMP::Mesh<DIM,CELL_NODES> &cath_mesh,
        const CLOUDEA::FemMats<DIM,CELL_NODES> &tis_fem, const CLOUDEA::FemMats<DIM,CELL_NODES> &cath_fem,
        const std::shared_ptr<Constitutive> &tis_elastic, const std::shared_ptr<Constitutive> &cath_elastic,
        const Eigen::MatrixXd &disp, Eigen::MatrixXd &forces);


    inline void UpdateInternalForces(const IMP::Voronoi<DIM> &tis_voro, const IMP::Voronoi<DIM> &cath_voro,
        const CLOUDEA::Fpm<DIM> &tis_fpm, const CLOUDEA::Fpm<DIM> &cath_fpm,
        const std::shared_ptr<Constitutive> &tis_elastic, const std::shared_ptr<Constitutive> &cath_elastic,
        const Eigen::MatrixXd &disp, Eigen::MatrixXd &forces);


    inline void UpdateExternalForces(const IMP::Mesh<DIM,CELL_NODES> &tis_mesh, const IMP::Mesh<DIM,CELL_NODES> &cath_mesh,
        const CLOUDEA::FemMats<DIM,CELL_NODES> &tis_fem, const CLOUDEA::FemMats<DIM,CELL_NODES> &cath_fem,
        const BoundConds<DIM> &tis_deform_bc, const BoundConds<DIM> &cath_deform_bc, int step, Eigen::MatrixXd &forces);


    inline void UpdateExternalForces(const IMP::Voronoi<DIM> &tis_voro, const IMP::Voronoi<DIM> &cath_voro,
        const CLOUDEA::Fpm<DIM> &tis_fpm, const CLOUDEA::Fpm<DIM> &cath_fpm,
        const BoundConds<DIM> &tis_deform_bc, const BoundConds<DIM> &cath_deform_bc, int step, Eigen::MatrixXd &forces);


    inline void UpdateDisplacements(const BoundConds<DIM> &tis_deform_bc, const BoundConds<DIM> &cath_deform_bc, int step,
        int tis_nodes_num, const Eigen::MatrixXd &ext_forces, const Eigen::MatrixXd &int_forces, const Eigen::MatrixXd &old_disp,
        const Eigen::MatrixXd &disp, Eigen::MatrixXd &new_disp, bool &disp_error_check);


    inline void ApplyContact(const ContactHandle<DIM,CELL_NODES> &contact_handle, const IMP::Mesh<DIM,CELL_NODES> &master_mesh,
        const IMP::Mesh<DIM,CELL_NODES> &slave_mesh, int master_pad, int slave_pad, Eigen::MatrixXd &disp);


    inline void ApplyContact(ContactHandle<DIM,CELL_NODES> &contact_handle, const IMP::Voronoi<DIM> &master_voro,
        const IMP::Voronoi<DIM> &slave_voro, int master_pad, int slave_pad, Eigen::MatrixXd &disp);


    inline void InitHistory(int instants_num, const Eigen::MatrixXd &disp);


    /**
     * \brief Add the displacements in the history container.
     * \param [in] pos The position of the container where the voltage is stored.
     * \return [void]
     */
    inline void AddInHistory(std::size_t pos, const Eigen::MatrixXd &disp);


    /**
     * \brief Get the critical time step of the explicit integration algorithms.
     * \return [double] The critical time step.
     */
    inline auto DtCritical() const { return this->dt_crit_; }


    /**
     * @brief
     * @return [double]
     */
    inline auto SimulationTime() const { return this->sim_time_; }


    /**
     * @brief
     * @return [double]
     */
    inline auto Dt() const { return this->dt_; }


    /**
     * @brief
     * @return const std::vector<Eigen::MatrixXd>&
     */
    inline auto & DispHistory() const { return this->disp_history_; }


    /**
     * @brief
     * @param id
     * @return const std::vector<Eigen::MatrixXd>&
     */
    inline auto & DispHistory(std::size_t id) const { return this->disp_history_[id]; }


protected:


    inline void CorrectSlavePos2D(const ContactHandle<DIM,CELL_NODES> &contact_handle, const IMP::Mesh<DIM,CELL_NODES> &master_mesh,
        const Eigen::MatrixXd &disp, int master_pad, const std::vector<IMP::Vec<DIM, double>> &master_pos, std::vector<IMP::Vec<DIM, double>> &slave_pos);


    inline void CorrectSlavePos3D(const ContactHandle<DIM,CELL_NODES> &contact_handle, const IMP::Mesh<DIM,CELL_NODES> &master_mesh,
        const Eigen::MatrixXd &disp, int master_pad, const std::vector<IMP::Vec<DIM, double>> &master_pos, std::vector<IMP::Vec<DIM, double>> &slave_pos);


    inline void CorrectSlavePos3D(ContactHandle<DIM,CELL_NODES> &contact_handle, const IMP::Voronoi<DIM> &master_voro,
        const Eigen::MatrixXd &disp, int master_pad, const std::vector<IMP::Vec<DIM, double>> &master_pos,
        std::vector<IMP::Vec<DIM, double>> &slave_pos);


    inline void CorrectSlavePos3D_v2(ContactHandle<DIM,CELL_NODES> &contact_handle, const IMP::Voronoi<DIM> &master_voro,
        const Eigen::MatrixXd &disp, int master_pad, const std::vector<IMP::Vec<DIM, double>> &master_pos,
        std::vector<IMP::Vec<DIM, double>> &slave_pos);


    /**
     * @brief
     * @param voro
     * @param fpm
     * @param hyperelastic
     * @param disp
     * @param forces
     * @return [void]
     */
    inline void ComputeFpmForceCorrection(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm,
        const std::shared_ptr<Constitutive> &hyperelastic, const Eigen::MatrixXd &disp, double penalty, Eigen::MatrixXd &forces);

};

/** \} End of Doxygen Groups */

} // End of namespace PNT

#endif //PHYNETOUCH_PHYSICS_DEFORMATION_HPP_

#include "PHYNETOUCH/engine/physics/deformation.tpp"