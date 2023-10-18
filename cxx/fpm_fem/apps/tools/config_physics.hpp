/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */

/**
   \file config_physics.hpp
   \brief ConfigPhysics class header file.
   \author Konstantinos A. Mountris
   \date 01/08/2021
*/

#pragma once
#ifndef PHYNETOUCH_APPS_TOOLS_CONFIG_PHYSICS_HPP_
#define PHYNETOUCH_APPS_TOOLS_CONFIG_PHYSICS_HPP_

#include "parser.hpp"

#include "PHYNETOUCH/engine/utilities/logger.hpp"
#include "PHYNETOUCH/engine/conditions/bound_conds.hpp"
#include "PHYNETOUCH/engine/materials/cardiac_tissue.hpp"
#include "PHYNETOUCH/engine/materials/catheter.hpp"
#include "PHYNETOUCH/engine/materials/constitutive.hpp"
#include "PHYNETOUCH/engine/physics/bioheat.hpp"
#include "PHYNETOUCH/engine/physics/contact_handle.hpp"
#include "PHYNETOUCH/engine/physics/deformation.hpp"
#include "PHYNETOUCH/engine/physics/electrical.hpp"

#include <IMP/IMP>
#include <CLOUDEA/CLOUDEA>
#include <termcolor/termcolor.hpp>
#include <Eigen/Eigen>

#include <string>
#include <filesystem>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <memory>

using namespace PNT;

namespace PNTSIM {

/** \addtogroup Application-Tools \{ */

/**
 * \class ConfigPhysics
 * \brief Class to configure the physics of a simulation.
 */
template <short DIM, short CELL_NODES>
class ConfigPhysics {


public:

    /**
     * \brief The constructor of the Config Physics.
     */
    ConfigPhysics();


    /**
     * \brief The destructor of the Config Physics.
     */
    virtual ~ConfigPhysics();


    /**
     * @brief Set the Contact Handle object
     * 
     * @param parser 
     * @param tis_mesh 
     * @param tis_voro 
     * @param cath_mesh 
     * @param cath_voro 
     * @param contact_handle 
     * @param stream 
     */
    void SetContactHandle(const Parser &parser, const MpiHandler &mpi_handler,const IMP::Mesh<DIM,CELL_NODES> &tis_mesh,
        const IMP::Voronoi<DIM> &tis_voro, const IMP::Mesh<DIM,CELL_NODES> &cath_mesh, const IMP::Voronoi<DIM> &cath_voro,
        ContactHandle<DIM,CELL_NODES> &contact_handle, std::ostream &stream) const;


    /**
     * \brief Solve an electrical problem.
     * \param [in] parser The parser of the simulation input file.
     * \param [in] tissue_mesh The mesh topology of the tissue.
     * \param [in] tissue_voro The voronoi tesselation of the tissue.
     * \param [in] catheter_mesh The mesh topology of the catheter.
     * \param [in] catheter_voro The voronoi tesselation of the catheter.
     * \param [in] tissue_mat The material of the tissue.
     * \param [in] catheter_mat The material of the catheter.
     * \param [in] tissue_fpm The Fragile Points Method (FPM) approximant of the tissue.
     * \param [in] catheter_fpm The Fragile Points Method (FPM) approximant of the catheter.
     * \param [in] tissue_bc The tissue boundary conditions for the electrical problem.
     * \param [in] catheter_bc The catheter boundary conditions for the electrical problem.
     * \param [in] units The measure units.
     * \param [out] electrical The physics of the electrical problem.
     * \param [out] stream The message output stream.
     * \return [void]
     */
    void SolveElectrical(const Parser &parser, const MpiHandler &mpi_handler, const IMP::Mesh<DIM,CELL_NODES> &tissue_mesh,
        const IMP::Voronoi<DIM> &tissue_voro, const IMP::Mesh<DIM,CELL_NODES> &catheter_mesh,
        const IMP::Voronoi<DIM> &catheter_voro, CardiacTissue &tissue_mat, const Catheter<DIM> &catheter_mat,
        const CLOUDEA::Fpm<DIM> &tissue_fpm, const CLOUDEA::Fpm<DIM> &catheter_fpm,
        const BoundConds<1> &tissue_bc, const BoundConds<1> &catheter_bc,
        Electrical<DIM,CELL_NODES> &electrical, std::ostream &stream) const;


    void SolveBioheat(const Parser &parser, const MpiHandler &mpi_handler,
        const IMP::Mesh<DIM,CELL_NODES> &mesh_tis, const IMP::Voronoi<DIM> &voro_tis,
        const IMP::Mesh<DIM,CELL_NODES> &mesh_cath, const IMP::Voronoi<DIM> &voro_cath,
        CardiacTissue &mat_tis, Catheter<DIM> &mat_cath, const CLOUDEA::Fpm<DIM> &fpm_tis,
        const CLOUDEA::Fpm<DIM> &fpm_cath, BoundConds<1> &bioheat_bc_tis, BoundConds<1> &bioheat_bc_cath,
        const MeasureUnits &units, Bioheat<DIM,CELL_NODES> &bioheat_tis,
        Bioheat<DIM,CELL_NODES> &bioheat_cath, std::ostream &stream) const;


    /**
     * @brief
     * @param tissue_mesh
     */
    void SolveDeformation(const Parser &parser, const MpiHandler &mpi_handler, const IMP::Mesh<DIM,CELL_NODES> &tis_mesh,
        const IMP::Voronoi<DIM> &tis_voro, const IMP::Mesh<DIM,CELL_NODES> &cath_mesh,
        const IMP::Voronoi<DIM> &cath_voro, const std::shared_ptr<Constitutive> &tis_elastic,
        const std::shared_ptr<Constitutive> &cath_elastic, const CLOUDEA::FemMats<DIM,CELL_NODES> &tis_fem,
        const CLOUDEA::FemMats<DIM,CELL_NODES> &cath_fem, const CLOUDEA::Fpm<DIM> &tis_fpm, const CLOUDEA::Fpm<DIM> &cath_fpm,
        const MeasureUnits &units, BoundConds<DIM> &tis_deform_bc, BoundConds<DIM> &cath_deform_bc,
        Deformation<DIM,CELL_NODES> &deformation, std::ostream &stream) const;


    /**
     * \brief Solve a multiphysics problem.
     * \param [in] parser The parser of the simulation input file.
     * \param [in] tissue_mesh The mesh topology of the tissue.
     * \param [in] tissue_voro The voronoi tesselation of the tissue.
     * \param [in] catheter_mesh The mesh topology of the catheter.
     * \param [in] catheter_voro The voronoi tesselation of the catheter.
     * \param [in] tissue_material The material of the tissue.
     * \param [in] catheter_material The material of the catheter.
     * \param [in] tissue_fpm The Fragile Points Method (FPM) approximant of the tissue.
     * \param [in] catheter_fpm The Fragile Points Method (FPM) approximant of the catheter.
     * \param [in] tissue_electrical_bc The tissue boundary conditions for the electrical problem.
     * \param [in] catheter_electrical_bc The catheter boundary conditioons for the electrical problem.
     * \param [in] tissue_bioheat_bc The tissue boundary conditions for the bioheat problem.
     * \param [in] catheter_bioheat_bc The catheter boundary conditions for the bioheat problem.
     * \param [in] units The measure units.
     * \param [out] electrical The physics of the electrical problem.
     * \param [out] bioheat The physics of the bioheat problem.
     * \param [out] stream The message output stream.
     * \return [void]
     */
    void SolveMultiphysics(const Parser &parser, const MpiHandler &mpi_handler, const IMP::Mesh<DIM,CELL_NODES> &tissue_mesh,
        const IMP::Voronoi<DIM> &tissue_voro, const IMP::Mesh<DIM,CELL_NODES> &catheter_mesh,
        const IMP::Voronoi<DIM> &catheter_voro, CardiacTissue &tissue_material, Catheter<DIM> &catheter_material,
        const CLOUDEA::Fpm<DIM> &tissue_fpm, const CLOUDEA::Fpm<DIM> &catheter_fpm, BoundConds<1> &tissue_electrical_bc,
        BoundConds<1> &catheter_electrical_bc, BoundConds<1> &tissue_bioheat_bc, BoundConds<1> &catheter_bioheat_bc,
        const MeasureUnits &units, Electrical<DIM,CELL_NODES> &electrical_tis, Electrical<DIM,CELL_NODES> &electrical_cath,
        Bioheat<DIM,CELL_NODES> &bioheat_tis, Bioheat<DIM,CELL_NODES> &bioheat_cath, std::ostream &stream) const;


protected:

    /**
     * \brief Solve a coupled bioheat problem with an electrical problem.
     * \param [in] parser The parser of the simulation input file.
     * \param [in] tissue_mesh The mesh topology of the tissue.
     * \param [in] tissue_voro The voronoi tesselation of the tissue.
     * \param [in] catheter_mesh The mesh topology of the catheter.
     * \param [in] catheter_voro The voronoi tesselation of the catheter.
     * \param [in] tissue_material The material of the tissue.
     * \param [in] catheter_material The material of the catheter.
     * \param [in] tissue_fpm The Fragile Points Method (FPM) approximant of the tissue.
     * \param [in] catheter_fpm The Fragile Points Method (FPM) approximant of the catheter.
     * \param [in] tissue_electrical_bc The tissue boundary conditions for the electrical problem.
     * \param [in] catheter_electrical_bc The catheter boundary conditioons for the electrical problem.
     * \param [in] tissue_bioheat_bc The tissue boundary conditions for the bioheat problem.
     * \param [in] catheter_bioheat_bc The catheter boundary conditions for the bioheat problem.
     * \param [in] units The measure units.
     * \param [out] electrical The physics of the electrical problem.
     * \param [out] bioheat The physics of the bioheat problem.
     * \param [out] stream The message output stream.
     * \return [void]
     */
    void SolveCoupledBioheat(const Parser &parser, const MpiHandler &mpi_handler, const IMP::Mesh<DIM,CELL_NODES> &tissue_mesh,
        const IMP::Voronoi<DIM> &tissue_voro, const IMP::Mesh<DIM,CELL_NODES> &catheter_mesh,
        const IMP::Voronoi<DIM> &catheter_voro, CardiacTissue &tissue_material, Catheter<DIM> &catheter_material,
        const CLOUDEA::Fpm<DIM> &tissue_fpm, const CLOUDEA::Fpm<DIM> &catheter_fpm,
        BoundConds<1> &tissue_electrical_bc, BoundConds<1> &catheter_electrical_bc,
        BoundConds<1> &tissue_bioheat_bc, BoundConds<1> &catheter_bioheat_bc,
        const MeasureUnits &units, Electrical<DIM,CELL_NODES> &electrical_tis, Electrical<DIM,CELL_NODES> &electrical_cath,
        Bioheat<DIM,CELL_NODES> &bioheat_tis, Bioheat<DIM,CELL_NODES> &bioheat_cath, std::ostream &stream) const;


};

/** \} End of Doxygen Groups*/

} //end of namespace PNTSIM

#endif //PHYNETOUCH_APPS_TOOLS_CONFIG_PHYSICS_HPP_

#include "config_physics.tpp"