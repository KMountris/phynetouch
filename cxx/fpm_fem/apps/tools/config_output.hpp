/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */

/**
   \file config_output.hpp
   \brief ConfigOutput class header file.
   \author Konstantinos A. Mountris
   \date 04/08/2021
*/

#pragma once
#ifndef PHYNETOUCH_APPS_TOOLS_CONFIG_OUTPUT_HPP_
#define PHYNETOUCH_APPS_TOOLS_CONFIG_OUTPUT_HPP_

#include "parser.hpp"

#include "PHYNETOUCH/engine/utilities/logger.hpp"
#include "PHYNETOUCH/engine/physics/electrical.hpp"
#include "PHYNETOUCH/engine/physics/bioheat.hpp"
#include "PHYNETOUCH/engine/physics/deformation.hpp"
#include "PHYNETOUCH/engine/exporters/ensight_exporter.hpp"

#include <Eigen/Eigen>
#include <IMP/IMP>
#include <termcolor/termcolor.hpp>

#include <string>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <memory>
#include <set>

using namespace PNT;

namespace PNTSIM {

/** \addtogroup Application-Tools \{ */

/**
 * \class ConfigOutput
 * \brief Class to configure the output of simulation results.
 */
template<short DIM, short CELL_NODES>
class ConfigOutput {

private:

    std::vector<std::string> states_types_;

    std::vector<std::string> field_names_;

    std::vector<std::string> field_types_;

    std::vector<std::size_t> field_steps_;


public:

    /**
     * \brief ConfigOutput object constructor.
     */
    ConfigOutput();


    /**
     * \brief ConfigOutput object destructor.
     */
    virtual ~ConfigOutput();


    /**
     * @brief 
     * 
     * @param parser 
     * @param tissue_mesh 
     * @param catheter_mesh 
     * @param electrical 
     * @param stream 
     */
    // void OutputElectrical(const Parser &parser, const MpiHandler &mpi_handler, const IMP::Mesh<DIM,CELL_NODES> &mesh_tis,
    //     const IMP::Mesh<DIM,CELL_NODES> &mesh_cath, const Electrical<DIM, CELL_NODES> &electrical, std::ostream &stream) const;


    /**
     * @brief 
     * 
     * @param parser 
     * @param tissue_mesh 
     * @param catheter_mesh 
     * @param electrical 
     * @param stream 
     */
    // void OutputBioheat(const Parser &parser, const MpiHandler &mpi_handler,
    //     const IMP::Mesh<DIM,CELL_NODES> &mesh_tis, const IMP::Mesh<DIM,CELL_NODES> &mesh_cath,
    //     const Bioheat<DIM, CELL_NODES> &bioheat, std::ostream &stream) const;


    /**
     * @brief 
     * 
     * @param parser 
     * @param tissue_mesh 
     * @param deformation 
     * @param stream 
     */
    void OutputDeformation(const Parser &parser, const MpiHandler &mpi_handler, const IMP::Mesh<DIM,CELL_NODES> &mesh_tis,
        const IMP::Mesh<DIM,CELL_NODES> &mesh_cath, const Deformation<DIM, CELL_NODES> &deformation, std::ostream &stream);


    void OutputMultiphysics(const Parser &parser, const MpiHandler &mpi_handler,
        const IMP::Mesh<DIM,CELL_NODES> &mesh_tis, const IMP::Mesh<DIM,CELL_NODES> &mesh_cath,
        const IMP::Voronoi<DIM> &voro_tis, const IMP::Voronoi<DIM> &voro_cath,
        const CardiacTissue &mat_tis, const Catheter<DIM> &mat_cath,
        const CLOUDEA::Fpm<DIM> &fpm_tis, const CLOUDEA::Fpm<DIM> &fpm_cath,
        const Electrical<DIM, CELL_NODES> &electrical_tis, const Electrical<DIM, CELL_NODES> &electrical_cath,
        const Bioheat<DIM, CELL_NODES> &bioheat_tis, const Bioheat<DIM, CELL_NODES> &bioheat_cath, std::ostream &stream);



protected:

    void Reset();

    void OutputMultiphysicsToEnsight(const Parser &parser, const MpiHandler &mpi_handler,
        const std::string &body_type, const IMP::Mesh<DIM,CELL_NODES> &mesh, const IMP::Voronoi<DIM> &voro,
        const CLOUDEA::Fpm<DIM> &fpm, const std::vector<double> &electrical_conductivity,
        const Electrical<DIM, CELL_NODES> &electrical, const Bioheat<DIM, CELL_NODES> &bioheat, std::ostream &stream);


    void OutputEnsightGeo(const Parser &parser, const MpiHandler &mpi_handler,
        const std::string &body_type, const IMP::Mesh<DIM,CELL_NODES> &mesh, std::ostream &stream);


    void OutputEnsightTemperature(const Parser &parser, const MpiHandler &mpi_handler,
        const std::string &body_type, const std::vector<Eigen::VectorXd> &temperatures, std::ostream &stream);


    void OutputEnsightVoltage(const Parser &parser, const MpiHandler &mpi_handler,
        const std::string &body_type, const std::vector<Eigen::VectorXd> &voltages, std::ostream &stream);


    void OutputEnsightHeatSource(const Parser &parser, const MpiHandler &mpi_handler,
    const std::string &body_type, const IMP::Mesh<DIM,CELL_NODES> &mesh, const IMP::Voronoi<DIM> &voro,
    const CLOUDEA::Fpm<DIM> &fpm, const std::vector<double> &electrical_conductivity,
    const std::vector<Eigen::VectorXd> &voltages, std::ostream &stream);


    void OutputEnsightElectricField(const Parser &parser, const MpiHandler &mpi_handler,
    const std::string &body_type, const IMP::Mesh<DIM,CELL_NODES> &mesh, const IMP::Voronoi<DIM> &voro,
    const CLOUDEA::Fpm<DIM> &fpm, const std::vector<Eigen::VectorXd> &voltages, std::ostream &stream);


    void OutputEnsightScalarStates(const Parser &parser, const MpiHandler &mpi_handler, const std::string &body_type,
        const std::string &states_type, const std::vector<Eigen::VectorXd> &scalars, std::ostream &stream);


    void OutputEnsightVectorStates(const Parser &parser, const MpiHandler &mpi_handler, const std::string &body_type,
        const std::string &states_type, const std::vector<Eigen::MatrixXd> &vectors, std::ostream &stream);

    void OutputEnsightVectorStatesDeprecated(const Parser &parser, const MpiHandler &mpi_handler, const std::string &body_type,
        const std::string &states_type, const std::vector<Eigen::MatrixXd> &vectors, int nodes_num, std::ostream &stream);

    void OutputEnsightAnimation(const Parser &parser, const MpiHandler &mpi_handler,
        const std::string &body_type, double time_inc, std::ostream &stream);

};

/** \} End of Doxygen Groups*/

} //end of namespace PNTSIM

#endif //PHYNETOUCH_APPS_TOOLS_CONFIG_OUTPUT_HPP_

#include "config_output.tpp"