/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */

/**
   \file config_approximation.hpp
   \brief ConfigApproximation class header file.
   \author Konstantinos A. Mountris
   \date 11/07/2021
*/

#pragma once
#ifndef PHYNETOUCH_APPS_TOOLS_CONFIG_APPROXIMATION_HPP_
#define PHYNETOUCH_APPS_TOOLS_CONFIG_APPROXIMATION_HPP_

#include "parser.hpp"

#include "PHYNETOUCH/engine/utilities/logger.hpp"

#include <IMP/IMP>
#include <CLOUDEA/CLOUDEA>

#include <termcolor/termcolor.hpp>
#include <Eigen/Dense>

#include <string>
#include <filesystem>
#include <utility>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <functional>

using namespace PNT;

namespace PNTSIM {

/** \addtogroup Application-Tools \{ */

/**
 * \class ConfigApproximation
 * \brief Class to configure the numerical approximation method of a simulation.
 */
template <short DIM, short CELL_NODES>
class ConfigApproximation {

public:

    /**
     * \brief Constructor of the Config Approximation.
     */
    ConfigApproximation();


    /**
     * \brief Destructor the Config Approximation.
     */
    virtual ~ConfigApproximation();


    /**
     * \brief Set the approximation using the Finite Element Method (FEM).
     * \param [in] mesh The mesh topology for which we will set the approximation.
     * \param [out] fem The FEM matrices for all the gauss points of the mesh.
     * \param [out] stream The message output stream.
     * \return [void]
     */
    void SetFemApproximation(const MpiHandler &mpi_handler, const IMP::Mesh<DIM,CELL_NODES> &mesh,
        CLOUDEA::FemMats<DIM,CELL_NODES> &fem, std::ostream &stream) const;


    /**
     * \brief Set the approximation using the Fragile Points Method (FPM).
     * \param [in] parser The parser of the simulation input file.
     * \param [in] voro The voronoi tesselation for which we will set the approximation.
     * \param [out] fpm_approx The FPM approximation for all the subdomains of the voronoi tesselation.
     * \param [out] stream The message output stream.
     * \return [void]
     */
    void SetFpmApproximation(const Parser &parser, const MpiHandler &mpi_handler, const IMP::Voronoi<DIM> &voro,
        CLOUDEA::Fpm<DIM> &fpm_approx, std::ostream &stream) const;


    void SetMcmApproximation(const Parser &parser, const MpiHandler &mpi_handler, const IMP::Grid<DIM, CELL_NODES> &grid,
        std::unique_ptr<CLOUDEA::Mfree<DIM>> &mcm_approx, IMP::NodeSet &neumann_nset, std::ostream &stream) const;


};

/** \} End of Doxygen Groups*/

} //end of namespace PNTSIM

#endif //PHYNETOUCH_APPS_TOOLS_CONFIG_APPROXIMATION_HPP_

#include "config_approximation.tpp"