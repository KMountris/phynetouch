/*
 * PHYNETOUCH. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/**
   \file ensight_exporter.hpp
   \brief EnsightExporter class header file.
   \author Konstantinos A. Mountris
   \date 25/09/2019
*/

#ifndef PHYNETOUCH_EXPORTERS_ENSIGHT_EXPORTER_HPP_
#define PHYNETOUCH_EXPORTERS_ENSIGHT_EXPORTER_HPP_

#include "PHYNETOUCH/engine/utilities/logger.hpp"

#include <boost/filesystem.hpp>
#include <boost/range.hpp>

#include <IMP/Vectors>
#include <IMP/Tesselations>

#include <Eigen/Dense>

#include <string>
#include <iostream>
#include <cstddef>
#include <vector>
#include <map>
#include <iterator>
#include <stdexcept>
#include <exception>
#include <filesystem>


namespace PNT {

/** \addtogroup Exporters \{ */

/**
 * \class EnsightExporter
 * \brief Class implemmenting the output of models and their solution fields for post-processing in ParaView.
 * \tparam DIM The dimensions of the domain. Supported: [1 | 2 | 3].
 * \tparam CELL_NODES The number of nodes composing the cells of the domain.
 */

template <short DIM, short CELL_NODES=1>
class EnsightExporter {

public:
    /**
     * \brief EnsightExporter constructor.
     */
    EnsightExporter();
    
    
    /**
     * \brief EnsightExporter destructor.
     */
    virtual ~EnsightExporter();
    
    
    /**
     * \brief  Save a Ensight geometry (.geo) file for the given model.
     * \param [in] mesh The model's mesh.
     * \return [void]
     */
    void SaveGeo(const std::vector<IMP::Vec<DIM, double>> &nodes, const std::vector<IMP::Cell<DIM, CELL_NODES>> &cells, const std::string &out_filename);
    
    
    //void AddPointNormals();

    
    /**
     * \brief 
     * \param scalar_field 
     * \param scalar_field_name 
     */
    void SaveScalarField(const Eigen::VectorXd &scalar_field, const std::string &out_filename);


    /**
     * \brief 
     * \param scalar_field 
     * \param scalar_field_name 
     */
    void SaveVectorField(const Eigen::MatrixXd &vector_field, const std::string &out_filename);

    /**
     * \brief Create a ParaView animation (.pvd) xml file.
     *
     * If a non-existing path is given for the animation file it is generated automatically by the exporter.
     * Similarly if the file's name has not .pvd extension it is added automatically.
     *
     * \param [in] vtu_files_dir The name of the directory where the vtu files to included in the animation are located.
     * \param [in] files_number The number of .vtu files to be included in the animation's collection.
     * \param [in] pvd_filename The filename of the ParaView animation (.pvd) xml file.
     * \return [void]
     */
    void SaveAnimation(const std::string &out_filename, const std::string &geo_filename,
            const std::vector<std::string> &states_files, const std::vector<std::string> &fields_names, const std::vector<std::string> &fields_types,
            const std::vector<std::size_t> &fields_steps, double time_inc);



};


/** \} End of Doxygen Groups */

} // End of namespace PNT.

#endif  //PHYNETOUCH_PARAVIEW_EXPORTER_PARAVIEW_EXPORTER_HPP_

#include "PHYNETOUCH/engine/exporters/ensight_exporter.tpp"
