/*
 * ELECTRA. Electrophysiology Simulation Software.
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
   \file bin_exporter.hpp
   \brief BinaryExporter class header file.
   \author Konstantinos A. Mountris
   \date 20/06/2020
*/

#ifndef ELECTRA_EXPORTERS_BINARY_EXPORTER_HPP_
#define ELECTRA_EXPORTERS_BINARY_EXPORTER_HPP_

#include "ELECTRA/engine/electrophysiology/ep_factory.hpp"
#include "ELECTRA/engine/utilities/logger.hpp"

#include <boost/filesystem.hpp>

#include <string>
#include <vector>
#include <iterator>
#include <stdexcept>
#include <exception>
#include <iostream>
#include <sstream>
#include <fstream>
#include <memory>


namespace ELECTRA {

/** \addtogroup Exporters \{ */

/**
 * \class BinaryExporter
 * \brief Class implemmenting output in binary format.
 */
class BinaryExporter {

public:
    /**
     * \brief BinaryExporter constructor.
     */
    BinaryExporter();


    /**
     * \brief BinaryExporter destructor.
     */
    virtual ~BinaryExporter();


    void WriteCellsState(const std::vector<std::unique_ptr<BasicEp>> &cells, const std::string &filename);

};


/** \} End of Doxygen Groups*/

} // End of namespace ELECTRA.

#endif  //ELECTRA_EXPORTERS_BINARY_EXPORTER_HPP_
