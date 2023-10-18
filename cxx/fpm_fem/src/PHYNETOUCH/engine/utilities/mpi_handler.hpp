/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */

/**
   \file mpi_handler.hpp
   \brief MpiHandler struct header file.
   \author Konstantinos A. Mountris
   \date 18/05/2022
*/

#ifndef PHYNETOUCH_UTILITIES_MPI_HANDLER_HPP_
#define PHYNETOUCH_UTILITIES_MPI_HANDLER_HPP_


namespace PNT {

/** \addtogroup Utilities \{ */


/**
 * \struct MpiHandler
 * \brief Struct implementing a handler for MPI related variables.
 */

struct MpiHandler
{
    int rank_id;                /*!< The index of the MPI rank. */

    int procs_num;              /*!< The number of the MPI processors. */
};

/** @} End of Doxygen Groups*/
} //end of namespace PNT

#endif //PHYNETOUCH_UTILITIES_MPI_HANDLER_HPP_
