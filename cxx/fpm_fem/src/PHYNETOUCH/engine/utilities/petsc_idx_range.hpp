/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */

/**
   \file mpi_index_range.hpp
   \brief MpiIndexRange struct header file.
   \author Konstantinos A. Mountris
   \date 28/05/2022
*/

#ifndef PHYNETOUCH_UTILITIES_PETSC_IDX_RANGE_HPP_
#define PHYNETOUCH_UTILITIES_PETSC_IDX_RANGE_HPP_

#include <petscsys.h>
#include <algorithm>

namespace PNT {

/** \addtogroup Utilities \{ */


/**
 * \class PetscIdxRange
 * \brief Class implementing an index range for Petsc parallel ranks.
 */

class PetscIdxRange
{
private:

    PetscInt id_start_;      /*!< The first index in the range */

    PetscInt id_end_;        /*!< The last index in the range */

    PetscInt ids_num_;       /*!< The number of indices in the range */

public:

    /**
     * \brief The object constructor.
    */
    PetscIdxRange();


    /**
     * \brief The object destructor.
    */
    virtual ~PetscIdxRange();


    /**
     * \brief Compute the range for a Petsc rank.
     * \param [in] mpi_handler The MPI handler to retrive the MPI rank.
     * \param [in] total_ids The total number of indices from where the range will be obtained.
     * \return [void]
     */
    void Compute(PetscInt ids, PetscInt rank, PetscInt ncpu);


    inline auto IdStart() const { return this->id_start_; }


    inline auto IdEnd() const { return this->id_end_; }


    inline auto IdsNum() const { return this->ids_num_; }

};

/** @} End of Doxygen Groups*/
} //end of namespace PNT

#endif //PHYNETOUCH_UTILITIES_PETSC_IDX_RANGE_HPP_
