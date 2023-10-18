/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */


#include "PHYNETOUCH/engine/utilities/petsc_idx_range.hpp"

namespace PNT {


PetscIdxRange::PetscIdxRange() : id_start_(0), id_end_(0), ids_num_(0)
{}


PetscIdxRange::~PetscIdxRange()
{}


void PetscIdxRange::Compute(PetscInt ids, PetscInt rank, PetscInt ncpu)
{
    this->id_start_ = rank * (ids / ncpu) + ((ids % ncpu) < rank ? (ids % ncpu) : rank);

    this->id_end_ = this->id_start_ + (ids / ncpu) + ((ids % ncpu) > rank);

    this->ids_num_ = this->id_end_ - this->id_start_;
}



} //end of namespace PNT
