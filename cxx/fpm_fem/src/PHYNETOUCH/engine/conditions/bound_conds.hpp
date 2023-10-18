/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
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
   \file bound_conds.hpp
   \brief Header file for boundary conditions.
   \author Konstantinos A. Mountris
   \date 19/09/2021
*/

#pragma once
#ifndef PHYNETOUCH_BOUND_CONDS_HPP_
#define PHYNETOUCH_BOUND_CONDS_HPP_


#include "PHYNETOUCH/engine/conditions/body_load_bc.hpp"
#include "PHYNETOUCH/engine/conditions/dirichlet_bc.hpp"
#include "PHYNETOUCH/engine/conditions/initial_bc.hpp"
#include "PHYNETOUCH/engine/conditions/prescribed_bc.hpp"


namespace PNT {

/** \addtogroup Conditions \{ */

/**
 * \class BoundConds
 * \author Konstantinos A. Mountris
 * \brief Boundary conditions.
 */
template<short DIM>
class BoundConds {

private:

    std::vector<DirichletBc<DIM>> dirichlet_;           /**< The Dirichlet boundary conditions */

    BodyLoadBc<DIM> body_load_;                         /**< The body loads conditions */

    InitialBc initial_;                                 /**< The initial boundary condition */

    PrescribedBc prescribed_;                   /**< The prescribed values boundary condition */

    bool has_body_load_;


public:

    /**
     * \brief The default constructor.
     */
    BoundConds() : dirichlet_(), body_load_(), initial_(), prescribed_(), has_body_load_(false) {}


    /**
     * \brief The destructor.
     */
    virtual ~BoundConds() {}


    /**
     * \brief Set the Dirichlet boundary conditions.
     * \param [in] dirichlet The Dirichlet boundary conditions.
     * \return [void]
     */
    inline void SetDirichlet(const std::vector<DirichletBc<DIM>> &dirichlet) { this->dirichlet_ = dirichlet; }


    /**
     * \brief Set the body loads conditions.
     * \param [in] body_loads The body loads conditions.
     * \return [void]
     */
    inline void SetBodyLoad(const BodyLoadBc<DIM> &body_load) {
        this->body_load_.SetValue(body_load.Value());
        this->body_load_.SetDirection(body_load.Direction());
        this->body_load_.SetLoadingCurve(body_load.LoadType(), body_load.LoadStart(), body_load.LoadDuration());

        this->has_body_load_ = true;
    }


    /**
     * \brief Set the initial boundary condition.
     * \param initial The initial boundary condition.
     * \return [void]
     */
    inline void SetInitial(const InitialBc &initial) { this->initial_ = initial; }


    /**
     * @brief Set the Prescribed Bc object
     * 
     * @param prescribed_bc 
     */
    inline void SetPrescribed(const PrescribedBc &prescribed) { this->prescribed_ = prescribed; }


    inline void UpdateDirichlet(int step, double dt) {
        for (auto &dir_bc : this->dirichlet_)
            dir_bc.UpdateValue(step, dt);
    }


    /**
     * \brief Get the Dirichlet boundary conditions.
     * \return [const std::vector<DirichletBc>&] The Dirichlet boundary conditions.
     */
    inline auto & Dirichlet() const { return this->dirichlet_; }


    /**
     * \brief Get the Dirichlet boundary conditions with editing access.
     * \return [std::vector<DirichletBc<DIM>>&]
     */
    inline auto & Dirichlet() { return this->dirichlet_; }


    /**
     * \brief Get the number of Dirichlet boundary conditions.
     * \return [int] The number of Dirichlet boundary conditions. 
     */
    inline auto DirichletNum() const { return static_cast<int>(this->dirichlet_.size()); }


    /**
     * \brief Get the body loads conditions.
     * \return [const std::vector<BodyLoadBc>&] The body loads boundary conditions.
     */
    inline auto & BodyLoad() const { return this->body_load_; }


    /**
     * \brief Get the body loads conditions with editing access.
     * \return [std::vector<BodyLoadBc<DIM>>&] The body loads boundary conditions.
     */
    inline auto & BodyLoad() { return this->body_load_; }


    /**
     * @brief 
     * 
     * @return true 
     * @return false 
     */
    inline bool HasBodyLoad() const { return this->has_body_load_; }


    /**
     * \brief Get the Initial boundary condition.
     * \return [const InitialBc&] The initial boundary condition.
     */
    inline auto & Initial() const { return this->initial_; }


    /**
     * @brief 
     * 
     * @return auto& 
     */
    inline auto & Prescribed() const { return this->prescribed_; }
    inline auto & Prescribed() { return this->prescribed_; }
};


/** \} End of Doxygen Groups */

} // End of namespace PNT

#endif //PHYNETOUCH_BOUND_CONDS_HPP_