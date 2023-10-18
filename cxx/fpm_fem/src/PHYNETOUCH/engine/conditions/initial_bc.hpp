/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */


/**
   \file initial_bc.hpp
   \brief Header file for Initial boundary condition.
   \author Konstantinos A. Mountris
   \date 19/09/2021
*/

#ifndef PHYNETOUCH_INITIAL_BC_HPP_
#define PHYNETOUCH_INITIAL_BC_HPP_


#include <Eigen/Eigen>

namespace PNT {

/** \addtogroup Conditions \{ */

/**
 * \class InitialBc
 * \author Konstantinos A. Mountris
 * \brief Initial boundary condition.
 */
class InitialBc {

private:

    Eigen::VectorXd values_;     /**< The initial values of the nodes */


public:

    /**
     * \brief The default constructor.
     */
    InitialBc() : values_() {}


    /**
     * \brief The destructor.
     */
    virtual ~InitialBc() {}


    /**
     * \brief Set the value of the initial condition.
     * \param [in] value The value of the initial condition.
     * \return [void]
     */
    inline void SetValues(int nodes_num, double value) { value * Eigen::VectorXd::Ones(nodes_num); }


    inline void SetValues(const Eigen::VectorXd &values) { this->values_ = values; }


    /**
     * \brief Get the value of the initial condition.
     * \return [double] The value of the initial condition.
     */
    inline auto & Values() const { return this->values_; }
    
    inline double Values(std::size_t id) const { return this->values_.coeffRef(id); }




};


/** \} End of Doxygen Groups */

} // End of namespace PNT

#endif //PHYNETOUCH_INITIAL_BC_HPP_