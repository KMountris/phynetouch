/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */

/**
   \file prescribed_bc.hpp
   \brief Header file for Prescribed boundary condition.
   \author Konstantinos A. Mountris
   \date 25/08/2022
*/

#ifndef PHYNETOUCH_CONDITIONS_PRESCRIBED_BC_HPP_
#define PHYNETOUCH_CONDITIONS_PRESCRIBED_BC_HPP_

#include <Eigen/Dense>

#include <vector>



namespace PNT {

/** \addtogroup Conditions \{ */

/**
 * \class PrescribedBc
 * \author Konstantinos A. Mountris
 * \brief Prescribed boundary condition. If the boundary condition is applied on a scalar field its dimension is 1.
 *        Otherwise it is the same with the dimension of the geometry.
 */
class PrescribedBc {

private:

    std::vector<int> node_ids_;                     /**< The indices of the nodes where the values are prescribed */

    std::vector<double> values_;                    /**< The values to be prescribed */

public:

    /**
     * \brief The default constructor.
     */
    PrescribedBc();


    /**
     * \brief The destructor.
     */
    virtual ~PrescribedBc();


    /**
     * @brief 
     * 
     * @param mat_cath 
     * @param values_cath 
     * @param nnum_tis 
     */
    void ConnectBodies(const std::vector<int> &connector_node_ids,
        const std::vector<std::vector<int>> &connected_node_ids,
        const Eigen::VectorXd &connector_values, int connected_body_nodes_num);


    /**
     * \brief Get the indices of the nodes where the Prescribed condition is applied.
     * \return [const std::vector<int>&] The nodes where the Prescribed condition is applied.
     */
    inline auto & NodeIds() const { return this->node_ids_; }


    /**
     * \brief Get the values of the Prescribed condition.
     * \return [const std::vector<double>&] The values of the Prescribed condition.
     */
    inline auto & Values() const { return this->values_; }


};


/** \} End of Doxygen Groups */

} // End of namespace PNT

#endif //PHYNETOUCH_CONDITIONS_PRESCRIBED_BC_HPP_
