/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */

/**
   \file KDTreeIMPVecVectorAdaptor.hpp
   \brief KDTreeIMPVecVectorAdaptor class header file.
   \author Konstantinos A. Mountris
   \date 10/07/2021
*/

#pragma once
#ifndef PHYNETOUCH_UTILITIES_KDTREEIMPVECVECTORADAPTOR_HPP_
#define PHYNETOUCH_UTILITIES_KDTREEIMPVECVECTORADAPTOR_HPP_

#include "PHYNETOUCH/engine/utilities/logger.hpp"

#include <nanoflann.hpp>

#include <IMP/Vectors>

#include <vector>
#include <stdexcept>
#include <exception>

// ===== This example shows how to use nanoflann with these types of containers:
// using my_vector_of_vectors_t = std::vector<std::vector<double> > ;
//
// The next one requires #include <Eigen/Dense>
// using my_vector_of_vectors_t = std::vector<Eigen::VectorXd> ;
// =============================================================================

namespace PNT {

/** \addtogroup Utilities \{ */

/** A simple vector-of-vectors adaptor for nanoflann, without duplicating the
 * storage. The i'th vector represents a point in the state space.
 *
 *  \tparam DIM If set to >0, it specifies a compile-time fixed dimensionality
 *      for the points in the data set, allowing more compiler optimizations.
 *  \tparam num_t The type of the point coordinates (typ. double or float).
 *  \tparam Distance The distance metric to use: nanoflann::metric_L1,
 *          nanoflann::metric_L2, nanoflann::metric_L2_Simple, etc.
 *  \tparam IndexType The type for indices in the KD-tree index
 *         (typically, size_t of int)
 */
template <
    class IMPVecVectorType, typename num_t = double, int DIM = -1,
    class Distance = nanoflann::metric_L2, typename IndexType = size_t>
struct KDTreeIMPVecVectorAdaptor
{
    using self_t =
        KDTreeIMPVecVectorAdaptor<IMPVecVectorType, num_t, DIM, Distance>;
    using metric_t =
        typename Distance::template traits<num_t, self_t>::distance_t;
    using index_t =
        nanoflann::KDTreeSingleIndexAdaptor<metric_t, self_t, DIM, IndexType>;

    /** The kd-tree index for the user to call its methods as usual with any
     * other FLANN index */
    index_t* index = nullptr;

    /// Constructor: takes a const ref to the vector of vectors object with the
    /// data points
    KDTreeIMPVecVectorAdaptor(const size_t , const IMPVecVectorType& mat, const int leaf_max_size = 10) : data_(mat)
    {
        assert(mat.size() != 0 && mat[0].Data().size() != 0);
        const size_t dims = mat[0].Data().size();
        if (DIM > 0 && static_cast<int>(dims) != DIM) {
            auto err_msg = "Data set dimensionality does not match the 'DIM' template argument";
            throw std::runtime_error(Logger::Error(err_msg));
        }
        index = new index_t(static_cast<int>(dims), *this, nanoflann::KDTreeSingleIndexAdaptorParams(leaf_max_size));
        index->buildIndex();
    }

    ~KDTreeIMPVecVectorAdaptor() { delete index; }

    const IMPVecVectorType& data_;

    /** Query for the \a num_closest closest points to a given point
     *  (entered as query_point[0:dim-1]).
     *  Note that this is a short-cut method for index->findNeighbors().
     *  The user can also call index->... methods as desired.
     *
     * \note nChecks_IGNORED is ignored but kept for compatibility with
     * the original FLANN interface.
     */
    inline void query(
        const num_t* query_point, const size_t num_closest,
        IndexType* out_indices, num_t* out_distances_sq,
        const int nChecks_IGNORED = 10) const
    {
        nanoflann::KNNResultSet<num_t, IndexType> resultSet(num_closest);
        resultSet.init(out_indices, out_distances_sq);
        index->findNeighbors(resultSet, query_point, nanoflann::SearchParams());
    }

    // Interface expected by KDTreeSingleIndexAdaptor
    const self_t& derived() const { return *this; }
    self_t&       derived() { return *this; }

    // Must return the number of data points
    inline size_t kdtree_get_point_count() const { return data_.size(); }

    // Returns the dim'th component of the idx'th point in the class:
    inline num_t kdtree_get_pt(const size_t idx, const size_t dim) const
    {
        return data_[idx][dim];
    }

    // Optional bounding-box computation: return false to default to a standard
    // bbox computation loop.
    // Return true if the BBOX was already computed by the class and returned
    // in "bb" so it can be avoided to redo it again. Look at bb.size() to
    // find out the expected dimensionality (e.g. 2 or 3 for point clouds)
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /*bb*/) const
    {
        return false;
    }

};

/** @} End of Doxygen Groups*/

} //end of namespace PNT

#endif // PHYNETOUCH_UTILITIES_KDTREEIMPVECVECTORADAPTOR_HPP_
