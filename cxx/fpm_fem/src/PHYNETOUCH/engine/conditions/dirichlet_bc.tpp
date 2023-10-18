/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */

#ifndef PHYNETOUCH_DIRICHLET_BC_TPP_
#define PHYNETOUCH_DIRICHLET_BC_TPP_

#include "PHYNETOUCH/engine/conditions/dirichlet_bc.hpp"


namespace PNT {


template<short DIM>
DirichletBc<DIM>::DirichletBc() : node_ids_(), dir_(), value_(0.), load_curve_()
{}


template<short DIM>
DirichletBc<DIM>::~DirichletBc() {}


template<short DIM>
void DirichletBc<DIM>::SetNodeIds(const std::vector<int> &node_ids)
{
    this->node_ids_ = node_ids;
}


template<short DIM>
void DirichletBc<DIM>::SetValue(double value)
{
    this->value_ = value;
}


template<short DIM>
void DirichletBc<DIM>::SetDirection(const IMP::Vec<DIM,double> &dir)
{
    this->dir_ = dir;
}


template<short DIM>
void DirichletBc<DIM>::SetLoadingCurve(LoadCurveType load_curve_type, double load_start, double load_duration)
{
    this->load_curve_ = LoadCurveFactory::Create(load_curve_type);
    this->load_curve_->SetStart(load_start);
    this->load_curve_->SetDuration(load_duration);
}


template<short DIM>
void DirichletBc<DIM>::LoadingDefinition(double dt)
{
    this->load_curve_->DefineLoad(this->value_, dt);
}


template<short DIM>
void DirichletBc<DIM>::UpdateValue(int step, double dt)
{
    auto time = step*dt;
    if (this->load_curve_->Data().size() > 0) {
        auto id = static_cast<std::size_t>(std::round(time/dt));
    
        // Load value for given time instance.
        auto new_val = this->load_curve_->Data().back();
        if (id < this->load_curve_->Data().size()) {
            new_val = this->load_curve_->Data()[id];
        }
        this->SetValue(new_val);
    }
}

template<short DIM>
void DirichletBc<DIM>::Apply(int padding, int step, double dt, Eigen::MatrixXd &displacement) const
{
    // Apply boundary condition if time has passed the starting time.
    double time = step*dt;
    if (time >= this->load_curve_->Start()) {
        // Compute the load curve data entry index for the current time instance.
        time -= this->load_curve_->Start();
        auto id = static_cast<std::size_t>(std::round(time/dt));

        // Load value for given time instance.
        auto val = this->load_curve_->Data().back();
        if (id < this->load_curve_->Data().size()) { val = this->load_curve_->Data()[id]; }

        // Apply load value.
        for (const auto &nid : this->node_ids_) {
            for (short d = 0; d != DIM; ++d) {
                if (std::fabs(this->dir_[d]) > 0.) {
                    displacement.coeffRef(nid+padding, d) = val*this->dir_[d];
                }
            }
        }
    }
}



} // End of namespace PNT

#endif //PHYNETOUCH_DIRICHLET_BC_TPP_