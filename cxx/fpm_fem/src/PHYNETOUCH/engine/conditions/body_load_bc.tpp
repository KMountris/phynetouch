/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */


#ifndef PHYNETOUCH_BODY_LOAD_BC_TPP_
#define PHYNETOUCH_BODY_LOAD_BC_TPP_

#include "PHYNETOUCH/engine/conditions/body_load_bc.hpp"


namespace PNT {


template<short DIM>
BodyLoadBc<DIM>::BodyLoadBc() : dir_(), value_(0.), load_curve_()
{}


template<short DIM>
BodyLoadBc<DIM>::~BodyLoadBc() {}


template<short DIM>
void BodyLoadBc<DIM>::SetValue(double value)
{
    this->value_ = value;
}


template<short DIM>
void BodyLoadBc<DIM>::SetDirection(const IMP::Vec<DIM,double> &dir)
{
    this->dir_ = dir;
}


template<short DIM>
void BodyLoadBc<DIM>::SetLoadingCurve(LoadCurveType load_curve_type, double load_start, double load_duration)
{
    this->load_curve_ = LoadCurveFactory::Create(load_curve_type);
    this->load_curve_->SetStart(load_start);
    this->load_curve_->SetDuration(load_duration);

    this->load_type_ = load_curve_type;
    this->load_start_ = load_start;
    this->load_duration_ = load_duration;
}


template<short DIM>
void BodyLoadBc<DIM>::LoadingDefinition(double dt)
{
    this->load_curve_->DefineLoad(this->value_, dt);
}


template<short DIM>
void BodyLoadBc<DIM>::Apply(int step, double dt, Eigen::RowVectorXd &load) const
{
    // Apply boundary condition if time has passed the starting time.
    auto time = static_cast<double>(step)*dt;
    if (time >= this->load_curve_->Start()) {
        // Compute the load curve data entry index for the current time instance.
        time -= this->load_curve_->Start();
        auto id = static_cast<std::size_t>(std::round(time/dt));

        // Load value for given time instance.
        auto val = this->load_curve_->Data().back();
        if (id < this->load_curve_->Data().size()) { val = this->load_curve_->Data()[id]; }

        // Apply load value.
        load.setZero();
        for (short d = 0; d != DIM; ++d) {
            if (std::fabs(this->dir_[d]) > 2.*std::numeric_limits<double>::epsilon()) {
                load.coeffRef(d) = val*this->dir_[d];
            }
        }
    }
}



} // End of namespace PNT

#endif //PHYNETOUCH_BODY_LOAD_BC_TPP_