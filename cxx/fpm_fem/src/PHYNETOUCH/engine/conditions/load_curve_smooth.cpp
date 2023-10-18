/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 */


#include "PHYNETOUCH/engine/conditions/load_curve_smooth.hpp"

namespace PNT {


LoadCurveSmooth::LoadCurveSmooth()
{}


LoadCurveSmooth::~LoadCurveSmooth() {}


void LoadCurveSmooth::DefineLoad(double load, double dt)
{
    int steps_num = 1;
    if (this->Duration() > 0.) {
        steps_num = static_cast<int>(std::round(this->Duration()/dt));
        this->SetDuration(dt*steps_num);
    }
    // Set the load curve data.
    double time_rel = 0.;
    this->data_.assign(steps_num,0.);
    for (int i = 0; i != steps_num; ++i) {
        time_rel = (i+1) * (dt/this->Duration());
        this->data_[i] = load * (10.*std::pow(time_rel,3) - 15.*std::pow(time_rel,4) + 6.*std::pow(time_rel,5));
    }

    // Ensure that the last entry of the data is equal to the load magnitude.
    if (std::fabs(this->data_.back()) < std::fabs(load) ||
            std::fabs(this->data_.back()) > std::fabs(load)) { 
        this->data_.back() = load;
    }

}



} // End of namespace PNT