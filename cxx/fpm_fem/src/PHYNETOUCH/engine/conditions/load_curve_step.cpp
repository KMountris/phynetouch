/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 */


#include "PHYNETOUCH/engine/conditions/load_curve_step.hpp"

namespace PNT {


LoadCurveStep::LoadCurveStep()
{}


LoadCurveStep::~LoadCurveStep() {}


void LoadCurveStep::DefineLoad(double load, double dt)
{
    int steps_num = 1;
    if (this->Duration() > 0.) {
        steps_num = static_cast<int>(std::round(this->Duration()/dt));
        this->SetDuration(dt*steps_num);

        // Set the load curve data.
        this->data_.assign(steps_num,0.);
        for (int i = 0; i < steps_num; ++i) {
            if (i*dt >= this->Start()) {
                this->data_[i] = load;
            }
        }
    } else {
        this->data_.assign(steps_num,load);
    }
}



} // End of namespace PNT