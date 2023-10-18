/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */


#include "PHYNETOUCH/engine/utilities/measure_units.hpp"


namespace PNT {

MeasureUnits::MeasureUnits() : ref_time_si_value_(1.), ref_length_si_value_(1.), ref_mass_si_value_(1000.),
    ref_temperature_si_value_(1.), ref_current_si_value_(1.), ref_voltage_si_value_(1.), units_()
{
    // Set the default values of the units.
    this->SetTimeUnits();
    this->SetLengthUnits();
    this->SetMassUnits();
    this->SetTemperatureUnits();
    this->SetCurrentUnits();
    this->SetVoltageUnits();

    // Dependent units.
    this->SetEnergyUnits();
    this->SetPressureUnits();
    this->SetConductanceUnits();
    this->SetCapacitanceUnits();
    this->SetPowerUnits();
    this->SetForceUnits();
}


MeasureUnits::~MeasureUnits()
{}


void MeasureUnits::SetRefTimeSIValue(double val) noexcept
{
    this->ref_time_si_value_ = val;

    // Reset the values of the time units to account for the new reference unit.
    this->SetTimeUnits();
}


void MeasureUnits::SetRefLengthSIValue(double val) noexcept
{
    this->ref_length_si_value_ = val;

    // Reset the values of the length units to account for the new reference unit.
    this->SetLengthUnits();
}


void MeasureUnits::SetRefMassSIValue(double val) noexcept
{
    this->ref_mass_si_value_ = val;

    // Reset the values of the mass units to account for the new reference unit.
    this->SetMassUnits();
}



void MeasureUnits::SetRefTemperatureSIValue(double val) noexcept
{
    this->ref_temperature_si_value_ = val;

    // Reset the values of the temperature units to account for the new reference unit.
    this->SetTemperatureUnits();
}


void MeasureUnits::SetRefCurrentSIValue(double val) noexcept
{
    this->ref_current_si_value_ = val;

    // Reset the values of the current units to account for the new reference unit.
    this->SetCurrentUnits();
}



void MeasureUnits::SetRefVoltageSIValue(double val) noexcept
{
    this->ref_voltage_si_value_ = val;

    // Reset the values of the voltage units to account for the new reference unit.
    this->SetVoltageUnits();
}


const double & MeasureUnits::operator[] (const std::string &key) const
{
    return this->units_.at(key);
}


void MeasureUnits::SetTimeUnits() noexcept
{

    this->units_["h"]   = 3600. / this->ref_time_si_value_;
    this->units_["min"] = 60. / this->ref_time_si_value_;
    this->units_["s"]   = 1. / this->ref_time_si_value_;
    this->units_["ms"]  = 1.E-3 / this->ref_time_si_value_;
    this->units_["us"]  = 1.E-6 / this->ref_time_si_value_;
    this->units_["ns"]  = 1.E-9 / this->ref_time_si_value_;
    this->units_["ps"]  = 1.E-12 / this->ref_time_si_value_;
    this->units_["fs"]  = 1.E-15 / this->ref_time_si_value_;
}


void MeasureUnits::SetLengthUnits() noexcept
{
    // Length units.
    this->units_["Mm"]  = 1.E6 / this->ref_length_si_value_;
    this->units_["km"]  = 1.E3 / this->ref_length_si_value_;
    this->units_["hm"]  = 1.E2 / this->ref_length_si_value_;
    this->units_["dam"] = 1.E1 / this->ref_length_si_value_;
    this->units_["m"]   = 1. / this->ref_length_si_value_;
    this->units_["dm"]  = 1.E-1 / this->ref_length_si_value_;
    this->units_["cm"]  = 1.E-2 / this->ref_length_si_value_;
    this->units_["mm"]  = 1.E-3 / this->ref_length_si_value_;
    this->units_["um"]  = 1.E-6 / this->ref_length_si_value_;
    this->units_["nm"]  = 1.E-9 / this->ref_length_si_value_;
    this->units_["pm"]  = 1.E-12 / this->ref_length_si_value_;
    this->units_["fm"]  = 1.E-15 / this->ref_length_si_value_;

    // Surface units.
    this->units_["Mm2"]  = this->units_.at("Mm")*this->units_.at("Mm");
    this->units_["km2"]  = this->units_.at("km")*this->units_.at("km");
    this->units_["hm2"]  = this->units_.at("hm")*this->units_.at("hm");
    this->units_["dam2"] = this->units_.at("dam")*this->units_.at("dam");
    this->units_["m2"]   = this->units_.at("m")*this->units_.at("m");
    this->units_["dm2"]  = this->units_.at("dm")*this->units_.at("dm");
    this->units_["cm2"]  = this->units_.at("cm")*this->units_.at("cm");
    this->units_["mm2"]  = this->units_.at("mm")*this->units_.at("mm");
    this->units_["um2"]  = this->units_.at("um")*this->units_.at("um");
    this->units_["nm2"]  = this->units_.at("nm")*this->units_.at("nm");
    this->units_["pm2"]  = this->units_.at("pm")*this->units_.at("pm");
    this->units_["fm2"]  = this->units_.at("fm")*this->units_.at("fm");

    // Volume units.
    this->units_["Mm3"]  = this->units_.at("Mm")*this->units_.at("Mm")*this->units_.at("Mm");
    this->units_["km3"]  = this->units_.at("km")*this->units_.at("km")*this->units_.at("km");
    this->units_["hm3"]  = this->units_.at("hm")*this->units_.at("hm")*this->units_.at("hm");
    this->units_["dam3"] = this->units_.at("dam")*this->units_.at("dam")*this->units_.at("dam");
    this->units_["m3"]   = this->units_.at("m")*this->units_.at("m")*this->units_.at("m");
    this->units_["dm3"]  = this->units_.at("dm")*this->units_.at("dm")*this->units_.at("dm");
    this->units_["cm3"]  = this->units_.at("cm")*this->units_.at("cm")*this->units_.at("cm");
    this->units_["mm3"]  = this->units_.at("mm")*this->units_.at("mm")*this->units_.at("mm");
    this->units_["um3"]  = this->units_.at("um")*this->units_.at("um")*this->units_.at("um");
    this->units_["nm3"]  = this->units_.at("nm")*this->units_.at("nm")*this->units_.at("nm");
    this->units_["pm3"]  = this->units_.at("pm")*this->units_.at("pm")*this->units_.at("pm");
    this->units_["fm3"]  = this->units_.at("fm")*this->units_.at("fm")*this->units_.at("fm");

}


void MeasureUnits::SetMassUnits() noexcept
{
    this->units_["t"]   = 1.E6 / this->ref_mass_si_value_;
    this->units_["kg"] = 1.E3 / this->ref_mass_si_value_;
    this->units_["g"]   = 1. / this->ref_mass_si_value_;
    this->units_["mg"]  = 1.E-3 / this->ref_mass_si_value_;
    this->units_["mug"]  = 1.E-6 / this->ref_mass_si_value_;
    this->units_["ng"]  = 1.E-9 / this->ref_mass_si_value_;
}



void MeasureUnits::SetCurrentUnits() noexcept
{
    this->units_["GA"] = 1.E9 / this->ref_current_si_value_;
    this->units_["MA"] = 1.E6 / this->ref_current_si_value_;
    this->units_["kA"] = 1.E3 / this->ref_current_si_value_;
    this->units_["A"]  = 1. / this->ref_current_si_value_;
    this->units_["dA"] = 1.E-1 / this->ref_current_si_value_;
    this->units_["cA"] = 1.E-2 / this->ref_current_si_value_;
    this->units_["mA"] = 1.E-3 / this->ref_current_si_value_;
    this->units_["uA"] = 1.E-6 / this->ref_current_si_value_;
    this->units_["nA"] = 1.E-9 / this->ref_current_si_value_;
    this->units_["pA"] = 1.E-12 / this->ref_current_si_value_;
    this->units_["fA"] = 1.E-15 / this->ref_current_si_value_;
}


void MeasureUnits::SetVoltageUnits() noexcept
{
    this->units_["GV"] = 1.E9 / this->ref_voltage_si_value_;
    this->units_["MV"] = 1.E6 / this->ref_voltage_si_value_;
    this->units_["kV"] = 1.E3 / this->ref_voltage_si_value_;
    this->units_["V"]  = 1. / this->ref_voltage_si_value_;
    this->units_["dV"] = 1.E-1 / this->ref_voltage_si_value_;
    this->units_["cV"] = 1.E-2 / this->ref_voltage_si_value_;
    this->units_["mV"] = 1.E-3 / this->ref_voltage_si_value_;
    this->units_["uV"] = 1.E-6 / this->ref_voltage_si_value_;
    this->units_["nV"] = 1.E-9 / this->ref_voltage_si_value_;
    this->units_["pV"] = 1.E-12 / this->ref_voltage_si_value_;
    this->units_["fV"] = 1.E-15 / this->ref_voltage_si_value_;
}


void MeasureUnits::SetTemperatureUnits() noexcept
{
    this->units_["K"]  = 1. / this->ref_temperature_si_value_;
}


void MeasureUnits::SetEnergyUnits() noexcept
{
    auto joule = this->units_.at("kg") * this->units_.at("m2") / std::pow(this->units_.at("s"), 2.);

    this->units_["GJ"] = 1.E9 * joule;
    this->units_["MJ"] = 1.E6 * joule;
    this->units_["kJ"] = 1.E3 * joule;
    this->units_["J"]  = 1.;//joule;
    this->units_["dJ"] = 1.E-1 * joule;
    this->units_["cJ"] = 1.E-2 * joule;
    this->units_["mJ"] = 1.E-3 * joule;
    this->units_["uJ"] = 1.E-6 * joule;
    this->units_["nJ"] = 1.E-9 * joule;
    this->units_["pJ"] = 1.E-12 * joule;
    this->units_["fJ"] = 1.E-15 * joule;
}


void MeasureUnits::SetPressureUnits() noexcept
{
    auto pascal = this->units_.at("kg") / (this->units_.at("m") * std::pow(this->units_.at("s"), 2.));

    this->units_["GPa"] = 1.E9 * pascal;
    this->units_["MPa"] = 1.E6 * pascal;
    this->units_["kPa"] = 1.E3 * pascal;
    this->units_["Pa"]  = pascal;
    this->units_["dPa"] = 1.E-1 * pascal;
    this->units_["cPa"] = 1.E-2 * pascal;
    this->units_["mPa"] = 1.E-3 * pascal;
    this->units_["uPa"] = 1.E-6 * pascal;
    this->units_["nPa"] = 1.E-9 * pascal;
    this->units_["pPa"] = 1.E-12 * pascal;
    this->units_["fPa"] = 1.E-15 * pascal;
}


void MeasureUnits::SetConductanceUnits() noexcept
{
    auto siemens = this->units_.at("A") / this->units_.at("V");

    this->units_["GS"] = 1.E9 * siemens;
    this->units_["MS"] = 1.E6 * siemens;
    this->units_["kS"] = 1.E3 * siemens;
    this->units_["S"]  = siemens;
    this->units_["dS"] = 1.E-1 * siemens;
    this->units_["cS"] = 1.E-2 * siemens;
    this->units_["mS"] = 1.E-3 * siemens;
    this->units_["uS"] = 1.E-6 * siemens;
    this->units_["nS"] = 1.E-9 * siemens;
    this->units_["pS"] = 1.E-12 * siemens;
    this->units_["fS"] = 1.E-15 * siemens;
}


void MeasureUnits::SetCapacitanceUnits() noexcept
{
    auto farad = std::pow(this->units_.at("s"),4.)*std::pow(this->units_.at("A"),2.) / (this->units_.at("kg")*this->units_.at("m2"));

    this->units_["GF"] = 1.E9 * farad;
    this->units_["MF"] = 1.E6 * farad;
    this->units_["kF"] = 1.E3 * farad;
    this->units_["F"]  = farad;
    this->units_["dF"] = 1.E-1 * farad;
    this->units_["cF"] = 1.E-2 * farad;
    this->units_["mF"] = 1.E-3 * farad;
    this->units_["uF"] = 1.E-6 * farad;
    this->units_["nF"] = 1.E-9 * farad;
    this->units_["pF"] = 1.E-12 * farad;
    this->units_["fF"] = 1.E-15 * farad;
}


void MeasureUnits::SetPowerUnits() noexcept
{
    auto watt = this->units_.at("kg")*this->units_.at("m2") / std::pow(this->units_.at("s"),3.);

    this->units_["GW"] = 1.E9 * watt;
    this->units_["MW"] = 1.E6 * watt;
    this->units_["kW"] = 1.E3 * watt;
    this->units_["W"]  = 1.;//watt;
    this->units_["dW"] = 1.E-1 * watt;
    this->units_["cW"] = 1.E-2 * watt;
    this->units_["mW"] = 1.E-3 * watt;
    this->units_["uW"] = 1.E-6 * watt;
    this->units_["nW"] = 1.E-9 * watt;
    this->units_["pW"] = 1.E-12 * watt;
    this->units_["fW"] = 1.E-15 * watt;
}


void MeasureUnits::SetForceUnits() noexcept
{
    auto newton = this->units_.at("kg") * this->units_.at("m") / std::pow(this->units_.at("s"), 2.);

    this->units_["GN"] = 1.E9 * newton;
    this->units_["MN"] = 1.E6 * newton;
    this->units_["kN"] = 1.E3 * newton;
    this->units_["N"]  = newton;
    this->units_["dN"] = 1.E-1 * newton;
    this->units_["cN"] = 1.E-2 * newton;
    this->units_["mN"] = 1.E-3 * newton;
    this->units_["uN"] = 1.E-6 * newton;
    this->units_["nN"] = 1.E-9 * newton;
    this->units_["pN"] = 1.E-12 * newton;
    this->units_["fN"] = 1.E-15 * newton;
}


} // End of namespace PNT