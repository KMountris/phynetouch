/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */

#include "config_units.hpp"

namespace PNTSIM
{


ConfigUnits::ConfigUnits()
{}


ConfigUnits::~ConfigUnits()
{}


void ConfigUnits::SetReferenceScale(const Parser &parser, const MpiHandler &mpi_handler,
    MeasureUnits &units, std::ostream& stream) const
{
    // Create the SI units.
    MeasureUnits si_units;

    // Set the reference values to the scaled measure units.
    auto ref_time_unit = parser.GetValue<std::string>("reference units.time");
    auto ref_length_unit = parser.GetValue<std::string>("reference units.length");
    auto ref_mass_unit = parser.GetValue<std::string>("reference units.mass");
    auto ref_temperature_unit = parser.GetValue<std::string>("reference units.temperature");
    auto ref_current_unit = parser.GetValue<std::string>("reference units.current");
    auto ref_voltage_unit = parser.GetValue<std::string>("reference units.voltage");

    // Set reference units.
    units.SetRefTimeSIValue(si_units[ref_time_unit]);
    units.SetRefLengthSIValue(si_units[ref_length_unit]);
    units.SetRefMassSIValue(si_units[ref_mass_unit]);
    units.SetRefTemperatureSIValue(si_units[ref_temperature_unit]);
    units.SetRefCurrentSIValue(si_units[ref_current_unit]);
    units.SetRefVoltageSIValue(si_units[ref_voltage_unit]);

    // Set dependent units.
    units.SetEnergyUnits();
    units.SetPressureUnits();
    units.SetConductanceUnits();
    units.SetCapacitanceUnits();
    units.SetPowerUnits();
    units.SetForceUnits();

    if (mpi_handler.rank_id == 0) {
        stream << Logger::Message("Reference units\n");
        stream << "                  - time: " << ref_time_unit << "\n";
        stream << "                  - length: " << ref_length_unit << "\n";
        stream << "                  - mass: "  << ref_mass_unit << "\n";
        stream << "                  - temperature: "  << ref_temperature_unit << "\n";
        stream << "                  - voltage: " << ref_voltage_unit << "\n";
    }
}


} // end of namespace PNTSIM