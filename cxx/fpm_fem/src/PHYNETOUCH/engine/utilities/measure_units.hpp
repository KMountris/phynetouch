/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */

/**
   \file measure_units.hpp
   \brief MeasureUnits class header file.
   \author Konstantinos A. Mountris
   \date 10/07/2021
*/

#ifndef PHYNETOUCH_UTILITIES_MEASURE_UNITS_HPP_
#define PHYNETOUCH_UTILITIES_MEASURE_UNITS_HPP_

#include <string>
#include <unordered_map>
#include <iostream>
#include <stdexcept>
#include <exception>
#include <cmath>


namespace PNT {

/** \addtogroup Utilities \{ */


/**
 * \class MeasureUnits
 * \brief Class implemmenting measure units conversion on pre-defined reference units in the SI system.
 * \author Konstantinos A. Mountris
 */
class MeasureUnits
{

private:

    double ref_time_si_value_;                          /**< The reference value of the time unit with respect to the SI system */

    double ref_length_si_value_;                        /**< The reference value of the length unit with respect to the SI system */

    double ref_mass_si_value_;                          /**< The reference value of the mass unit with respect to the SI system */

    double ref_temperature_si_value_;                   /**< The reference value of the temperature unit with respect to the SI system */

    double ref_current_si_value_;                       /**< The reference value of the current unit with respect to the SI system */

    double ref_voltage_si_value_;                       /**< The reference value of the voltage unit with respect to the SI system */

    std::unordered_map<std::string, double> units_;     /**< The measure units values with respect to the available reference unit values */


protected:

    /**
     * \brief Set the corresponding values of the time measure units considering the given reference values.
     * \return [void]
     */
    void SetTimeUnits() noexcept;


    /**
     * \brief Set the corresponding values of the length measure units considering the given reference values.
     * \return [void]
     */
    void SetLengthUnits() noexcept;


    /**
     * \brief Set the corresponding values of the mass measure units considering the given reference values.
     * \return [void]
     */
    void SetMassUnits() noexcept;


    /**
     * \brief Set the corresponding values of the temperature measure units considering the given reference values.
     * \return [void]
     */
    void SetTemperatureUnits() noexcept;


    /**
     * \brief Set the corresponding values of the current measure units considering the given reference values.
     * \return [void]
     */
    void SetCurrentUnits() noexcept;


    /**
     * \brief Set the corresponding values of the voltage measure units considering the given reference values.
     * \return [void]
     */
    void SetVoltageUnits() noexcept;


public:

    /**
     * \brief The MeasureUnits default constructor.
     */
    MeasureUnits();


    /**
     * \brief The MeasureUnits destructor.
     */
    virtual ~MeasureUnits();


    /**
     * \brief Set the reference time unit value with respect to the SI system (e.g., s=1., ds=0.1, ms=0.01, etc.).
     * \param [in] val The value of the reference time unit with respect to the SI system.
     * \return [void]
     */
    void SetRefTimeSIValue(double val) noexcept;


    /**
     * \brief Set the reference length unit value with respect to the SI system (e.g., m=1., dm=0.1, cm=0.01, etc.).
     * \param [in] val The value of the reference length unit with respect to the SI system.
     * \return [void]
     */
    void SetRefLengthSIValue(double val) noexcept;


    /**
     * \brief Set the reference mass unit value with respect to the SI system (e.g., kg=1000., g=1, mg=0.001, etc.).
     * \param [in] val The value of the reference mass unit with respect to the SI system.
     * \return [void]
     */
    void SetRefMassSIValue(double val) noexcept;


    /**
     * \brief Set the reference temperature unit value with respect to the SI system (e.g., kK=1000., K=1, dK=0.1, etc.).
     * \param [in] val The value of the reference temperature unit with respect to the SI system.
     * \return [void]
     */
    void SetRefTemperatureSIValue(double val) noexcept;


    /**
     * \brief Set the reference current unit value with respect to the SI system (e.g., A=1., dA=0.1, cA=0.01, etc.).
     * \param [in] val The value of the reference current unit with respect to the SI system.
     * \return [void]
     */
    void SetRefCurrentSIValue(double val) noexcept;


    /**
     * \brief Set the reference voltage unit value with respect to the SI system (e.g., V=1., dV=0.1, cV=0.01, mV=0.001, etc.).
     * \param [in] val The value of the reference voltage unit with respect to the SI system.
     * \return [void]
     */
    void SetRefVoltageSIValue(double val) noexcept;


    /**
     * \brief Set the corresponding values of the pressure measure units considering the given reference values.
     * \return [void]
     */
    void SetPressureUnits() noexcept;


    /**
     * \brief Set the corresponding values of the conductance measure units considering the given reference values.
     * \return [void]
     */
    void SetConductanceUnits() noexcept;


    /**
     * \brief Set the corresponding values of the energy measure units considering the given reference values.
     * \return [void]
     */
    void SetEnergyUnits() noexcept;


    /**
     * \brief Set the corresponding values of the capacitance measure units considering the given reference values.
     * \return [void]
     */
    void SetCapacitanceUnits() noexcept;



    /**
     * \brief Set the corresponding values of the power measure units considering the given reference values.
     * \return [void]
     */
    void SetPowerUnits() noexcept;


    /**
     * \brief Set the corresponding values of the force measure units considering the given reference values.
     * \return [void]
     */
    void SetForceUnits() noexcept;


    /**
     * \brief Operator to access values from the units container by key.
     * \note The returned value is expressed with respect to the reference unit values.
     *       For example: MeasureUnits["cm"] = 0.01 if the reference length unit value is in meters
     *       and MeasureUnits["cm"] = 1 if the reference length unit value is in centimeters.
     * \param [in] key The name of the unit to retrieve its value.
     * \return [const double&] The value of the unit corresponding to the given key.
     */
    const double & operator [] (const std::string &key) const;

};



/** \} End of Doxygen Groups */

} // End of namespace PNT


#endif // PHYNETOUCH_UTILITIES_MEASURE_UNITS_HPP_