/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */

/**
   \file timer.hpp
   \brief Timer class header file.
   \author Konstantinos A. Mountris
   \date 10/07/2021
*/

#ifndef PHYNETOUCH_UTILITIES_TIMER_HPP_
#define PHYNETOUCH_UTILITIES_TIMER_HPP_


#include <iostream>
#include <chrono>
#include <string>

namespace PNT {

/** \addtogroup Utilities \{ */


/**
 * \class Timer
 * \brief Class implemmenting a timer for profiling of PHYNETOUCH.
 */

class Timer
{
public:

    /*!
     * \brief Timer constructor.
     */
    Timer();


    /*!
     * \brief Timer destructor.
     */
    virtual ~Timer();


    /*!
     * \brief Reset the timer.
     * \return [void]
     */
    void Reset();


    /*!
     * \brief Get the elapsed time in milliseconds.
     * \return [double] The elapsed time in milliseconds.
     */
    double ElapsedMilliSecs() const;


    /*!
     * \brief Get the elapsed time in seconds.
     * \return [double] The elapsed time in seconds.
     */
    double ElapsedSecs() const;


    /*!
     * \brief Get the elapsed time in minutes.
     * \return [double] The elapsed time in minutes.
     */
    double ElapsedMinutes() const;


    /*!
     * \brief Get the elapsed time in hours.
     * \return [double] The elapsed time in hours.
     */
    double ElapsedHours() const;


    /*!
     * \brief Print the elapsed time either in milliseconds, seconds, or minutes in string format.
     *        The result depends on the total elapsed time.
     * \return [std::string] The elapsed time.
     */
    std::string PrintElapsedTime() const;

private:
    typedef std::chrono::high_resolution_clock timer_clock_;                /*!< The timer's clock. */

    typedef std::chrono::duration<double, std::ratio<1> > timer_second_;    /*!< The time in seconds. */

    typedef std::chrono::duration<double, std::milli > timer_millisec_;     /*!< The time in milliseconds. */

    std::chrono::time_point<timer_clock_> beg_;                             /*!< The beginning of time. */
};

/** @} End of Doxygen Groups*/
} //end of namespace PNT

#endif //PHYNETOUCH_UTILITIES_EXPLICIT_SIM_TIMER_HPP_
