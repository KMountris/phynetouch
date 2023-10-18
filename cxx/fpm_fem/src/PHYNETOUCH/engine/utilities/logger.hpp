/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */

/**
   \file logger.hpp
   \brief Logger class header file.
   \author Konstantinos A. Mountris
   \date 10/07/2021
*/

#ifndef PHYNETOUCH_UTILITIES_LOGGER_HPP_
#define PHYNETOUCH_UTILITIES_LOGGER_HPP_

#include <string>

namespace PNT {

/** \addtogroup Utilities \{ */

/**
 * \class Logger
 * \brief Class implemmenting output messages for logging of PHYNETOUCH.
 */

class Logger
{
public:

    /**
     * \brief Logger default constructor.
     */
    Logger() {}


    /**
     * \brief Logger default destructor.
     */
    virtual ~Logger() {}


    /**
     * \brief Log an PHYNETOUCH message.
     * \param [in] msg The message to be logged.
     * \return [std::string] The logged message.
     */
    inline static std::string Message(const std::string &msg) {
        return "[PHYNETOUCH] " + msg;
    }


    /**
     * \brief Log an PHYNETOUCH error.
     * \param [in] err The error to be logged.
     * \return [std::string] The logged error.
     */
    inline static std::string Error(const std::string &err) {
        return "[PHYNETOUCH ERROR] " + err;
    }


    /**
     * \brief Log an PHYNETOUCH warning.
     * \param [in] wrng The warning to be logged.
     * \return [std::string] The logged warning.
     */
    inline static std::string Warning(const std::string &wrng) {
        return "\n[PHYNETOUCH WARNING] " + wrng + "\n";
    }


    /**
     * \brief Log an PHYNETOUCH ToDO message.
     * \param [in] todo The ToDo message to be logged.
     * \return [std::string] The logged ToDo message.
     */
    inline static std::string ToDo(const std::string &todo) {
        return "[PHYNETOUCH TODO] " + todo;
    }


};

/** @} End of Doxygen Groups*/

} //end of namespace PNT

#endif //PHYNETOUCH_UTILITIES_LOGGER_HPP_
