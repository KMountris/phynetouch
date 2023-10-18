/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */


/**
 * \file parser.hpp
 * \brief Parser class header file
 * \author Konstantinos A. Mountris
 * \date 10/05/2021
 */

#pragma once
#ifndef PHYNETOUCH_APPS_TOOLS_PARSER_HPP_
#define PHYNETOUCH_APPS_TOOLS_PARSER_HPP_

#include "PHYNETOUCH/engine/utilities/logger.hpp"

#include <nlohmann/json.hpp>

#include <string>
#include <sstream>
#include <fstream>
#include <streambuf>
#include <iostream>
#include <stdexcept>
#include <exception>
#include <filesystem>

using namespace PNT;

/** \addtogroup Application-Tools \{ */

/**
 * \namespace PNTSIM
 * \brief Collects tools for configuring and executing PHYNETOUCH applications.
 */
namespace PNTSIM {

class Parser
{
private:

    nlohmann::json json_;           /**< The parser of the PHYNETOUCH application JSON file  */

    std::string parent_path_;       /**< The absolute path of the PHYNETOUCH application JSON file */


protected:

    /**
     * \brief Removes comment lines from a string.
     * Comments are expected in the style of C++ and JavaScript.
     * \param [in] str The string containing comment lines.
     * \return [std::string] A new string whith the comment lines removed.
     */
    std::string RemoveCommentLines(const std::string &str);

public:

    /**
     * \brief Construct a new Parser object.
     * \param [in] json_filename The filename of the configuration file of the PHYNETOUCH application including the path.
     */
    Parser(const std::string &json_filename);


    /**
     * \brief Destroy the Parser object.
     */
    virtual ~Parser();


    /**
     * \brief Get the value of an attribute in the configuration file.
     * \tparam T The type of the value to be retrieved.
     * \param [in] attribute The name of the requested attribute including the complete path hierarchy.
     * \return [T] The value of an attribute in the configuration file.
     */
    template <class T>
    T GetValue(const std::string &attribute) const;


    /**
     * @brief Get the Object object
     * 
     * @param attribute 
     * @return nlohmann::json 
     */
    nlohmann::json GetObject(const std::string &attribute) const;


    /**
     * \brief Check if an attribute exists in the configuration file.
     * \param attribute The name of the requested attribute including the complete path hierarchy.
     * \return [true] The attribute exists in the configuration file.
     * \return [false] The attribute does not exist in the configuration file.
     */
    bool HasAttribute(const std::string &attribute) const;


    /**
     * @brief Check if the data of the attribute is a single value.
     * @param attribute The name of the requested attribute including the complete path hierarchy.
     * @return [true] The data of the attribute is a single value.
     * @return [false] The data of the attribute is not a single value.
     */
    bool IsSingleValue(const std::string& attribute) const;


    /**
     * \brief Check if the data of the attribute is an array.
     * \param attribute The name of the requested attribute including the complete path hierarchy.
     * \return [true] The data of the attribute is an array.
     * \return [false] The data of the attribute is not an array.
     */
    bool IsArray(const std::string& attribute) const;


    /**
     * \brief Check if the data of the attribute is an array of arrays.
     * \param attribute The name of the requested attribute including the complete path hierarchy.
     * \return [true] The data of the attribute is an array of arrays.
     * \return [false] The data of the attribute is not an array of arrays.
     */
    bool IsMultiArray(const std::string& attribute) const;


    /**
     * \brief Resolve the path of a given file. If the file's path is not valid then the parent path of the parser
     * is used.
     * \param [out] filename The name of the file for which the path will be resolved.
     * \return [void]
     */
    void ResolveFilePath(std::string &filename) const;


    /**
     * \brief Get the parent path of the parsed simulation configuration file.
     * \return [const std::string&] The parent path of the parsed simulation configuration file.
     */
    inline const std::string & ParentPath() const { return this->parent_path_; }


};


} // end of namespace PNTSIM
/** \} End of Doxygen Groups */

#endif //PHYNETOUCH_APPS_TOOLS_PARSER_HPP_

#include "parser.tpp"