/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file constitutive_factory.hpp
   \brief ConstitutiveFactory class header file.
   \author Konstantinos A. Mountris
   \date 23/10/2019
*/

#pragma once
#ifndef PHYNETOUCH_MATERIALS_CONSTITUTIVE_FACTORY_HPP_
#define PHYNETOUCH_MATERIALS_CONSTITUTIVE_FACTORY_HPP_

#include "PHYNETOUCH/engine/materials/constitutive.hpp"
#include "PHYNETOUCH/engine/materials/neohookean.hpp"
#include "PHYNETOUCH/engine/materials/rigid.hpp"


#include <memory>
#include <string>


namespace PNT {

/** \addtogroup Materials \{ */

/**
 * \class ConstitutiveFactory
 * \brief Class implementing a factory for constitutive law generation.
 */
class ConstitutiveFactory
{
public:

    /**
     * \brief Creates a constitutive law material according to the provided type.
     * \param [in] constitutive_type The type of the desired constitutive law.
     * \return [std::shared_ptr<Constitutive>] Shared pointer to the created constitutive law.
     */
    inline static std::shared_ptr<Constitutive> Create(ConstitutiveType constitutive_type) {

        std::shared_ptr<Constitutive> ct_ptr;

        switch (constitutive_type) {
        case ConstitutiveType::neohookean :
            ct_ptr = std::make_shared<Neohookean>();
            break;
        case ConstitutiveType::rigid :
            ct_ptr = std::make_shared<Rigid>();
            break;
        default:
            auto err_msg = Logger::Error("Could not create constitutive model of unsupported type.");
            throw std::invalid_argument(Logger::Error(err_msg));
            break;
        }
        return ct_ptr;
    }
};


/** \} End of Doxygen Groups */

} // End of namespace PNT

#endif // PHYNETOUCH_MATERIALS_CONSTITUTIVE_FACTORY_HPP_