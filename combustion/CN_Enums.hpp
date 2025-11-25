//
// Created by Christian Messe on 25.09.19.
//

#ifndef BELFEM_CN_ENUMS_HPP
#define BELFEM_CN_ENUMS_HPP

#include "typedefs.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    enum class Oxidizer
    {
        LOX,
        AIR,
        H2O2,
        UNDEFINED
    };

//------------------------------------------------------------------------------

    enum class Fuel
    {
        LH2,
        LCH4,
        LNG,
        LC3H6,
        C2H5OH,
        UNDEFINED
    };

//------------------------------------------------------------------------------

    string
    fuel_to_string( const Fuel & aFuel );

//------------------------------------------------------------------------------

    Fuel
    string_to_fuel( const string & aFuel );

//------------------------------------------------------------------------------

    string
    oxidizer_to_string( const Oxidizer & aOxidizer );

//------------------------------------------------------------------------------

    Oxidizer
    string_to_oxidizer( const string & aOxidizer );

//------------------------------------------------------------------------------
}

#endif //BELFEM_CN_ENUMS_HPP
