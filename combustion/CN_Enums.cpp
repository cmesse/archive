//
// Created by Christian Messe on 25.09.19.
//

#include "CN_Enums.hpp"
#include "assert.hpp"
#include "stringtools.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    string
    fuel_to_string( const Fuel & aFuel )
    {
        switch( aFuel )
        {
            case( Fuel::LH2 ):
            {
                return "H2";
            }
            case( Fuel::LCH4 ):
            {
                return "CH4";
            }
            case( Fuel::LNG ):
            {
                return "LNG";
            }
            case( Fuel::LC3H6 ):
            {
                return "C3H6";
            }
            case( Fuel::C2H5OH ):
            {
                return "C2H5OH";
            }
            default :
            {
                BELFEM_ERROR( false, "Unknown Fuel." );
                return "unknown";
            }
        }
    }

//------------------------------------------------------------------------------

    Fuel
    string_to_fuel( const string & aFuel )
    {
        if( aFuel == "H2" || aFuel == "LH2" || aFuel == "HYDROGEN" )
        {
            return Fuel::LH2 ;
        }
        else if( aFuel == "CH4" || aFuel == "LCH4" || aFuel == "METHANE" )
        {
            return Fuel::LCH4 ;
        }
        else if( aFuel == "C3H6" || aFuel == "LC3H6" || aFuel == "PROPENE" )
        {
            return Fuel::LC3H6 ;
        }
        else if( aFuel == "LNG" )
        {
            return Fuel::LNG ;
        }
        else
        {
            BELFEM_ERROR( false, "Unknown Fuel : %s.", aFuel.c_str() );
            return Fuel::UNDEFINED;
        }
    }

//------------------------------------------------------------------------------

    string
    oxidizer_to_string( const Oxidizer & aOxidizer )
    {
        switch( aOxidizer )
        {
            case( Oxidizer::LOX ):
            {
                return "O2";
            }
            case( Oxidizer::H2O2 ):
            {
                return "H2O2";
            }
            case( Oxidizer::AIR  ):
            {
                return "Air";
            }
            default :
            {
                BELFEM_ERROR( false, "Unknown Oxidizer." );
                return "unknown";
            }
        }
    }

//------------------------------------------------------------------------------

    Oxidizer
    string_to_oxidizer( const string & aOxidizer )
    {
        if( aOxidizer == "O2" || aOxidizer == "LOX" || aOxidizer == "LO2" || aOxidizer == "OXYGEN" )
        {
            return Oxidizer::LOX;
        }
        if( aOxidizer == "H2O2" || aOxidizer == "HTP" )
        {
            return Oxidizer::H2O2;
        }
        else if( string_to_upper( aOxidizer ) == "AIR" )
        {
            return Oxidizer::AIR;
        }
        else
        {
            BELFEM_ERROR( false, "Unknown Oxidizer : %s.", aOxidizer.c_str() );
            return Oxidizer::UNDEFINED;
        }
    }

//------------------------------------------------------------------------------
}