//
// Created by Christian Messe on 04.01.21.
//
#include "assert.hpp"
#include "CH_Enums.hpp"
#include "stringtools.hpp"

namespace belfem
{
    string
    to_string( const channel::BoundaryLayerMethod & aMethod )
    {
        switch ( aMethod )
        {
            case ( channel::BoundaryLayerMethod::Messe ) :
            {
                return "Messe";
            }
            case ( channel::BoundaryLayerMethod::Bartz ) :
            {
                return "Bartz";
            }
            case ( channel::BoundaryLayerMethod::Eckert ) :
            {
                return "Eckert";
            }
            case( channel::BoundaryLayerMethod::Pizzarelli ) :
            {
                return "Pizzarelli" ;
            }
            case( channel::BoundaryLayerMethod::LebedinskyKalmykov ) :
            {
                return "Lebedinsky/Kalmykov" ;
            }
            /*case( channel::BoundaryLayerMethod::vanDriest ) :
            {
                return "van Driest" ;
            }*/
            default:
            {
                BELFEM_ERROR( false, "Invalid Boundary layer Method" );
                return "unknown";
            }
        }
    }

//------------------------------------------------------------------------------
    namespace channel
    {
        BoundaryLayerMethod
        boundary_layer_method( const string & aString )
        {
            // convert string to lower case
            string tString = string_to_lower( aString );

            if ( tString == "messe" )
            {
                return BoundaryLayerMethod::Messe;
            }
            else if ( tString == "bartz" )
            {
                return BoundaryLayerMethod::Bartz;
            }
            else if ( tString == "eckert" )
            {
                return BoundaryLayerMethod::Eckert;
            }
            else if ( tString == "pizzarelli" )
            {
                return BoundaryLayerMethod::Pizzarelli ;
            }
            else if ( tString == "lebedinskykalmykov" ||
                      tString == "lebedinsky kalmykov" ||
                      tString == "lebedinsky-kalmykov" ||
                      tString == "lebedinsky/kalmykov" )
            {
                return BoundaryLayerMethod::LebedinskyKalmykov ;
            }
            else if ( tString == "van driest" || "vandriest" )
            {
                return BoundaryLayerMethod::vanDriest ;
            }
            else
            {
                BELFEM_ERROR( false, "Unknown boundary layer method" );
                return BoundaryLayerMethod::UNDEFINED;
            }
        }
//------------------------------------------------------------------------------
    }
}