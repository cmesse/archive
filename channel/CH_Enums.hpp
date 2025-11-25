//
// Created by Christian Messe on 21.11.20.
//

#ifndef BELFEM_CH_ENUMS_HPP
#define BELFEM_CH_ENUMS_HPP

#include "typedefs.hpp"

namespace belfem
{
    enum class ChannelType
    {
        CoolingChannel,
        CombustionChamber
    };

//------------------------------------------------------------------------------

    enum class EngineType
    {
        Cylinder,
        Hyperboloid,
        UNDEFINED
    };

//------------------------------------------------------------------------------
    namespace channel
    {
//------------------------------------------------------------------------------

        enum class BoundaryLayerMethod
        {
            Messe,
            Eckert,
            Bartz,              // hotgas only
            Pizzarelli,         // methane only
            LebedinskyKalmykov, // methane only
            vanDriest,          // Nozzle only
            UNDEFINED
        };

        BoundaryLayerMethod
        boundary_layer_method( const string & aString );

    }

//------------------------------------------------------------------------------

        string
        to_string( const channel::BoundaryLayerMethod & aMethod );

//------------------------------------------------------------------------------
}

#endif //BELFEM_CH_ENUMS_HPP
