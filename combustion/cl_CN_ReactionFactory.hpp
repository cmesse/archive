//
// Created by Christian Messe on 26.09.19.
//

#ifndef BELFEM_CL_CN_REACTIONFACTORY_HPP
#define BELFEM_CL_CN_REACTIONFACTORY_HPP

#include "cl_CN_Entry.hpp"
#include "cl_CN_Reaction.hpp"
#include "cl_Map.hpp"

namespace belfem
{
    namespace combustion
    {
        class Scheme;

        class ReactionFactory
        {
//------------------------------------------------------------------------------

            Scheme             * mScheme;
            Map< string, uint >  mSpeciesMap;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            ReactionFactory( Scheme * aScheme );

//------------------------------------------------------------------------------

            ~ReactionFactory() = default;

//------------------------------------------------------------------------------

            Reaction *
            create_reaction( const Entry * aEntry );

//------------------------------------------------------------------------------
        };
    }
}
#endif //BELFEM_CL_CN_REACTIONFACTORY_HPP
