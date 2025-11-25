//
// Created by Christian Messe on 26.09.19.
//

#include "cl_CN_ReactionFactory.hpp"
#include "cl_CN_Scheme.hpp"
#include "fn_CN_molars_from_string.hpp"
#include "cl_Gas.hpp"

#include "cl_CN_Reaction_Arrhenius.hpp"
#include "cl_CN_Reaction_Duplicate.hpp"
#include "cl_CN_Reaction_Lindemann.hpp"
#include "cl_CN_Reaction_Troe.hpp"

namespace belfem
{
    namespace combustion
    {
//------------------------------------------------------------------------------

        ReactionFactory::ReactionFactory( Scheme * aScheme ) :
            mScheme( aScheme )
        {
            // get reference to combustion gas
            Gas * tGas = aScheme->combgas();

            // create the species map
            for ( uint k = 0; k < aScheme->number_of_reacting_species(); ++k )
            {
                mSpeciesMap[ tGas->data( k )->label() ] = k;
            }

        }

//------------------------------------------------------------------------------

        Reaction *
        ReactionFactory::create_reaction( const Entry * aEntry )
        {
            Cell< string > tEductLabels;
            Vector< uint > tEductIndices;
            Vector< real > tEductMolars;
            Cell< string > tProductLabels;
            Vector< uint > tProductIndices;
            Vector< real > tProductMolars;
            Vector< uint > tThirdBodyIndices;

            bool tHasThirdBody = false;

            // interpret reaction string
            molars_from_string(
                    aEntry->reaction(),
                    tEductLabels,
                    tEductMolars,
                    tProductLabels,
                    tProductMolars,
                    tHasThirdBody );

            // create index list for educts
            tEductIndices.set_size( tEductLabels.size() );

            for( uint k=0; k<tEductLabels.size(); ++k )
            {
                tEductIndices( k ) = mSpeciesMap( tEductLabels( k ) );
            }

            // create index list for products
            tProductIndices.set_size( tProductLabels.size() );

            for( uint k=0; k<tProductLabels.size(); ++k )
            {
                tProductIndices( k ) = mSpeciesMap( tProductLabels( k ) );
            }

            Reaction * aReaction = nullptr;

            // make sure that this is an active reaction
            BELFEM_ERROR( aEntry->is_active(),
                "can't create Reaction object from deactivated entry ( duplicate of %s)",
                aEntry->reaction().c_str() );

            if( aEntry->is_duplicate() )
            {

                BELFEM_ERROR( aEntry->duplicate().length() > 0,
                    "Reaction %s is marked as duplicate but has no duplicate assigned",
                    aEntry->reaction().c_str() );

                // create duplicate reaction
                aReaction = new Reaction_Duplicate( *mScheme,
                                                    tEductIndices,
                                                    tEductMolars,
                                                    tProductIndices,
                                                    tProductMolars,
                                                    tHasThirdBody,
                                                    aEntry->coeffs(),
                                                    aEntry->duplicate() );
            }
            else if( aEntry->has_troe() )
            {
                BELFEM_ERROR( tHasThirdBody,
                             "Reaction %s is supposed to be a Troe type reaction, but does not seem to have any inert partners",
                             aEntry->reaction().c_str() );

                BELFEM_ERROR( aEntry->has_low(),
                        "Reaction %s has TROE entry but no LOW entry.",
                              aEntry->reaction().c_str() );

                aReaction = new Reaction_Troe( *mScheme,
                                               tEductIndices,
                                               tEductMolars,
                                               tProductIndices,
                                               tProductMolars,
                                               tHasThirdBody,
                                               aEntry->low(),
                                               aEntry->coeffs(),
                                               aEntry->troe() );
            }
            else if ( aEntry->has_low() )
            {
                BELFEM_ERROR( tHasThirdBody,
                    "Reaction %s is supposed to be a Lindemann type reaction, but does not seem to have any inert partners",
                    aEntry->reaction().c_str() );

                aReaction = new Reaction_Lindemann( *mScheme,
                                                    tEductIndices,
                                                    tEductMolars,
                                                    tProductIndices,
                                                    tProductMolars,
                                                    tHasThirdBody,
                                                    aEntry->low(),
                                                    aEntry->coeffs() );
            }
            else
            {
                // create Arrhenius reaction
                aReaction = new Reaction_Arrhenius(
                        *mScheme,
                        tEductIndices,
                        tEductMolars,
                        tProductIndices,
                        tProductMolars,
                        tHasThirdBody,
                        aEntry->coeffs() );
            }

            // create index list for third body indices
            const Cell< string > & tThirdBodyLabels = aEntry->third_body_species();

            if( tThirdBodyLabels.size() > 0 )
            {
                tThirdBodyIndices.set_size( tThirdBodyLabels.size() );

                for( uint k=0; k< tThirdBodyLabels.size(); ++k )
                {
                    tThirdBodyIndices( k ) = mSpeciesMap( tThirdBodyLabels( k ) );
                }

                aReaction->set_third_body(
                        tThirdBodyIndices,
                        aEntry->third_body_weights() );
            }

            return aReaction;
        }

//------------------------------------------------------------------------------
    }
}
