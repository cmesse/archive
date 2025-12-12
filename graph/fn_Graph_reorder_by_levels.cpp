// fn_Graph_reorder_by_levels.cpp

#include <algorithm>
#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Queue.hpp"
#include "cl_DynamicBitset.hpp"
#include "op_Graph_Vertex_Index.hpp"
#include "fn_Graph_symrcm.hpp"

namespace belfem
{
    namespace graph
    {
//------------------------------------------------------------------------------

        // BFS to compute distances from a set of source vertices
        real
        compute_distances(
            Graph & aGraph,
            Graph & aSources,
            Vector< real > & aDistance )
        {
            index_t tNumVertices = aGraph.size();
            aDistance.set_size( tNumVertices, BELFEM_QUIET_NAN );

            for ( Vertex * tVertex : aSources )
            {
                tVertex->unflag() ;
            }

            Queue< Vertex * > tQueue;

            // Initialize sources at distance 0
            for ( Vertex * tVertex : aSources )
            {
                aDistance( tVertex->index() ) = 0.0;
                tQueue.push( tVertex );
                tVertex->flag();
            }

            // BFS
            real aMaxDistance = 0.0 ;

            while ( !tQueue.empty() )
            {
                Vertex * tVertex = tQueue.pop();
                real tCurrentDist = aDistance( tVertex->index() );

                // Visit all neighbors
                for ( uint k = 0; k < tVertex->number_of_vertices(); ++k )
                {
                    Vertex * tNeighbor = tVertex->vertex( k );
                    if ( ! tNeighbor->is_flagged() )
                    {
                        real tNextDist = tCurrentDist + 1.0;
                        aMaxDistance = std::max( aMaxDistance, tNextDist );
                        aDistance( tNeighbor->index() ) = tNextDist;
                        tNeighbor->flag();
                        tQueue.push( tNeighbor );
                    }
                }
            }

            return aMaxDistance;
        }

//------------------------------------------------------------------------------

        // Compute pseudo-temperature from distances to both boundaries
        // Returns value in [0, 1] similar to Poisson solution
        void
        compute_pseudo_temperature(
            const Vector< real > & aDist0,
            const Vector< real > & aDist1,
            Vector< real > & aPseudoT )
        {
            index_t tNumVertices = aDist0.length();
            aPseudoT.set_size( tNumVertices, 0.5 );

            for ( index_t k =0 ; k<tNumVertices; ++k )
            {
                real d0 = aDist0( k );
                real d1 = aDist1( k );

                if ( d0 == BELFEM_QUIET_NAN && d1 == BELFEM_QUIET_NAN )
                {
                    // Disconnected vertex - place in middle
                    aPseudoT( k ) = 0.5;
                }
                else if ( d0 == BELFEM_QUIET_NAN )
                {
                    // Only reachable from Γ₁
                    aPseudoT( k ) = 1.0;
                }
                else if ( d1 == BELFEM_QUIET_NAN )
                {
                    // Only reachable from Γ₀
                    aPseudoT( k ) = 0.0;
                }
                else if ( d0 == 0.0 && d1 == 0.0 )
                {
                    aPseudoT( k ) = 0.0;
                }
                else
                {
                    // Reachable from both - interpolate
                    aPseudoT( k ) = d0 / ( d0 + d1 );
                }
            }
        }

//------------------------------------------------------------------------------

        // Build subgraph for vertices in a given index range
        void
        build_level_subgraph(
            Graph & aGraph,
            Cell< index_t > & aLevelIndices,
            Graph & aSubgraph,
            Map< index_t, index_t > & aGlobalToLocal )
        {
            index_t tNumLocal = aLevelIndices.size();
            aSubgraph.set_size( tNumLocal, nullptr );
            aGlobalToLocal.clear();

            // Create local vertices and mapping
            for ( index_t k = 0; k < tNumLocal; ++k )
            {
                index_t tGlobal = aLevelIndices( k );
                aGlobalToLocal[ tGlobal ] = k;

                aSubgraph( k ) = new Vertex();
                aSubgraph( k )->set_id( aGraph( tGlobal )->id() );
                aSubgraph( k )->set_index( k );
            }

            // Build edges (only between vertices in this level)
            for ( index_t k = 0; k < tNumLocal; ++k )
            {
                index_t tGlobal = aLevelIndices( k );
                Vertex * tOriginal = aGraph( tGlobal );

                // Count neighbors in this level
                uint tCount = 0;
                for ( uint n = 0; n < tOriginal->number_of_vertices(); ++n )
                {
                    index_t tNeighborGlobal = tOriginal->vertex( n )->index();
                    if ( aGlobalToLocal.key_exists( tNeighborGlobal ) )
                    {
                        ++tCount;
                    }
                }

                aSubgraph( k )->init_vertex_container( tCount );

                for ( uint n = 0; n < tOriginal->number_of_vertices(); ++n )
                {
                    index_t tNeighborGlobal = tOriginal->vertex( n )->index();
                    if ( aGlobalToLocal.key_exists( tNeighborGlobal ) )
                    {
                        index_t tNeighborLocal = aGlobalToLocal( tNeighborGlobal );
                        aSubgraph( k )->insert_vertex( aSubgraph( tNeighborLocal ) );
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        delete_subgraph( Graph & aSubgraph )
        {
            for ( Vertex * tVertex : aSubgraph )
            {
                delete tVertex;
            }
            aSubgraph.clear();
        }

//------------------------------------------------------------------------------

        // Refine ordering within a level using matching
        // Matched pairs get consecutive indices
        void
        refine_level_with_matching(
            Graph & aGraph,
            Cell< index_t > & aLevelIndices,
            index_t & aCurrentIndex )
        {
            index_t tNumInLevel = aLevelIndices.size();

            if ( tNumInLevel == 0 )
            {
                return;
            }

            if ( tNumInLevel == 1 )
            {
                aGraph( aLevelIndices( 0 ) )->set_index( aCurrentIndex++ );
                return;
            }

            // Build subgraph for this level
            Graph tSubgraph;
            Map< index_t, index_t > tGlobalToLocal;
            build_level_subgraph( aGraph, aLevelIndices, tSubgraph, tGlobalToLocal );

            // Find matching
            Cell< index_t > tMatch;
            max_cardinality_matching( tSubgraph, tMatch );

            // Assign indices: matched pairs get consecutive indices
            DynamicBitset tAssigned( tNumInLevel );

            for ( index_t k = 0; k < tNumInLevel; ++k )
            {
                if ( !tAssigned.test( k ) )
                {
                    index_t tGlobal = aLevelIndices( k );
                    aGraph( tGlobal )->set_index( aCurrentIndex++ );
                    tAssigned.set( k );

                    // If matched, assign partner next
                    if ( tMatch( k ) != gNoIndex && !tAssigned.test( tMatch( k ) ) )
                    {
                        index_t tPartnerGlobal = aLevelIndices( tMatch( k ) );
                        aGraph( tPartnerGlobal )->set_index( aCurrentIndex++ );
                        tAssigned.set( tMatch( k ) );
                    }
                }
            }

            // Cleanup
            delete_subgraph( tSubgraph );
        }

//------------------------------------------------------------------------------

        void
        reorder_by_levels(
            Graph & aGraph,
            Graph & aSinks,
            Graph & aSources,
            Map< id_t, real > * aField,
            const bool aSort )
        {
            if ( aGraph.size() == 0 )
            {
                return;
            }

            index_t tNumVertices = aGraph.size();

            // Ensure consistent indexing
            for ( index_t k = 0; k < tNumVertices; ++k )
            {
                aGraph( k )->set_index( k );
            }

            // Step 1: Compute distances from both boundaries
            Vector< real > tDist0;
            Vector< real > tDist1;
            real tMaxA = compute_distances( aGraph, aSinks, tDist0 );
            real tMaxB = compute_distances( aGraph, aSources, tDist1 );

            real tMaxDist = tMaxA > tMaxB ? tMaxA : tMaxB;

            // Step 2: Compute pseudo-temperature T ∈ [0, 1]
            Vector< real > tPseudoT;
            compute_pseudo_temperature( tDist0, tDist1, tPseudoT );

            if ( aField != nullptr )
            {
                Map< id_t, real > & tField = *aField;
                tField.clear();

                for ( index_t k = 0; k < tNumVertices; ++k )
                {
                    tField[ aGraph( k )->id() ] = tPseudoT( k );
                }
            }

            // Use finer binning for smoother result
            index_t tNumLevels = std::max( static_cast< index_t>(tMaxDist + 1.),  static_cast< index_t>(10) );

            // Step 4: Bin vertices by pseudo-temperature
            Cell< Cell< index_t > > tLevelBins( tNumLevels, {} );

            for ( index_t k = 0; k < tNumVertices; ++k )
            {
                // Map T ∈ [0,1] to level ∈ [0, tNumLevels-1]
                index_t tLevel = static_cast< index_t >( tPseudoT( k ) * ( tNumLevels - 1 ) + 0.5 );
                tLevel = std::min( tLevel, tNumLevels - 1 );
                tLevelBins( tLevel ).push( k );
            }

            // Step 5: Process each level with matching refinement
            index_t tCurrentIndex = 0;

            for ( index_t L = 0; L < tNumLevels; ++L )
            {
                refine_level_with_matching( aGraph, tLevelBins( L ), tCurrentIndex );
            }

            // Step 6: Sort graph by new indices
            if ( aSort ) sort( aGraph, opVertexIndex );
        }

//------------------------------------------------------------------------------
    }
}