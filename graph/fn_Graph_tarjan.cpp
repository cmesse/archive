//
// Pocket Detection for Mesh Cleanup
// Implements Tarjan's algorithm with extensions for unicyclic component detection
//


#include <algorithm>
#include <utility>       // For std::pair
#include <limits>


#include "fn_Graph_tarjan.hpp"
#include "cl_Queue.hpp"
#include "cl_Map.hpp"    // BELFEM's Map class


namespace belfem
{
    namespace graph
    {
        ArticulationPointInfo::ArticulationPointInfo(
                Vertex * aVertex,
                Cell< Vertex * > & aGraph,
                Cell< Cell< Vertex * > > aSeparatedComponents,
                ArticulationPointWork & aWork ) :
                mVertex( aVertex ),
                mGraph( aGraph ),
                mIsUnicyclic( new DynamicBitset( aSeparatedComponents.size()) ),
                mVisited( aWork.mVisited ),
                mInComponent( aWork.mInComponent ),
                mParent( aWork.mParent ),
                mDepth( aWork.mDepth )
            {
                index_t tSize = aSeparatedComponents.size();
                mComponentSizes.set_size( tSize, 0 );
                mComponentEdgeCounts.set_size( tSize, 0 );
                mCycleLength.set_size( tSize, 0 );
                mSeparatedComponents = std::move( aSeparatedComponents );
                DynamicBitset * tContains = aWork.mContains;


                for ( index_t k=0 ; k<tSize; ++k )
                {
                    Cell< Vertex * > & tComp = mSeparatedComponents( k );

                    // copy the compinent size
                    mComponentSizes( k ) = tComp.size();

                    // set contain checks
                    tContains->reset();
                    for ( Vertex * tVertex : tComp )
                    {
                        tContains->set( tVertex->index() );
                    }

                    // count edges
                    index_t tCount = 0 ;
                    for ( Vertex * tVertex : tComp )
                    {
                        for ( uint v=0; v<tVertex->number_of_vertices(); ++v )
                        {
                            if ( tContains->test( tVertex->vertex( v )->index() ) )
                            {
                                ++tCount;
                            }
                        }
                    }

                    mComponentEdgeCounts( k ) = tCount / 2 ;

                    // Check for unicyclic: exactly one cycle if edges == vertices (connected component)
                    bool tIsUnicyclic = ( mComponentEdgeCounts( k ) == tComp.size() );

                    if ( tIsUnicyclic )
                    {
                        mIsUnicyclic->set( k );

                        auto [ tHasCycle, tLen ] = this->detect_cycle_properties( k );
                        BELFEM_ASSERT( tHasCycle, "Expected cycle in unicyclic component" );
                        mCycleLength( k ) = tLen;
                    }
                }

            }

        /**
         * Find articulation points using Tarjan's algorithm
         * @param aGraph The graph represented as Cell of Vertex pointers
         * @return Cell containing all articulation points
         */
        /*Cell< Vertex * >
        find_articulation_points( Cell< Vertex * > & aGraph )
        {
            if ( aGraph.size() == 0 )
            {
                return Cell< Vertex * >();
            }

            const index_t tN = aGraph.size();

            // Ensure all vertices have proper indices
            index_t tIndex = 0 ;
            for ( Vertex * tVertex : aGraph )
            {
                tVertex->set_index( tIndex++ );
            }

            // Use DynamicBitset for efficient boolean tracking
            DynamicBitset tVisited( tN );
            DynamicBitset tIsArticulation( tN );

            // Tarjan's algorithm data structures
            Vector< index_t > tDiscoveryTime( tN, gNoIndex );
            Vector< index_t > tLowLink( tN, gNoIndex );
            Vector< index_t > tParent( tN, gNoIndex );

            index_t tTimer = 0;

            // DFS lambda with Tarjan's logic
            auto tDFS = [ & ]( auto&& aSelf, index_t aV ) -> void
            {
                index_t tChildren = 0;
                tVisited.set( aV );
                tDiscoveryTime( aV ) = tLowLink( aV ) = tTimer++;

                Vertex * tVertex = aGraph( aV );

                // Process all adjacent vertices
                for ( uint k = 0; k < tVertex->number_of_vertices(); ++k )
                {
                    Vertex * tNeighbor = tVertex->vertex( k );
                    index_t tU = tNeighbor->index();

                    if ( !tVisited.test( tU ) )
                    {
                        ++tChildren;
                        tParent( tU ) = aV;

                        // Recursive DFS
                        aSelf( aSelf, tU );

                        // Update low link value
                        tLowLink( aV ) = std::min( tLowLink( aV ), tLowLink( tU ) );

                        // Non-root articulation point check
                        if ( tParent( aV ) != gNoIndex && tLowLink( tU ) >= tDiscoveryTime( aV ) )
                        {
                            tIsArticulation.set( aV );
                        }
                    }
                    else if ( tU != tParent( aV ) )
                    {
                        // Back edge - update low link
                        tLowLink( aV ) = std::min( tLowLink( aV ), tDiscoveryTime( tU ) );
                    }
                }

                // Root articulation point check (after processing all children)
                if ( tParent( aV ) == gNoIndex && tChildren > 1 )
                {
                    tIsArticulation.set( aV );
                }
            };

            // Run DFS from each unvisited vertex (handles disconnected graphs)
            for ( index_t i = 0; i < tN; ++i )
            {
                if ( !tVisited.test( i ) )
                {
                    tDFS( tDFS, i );
                }
            }

            // Extract articulation point indices efficiently
            Cell< index_t > tArticulationIndices;
            tIsArticulation.where( tArticulationIndices, false ); // false = dense allocation

            // Convert indices to vertex pointers
            Cell< Vertex * > tArticulationPoints;
            tArticulationPoints.reserve( tArticulationIndices.size() );

            for ( index_t idx : tArticulationIndices )
            {
                tArticulationPoints.push( aGraph( idx ) );
            }

            return tArticulationPoints;
        }*/



        /**
         * Extended articulation point finder with component analysis
         * @param aGraph The graph to analyze
         * @return Detailed information about each articulation point
         */
        Cell< ArticulationPointInfo * >
        find_articulation_points_with_components( Cell< Vertex * > & aGraph )
        {
            Cell< ArticulationPointInfo * > aAPInfo;

            if ( aGraph.size() == 0 )
            {
                return aAPInfo;
            }

            // Ensure all vertices have proper indices
            index_t tIndex = 0;
            for ( Vertex * tVertex : aGraph )
            {
                tVertex->set_index( tIndex++ );
            }

            const index_t tN = aGraph.size();

            // Use DynamicBitset for efficient boolean tracking
            DynamicBitset tVisited( tN );

            // Tarjan's algorithm data structures
            Vector< index_t > tDiscoveryTime( tN, gNoIndex );
            Vector< index_t > tLowLink( tN, gNoIndex );
            Vector< index_t > tParent( tN, gNoIndex );
            index_t tTimer = 0;

            // the work object
            ArticulationPointWork tWork( tN );

            // DFS lambda that returns the subtree vertices
            auto tDFS = [ & ]( auto&& aSelf, Cell< Vertex * > & aGraph, index_t aV ) -> Cell< Vertex * >
            {
                Vertex * tVertex = aGraph( aV );
                tVisited.set( aV );
                tDiscoveryTime( aV ) = tLowLink( aV ) = tTimer++;


                Cell< Cell< Vertex * > > tChildSubtrees; // Collect all direct child subtrees
                Cell< Cell< Vertex * > > tSeparated; // Separated components for this vertex
                index_t tChildren = 0;

                Cell< Vertex * > aSubtree;
                aSubtree.push( tVertex ); // Start with self

                // Process all adjacent vertices
                for ( uint k = 0; k < tVertex->number_of_vertices(); ++k )
                {
                    Vertex * tNeighbor = tVertex->vertex( k );
                    index_t tU = tNeighbor->index();

                    if ( tU == tParent( aV ) ) continue; // Skip parent

                    if ( !tVisited.test( tU ) )
                    {
                        ++tChildren;
                        tParent( tU ) = aV;

                        // Recursive DFS, get child subtree
                        Cell< Vertex * > tChildSubtree = aSelf( aSelf, aGraph, tU );

                        // Update low link
                        tLowLink( aV ) = std::min( tLowLink( aV ), tLowLink( tU ) );

                        // Check for separated component (non-root)
                        if ( tParent( aV ) != gNoIndex && tLowLink( tU ) >= tDiscoveryTime( aV ) )
                        {
                            tSeparated.push( tChildSubtree ); // Copy (cheap for pointers)
                        }

                        // Always add to full subtree and child subtrees
                        append( aSubtree, tChildSubtree ); // Append
                        tChildSubtrees.push( tChildSubtree ); // Copy
                    }
                    else
                    {
                        // Back edge - update low link
                        tLowLink( aV ) = std::min( tLowLink( aV ), tDiscoveryTime( tU ) );
                    }
                }

                // Handle root case separately after processing all children
                if ( tParent( aV ) == gNoIndex && tChildren > 1 )
                {
                    tSeparated = tChildSubtrees; // Copy all child subtrees as separated
                }

                if ( !tSeparated.empty() )
                {
                    // create a new info point
                    aAPInfo.push( new ArticulationPointInfo( tVertex, aGraph, tSeparated, tWork ) );
                }

                return aSubtree ;
            };

            for ( index_t i = 0; i < tN; ++i )
            {
                if ( ! tVisited.test( i ) )
                {
                    tDFS( tDFS, aGraph, i );
                }
            }

            return aAPInfo;
        }

        std::pair< bool, index_t >
        ArticulationPointInfo::detect_cycle_properties( const index_t aIndex )
        {
            const Cell< Vertex * > & tComp = mSeparatedComponents( aIndex );
            if ( tComp.size() < 3 ) return { false, 0 };

            // Reset work for this component
            mVisited->reset();
            mParent.fill( gNoIndex );
            mDepth.fill( 0 );

            // Set in_component for membership (global)
            mInComponent->reset();
            for ( Vertex * tV : tComp )
            {
                mInComponent->set( tV->index() );
            }

            bool tCycleFound = false;
            index_t tCycleLength = 0;

            auto tDFS = [ & ]( auto&& aSelf, index_t aV, index_t aCurrentDepth ) -> bool
            {
                mVisited->set( aV );
                mDepth( aV ) = aCurrentDepth;
                Vertex * tVertex = mGraph( aV );  // Global lookup
                for ( uint k = 0; k < tVertex->number_of_vertices(); ++k )
                {
                    Vertex * tNeighbor = tVertex->vertex( k );
                    index_t tNeighborIdx = tNeighbor->index();
                    if ( ! mInComponent->test( tNeighborIdx ) ) continue;
                    if ( ! mVisited->test( tNeighborIdx ) ) {
                        mParent( tNeighborIdx ) = aV;
                        if ( aSelf( aSelf, tNeighborIdx, aCurrentDepth + 1 ) )
                        {
                            return true;
                        }
                    }
                    else if ( mParent( aV ) != tNeighborIdx && ! tCycleFound )
                    {
                        tCycleFound = true;
                        tCycleLength = mDepth( aV ) - mDepth( tNeighborIdx ) + 1;
                        return true;
                    }
                }
                return false;
            };

            // Start from global indices of component vertices
            for ( Vertex * tV : tComp )
            {
                index_t g = tV->index();
                if ( ! mVisited->test( g ) )
                {
                    if ( tDFS( tDFS, g, 0 ) )
                    {
                        break;
                    }
                }
            }

            return { tCycleFound, tCycleLength };
        }

        /**
         * Detect all pockets in the graph
         * @param aGraph The graph to analyze
         * @param aMaxPocketSize Maximum size for a valid pocket
         * @param aRequireUnicyclic Whether to require exactly one cycle
         * @return Collection of detected pockets with metadata
         */
        Cell< PocketInfo >
        detect_pockets_with_info( Cell< Vertex * > & aGraph,
                                 const index_t aMaxPocketSize,
                                 const bool aRequireUnicyclic )
        {
            Cell< PocketInfo > tPockets;

            auto tAPInfo = find_articulation_points_with_components( aGraph );

            for ( ArticulationPointInfo * tInfo : tAPInfo )
            {
                for ( index_t i = 0; i < tInfo->size(); ++i )
                {
                    // Check pocket criteria
                    if ( tInfo->component_size( i ) <= aMaxPocketSize &&
                         ( !aRequireUnicyclic || tInfo->unicyclic_test( i ) ) )
                    {
                        PocketInfo tPocket;
                        tPocket.mNeckVertex = tInfo->vertex() ;
                        tPocket.mPocketVertices = tInfo->separated_components( i );
                        tPocket.mSize = tInfo->component_size( i );
                        tPocket.mEdgeCount = tInfo->component_edge_count( i );
                        tPocket.mIsUnicyclic = tInfo->unicyclic_test( i );
                        tPocket.mCycleLength = tInfo->cycle_length( i );

                        tPockets.push( std::move( tPocket ) );
                    }
                }

                // tidy up
                delete tInfo;
            }

            // Sort pockets by quality score for prioritized removal
            std::sort( tPockets.vector_data().begin(),
                      tPockets.vector_data().end(),
                      []( const PocketInfo & a, const PocketInfo & b )
                      {
                          return a.quality_score() > b.quality_score();
                      });

            return tPockets;
        }

        /**
         * Remove detected pockets from the graph
         * @param aGraph The graph to modify
         * @param aPockets The pockets to remove
         * @return Number of vertices removed
         */
        index_t
        remove_pockets( Cell< Vertex * > & aGraph,
                       const Cell< PocketInfo > & aPockets )
        {
            if ( aPockets.size() == 0 )
            {
                return 0;
            }

            index_t tRemovedCount = 0;
            DynamicBitset tToRemove( aGraph.size() );

            // Mark all pocket vertices for removal
            for ( const auto & tPocket : aPockets )
            {
                for ( Vertex * tV : tPocket.mPocketVertices )
                {
                    tToRemove.set( tV->index() );
                    ++tRemovedCount;
                }
            }

            // Create new graph without pocket vertices
            Cell< Vertex * > tNewGraph;
            tNewGraph.reserve( aGraph.size() - tRemovedCount );

            for ( index_t i = 0; i < aGraph.size(); ++i )
            {
                if ( !tToRemove.test( i ) )
                {
                    tNewGraph.push( aGraph( i ) );
                }
            }

            // Update the original graph
            aGraph = std::move( tNewGraph );

            // Re-index vertices
            for ( index_t i = 0; i < aGraph.size(); ++i )
            {
                aGraph( i )->set_index( i );
            }

            return tRemovedCount;
        }

    } // namespace graph
} // namespace belfem