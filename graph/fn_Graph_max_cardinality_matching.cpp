/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * Implementation of the Micali-Vazirani maximum cardinality matching algorithm
 * for general graphs.
 *
 * ALGORITHM OVERVIEW:
 * The Micali-Vazirani algorithm finds a maximum cardinality matching in an
 * undirected graph in O(√V × E) time. It works by iteratively finding
 * augmenting paths of increasing length, using a level-based BFS structure
 * and handling odd-length cycles (blossoms) through contraction.
 *
 * KEY DATA STRUCTURES:
 * - CSR (Compressed Sparse Row) format for cache-efficient adjacency access
 * - Edge status array aligned with CSR for O(1) lookups
 * - Level arrays (mEvenLevel, mOddLevel) for alternating forest structure
 * - Free-list based stacks to avoid dynamic allocation in inner loops
 * - DynamicBitset for efficient set membership tests
 *
 * IMPLEMENTATION NOTES:
 * - Consistent 0-based indexing throughout
 * - gNoIndex constant used for undefined/uninitialized values
 * - Manual stack in find_path prevents overflow on deep recursion
 * - Initial greedy matching provides better starting point
 */

#include <algorithm>
#include <stack>
#include <cstdio>
#include <unistd.h>

#include "fn_Graph_max_cardinality_matching.hpp"
#include "cl_DynamicBitset.hpp"
#include "cl_Cell.hpp"
#include "cl_Map.hpp"
#include "cl_Queue.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    namespace graph
    {
        // Graph typedef for cleaner signatures
        typedef Cell< Vertex * > Graph;

//------------------------------------------------------------------------------

        /**
         * @brief Implementation of the Micali-Vazirani maximum matching algorithm
         *
         * This class encapsulates all data structures and methods needed to find
         * a maximum cardinality matching in an undirected graph.
         */
        class MicaliVazirani
        {
            // Graph properties
            index_t mNumVertices;       // Number of vertices in the graph
            index_t mNumEdges;           // Number of edges in the graph
            index_t mInfinity;           // Sentinel value for "infinity" (max(V,E)+1)
            index_t mSearchLevelLimit;   // Maximum search depth (V/2 + 1)

            Graph & mGraph;              // Reference to input graph

            // CSR (Compressed Sparse Row) adjacency storage
            // Provides O(1) neighbor access with better cache locality than adjacency lists
            Vector< index_t > mAdjStart;   // Start index in mAdjList for each vertex
            Vector< index_t > mAdjList;    // Flattened array of all neighbors
            Map< key_t, index_t > mEdgeMap ;

            // Edge status tracking (aligned with CSR for O(1) access)
            // Stores usage flags for each edge during path search
            Vector< index_t > mEdgeStatus; // Status codes for edges

            // Core algorithm arrays (all 0-based indexing)
            Vector< index_t > mMatch;      // Current matching: mMatch[v] = matched vertex or gNoIndex
            Vector< index_t > mEvenLevel;  // BFS level when vertex reached as even vertex
            Vector< index_t > mOddLevel;   // BFS level when vertex reached as odd vertex
            Vector< index_t > mBlossom;    // Blossom ID containing this vertex (gNoIndex if none)
            Vector< index_t > mForest;     // Parent in alternating forest for path reconstruction
            Vector< index_t > mPath;       // Temporary path storage during augmentation
            Vector< int >     mLR;         // Left/Right marker (+/-mNumCalls) for DCV detection
            Vector< index_t > mPrdctr;     // Predecessor count for topological erase

            // Blossom data structures (indexed by blossom ID)
            Vector< index_t > mBase;       // Base vertex of each blossom (DCV)
            Vector< index_t > mBstar;      // Current tip of blossom during expansion
            Vector< index_t > mPeakLeft;   // Left peak vertex (bridge endpoint)
            Vector< index_t > mPeakRight;  // Right peak vertex (bridge endpoint)

            // Free-list based stacks (avoid allocation in inner loops)
            // Each stack uses three arrays: Stack0 (next pointer), Stack1 (data), Index (head)
            Vector< index_t > mPred0, mPred1, mPIndex;     // Predecessor stack for path finding
            Vector< index_t > mDerp0, mDerp1, mDIndex;     // Dependency stack for topological erase
            Vector< index_t > mAnom0, mAnom1, mAIndex;     // Anomaly stack for cross-level edges
            Vector< index_t > mBridge0, mBridge1, mBridge2, mBIndex; // Bridge triple stack

            // Bitsets for efficient set operations
            DynamicBitset mMark;         // Marks erased vertices during topological cleanup
            DynamicBitset mVisited;      // Tracks visited vertices during path search

            // Free list management and call counting
            index_t mPFree, mDFree, mAFree, mBFree; // Free list heads for each stack
            index_t mNumCalls;           // Call counter for unique LR marking per blsaug call

        public:
            /**
             * @brief Constructs the algorithm with a graph reference
             * @param aGraph Reference to graph (will set vertex indices internally)
             */
            MicaliVazirani( Graph & aGraph );

            /**
             * @brief Runs the matching algorithm and returns the result
             * @param aMatch Output array: aMatch[v] = matched vertex or gNoIndex
             * @return Cardinality of the maximum matching
             */
            index_t run( Cell< index_t > & aMatch );

        private:
            // CSR construction and edge status management
            void build_csr();

            inline index_t find_edge_index( index_t aU, index_t aV ) const;
            inline index_t get_edge_status( index_t aU, index_t aV ) const;
            inline void add_edge_status( index_t aU, index_t aV, index_t aCode );

            // Free-list stack operations
            void add_to_stack( Vector< index_t > & aStack0,
                           Vector< index_t > & aStack1,
                           Vector< index_t > & aIndex,
                           index_t & aFree,
                           index_t aU, index_t aV );

            void add_bridge( index_t aU, index_t aV, index_t aBr );

            // Core algorithm components
            index_t ancest( index_t aV, bool aCheckUnused, index_t & aIndex );
            void bastar( index_t & aV, index_t & aU );
            void find_path( index_t aHigh, index_t aLow, index_t aB, int aJob );
            bool blsaug( index_t aW1, index_t aW2,
                        index_t aSearchLevel, index_t & aBlossomCounter,
                        index_t & aCardinality );
            void search( index_t & aCardinality );
            void compute_initial_matching( index_t & aCardinality );

            // Phase management and level processing
            void reset_phase_arrays();
            void initialize_free_lists();
            void process_even_level( Cell< index_t > & aLevelVertices, index_t aSearchLevel );
            void process_odd_level( Cell< index_t > & aLevelVertices, index_t aSearchLevel );
            bool process_bridges( index_t aSearchLevel, index_t & aBlossomCounter, index_t & aCardinality );
        };

//------------------------------------------------------------------------------

        MicaliVazirani::MicaliVazirani( Graph & aGraph )
            : mGraph( aGraph ),
              mMark( aGraph.size() > 0 ? aGraph.size() : 1 ),
              mVisited( aGraph.size() > 0 ? aGraph.size() : 1 )
        {
            mNumVertices = aGraph.size();
            if ( mNumVertices == 0 ) return;

            // Set consecutive 0-based indices and cache degree in level()
            mNumEdges = 0;
            for ( index_t k = 0; k < mNumVertices; ++k )
            {
                mGraph( k )->set_index( k );
                mGraph( k )->set_level( mGraph( k )->number_of_vertices() );
                mNumEdges += mGraph( k )->number_of_vertices();
            }
            mNumEdges /= 2;

            mInfinity = std::max( mNumVertices, mNumEdges ) + 1;
            mSearchLevelLimit = mNumVertices / 2 + 1;

            // Build CSR adjacency structure
            build_csr();

            // Allocate edge status array (aligned with CSR)
            mEdgeStatus.set_size( mNumEdges * 2, 0 );

            // Allocate core arrays
            mMatch.set_size( mNumVertices, gNoIndex );
            mEvenLevel.set_size( mNumVertices, mInfinity );
            mOddLevel.set_size( mNumVertices, mInfinity );
            mBlossom.set_size( mNumVertices, gNoIndex );
            mForest.set_size( mNumVertices, gNoIndex );
            mPath.set_size( mNumVertices, gNoIndex );
            mLR.set_size( mNumVertices, 0 );
            mPrdctr.set_size( mNumVertices, 0 );

            index_t tBlossomSize = mNumVertices / 2 + 1;
            mBase.set_size( tBlossomSize, gNoIndex );
            mBstar.set_size( tBlossomSize, gNoIndex );
            mPeakLeft.set_size( tBlossomSize, gNoIndex );
            mPeakRight.set_size( tBlossomSize, gNoIndex );
            mBIndex.set_size( tBlossomSize, 0 );

            // Stack arrays
            index_t tStackSize = mInfinity;
            mPred0.set_size( tStackSize, 0 );
            mPred1.set_size( tStackSize, 0 );
            mPIndex.set_size( mNumVertices, 0 );
            mDerp0.set_size( tStackSize, 0 );
            mDerp1.set_size( tStackSize, 0 );
            mDIndex.set_size( mNumVertices, 0 );
            mAnom0.set_size( mNumVertices, 0 );
            mAnom1.set_size( mNumVertices, 0 );
            mAIndex.set_size( mNumVertices, 0 );
            mBridge0.set_size( mNumEdges + 1, 0 );
            mBridge1.set_size( mNumEdges + 1, 0 );
            mBridge2.set_size( mNumEdges + 1, 0 );
        }

//------------------------------------------------------------------------------

        /**
         * @brief Builds CSR (Compressed Sparse Row) adjacency structure
         *
         * Converts the adjacency list representation into CSR format for cache-
         * efficient neighbor access. The CSR format stores all neighbors in a
         * contiguous array with an offset array for each vertex.
         *
         * Time: O(V + E)
         * Space: O(V + E)
         */
        void
        MicaliVazirani::build_csr()
        {
            // Build offset array using prefix sum of degrees
            mAdjStart.set_size( mNumVertices + 1, 0 );

            for ( index_t k = 0; k < mNumVertices; ++k )
            {
                mAdjStart( k + 1 ) = mAdjStart( k ) + mGraph( k )->level();
            }

            // Flatten all adjacency lists into a single contiguous array
            mAdjList.set_size( mNumEdges * 2, 0 );

            for ( index_t k = 0; k < mNumVertices; ++k )
            {
                Vertex * tVertex = mGraph( k );
                index_t tStart = mAdjStart( k );
                uint tDegree = tVertex->level();

                for ( uint j = 0; j < tDegree; ++j )
                {
                    mAdjList( tStart + j ) = tVertex->vertex( j )->index();
                }
            }

            for ( index_t u = 0; u < mNumVertices; ++u )
            {

                for ( index_t k = mAdjStart( u ); k < mAdjStart( u+1 ); ++k )
                {
                    index_t v = mAdjList( k );

                    key_t tKey = u > v ? u * mNumVertices + v : v * mNumVertices + u;
                    mEdgeMap[ tKey ] = k;
                }
            }
        }

//------------------------------------------------------------------------------

        /**
         * @brief Finds the CSR index for an edge (u,v)
         *
         * Searches the adjacency list of the lower-indexed vertex to find the
         * position of the edge in the CSR structure. This index can be used to
         * access both the neighbor in mAdjList and the edge status in mEdgeStatus.
         *
         * @param aU First vertex
         * @param aV Second vertex
         * @return Index into mAdjList and mEdgeStatus arrays
         *
         * Time: O(degree(min(u,v)))
         */
        inline index_t
        MicaliVazirani::find_edge_index( index_t aU, index_t aV ) const
        {
            key_t tMin = std::min( aU, aV );
            key_t tMax = std::max( aU, aV );

            // not sure if this is faster. Maybe?
            return mEdgeMap( tMax * mNumVertices + tMin );
            // Search in adjacency list of vertex with smaller index
            /*for ( index_t k = mAdjStart( tMin ); k < mAdjStart( tMin + 1 ); ++k )
            {
                if ( mAdjList( k ) == tMax )
                {
                    return k;
                }
            }*/

            BELFEM_ERROR( false, "Edge not found: %lu - %lu",
                        (long unsigned int) aU, (long unsigned int) aV );
            return 0;
        }

//------------------------------------------------------------------------------

        /**
         * @brief Gets the status code for an edge
         *
         * Edge status codes track usage during path search:
         * - Even values: edge not yet used
         * - Odd values: edge has been used
         *
         * @return Edge status (combination of usage flags)
         */
        inline index_t
        MicaliVazirani::get_edge_status( index_t aU, index_t aV ) const
        {
            return mEdgeStatus( find_edge_index( aU, aV ) );
        }

//------------------------------------------------------------------------------

        /**
         * @brief Adds a status code to an edge
         * @param aCode Status flag to add (typically 1 or 2)
         */
        inline void
        MicaliVazirani::add_edge_status( index_t aU, index_t aV, index_t aCode )
        {
            mEdgeStatus( find_edge_index( aU, aV ) ) += aCode;
        }

//------------------------------------------------------------------------------

        /**
         * @brief Adds an entry to a free-list based linked-list stack
         *
         * This implements a linked-list stack using free-list allocation to avoid
         * dynamic memory allocation in inner loops. Each stack uses three arrays:
         * - Stack0: next pointer (linked list)
         * - Stack1: data (the value being stored)
         * - Index: head pointer for each vertex
         *
         * @param aStack0 Next pointer array
         * @param aStack1 Data array
         * @param aIndex Head pointer array (indexed by vertex)
         * @param aFree Free list head (modified by this call)
         * @param aU Vertex index (key)
         * @param aV Value to push
         *
         * Time: O(1)
         */
        void
        MicaliVazirani::add_to_stack( Vector< index_t > & aStack0,
                                   Vector< index_t > & aStack1,
                                   Vector< index_t > & aIndex,
                                   index_t & aFree,
                                   index_t aU, index_t aV )
        {
            BELFEM_ASSERT( aFree != 0, "Stack overflow" );

            // Allocate node from free list
            index_t tNext = aFree;
            aFree = aStack0( tNext - 1 );

            // Link into vertex's stack
            aStack0( tNext - 1 ) = aIndex( aU );
            aStack1( tNext - 1 ) = aV;
            aIndex( aU ) = tNext;
        }

//------------------------------------------------------------------------------

        /**
         * @brief Adds a bridge edge to the bridge stack for a specific level
         *
         * Bridges are edges connecting two even-level vertices at the same or
         * different levels. They are processed to find augmenting paths or form
         * blossoms. This function stores bridge triples (u, v, level) organized
         * by their associated level.
         *
         * @param aU First vertex of bridge
         * @param aV Second vertex of bridge
         * @param aBr Bridge level (search level where bridge was discovered)
         *
         * Time: O(1)
         */
        void
        MicaliVazirani::add_bridge( index_t aU, index_t aV, index_t aBr )
        {
            // Defensive: ignore bridges with undefined vertices
            if ( aU == gNoIndex || aV == gNoIndex || aBr == gNoIndex )
            {
                return;
            }

            BELFEM_ASSERT( mBFree != 0, "Bridge stack overflow" );

            // Allocate node from free list
            index_t tNext = mBFree;
            mBFree = mBridge0( tNext - 1 );

            // Store bridge triple and link to level's bridge list
            mBridge0( tNext - 1 ) = mBIndex( aBr );
            mBridge1( tNext - 1 ) = aU;
            mBridge2( tNext - 1 ) = aV;
            mBIndex( aBr ) = tNext;
        }

//------------------------------------------------------------------------------

        /**
         * @brief Finds an ancestor vertex in the predecessor stack
         *
         * Walks the predecessor stack for vertex aV looking for the first ancestor
         * that satisfies the edge status condition. The predecessor stack stores
         * potential path vertices discovered during the search.
         *
         * @param aV Vertex whose ancestors to search
         * @param aCheckUnused If true, return first ancestor with unused edge (even status code)
         *                     If false, return first ancestor with status < 2
         * @param aIndex Head of predecessor stack (modified during traversal)
         * @return First qualifying ancestor, or gNoIndex if none found
         *
         * Time: O(stack depth)
         */
        index_t
        MicaliVazirani::ancest( index_t aV, bool aCheckUnused, index_t & aIndex )
        {
            // Walk the predecessor linked list
            while ( aIndex != 0 )
            {
                index_t tW = mPred1( aIndex - 1 );  // Get vertex from stack
                aIndex = mPred0( aIndex - 1 );       // Move to next in list

                // Skip erased vertices
                if ( mMark.test( tW ) ) continue;

                // Check edge status for (tW, aV)
                index_t tCode = get_edge_status( tW, aV );
                if ( aCheckUnused )
                {
                    // Return if edge is unused (even status)
                    if ( tCode % 2 == 0 ) return tW;
                }
                else
                {
                    // Return if edge has low usage
                    if ( tCode < 2 ) return tW;
                }
            }
            return gNoIndex;
        }

//------------------------------------------------------------------------------

        /**
         * @brief Expands a blossom along a path (blossom-star operation)
         *
         * This function "unwinds" a blossom structure to find the base vertex.
         * When a vertex is part of a blossom, we need to traverse the blossom
         * structure to find the actual vertex to work with. This updates mForest
         * pointers and mBstar arrays as it goes.
         *
         * @param aV Input: forest parent vertex (or gNoIndex to start fresh)
         *           Output: last vertex before base
         * @param aU Input: vertex to expand from
         *           Output: base vertex (not in a blossom)
         *
         * The algorithm follows mBstar pointers through nested blossoms until
         * finding a vertex not in any blossom.
         *
         * Time: O(blossom nesting depth)
         */
        void
        MicaliVazirani::bastar( index_t & aV, index_t & aU )
        {
            index_t tVOld = aV;

            // Traverse blossom structure until reaching a vertex not in a blossom
            while ( mBlossom( aU ) != gNoIndex )
            {
                mForest( aU ) = aV;
                aV = aU;
                aU = mBstar( mBlossom( aU ) );
            }

            // Update forest pointers along the path
            index_t tW = mForest( aV );
            if ( tVOld == gNoIndex )
            {
                mForest( aU ) = aV;
                aV = gNoIndex;
            }

            // Update blossom star pointers back to base
            while ( tW != tVOld && tW != gNoIndex )
            {
                if ( mBlossom( tW ) != gNoIndex )
                {
                    mBstar( mBlossom( tW ) ) = aU;
                }
                tW = mForest( tW );
            }
        }

//------------------------------------------------------------------------------

        /**
         * @brief Constructs an alternating path through blossoms and forest
         *
         * This function builds a path from aHigh to aLow through the alternating
         * forest, handling blossom expansion as needed. The path is stored in mPath
         * and oriented according to aJob.
         *
         * Uses manual stack management to prevent stack overflow on deep graphs
         * that would occur with recursive implementation.
         *
         * @param aHigh Starting vertex (root side)
         * @param aLow  Ending vertex (destination)
         * @param aB    Blossom ID context (gNoIndex if not in a blossom)
         * @param aJob  Path orientation:
         *               1 = direct (follow mForest direction)
         *              -1 = inverted (reverse at end)
         *               2 = through blossom
         *
         * On completion, mPath contains the path: mPath[v] = next vertex on path
         *
         * Time: O(path length × blossom depth)
         */
        void
        MicaliVazirani::find_path( index_t aHigh, index_t aLow, index_t aB, int aJob )
        {
            // Stack entry for simulating recursion
            struct StackEntry
            {
                index_t high, low, b;
                int job;
                index_t entrance, bass, lastB;
                int retadd;  // Return address code
            };

            std::stack< StackEntry > tStack;
            index_t tHi = aHigh, tLo = aLow, tB = aB;
            int tJ = aJob;

            // Variables used across goto labels (must be in outer scope)
            index_t tV, tVIndex, tU;
            index_t tEntrance, tBass, tLastB;
            int tRetadd;

            // Clear visited marks from previous find_path calls
            mVisited.reset();

            index_t tCountA = 0 ;
            index_t tCountB = 0 ;
            while ( true )
            {
                // #guard
                BELFEM_ERROR( tCountA++ < mNumVertices, "could not find path" );

            start_outer:
                // Base case: already at destination
                if ( tHi == tLo )
                {
                    if ( tStack.empty() ) return;

                    // Restore from stack (simulated return from recursion)
                    StackEntry tE = tStack.top(); tStack.pop();
                    tHi = tE.high; tLo = tE.low; tB = tE.b; tJ = tE.job;
                    tEntrance = tE.entrance; tBass = tE.bass; tLastB = tE.lastB;
                    tRetadd = tE.retadd;

                    // Jump to appropriate continuation point
                    if ( tRetadd == 777 ) { tEntrance = tBass; goto step7; }
                    else if ( tRetadd == 902 )  // Left peak case
                    {
                        mPath( mPeakLeft( tLastB ) ) = mPeakRight( tLastB );
                        tStack.push( { tHi, tLo, tB, tJ, tEntrance, tBass, tLastB, 777 } );
                        tHi = mPeakRight( tLastB ); tLo = tBass; tJ = 1; tB = tLastB;
                        continue;
                    }
                    else  // Right peak case (904)
                    {
                        mPath( mPeakRight( tLastB ) ) = mPeakLeft( tLastB );
                        tStack.push( { tHi, tLo, tB, tJ, tEntrance, tBass, tLastB, 777 } );
                        tHi = mPeakLeft( tLastB ); tLo = tEntrance; tJ = 1; tB = tLastB;
                        continue;
                    }
                }

                // Initialize search from high vertex
                tV = tHi;
                if ( tV == gNoIndex )
                {
                    // Invalid start vertex - abort this search branch
                    if ( tStack.empty() ) return;
                    StackEntry tE = tStack.top(); tStack.pop();
                    tHi = tE.high; tLo = tE.low; tB = tE.b; tJ = tE.job;
                    tEntrance = tE.entrance; tBass = tE.bass;
                    tLastB = tE.lastB; tRetadd = tE.retadd;
                    continue;
                }
                tVIndex = mPIndex( tV );

                // Main search loop: walk predecessors to find path to destination
                tCountB = 0 ;
                while ( true )
                {
                    // #guard
                    BELFEM_ERROR( tCountB++ < mNumVertices, "could not find path" );

                    // Find next ancestor in predecessor chain
                    tU = ( tVIndex != 0 && tV != gNoIndex )
                         ? ancest( tV, false, tVIndex )
                         : gNoIndex;

                    // No more ancestors - need to backtrack via forest
                    if ( tU == gNoIndex )
                    {
                        if ( mForest( tV ) == gNoIndex )
                        {
                            // Fallback: if we've climbed to the intended destination (tLo)
                            // but ancest() can't yield it due to edge-status filtering,
                            // accept the forest chain tHi -> ... -> tLo as the path.
                            if ( tV == tLo )
                            {
                                // Reconstruct path by walking forward along forest from tHi to tLo
                                // Guard against accidental cycles with iteration cap
                                index_t cur = tHi;
                                int cap = 0;
                                while ( cur != gNoIndex && cur != tLo && cap < (int)mNumVertices )
                                {
                                    index_t nxt = mForest( cur );
                                    if ( nxt == gNoIndex ) break;
                                    mPath( cur ) = nxt;
                                    cur = nxt;
                                    ++cap;
                                }

                                if ( cur == tLo )
                                {
                                    // Successful reconstruction; proceed with standard post-processing
                                    tEntrance = tHi;
                                    goto step7;
                                }
                            }

                            BELFEM_ERROR( false,
                                "Cannot find path in find_path | side=%s job=%d numCalls=%lu aHigh=%lu aLow=%lu tV=%lu",
                                ( aJob == 1 ? "LEFT" : ( aJob == -1 ? "RIGHT" : "UNKNOWN" ) ),
                                (int)aJob, (unsigned long)mNumCalls,
                                (unsigned long)aHigh, (unsigned long)aLow, (unsigned long)tV );
                        }

                        // Backtrack via forest
                        index_t tNextV = mForest( tV );

                        // If the forest points to itself or directly to the destination,
                        // we have reached the end of the path.
                        if ( tNextV == tV || tNextV == tLo )
                        {
                            // We've reached the destination - reconstruct path
                            index_t cur = tHi;
                            int cap = 0;
                            while ( cur != gNoIndex && cur != tLo && cap < (int)mNumVertices )
                            {
                                index_t nxt = mForest( cur );
                                if ( nxt == gNoIndex || nxt == cur ) break;
                                mPath( cur ) = nxt;
                                cur = nxt;
                                ++cap;
                            }
                            // Set final link to destination
                            if ( cur != tLo && cur != gNoIndex )
                            {
                                mPath( cur ) = tLo;
                            }
                            tEntrance = tHi;
                            goto step7;
                        }

                        // If the next vertex is gNoIndex we cannot continue.
                        BELFEM_ERROR( tNextV != gNoIndex,
                                     "find_path: reached gNoIndex vertex while back‑tracking "
                                     "(high=%lu low=%lu)", (unsigned long)aHigh,
                                     (unsigned long)aLow );

                        tV = tNextV;

                        tVIndex = mPIndex( tV );
                        tVIndex = mPIndex( tV );
                        continue;
                    }

                    // Handle blossom membership
                    if ( mBlossom( tV ) == tB ) { add_edge_status( tU, tV, 2 ); }
                    else if ( mBlossom( tV ) != gNoIndex ) { tU = mBase( mBlossom( tV ) ); }

                    // Haven't reached destination yet - continue search
                    if ( tU != tLo )
                    {
                        if ( !mVisited.test( tU ) )
                        {
                            index_t tMinU = std::min( mEvenLevel( tU ), mOddLevel( tU ) );
                            index_t tMinLo = std::min( mEvenLevel( tLo ), mOddLevel( tLo ) );
                            if ( tMinU > tMinLo )
                            {
                                if ( tJ == 2 || !( mBlossom( tU ) == tB && mLR( tU ) == -mLR( tHi ) ) )
                                {
                                    mVisited.set( tU );
                                    mForest( tU ) = tV;
                                    tV = tU;
                                    tVIndex = mPIndex( tV );
                                }
                            }
                        }
                        continue;
                    }

                    // Destination reached! Reconstruct path from tHi to tLo
                    mPath( tV ) = tLo;
                    while ( tV != tHi )
                    {
                        index_t tTemp = tV;
                        tV = mForest( tV );
                        mPath( tV ) = tTemp;
                    }

                    tEntrance = tHi;
                step7:
                    // Post-process path to handle nested blossoms
                    while ( tEntrance != tLo )
                    {
                        tBass = mPath( tEntrance );
                        if ( mBlossom( tEntrance ) == tB )
                        {
                            tEntrance = tBass;
                            continue;
                        }

                        tLastB = mBlossom( tEntrance );
                        if ( tLastB == gNoIndex )
                        {
                            tEntrance = tBass;
                            continue;
                        }

                        if ( mEvenLevel( tEntrance ) <= mOddLevel( tEntrance ) )
                        {
                            tStack.push( { tHi, tLo, tB, tJ, tEntrance, tBass, tLastB, 777 } );
                            tHi = tEntrance; tLo = tBass; tJ = 2; tB = tLastB;
                            goto start_outer;
                        }
                        else
                        {
                            tRetadd = ( mLR( tEntrance ) > 0 ) ? 902 : 904;
                            index_t tLastHigh = ( mLR( tEntrance ) > 0 ) ?
                                mPeakLeft( tLastB ) : mPeakRight( tLastB );
                            tStack.push( { tHi, tLo, tB, tJ, tEntrance, tBass, tLastB, tRetadd } );
                            tHi = tLastHigh; tLo = tEntrance; tJ = -1; tB = tLastB;
                            goto start_outer;
                        }
                    }

                    if ( tJ == -1 )
                    {
                        // Invert path
                        index_t tPree = gNoIndex, tPntr = tHi, tSucc = mPath( tHi );
                        while ( tPntr != tLo )
                        {
                            mPath( tPntr ) = tPree;
                            tPree = tPntr;
                            tPntr = tSucc;

                            // #guard
                            if ( tPntr == gNoIndex ) break;
                            tSucc = mPath( tPntr );
                        }
                        // #guard
                        if ( tPntr == gNoIndex ) break;

                        mPath( tPntr ) = tPree;
                    }

                    if ( tStack.empty() ) return;
                    {
                        StackEntry tE = tStack.top(); tStack.pop();
                        tHi = tE.high; tLo = tE.low; tB = tE.b; tJ = tE.job;
                        tEntrance = tE.entrance;
                        tBass = tE.bass; tLastB = tE.lastB;
                        tRetadd = tE.retadd;
                    }

                    if ( tRetadd == 777 ) { tEntrance = tBass; goto step7; }
                    else if ( tRetadd == 902 )
                    {
                        mPath( mPeakLeft( tLastB ) ) = mPeakRight( tLastB );
                        tStack.push( { tHi, tLo, tB, tJ, tEntrance, tBass, tLastB, 777 } );
                        tHi = mPeakRight( tLastB ); tLo = tBass; tJ = 1; tB = tLastB;
                        goto start_outer;
                    }
                    else
                    {
                        mPath( mPeakRight( tLastB ) ) = mPeakLeft( tLastB );
                        tStack.push( { tHi, tLo, tB, tJ, tEntrance, tBass, tLastB, 777 } );
                        tHi = mPeakLeft( tLastB ); tLo = tBass; tJ = 1; tB = tLastB;
                        goto start_outer;
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        /**
         * @brief Blossom augmentation - finds augmenting paths or creates blossoms
         *
         * This is the core subroutine of the Micali-Vazirani algorithm. Given a
         * bridge edge (aW1, aW2) connecting two even-level vertices, it performs
         * a bidirectional search to either:
         * 1. Find an augmenting path (if both sides reach unmatched vertices)
         * 2. Create a new blossom (if the paths meet at a common vertex - DCV)
         *
         * The algorithm maintains two frontiers (left from aW1, right from aW2)
         * and alternates processing them based on their levels. When the frontiers
         * meet, a Deepest Common Vertex (DCV) is identified as the blossom base.
         *
         * @param aW1 Left bridge endpoint
         * @param aW2 Right bridge endpoint
         * @param aSearchLevel Current BFS level
         * @param aBlossomCounter Counter for blossom IDs (incremented if blossom formed)
         * @param aCardinality Current matching cardinality (incremented if augmented)
         * @return true if augmenting path found, false otherwise
         *
         * Time: O(V) per call
         */
        bool
        MicaliVazirani::blsaug(
            index_t aW1, index_t aW2,
            index_t aSearchLevel, index_t & aBlossomCounter, index_t & aCardinality )
        {
            ++mNumCalls;  // Unique ID for this blsaug call (for LR marking)

            Cell< index_t > tEstk;   // Erase stack for topological cleanup
            tEstk.reserve( mNumVertices );
            Cell< index_t > tMember; // Blossom member vertices
            tMember.reserve( mNumVertices );

            // Expand blossoms at bridge endpoints to get base vertices
            index_t tLeftVertex = aW1, tLeftJump = gNoIndex;
            if ( mBlossom( aW1 ) != gNoIndex )
            {
                index_t tZero = gNoIndex;
                bastar( tZero, tLeftVertex );
            }

            index_t tRightVertex = aW2, tRightJump = gNoIndex;
            if ( mBlossom( aW2 ) != gNoIndex )
            {
                index_t tZero = gNoIndex;
                bastar( tZero, tRightVertex );
            }

            // If both sides are the same vertex, no augmentation possible
            if ( tRightVertex == tLeftVertex ) return false;

            // Initialize bidirectional search
            index_t tLeftIndex = mPIndex( tLeftVertex );
            index_t tRightIndex = mPIndex( tRightVertex );

            mLR( tLeftVertex ) = mNumCalls;      // Mark left with +mNumCalls
            mLR( tRightVertex ) = -( int ) mNumCalls;  // Mark right with -mNumCalls
            tMember.push( tLeftVertex );
            tMember.push( tRightVertex );
            mForest( aW1 ) = gNoIndex;
            index_t tDcv = gNoIndex;  // Deepest Common Vertex (blossom base)
            index_t tBarier = aW2;    // Barrier for right-side search

            bool tFormBlossom = false;

            // Bidirectional search: alternate between left and right frontiers
            index_t tCount = 0 ;

            while ( true )
            {
                // #guard
                if ( tCount++ > mNumVertices ) return false ;

                // Check for augmenting path (both endpoints unmatched)
                if ( mMatch( tLeftVertex ) == gNoIndex && mMatch( tRightVertex ) == gNoIndex )
                {
                    // Anchor bridge vertices as forest roots for path reconstruction
                    if ( aW1 != gNoIndex ) mForest( aW1 ) = aW1;
                    if ( aW2 != gNoIndex ) mForest( aW2 ) = aW2;

                    // Build paths from roots to bridge endpoints
                    // Left side: direct path (Root -> ... -> Bridge)
                    find_path( tLeftVertex, aW1, gNoIndex, 1 );

                    // Right side: inverted path (Bridge -> ... -> Root)
                    find_path( tRightVertex, aW2, gNoIndex, -1 );

                    // Connect the two paths at the bridge
                    mPath( aW1 ) = aW2;

                    // Validate path integrity before applying matching
                    {
                        bool valid = true;
                        index_t cur = tLeftVertex;
                        index_t steps = 0;
                        while ( cur != tRightVertex )
                        {
                            if ( steps > mNumVertices + 1 ) { valid = false; break; }
                            index_t nxt = mPath( cur );
                            if ( nxt == gNoIndex ) { valid = false; break; }
                            cur = nxt;
                            ++steps;
                        }

                        if ( !valid )
                        {
                            // Abort this augmentation attempt cleanly
                            return false;
                        }
                    }

                    // Apply matching along augmenting path
                    index_t tP1 = tLeftVertex;
                    while ( true )
                    {
                        index_t tP2 = mPath( tP1 );
                        if ( tP2 == gNoIndex )
                        {
                            // Path broken - abort safely
                            return false;
                        }
                        mMatch( tP1 ) = tP2;
                        mMatch( tP2 ) = tP1;
                        tEstk.push( tP1 );
                        tEstk.push( tP2 );
                        tP1 = mPath( tP2 );
                        if ( tP2 == tRightVertex ) break;
                    }
                    ++aCardinality;

                    // Topological erase using bitset
                    while ( !tEstk.empty() )
                    {
                        tP1 = tEstk.pop();
                        if ( !mMark.test( tP1 ) )
                        {
                            mMark.set( tP1 );
                            index_t tNext = mDIndex( tP1 );
                            while ( tNext != 0 )
                            {
                                index_t tP2 = mDerp1( tNext - 1 );
                                tNext = mDerp0( tNext - 1 );
                                if ( --mPrdctr( tP2 ) == 0 ) tEstk.push( tP2 );
                            }
                        }
                    }

                    return true;  // Augmentation found
                }

                index_t tMinLeft = std::min( mEvenLevel( tLeftVertex ), mOddLevel( tLeftVertex ) );
                index_t tMinRight = std::min( mEvenLevel( tRightVertex ), mOddLevel( tRightVertex ) );

                if ( tMinLeft >= tMinRight )
                {
                    // Process left side
                    index_t tU = ( tLeftIndex != 0 ) ? ancest( tLeftVertex, true, tLeftIndex ) : gNoIndex;
                    if ( tU == gNoIndex )
                    {
                        if ( mForest( tLeftVertex ) == gNoIndex )
                        {
                            if ( tDcv != gNoIndex ) tFormBlossom = true;
                            break;
                        }
                        tLeftVertex = mForest( tLeftVertex );
                        tLeftIndex = mPIndex( tLeftVertex );
                    }
                    else
                    {
                        add_edge_status( tLeftVertex, tU, 1 );
                        tLeftJump = tLeftVertex;
                        if ( mBlossom( tU ) != gNoIndex )
                        {
                            bastar( tLeftVertex, tU );
                            tLeftIndex = mPIndex( tLeftVertex );
                        }

                        if ( std::abs( mLR( tU ) ) != ( int ) mNumCalls )
                        {
                            mLR( tU ) = mNumCalls;
                            tMember.push( tU );
                            mForest( tU ) = tLeftVertex;
                            tLeftVertex = tU;
                            tLeftIndex = mPIndex( tLeftVertex );
                        }
                        else if ( mLR( tU ) == -( int ) mNumCalls )
                        {
                            // Found Deepest Common Vertex (DCV) - collision between sides
                            mLR( tU ) = mNumCalls;
                            tMember.push( tU );
                            tRightVertex = mForest( tRightVertex );
                            if ( tRightJump != gNoIndex ) tRightVertex = tRightJump;
                            tRightIndex = mPIndex( tRightVertex );
                            mForest( tU ) = tLeftVertex;
                            tLeftVertex = tU;
                            tLeftIndex = mPIndex( tLeftVertex );
                            tDcv = tU;
                        }
                    }
                }
                else
                {
                    // Process right side
                    index_t tU = ( tRightIndex != 0 ) ? ancest( tRightVertex, true, tRightIndex ) : gNoIndex;
                    if ( tU == gNoIndex )
                    {
                        if ( tRightVertex == tBarier )
                        {
                            if ( tDcv == gNoIndex )
                            {
                                break;
                            }
                            // Switch barrier
                            tRightVertex = tDcv;
                            tRightIndex = mPIndex( tRightVertex );
                            tBarier = tDcv;
                            mLR( tRightVertex ) = -( int ) mNumCalls;
                            tMember.push( tRightVertex );
                            if ( mForest( tLeftVertex ) == gNoIndex )
                            {
                                tFormBlossom = true;
                                break;
                            }
                            tLeftVertex = mForest( tLeftVertex );
                            if ( tLeftJump != gNoIndex ) tLeftVertex = tLeftJump;
                            tLeftIndex = mPIndex( tLeftVertex );
                        }
                        else
                        {
                            // #guard
                            if ( mForest( tRightVertex ) == gNoIndex ) return false;
                            tRightVertex = mForest( tRightVertex );
                            tRightIndex = mPIndex( tRightVertex );
                        }
                    }
                    else
                    {
                        add_edge_status( tRightVertex, tU, 1 );
                        tRightJump = tRightVertex;
                        if ( mBlossom( tU ) != gNoIndex )
                        {
                            bastar( tRightVertex, tU );
                            tRightIndex = mPIndex( tRightVertex );
                        }

                        if ( std::abs( mLR( tU ) ) != ( int ) mNumCalls )
                        {
                            mLR( tU ) = -( int ) mNumCalls;
                            tMember.push( tU );
                            mForest( tU ) = tRightVertex;
                            tRightVertex = tU;
                            tRightIndex = mPIndex( tRightVertex );
                        }
                        else if ( mLR( tU ) == ( int ) mNumCalls )
                        {
                            // Found Deepest Common Vertex (DCV) - collision between sides
                            mLR( tU ) = -( int ) mNumCalls;
                            tMember.push( tU );
                            tLeftVertex = mForest( tLeftVertex );
                            if ( tLeftJump != gNoIndex ) tLeftVertex = tLeftJump;
                            tLeftIndex = mPIndex( tLeftVertex );
                            mForest( tU ) = tRightVertex;
                            tRightVertex = tU;
                            tRightIndex = mPIndex( tRightVertex );
                            tDcv = tU;
                        }
                    }
                }
            }

            // Form new blossom if needed
            if ( tFormBlossom && tDcv != gNoIndex )
            {
                mLR( tDcv ) = 0;
                ++aBlossomCounter;

                for ( index_t k = 0; k < tMember.size(); ++k )
                {
                    index_t tU = tMember( k );
                    if ( tU != tDcv && mBlossom( tU ) == gNoIndex )
                    {
                        mBlossom( tU ) = aBlossomCounter;
                        if ( mEvenLevel( tU ) < mOddLevel( tU ) )
                        {
                            mOddLevel( tU ) = 2 * aSearchLevel + 1 - mEvenLevel( tU );
                        }
                        else
                        {
                            mEvenLevel( tU ) = 2 * aSearchLevel + 1 - mOddLevel( tU );
                            index_t tIdx = mAIndex( tU );
                            while ( tIdx != 0 )
                            {
                                index_t tV = mAnom1( tIdx - 1 );
                                tIdx = mAnom0( tIdx - 1 );
                                add_bridge( tU, tV, ( mEvenLevel( tU ) + mEvenLevel( tV ) ) / 2 );
                                add_edge_status( tU, tV, 1 );
                            }
                        }
                    }
                }
                mPeakLeft( aBlossomCounter ) = aW1;
                mPeakRight( aBlossomCounter ) = aW2;
                mBase( aBlossomCounter ) = tDcv;
                mBstar( aBlossomCounter ) = tDcv;
            }

            return false;
        }

//------------------------------------------------------------------------------

        /**
         * @brief Computes an initial greedy matching
         *
         * Provides a good starting point for the algorithm by finding a maximal
         * (though not necessarily maximum) matching using a simple greedy approach.
         * Vertices are processed in order of increasing degree to maximize chances
         * of matching low-degree vertices.
         *
         * @param aCardinality Updated with the size of the initial matching
         *
         * Time: O(V log V + E)
         */
        void
        MicaliVazirani::compute_initial_matching( index_t & aCardinality )
        {
            // Sort vertices by degree (ascending) for better greedy matching
            Cell< std::pair< uint, index_t > > tDegreeOrder;
            tDegreeOrder.reserve( mNumVertices );

            for ( index_t k = 0; k < mNumVertices; ++k )
            {
                tDegreeOrder.push( { mGraph( k )->level(), k } );
            }

            std::sort( tDegreeOrder.begin(), tDegreeOrder.end() );

            // Track unmatched vertices with bitset
            DynamicBitset tUnmatched( mNumVertices );
            for ( index_t k = 0; k < mNumVertices; ++k )
            {
                tUnmatched.set( k );
            }

            // Greedily match vertices
            for ( auto & tPair : tDegreeOrder )
            {
                index_t tV = tPair.second;
                if ( !tUnmatched.test( tV ) ) continue;

                // Find first unmatched neighbor
                for ( index_t j = mAdjStart( tV ); j < mAdjStart( tV + 1 ); ++j )
                {
                    index_t tU = mAdjList( j );
                    if ( tUnmatched.test( tU ) )
                    {
                        mMatch( tU ) = tV;
                        mMatch( tV ) = tU;
                        tUnmatched.reset( tU );
                        tUnmatched.reset( tV );
                        ++aCardinality;
                        break;
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        /**
         * @brief Resets all phase-specific data structures
         *
         * Called at the start of each phase to clear level assignments, blossom
         * memberships, forest pointers, and other phase-specific state.
         *
         * Time: O(V + B + E) where B is the number of blossoms
         */
        void
        MicaliVazirani::reset_phase_arrays()
        {
            // Reset per-vertex arrays
            for ( index_t k = 0; k < mNumVertices; ++k )
            {
                mEvenLevel( k ) = mInfinity;
                mOddLevel( k ) = mInfinity;
                mBlossom( k ) = gNoIndex;
                mPIndex( k ) = 0;
                mDIndex( k ) = 0;
                mAIndex( k ) = 0;
                mForest( k ) = gNoIndex;
                mLR( k ) = 0;
                mPrdctr( k ) = 0;
            }

            // Reset bridge index array
            mBIndex.fill( 0.0 );

            // Reset edge status
            mEdgeStatus.fill( 0 );

            // Reset bitsets
            mMark.reset();
            mVisited.reset();
        }

//------------------------------------------------------------------------------

        /**
         * @brief Initializes free lists for stack allocation
         *
         * Sets up the free-list structure for O(1) stack push/pop operations
         * during the search. Each stack element is pre-linked to form a singly
         * linked list of available nodes.
         *
         * Time: O(E + V)
         */
        void
        MicaliVazirani::initialize_free_lists()
        {
            // Initialize free list heads
            mPFree = mDFree = mAFree = mBFree = 1;
            index_t tMaxStack = mPred0.length() < mNumEdges ? mPred0.length() : mNumEdges;

            // Link predecessor, dependency, and bridge stacks
            for ( index_t k = 0; k + 1 < tMaxStack; ++k )
            {
                mPred0( k ) = k + 2;    // Next free node
                mDerp0( k ) = k + 2;
                mBridge0( k ) = k + 2;
            }
            if ( tMaxStack > 0 )
            {
                mPred0( tMaxStack - 1 ) = 0;    // End of list
                mDerp0( tMaxStack - 1 ) = 0;
                mBridge0( tMaxStack - 1 ) = 0;
            }

            // Link anomaly stack
            for ( index_t k = 0; k + 1 < mNumVertices; ++k )
            {
                mAnom0( k ) = k + 2;
            }
            if ( mNumVertices > 0 )
            {
                mAnom0( mNumVertices - 1 ) = 0;
            }
        }

//------------------------------------------------------------------------------

        /**
         * @brief Processes vertices at an even level in the BFS
         *
         * For each even-level vertex, examines all incident edges. Depending on
         * the neighbor's level, either creates a bridge (if neighbor is also even),
         * or adds the neighbor to the odd level and updates predecessor stacks.
         *
         * @param aLevelVertices Vertices at this even level
         * @param aSearchLevel Current BFS level
         *
         * Time: O(E) total across all levels
         */
        void
        MicaliVazirani::process_even_level( Cell< index_t > & aLevelVertices, index_t aSearchLevel )
        {
            for ( index_t i = 0; i < aLevelVertices.size(); ++i )
            {
                index_t tV = aLevelVertices( i );
                // Process unmatched vertices at level 0, all vertices at higher levels
                if ( aSearchLevel != 0 || mMatch( tV ) == gNoIndex )
                {
                    // Examine all neighbors
                    for ( index_t j = mAdjStart( tV ); j < mAdjStart( tV + 1 ); ++j )
                    {
                        index_t tU = mAdjList( j );
                        index_t tCode = mEdgeStatus( j );

                        // Skip matched edges and used edges
                        if ( mMatch( tU ) != tV && tCode % 2 == 0 )
                        {
                            if ( mEvenLevel( tU ) != mInfinity )
                            {
                                // Both vertices even - create bridge
                                index_t tTemp = ( mEvenLevel( tU ) + aSearchLevel ) / 2;
                                if ( mEvenLevel( tU ) != aSearchLevel || tU < tV )
                                {
                                    add_bridge( tU, tV, tTemp );
                                }
                            }
                            else
                            {
                                // Add neighbor to odd level
                                if ( mOddLevel( tU ) == mInfinity )
                                {
                                    mOddLevel( tU ) = aSearchLevel + 1;
                                }
                                if ( mOddLevel( tU ) == aSearchLevel + 1 )
                                {
                                    add_to_stack( mPred0, mPred1, mPIndex, mPFree, tU, tV );
                                    add_to_stack( mDerp0, mDerp1, mDIndex, mDFree, tV, tU );
                                    ++mPrdctr( tU );
                                }
                                else if ( mOddLevel( tU ) < aSearchLevel )
                                {
                                    // Cross-level edge (anomaly)
                                    add_to_stack( mAnom0, mAnom1, mAIndex, mAFree, tU, tV );
                                }
                            }
                        }
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        /**
         * @brief Processes vertices at an odd level in the BFS
         *
         * For each odd-level vertex, follows the matching edge to add the matched
         * vertex to the next even level. This maintains the alternating forest
         * structure.
         *
         * @param aLevelVertices Vertices at this odd level
         * @param aSearchLevel Current BFS level
         *
         * Time: O(V) total across all levels
         */
        void
        MicaliVazirani::process_odd_level( Cell< index_t > & aLevelVertices, index_t aSearchLevel )
        {
            for ( index_t i = 0; i < aLevelVertices.size(); ++i )
            {
                index_t tV = aLevelVertices( i );
                if ( mBlossom( tV ) != gNoIndex ) continue;

                index_t tU = mMatch( tV );
                if ( tU == gNoIndex ) continue;

                if ( mOddLevel( tU ) == aSearchLevel )
                {
                    // Both endpoints odd at same level - create bridge
                    if ( tU < tV ) add_bridge( tU, tV, aSearchLevel );
                }
                else if ( mOddLevel( tU ) == mInfinity )
                {
                    // Add matched vertex to next even level
                    mEvenLevel( tU ) = aSearchLevel + 1;

                    // Free predecessor stack for tU
                    index_t tNext = mPIndex( tU );
                    while ( tNext != 0 )
                    {
                        index_t tOld = mPred0( tNext - 1 );
                        mPred0( tNext - 1 ) = mPFree;
                        mPFree = tNext;
                        tNext = tOld;
                    }
                    mPIndex( tU ) = 0;

                    // Add to stacks
                    add_to_stack( mPred0, mPred1, mPIndex, mPFree, tU, tV );
                    add_to_stack( mDerp0, mDerp1, mDIndex, mDFree, tV, tU );
                    ++mPrdctr( tU );
                }
            }
        }

//------------------------------------------------------------------------------

        /**
         * @brief Processes all bridges at a given search level
         *
         * Attempts to find augmenting paths or form blossoms for each bridge edge
         * discovered at this level. Returns true if an augmenting path is found.
         *
         * @return true if augmentation found or maximum matching reached
         */
        bool
        MicaliVazirani::process_bridges( index_t aSearchLevel, index_t & aBlossomCounter, index_t & aCardinality )
        {
            index_t tNext = ( aSearchLevel != 0 ) ? mBIndex( aSearchLevel ) : 0;

            index_t tCount = 0 ;

            while ( tNext != 0  )
            {
                // #guard
                if ( tCount++ > mNumVertices ) return false;

                index_t tU = mBridge1( tNext - 1 );
                index_t tV = mBridge2( tNext - 1 );

                tNext = mBridge0( tNext - 1 );

                // Skip erased vertices and vertices in the same blossom
                if ( mMark.test( tU ) || mMark.test( tV ) ) continue;
                if ( mBlossom( tU ) != gNoIndex && mBlossom( tV ) == mBlossom( tU ) ) continue;

                if ( blsaug( tU, tV, aSearchLevel, aBlossomCounter, aCardinality ) )
                {
                    return true;  // Found augmenting path
                }
                if ( aCardinality == mNumVertices / 2 )
                {
                    return true;  // Perfect matching found
                }
            }
            return false;
        }

//------------------------------------------------------------------------------

        /**
         * @brief Main search routine - finds maximum cardinality matching
         *
         * Implements the complete Micali-Vazirani algorithm:
         * 1. Compute initial greedy matching
         * 2. For each phase:
         *    a. Build alternating forest level-by-level
         *    b. Discover bridges between even vertices
         *    c. Process bridges to find augmenting paths or form blossoms
         * 3. Repeat until no augmenting paths found
         *
         * @param aCardinality Updated with final matching cardinality
         *
         * Time: O(√V × E) total
         */
        void
        MicaliVazirani::search( index_t & aCardinality )
        {
            // Start with empty matching
            for ( index_t k = 0; k < mNumVertices; ++k )
            {
                mMatch( k ) = gNoIndex;
            }

            // Compute initial greedy matching
            compute_initial_matching( aCardinality );

            if ( aCardinality == mNumVertices / 2 ) return;  // Perfect matching found

            // Main phase loop
            while ( true )
            {
                // Reset all phase-specific data
                reset_phase_arrays();
                initialize_free_lists();

                index_t tSearchLevel = 0;
                mNumCalls = 0;
                index_t tBlossomCounter = 0;

                // Initialize level 0 with all unmatched vertices
                for ( index_t k = 0; k < mNumVertices; ++k )
                {
                    if ( mMatch( k ) == gNoIndex )
                    {
                        mEvenLevel( k ) = 0;
                    }
                }

                bool tFoundAugmentation = false;

                // Level-by-level BFS
                index_t tCount = 0;
                while ( true )
                {
                    std::cout << "#search " << tCount << std::endl;
                    BELFEM_ERROR( tSearchLevel < mSearchLevelLimit, "Search level limit exceeded" );
                    BELFEM_ERROR( tCount++ < mNumVertices, "infinite loop" );

                    // Collect all vertices at current level
                    Cell< index_t > tLevelVertices;
                    tLevelVertices.reserve( mNumVertices );
                    for ( index_t k = 0; k < mNumVertices; ++k )
                    {
                        if ( mEvenLevel( k ) == tSearchLevel || mOddLevel( k ) == tSearchLevel )
                        {
                            tLevelVertices.push( k );
                        }
                    }

                    if ( tLevelVertices.empty() ) break;

                    // Process vertices based on parity of level
                    if ( tSearchLevel % 2 == 0 )
                    {
                        process_even_level( tLevelVertices, tSearchLevel );
                    }
                    else
                    {
                        process_odd_level( tLevelVertices, tSearchLevel );
                    }

                    // Try to find augmenting paths via bridges
                    if ( process_bridges( tSearchLevel, tBlossomCounter, aCardinality ) )
                    {
                        tFoundAugmentation = true;
                        if ( aCardinality == mNumVertices / 2 ) return;
                        break;
                    }

                    ++tSearchLevel;
                }

                // If no augmentation found, matching is maximum
                if ( !tFoundAugmentation ) break;
            }
        }

//------------------------------------------------------------------------------

        /**
         * @brief Public interface to run the matching algorithm
         *
         * @param aMatch Output: aMatch[v] = matched vertex or gNoIndex if unmatched
         * @return Cardinality of the maximum matching
         */
        index_t
        MicaliVazirani::run( Cell< index_t > & aMatch )
        {
            if ( mNumVertices == 0 )
            {
                aMatch.clear();
                return 0;
            }

            index_t tCardinality = 0;
            search( tCardinality );

            // Copy result to output array
            aMatch.set_size( mNumVertices, gNoIndex );
            for ( index_t k = 0; k < mNumVertices; ++k )
            {
                if ( mMatch( k ) != gNoIndex )
                {
                    aMatch( k ) = mMatch( k );
                }
            }
            return tCardinality;
        }

//------------------------------------------------------------------------------

        /**
         * @brief Computes a maximum cardinality matching in an undirected graph
         *
         * Public API function for computing maximum cardinality matching using the
         * Micali-Vazirani algorithm.
         *
         * @param aGraph Input graph (vertex indices will be set automatically)
         * @param aMatch Output: aMatch[v] = matched vertex or gNoIndex if unmatched
         * @return Cardinality of the maximum matching found
         *
         * Time: O(√V × E)
         * Space: O(V + E)
         */
        index_t
        max_cardinality_matching( Graph & aGraph, Cell< index_t > & aMatch )
        {
            MicaliVazirani tAlgorithm( aGraph );
            return tAlgorithm.run( aMatch );
        }

//------------------------------------------------------------------------------
    }
}
