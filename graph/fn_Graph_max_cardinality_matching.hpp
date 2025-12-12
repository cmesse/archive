/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */


#ifndef BELFEM_FN_GRAPH_MAX_CARDINALITY_MATCHING_HPP
#define BELFEM_FN_GRAPH_MAX_CARDINALITY_MATCHING_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_DynamicBitset.hpp"
#include "cl_Graph_Vertex.hpp"

namespace belfem
{
    namespace graph
    {
//------------------------------------------------------------------------------

        /**
         * Computes a maximum cardinality matching for a general graph using the
         * Micali-Vazirani algorithm. This is an O(sqrt(V) * E) algorithm.
         *
         * The algorithm finds the largest possible set of edges such that no two
         * edges share a common vertex (maximum matching).
         *
         * @param aGraph       Cell containing all vertices in the graph.
         *                     The vertex indices must be consecutive starting from 0.
         *                     Edges are defined by the vertex neighbor relationships.
         *
         * @param aMatch       Output: For each vertex index i, aMatch(i) contains the
         *                     index of its matched partner, or gNoIndex if unmatched.
         *
         * @return             The cardinality (size) of the maximum matching,
         *                     i.e., the number of matched pairs (edges in matching).
         *
         * @note               The input graph must be undirected, meaning if vertex A
         *                     is a neighbor of vertex B, then B must also be a neighbor
         *                     of A.
         *
         * @note               The vertex level() field is used internally to store
         *                     degree information during computation.
         *
         * Example usage:
         * @code
         *     Graph tGraph;
         *     // ... populate graph with vertices and edges ...
         *
         *     Cell< index_t > tMatch;
         *     index_t tCardinality = max_cardinality_matching( tGraph, tMatch );
         *
         *     // Now tMatch(i) contains the matched partner of vertex i
         *     // tCardinality is the number of matched pairs
         * @endcode
         */
        index_t
        max_cardinality_matching(
                Graph & aGraph,
                Cell< index_t >  & aMatch );

//------------------------------------------------------------------------------
    }
}

#endif // BELFEM_FN_GRAPH_MAX_CARDINALITY_MATCHING_HPP
