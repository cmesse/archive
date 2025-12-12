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

#ifndef BELFEM_FN_GRAPH_REORDER_BY_LEVELS_HPP
#define BELFEM_FN_GRAPH_REORDER_BY_LEVELS_HPP

#include "cl_Cell.hpp"
#include "cl_Graph_Vertex.hpp"
#include "cl_Map.hpp"

namespace belfem
{
    namespace graph
    {
        /**
         * Reorders vertices based on graph distance from two boundary sets.
         *
         * Creates an ordering similar to solving ∇²T = 0 with T=0 on Γ₀, T=1 on Γ₁,
         * but using purely combinatorial graph methods instead of a PDE solve.
         *
         * After calling:
         *   - vertex->index() contains the new (permuted) index
         *   - aGraph is sorted by the new indices
         *
         * @param aGraph     All vertices in the graph (will be reordered)
         * @param aBoundary0 Vertices on Γ₀ boundary (analogous to T=0)
         * @param aBoundary1 Vertices on Γ₁ boundary (analogous to T=1)
         */
        void
        reorder_by_levels(
            Graph & aGraph,
            Graph & aSinks,
            Graph & aSources,
            Map< id_t, real > * aField = nullptr,
            const bool aSort = true );
    }
}

#endif //BELFEM_FN_GRAPH_REORDER_BY_LEVELS_HPP