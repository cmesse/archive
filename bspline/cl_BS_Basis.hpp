/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_CL_BS_BASIS_HPP
#define BELFEM_CL_BS_BASIS_HPP

#include "typedefs.hpp"
#include "cl_Graph_Vertex.hpp"
#include "cl_Cell.hpp"
namespace belfem
{
    namespace bspline
    {
        // forward declaration of element
        class Element ;

        class Basis : public graph::Vertex
        {
            Cell< Element * > mElements ;

            uint mElementCounter = 0;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Basis( const id_t aID ) ;

//------------------------------------------------------------------------------

            ~Basis();

//------------------------------------------------------------------------------

            inline void
            increment_element_counter()
            {
                ++mElementCounter ;
            }

//------------------------------------------------------------------------------

            void
            init_element_container();

//------------------------------------------------------------------------------

            void
            insert_element( Element * aElement );

//------------------------------------------------------------------------------

            void
            link_basis();

//------------------------------------------------------------------------------
        };
    }
}
#endif //BELFEM_CL_BS_BASIS_HPP
