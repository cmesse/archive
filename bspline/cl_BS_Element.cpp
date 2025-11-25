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

#include "cl_BS_Element.hpp"

namespace belfem
{
    namespace bspline
    {
//------------------------------------------------------------------------------

        Element::Element( mesh::Element * aElement ) :
            mElement( aElement )
        {
            mBasis.set_size( aElement->number_of_nodes(), nullptr );
        }

//------------------------------------------------------------------------------

        void
        Element::set_basis(  Basis * aBasis, const uint aIndex )
        {
            mBasis( aIndex ) = aBasis;
        }

//------------------------------------------------------------------------------

        void
        Element::unflag_basis()
        {
            for( Basis * tBasis : mBasis )
            {
                tBasis->unflag();
            }
        }

//------------------------------------------------------------------------------
    }
}
