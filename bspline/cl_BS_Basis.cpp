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

#include "cl_BS_Basis.hpp"
#include "cl_BS_Element.hpp"

namespace belfem
{
    namespace bspline
    {
//------------------------------------------------------------------------------

        Basis::Basis( const id_t aID )
        {
            this->set_id( aID );
        }

//------------------------------------------------------------------------------

        Basis::~Basis()
        {
            this->reset_vertex_container();
        }

//------------------------------------------------------------------------------

        void
        Basis::init_element_container()
        {
            // allocate cell
            mElements.set_size( mElementCounter, nullptr );

            // reset conter
            mElementCounter = 0;
        }

//------------------------------------------------------------------------------

        void
        Basis::insert_element( Element * aElement )
        {
            mElements( mElementCounter++ ) = aElement ;
        }

//------------------------------------------------------------------------------

        void
        Basis::link_basis()
        {
            // step 1: unflag all basis
            for( Element * tElement : mElements )
            {
                tElement->unflag_basis();
            }

            // step 2: count connected basis
            for( Element * tElement : mElements )
            {
                for( Basis * tBasis: tElement->basis() )
                {
                    // check if basis is flagged
                    if( ! tBasis->is_flagged() )
                    {
                        // flag basis
                        tBasis->flag();

                        // increment counter
                        this->increment_vertex_counter();
                    }
                }
            }

            // step 3: init container
            this->init_vertex_container();

            // step 4: insert basis
            for( Element * tElement : mElements )
            {
                for( Basis * tBasis: tElement->basis() )
                {
                    // check if basis is flagged
                    if( tBasis->is_flagged() )
                    {
                        // flag basis
                        tBasis->unflag();

                        this->insert_vertex( tBasis );
                    }
                }
            }
        }

//------------------------------------------------------------------------------
    }
}