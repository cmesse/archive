//
// Created by christian on 9/14/24.
//
#include "cl_MantaSurface.hpp"

namespace belfem
{
    MantaSurface::MantaSurface( const index_t aIndex, mesh::Element * aElement, real & aTemperature ) :
            mElement( aElement ),
            mTemperature( aTemperature )
    {
        this->set_index( aIndex ) ;
        if( aElement != nullptr )
        {
            this->set_id( aElement->id() );
        }
    }

    MantaSurface::~MantaSurface()
    {
        this->reset_vertex_container();
    }


}