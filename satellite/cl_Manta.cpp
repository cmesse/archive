//
// Created by christian on 9/11/24.
//
#include "commtools.hpp"
#include "cl_Manta.hpp"
#include "HDF5_Tools.hpp"
#include "cl_Element_Factory.hpp"
#include "fn_sum.hpp"
#include "fn_unique.hpp"
#include "cl_Block.hpp"

namespace belfem
{


    Manta::Manta() :
        mCommRank( comm_rank() )
    {

    }

    Manta::~Manta()
    {
        if( mMesh != nullptr )
        {
            delete mMesh ;
        }

    }


}