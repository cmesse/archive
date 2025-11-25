//
// Created by christian on 9/14/24.
//

#ifndef BELFEM_FN_CREATE_MANTA_MESH_HPP
#define BELFEM_FN_CREATE_MANTA_MESH_HPP

#include "typedefs.hpp"
#include "cl_Mesh.hpp"
namespace belfem
{
    Mesh *
    create_manta_mesh( const string & aPath );

}

#endif //BELFEM_FN_CREATE_MANTA_MESH_HPP
