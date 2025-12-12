//
// Created by christian on 9/11/24.
//

#ifndef BELFEM_CL_MANTA_HPP
#define BELFEM_CL_MANTA_HPP
#include "typedefs.hpp"
#include "cl_HDF5.hpp"
#include "cl_Mesh.hpp"
#include "cl_SpMatrix.hpp"
#include "cl_Element.hpp"

namespace belfem
{


    class Manta
    {
        const proc_t  mCommRank ;
        Mesh     * mMesh = nullptr ;

        public:

            Manta();

            ~Manta() ;


        // private:

    };
}
#endif //BELFEM_CL_MANTA_HPP
