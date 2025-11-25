//
// Created by Christian Messe on 25.09.19.
//

#ifndef BELFEM_FN_CN_MOLARS_FROM_STRING_HPP
#define BELFEM_FN_CN_MOLARS_FROM_STRING_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"


namespace belfem
{
    namespace combustion
    {
//------------------------------------------------------------------------------

        void
        molars_from_string(
                const   string & aString,
                Cell< string > & aEductLabels,
                Vector< real > & aEductMolars,
                Cell< string > & aProductLabels,
                Vector< real > & aProductMolars,
                          bool & aHasThirdBody );

//------------------------------------------------------------------------------

        void
        split_to_molar_and_label(
                const string & aString,
                      string & aLabel,
                      real   & aMolar );

//------------------------------------------------------------------------------

    }
}
#endif //BELFEM_FN_CN_MOLARS_FROM_STRING_HPP
