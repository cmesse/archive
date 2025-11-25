//
// Created by Christian Messe on 08.09.20.
//

#ifndef BELFEM_EN_EN_ENUMS_HPP
#define BELFEM_EN_EN_ENUMS_HPP

namespace belfem
{
    namespace engine
    {
        enum class NozzleMode
        {
            ComputeCrossSection,
            ComputeExitPressure,
            UNDEFINED
        };

        enum class UnitSystem
        {
            Metric,
            Imperial,
            UNDEFINED
        };

    }
}
#endif //BELFEM_EN_EN_ENUMS_HPP
