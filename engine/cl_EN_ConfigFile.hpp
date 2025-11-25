//
// Created by Christian Messe on 15.04.21.
//

#ifndef BELFEM_CL_EN_CONFIGFILE_HPP
#define BELFEM_CL_EN_CONFIGFILE_HPP

#include "cl_EN_ConfigChamber.hpp"
namespace belfem
{
    namespace engine
    {
        class ConfigFile
        {
            ConfigChamber * mChamber = nullptr ;


        public:
            ConfigFile() = default ;

            ~ConfigFile() = default ;
        };
    }

}
#endif //BELFEM_CL_EN_CONFIG_HPP
