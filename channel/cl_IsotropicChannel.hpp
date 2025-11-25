//
// Created by Christian Messe on 07.12.20.
//

#ifndef BELFEM_CL_ISOTROPICCHANNEL_HPP
#define BELFEM_CL_ISOTROPICCHANNEL_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_HDF5.hpp"
#include "cl_Mesh.hpp"
#include "cl_Gas.hpp"
#include "CH_Enums.hpp"
#include "cl_CH_Segment.hpp"
#include "cl_CH_Geometry.hpp"
#include "cl_CH_Boundarylayer.hpp"
#include "cl_Spline.hpp"
#include "cl_CH_Element.hpp"

namespace belfem
{
    enum class IsotropicChannelType
    {
        CylindricChamber,
        Nozzle,
        UNDEFINED
    };

//------------------------------------------------------------------------------

    /**
     * a class for the rocket engine combustor that runs without
     * an ODE, assumes that heat loss is small compared to flow energy
     */
    class IsotropicChannel
    {
        const IsotropicChannelType mType ;

        Gas & mGas;

        Mesh & mMesh1;

        // geometry object
        channel::Geometry * mGeometry = nullptr ;

        Cell< channel::Segment * > mSegments ;

        // Boundary layer object
        channel::Boundarylayer * mBoundaryLayer ;

        Vector< real > mInitialMolarFractions ;

        Cell< Vector< real > > mMolarFractions ;

        Cell< Matrix< real > > mHeatData ;
        Cell< Matrix< real > > mViscosityData ;
        Cell< Matrix< real > > mConductivityData ;

        uint mNumberOfChannels = 1 ;

        // flag that tells if this channel reacts
        bool mIsReacting = true ;

        // Funcition pointer for computation object
        void
        ( IsotropicChannel:: * mStatesFuncition )
                (  const real          & aTt,
                   const real          & aPt ) ;

        // Funcition pointer for heatload object
        void
        ( IsotropicChannel::  * mHeatloadsFunction )();

        // tell if the heatloads are computed in reverse order
        bool mComputeInReverseOrder = false ;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        IsotropicChannel( const IsotropicChannelType aType,
                           const channel::BoundaryLayerMethod aMethod,
                           HDF5 & aDatabase,
                           Gas  & aGas,
                           Mesh * aMesh1,
                           Mesh * aMesh2 = nullptr );

//------------------------------------------------------------------------------

        ~IsotropicChannel();

//------------------------------------------------------------------------------

        void
        print();

//------------------------------------------------------------------------------

        void
        pull_temperatures();

//------------------------------------------------------------------------------

        void
        push_heatloads();

//------------------------------------------------------------------------------

        void
        push_flowdata();

//------------------------------------------------------------------------------

        void
        compute_states(
                const real & aTt,
                const real & aPt );

//------------------------------------------------------------------------------

        void
        compute_heatloads();

//------------------------------------------------------------------------------

        void
        set_surface_roughness( const real & aRa );

//------------------------------------------------------------------------------

        void
        set_bartz_geometry_params( const real & aHydraulicDiameter,
                                   const real & aNozzleCurvature ) ;

//------------------------------------------------------------------------------

        void
        set_friction_method( const channel::BoundaryLayerMethod aMethod );

//------------------------------------------------------------------------------

        /**
         * expose the segment container
         */
         Cell< channel::Segment * > &
         segments() ;

//------------------------------------------------------------------------------

        void
        set_reverse_order_flag( const bool aFlag );

//------------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------------

        void
        compute_states_chamber(
                const real & aTt,
                const real & aPt );

//------------------------------------------------------------------------------

        void
        compute_states_nozzle(
                const real & aTt,
                const real & aPt );

//------------------------------------------------------------------------------

        void
        compute_states_nozzle_from_characteristics(
                const real & aTt,
                const real & aPt );

//------------------------------------------------------------------------------

        void
        compute_pressure_derivatives();

//------------------------------------------------------------------------------

        void
        compute_massflow(
                const real & aTt,
                const real & aPt,
                real & aT,
                real & aP,
                real & aDotm );

//------------------------------------------------------------------------------

        void
        compute_heatloads_forward();

//------------------------------------------------------------------------------

        void
        compute_heatloads_backward();

//------------------------------------------------------------------------------

        /*
         * unline the x-coordinate, the surface coordinate
         * runs along the wetted surface of the nozzle
         */
        void
        compute_surface_coordinates();

//------------------------------------------------------------------------------

        /**
         * load the moc data from a file and project them onto the segments
         */
         void
         load_moc_data( const string & aFilePath );

//------------------------------------------------------------------------------

        /* real
        compute_first_segment( const real & aXoff );

        real
        compute_segment( const index_t & aIndex ) ; */

//------------------------------------------------------------------------------
    };
//------------------------------------------------------------------------------

    inline Cell< channel::Segment * > &
    IsotropicChannel::segments()
    {
        return mSegments ;
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_CL_ISOTROPICCHANNEL_HPP
