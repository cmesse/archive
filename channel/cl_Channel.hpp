//
// Created by Christian Messe on 23.11.20.
//

#ifndef BELFEM_CL_CHANNEL_HPP
#define BELFEM_CL_CHANNEL_HPP
#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_HDF5.hpp"
#include "cl_Mesh.hpp"
#include "cl_Gas.hpp"
#include "CH_Enums.hpp"
#include "cl_ODE_Integrator.hpp"
#include "fn_ODE_RK45.hpp"
#include "cl_CH_Segment.hpp"
#include "cl_CH_ChannelODE.hpp"
#include "cl_CH_Geometry.hpp"
#include "cl_CH_Boundarylayer.hpp"
namespace belfem
{
    namespace channel
    {
        class Element ;
    }

    /**
     * the main Channel object
     */
    class Channel
    {
        // const ChannelType mType ;

        Gas & mGas ;

        Mesh & mMesh1 ;
        // Mesh & mMesh2 ;

        // ODE Object
        channel::ChannelODE * mOde ;

        // Integrator
        ode::Integrator * mIntegrator ;

        // geometry object
        channel::Geometry * mGeometry = nullptr ;

        Cell< channel::Segment * > mSegments ;
        Cell< channel::Element * > mElements ;

        // Boundary layer object
        channel::Boundarylayer * mBoundaryLayer ;

        bool mIsReacting = false ;

        uint mNumberOfChannels = BELFEM_UINT_MAX ;

        // initial mixure
        Vector< real > mInitialMolarFractions ;

        // only needed for inverse step
        real mTt ;
        real mPt ;

        // pointer to last segment
        channel::Segment * mLastSegment = nullptr ;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        Channel( const ChannelType aType,
                 const channel::BoundaryLayerMethod aMethod,
                 HDF5 & aDatabase,
                 Gas  & aGas,
                 Mesh * aMesh1,
                 Mesh * aMesh2 = nullptr );

//------------------------------------------------------------------------------

        ~Channel();

//------------------------------------------------------------------------------

        void
        set_inflow_conditions_total(
                const real & aTt,
                const real & aPt,
                const real & aDotm );

//------------------------------------------------------------------------------

        void
        set_inflow_conditions(
                const real & aT,
                const real & aP,
                const real & aMa );

//------------------------------------------------------------------------------

        void
        set_surface_roughness( const real & aRa );

//------------------------------------------------------------------------------

        void
        save_data( const string & aPath );

//------------------------------------------------------------------------------

        void
        run();

//------------------------------------------------------------------------------

        void
        run_simple();

//------------------------------------------------------------------------------

        void
        run_inverse( const real & aTthroat, const real & aPthroat,
                     const real & aTt, const real & aPt );

//------------------------------------------------------------------------------

        void
        compute_inverse_step( const real & aTthroat, const real & aPthroat, Vector< real > & aF );

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
        set_reacting_flag();

//------------------------------------------------------------------------------

        void
        unset_reacting_flag();

//------------------------------------------------------------------------------

        /**
         * expose the gas
         */
         Gas *
         gas();

//------------------------------------------------------------------------------

        /**
         * return the exit temperature
         */
         const real &
         exit_temperature() const ;

//------------------------------------------------------------------------------

        /**
         * return the exit pressure
         */
        const real &
        exit_pressure() const ;

//------------------------------------------------------------------------------

        /**
         * return the exit velocity
         */
        const real &
        exit_velocity() const ;

//------------------------------------------------------------------------------

        /**
         * return the exit Mach number
         */
        const real &
        exit_mach() const ;

//------------------------------------------------------------------------------


        /**
         * compute the difference in total enthalpy
         */
         real
         compute_total_enthalpy_change();

//------------------------------------------------------------------------------

    };

//------------------------------------------------------------------------------

    inline Gas *
    Channel::gas()
    {
        return & mGas ;
    }

//------------------------------------------------------------------------------

    inline const real &
    Channel::exit_temperature() const
    {
        return mLastSegment->value( BELFEM_CHANNEL_TM );
    }

//------------------------------------------------------------------------------

    inline const real &
    Channel::exit_pressure() const
    {
        return mLastSegment->value( BELFEM_CHANNEL_PM );
    }

//------------------------------------------------------------------------------

    inline const real &
    Channel:: exit_velocity() const
    {
        return mLastSegment->value( BELFEM_CHANNEL_UM );
    }

//------------------------------------------------------------------------------

    inline const real &
    Channel:: exit_mach() const
    {
        return mLastSegment->value( BELFEM_CHANNEL_MAM );
    }

}
#endif //BELFEM_CL_CHANNEL_HPP
