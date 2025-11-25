//
// Created by Christian Messe on 23.11.20.
//
#include "assert.hpp"
#include "fn_gesv.hpp"
#include "fn_trans.hpp"
#include "cl_HDF5.hpp"

#include "cl_Channel.hpp"
#include "cl_CH_Factory.hpp"
#include "cl_CH_Boundarylayer.hpp"
#include "cl_CH_Element.hpp"
#include "fn_norm.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    Channel::Channel(
            const ChannelType aType,
            const channel::BoundaryLayerMethod aMethod,
            HDF5 & aDatabase,
            Gas & aGas,
            Mesh * aMesh1,
            Mesh * aMesh2 ) :
    //mType( aType ),
            mGas( aGas ),
            mMesh1( *aMesh1 )
    //mMesh2( aMesh2 == nullptr ? *aMesh1 : *aMesh2 )
    {
        // remember the molar fractions
        mInitialMolarFractions = mGas.molar_fractions();

        // create a factory
        channel::Factory tFactory( aDatabase );


        switch ( aType )
        {
            case ( ChannelType::CoolingChannel ) :
            {
                mOde = new channel::ChannelODE( aGas, channel::ChannelMode::Channel );
                mIntegrator = new ode::Integrator( *mOde, ode::Type::RK45 );

                tFactory.create_channels( "Liner", mMesh1, mSegments );
                aDatabase.load_data( "NumChannels", mNumberOfChannels );


                mBoundaryLayer = new channel::Boundarylayer( aGas,
                                                             aMethod,
                                                             channel::SigmaRecoveryMode::Petrukov,
                                                             50 );

                // activate autostep
                mIntegrator->set_auto_timestep( true );

                break;
            }
            case ( ChannelType::CombustionChamber ) :
            {
                mOde = new channel::ChannelODE( aGas, channel::ChannelMode::Channel );
                mIntegrator = new ode::Integrator( *mOde, ode::Type::RK45 );

                mBoundaryLayer = new channel::Boundarylayer( aGas,
                                                             aMethod,
                                                             channel::SigmaRecoveryMode::vanDriest,
                                                             60 );

                // create the geometry object
                mGeometry = tFactory.create_cylinder_geometry();

                // link the ODE with this geometry
                mOde->link_geometry( mGeometry, true );

                // deactivate autostep
                mIntegrator->set_auto_timestep( true );

                // create the cylinder object
                tFactory.create_cylinder_segments( mGeometry, mMesh1, mSegments, true );

                mNumberOfChannels = 1;
                break;
            }
            default:
            {
                BELFEM_ERROR( false, "unknown channel type" );
            }
        }



        // create the elements
        uint tNumChannels = ( mSegments.size() - 1 ) / 2;

        mElements.set_size( tNumChannels, nullptr );

        uint tCount = 0;
        for ( uint e = 0; e < tNumChannels; ++e )
        {
            mElements( e ) = new channel::Element(
                    *mSegments( tCount ),
                    *mSegments( tCount + 2 ),
                    *mSegments( tCount + 1 ));
            tCount += 2;
        }

        // link last segment
        mLastSegment = mSegments( mSegments.size() - 1 );
    }

//------------------------------------------------------------------------------

    Channel::~Channel()
    {
        delete mIntegrator;
        delete mOde;

        // delete the elements
        for ( channel::Element * tElement : mElements )
        {
            delete tElement;
        }

        // delete the segments
        for ( channel::Segment * tSegment : mSegments )
        {
            delete tSegment;
        }

        // check if geometry object exists and delete if so
        if ( mGeometry != nullptr )
        {
            delete mGeometry;
        }

        delete mBoundaryLayer;
    }

//------------------------------------------------------------------------------

    void
    Channel::set_inflow_conditions_total(
            const real & aTt,
            const real & aPt,
            const real & aDotm )
    {
        // relaxation factor
        real tOmega = 0.3;

        // energy state
        real tHt = mGas.h( aTt, aPt );

        // cross section
        real tA = mSegments( 0 )->cross_section() * mNumberOfChannels;

        // entropy state
        real tS = mGas.s( aTt, aPt );

        Vector< real > tRHS( 2 );

        real tT;
        real tP;
        real tU;
        real tV;

        // guess initial values
        tU = aDotm / ( mGas.rho( aTt, aPt ) * tA );
        real tMa = tU / mGas.c( aTt, aPt );

        real tGamma = mGas.gamma( aTt, aPt );

        tT = aTt / ( 1.0 + 0.5 * ( tGamma - 1 ) * tMa * tMa );
        tP = aPt * std::pow( tT / aTt, tGamma / ( tGamma - 1 ));


        real tError = 1e12;
        real tH;

        uint tCount = 0;

        while ( tError > 1e-4 )
        {
            // update inverse density
            tV = mGas.v( tT, tP );

            // compute right hand side
            tU = aDotm * tV / tA;

            tRHS( 0 ) = mGas.h( tT, tP ) + 0.5 * tU * tU - tHt;
            tRHS( 1 ) = mGas.s( tT, tP ) - tS;

            // compute error
            tError = std::sqrt( std::pow( tRHS( 0 ) / tHt, 2 )
                                + std::pow( tRHS( 1 ) / tS, 2 ));

            tH = tHt - 0.5 * tU * tU;

            tP = mGas.isen_p( aTt, aPt, tT );

            tT *= ( 1 - tOmega );
            tT += tOmega * mGas.T_from_h( tH, tP );


            if ( mIsReacting )
            {
                mGas.remix( mInitialMolarFractions, false, false );
                mGas.remix_to_equilibrium( tT, tP, true, false );
            }

            // Check for infintie loop
            BELFEM_ERROR( tCount++ < 100, "Too many iterations" );
        }


        if ( mIsReacting )
        {
            mGas.remix_to_equilibrium( tT, tP, true, true );
        }

        // compute the first segment
        mBoundaryLayer->use_input_from_parameters( false );
        mBoundaryLayer->set_hydraulic_diameter( mSegments( 0 )->value( BELFEM_CHANNEL_DH ));
        mBoundaryLayer->set_flow_conditions( tT, tP, tU );
        mBoundaryLayer->set_center_conditions( tT, tU );
        mBoundaryLayer->compute_initial_guesses();
        mBoundaryLayer->compute( mSegments( 0 )->data());

        // copy results into other segments for first run
        real tTauW = mSegments( 0 )->value( BELFEM_CHANNEL_TAUW );
        real tDotQ = mSegments( 0 )->value( BELFEM_CHANNEL_DOTQ );
        for ( channel::Segment * tSegment : mSegments )
        {
            tSegment->set_value( BELFEM_CHANNEL_TAUW, tTauW );
            tSegment->set_value( BELFEM_CHANNEL_DOTQ, tDotQ );
        }

    }

//------------------------------------------------------------------------------

    void
    Channel::set_inflow_conditions(
            const real & aT,
            const real & aP,
            const real & aMa )
    {
        if ( mIsReacting )
        {
            mGas.remix( mInitialMolarFractions, false, false );
            mGas.remix_to_equilibrium( aT, aP, true, true );
        }

        real tW = mGas.c( aT, aP );
        real tU = tW * aMa;

        // compute the first segment
        mBoundaryLayer->use_input_from_parameters( false );
        mBoundaryLayer->set_hydraulic_diameter( mSegments( 0 )->value( BELFEM_CHANNEL_DH ));
        mBoundaryLayer->set_flow_conditions( aT, aP, tU );
        mBoundaryLayer->set_center_conditions( aT, tU );
        mBoundaryLayer->compute_initial_guesses();
        mBoundaryLayer->compute( mSegments( 0 )->data());

        // copy results into other segments for first run
        real tTauW = mSegments( 0 )->value( BELFEM_CHANNEL_TAUW );
        real tDotQ = mSegments( 0 )->value( BELFEM_CHANNEL_DOTQ );
        for ( channel::Segment * tSegment : mSegments )
        {
            tSegment->set_value( BELFEM_CHANNEL_TAUW, tTauW );
            tSegment->set_value( BELFEM_CHANNEL_DOTQ, tDotQ );
        }

    }

//------------------------------------------------------------------------------

    void
    Channel::set_surface_roughness( const real & aRa )
    {
        mBoundaryLayer->set_surface_roughness( aRa );
    }

//------------------------------------------------------------------------------

    void
    Channel::run()
    {
        Vector< real > tY( 3 );

        real & tV = tY( 0 );
        real & tU = tY( 1 );
        real & tT = tY( 2 );

        tT = mSegments( 0 )->value( BELFEM_CHANNEL_TM );
        real tP = mSegments( 0 )->value( BELFEM_CHANNEL_PM );
        tU = mSegments( 0 )->value( BELFEM_CHANNEL_UM );
        tV = mGas.v( tT, tP );

        // initial state
        Vector< real > tY0( 3 );
        Vector< real > tY1( 3, 0.0 );

        mIntegrator->timestep() = 0.0001;

        real tTm;
        real tPm;
        real tUm;

        real tX;
        real tError;
        uint tCount;
        real tR0;
        real tR = mGas.R( tT, tP );
        mElements( 0 )->segment0()->set_value( BELFEM_CHANNEL_RM, tR );

        Vector< real > tdYdx( mGas.number_of_components());
        mBoundaryLayer->use_input_from_parameters( false );
        mBoundaryLayer->set_flow_conditions( tT, tP, tU );
        mBoundaryLayer->set_wall_temperature( mSegments( 0 )->value( BELFEM_CHANNEL_TW1 ));

        mBoundaryLayer->compute_initial_guesses();

        // compute first cell
        mBoundaryLayer->compute( mElements( 0 )->segment0()->data());

        // push data
        mElements( 0 )->segment0()->push_heatloads();

        // change read mode in boundary layer
        mBoundaryLayer->use_input_from_parameters( true );


        for ( channel::Element * tElement : mElements )
        {
            tY0 = tY;
            mOde->link_element( tElement );

            // reset counter
            tCount = 0;

            // reset error
            tError = BELFEM_REAL_MAX;

            tR0 = tR;
            tR = mGas.R( tT, tP );

            // forward copy
            tElement->segment1()->set_value( BELFEM_CHANNEL_TAUW,
                                             tElement->segment0()->value( BELFEM_CHANNEL_TAUW ));
            tElement->segment1()->set_value( BELFEM_CHANNEL_DOTQ,
                                             tElement->segment0()->value( BELFEM_CHANNEL_DOTQ ));
            tElement->segment1()->set_value( BELFEM_CHANNEL_RM,
                                             tR );

            tElement->segment2()->set_value( BELFEM_CHANNEL_TAUW,
                                             tElement->segment0()->value( BELFEM_CHANNEL_TAUW ));
            tElement->segment2()->set_value( BELFEM_CHANNEL_DOTQ,
                                             tElement->segment0()->value( BELFEM_CHANNEL_DOTQ ));
            tElement->segment2()->set_value( BELFEM_CHANNEL_RM,
                                             0.5 * ( tR0 + tR ));

            // iterate this element
            while ( tError > 1e-6 )
            {
                tX = tElement->x0();

                // shift state
                tY1 = tY;
                tY = tY0;

                mIntegrator->maxtime() = tElement->x2();
                //mIntegrator->timestep() = 0.00001;
                uint tCount2 = 0;

                // first half of channel
                while ( tX < tElement->x2() )
                {
                    mIntegrator->step( tX, tY );
                    tP = mGas.p( tT, tV );

                    BELFEM_ERROR( tCount2++ < 10000, "Too many iterations" );
                }

                // remember values
                tTm = tT;
                tPm = tP;
                tUm = tU;

                tCount2 = 0;

                // second half of channel
                mIntegrator->maxtime() = tElement->x1();
                //mIntegrator->timestep() = 0.00001;

                while ( tX < tElement->x1() )
                {
                    mIntegrator->step( tX, tY );

                    tP = mGas.p( tT, tV );
                    BELFEM_ERROR( tCount2++ < 10000, "Too many iterations" );
                }

                // compute segment 2
                tElement->segment2()->set_value( BELFEM_CHANNEL_TM, tTm );
                tElement->segment2()->set_value( BELFEM_CHANNEL_PM, tPm );
                tElement->segment2()->set_value( BELFEM_CHANNEL_UM, tUm );
                mBoundaryLayer->compute( tElement->segment2()->data() );



                // compute segment 1
                tElement->segment1()->set_value( BELFEM_CHANNEL_TM, tT );
                tElement->segment1()->set_value( BELFEM_CHANNEL_PM, tP );
                tElement->segment1()->set_value( BELFEM_CHANNEL_UM, tU );
                mBoundaryLayer->compute( tElement->segment1()->data() );

                // write data onto mesh ( for heatloads )
                tElement->segment2()->push_heatloads();

                // write data onto mesh ( for heatloads )
                tElement->segment1()->push_heatloads();

                // compute error
                tError = 0.0;
                for ( uint i = 0; i < 3; ++i )
                {
                    tError += std::pow(( tY( i ) - tY1( i )) / tY1( i ), 2 );
                }
                tError = std::sqrt( tError );

                // save boundary layer info
                BELFEM_ERROR( tCount++ < 500,
                             "Too many iterations. \n    Tm=%12.3f K, pm=%12.3f bar, um=%12.3f m/s, Tw=%12.3f K, Error=%8.3g",
                             ( double ) tTm, ( double ) tPm * 1e-5, ( double ) tUm, ( double ) mBoundaryLayer->Tw(), ( double ) tError );
            }

            // make alpha linear in middle segment
            // to avoid checkerboarding
            tElement->segment2()->set_value( BELFEM_CHANNEL_ALPHA1,
                                             0.5 * ( tElement->segment0()->value( BELFEM_CHANNEL_ALPHA1 )
                                                     + tElement->segment1()->value( BELFEM_CHANNEL_ALPHA1 ) ) );

            if( tElement->segment0()->num_walls() > 1 )
            {
                tElement->segment2()->set_value( BELFEM_CHANNEL_ALPHA2,
                                                 0.5 * ( tElement->segment0()->value( BELFEM_CHANNEL_ALPHA2 )
                                                         + tElement->segment1()->value( BELFEM_CHANNEL_ALPHA2 ) ) );
            }

            // write data onto mesh ( for alpha )
            tElement->segment2()->push_heatloads();

            // write data onto mesh ( for alpha )
            tElement->segment1()->push_heatloads();

            if ( mIsReacting )
            {
                // copy mass_fractions
                tdYdx = mGas.mass_fractions();

                mGas.remix( mInitialMolarFractions, false, false );
                mGas.remix_to_equilibrium( tT, tP, true, true );

                // mBoundaryLayer->compute( tElement->segment1()->data() );
                tdYdx -= mGas.mass_fractions();

                tdYdx /= -tElement->length();

                mOde->set_composition_change(
                        ( tR - tR0 ) / ( tR * tElement->length()), tdYdx );

                //mOde->set_composition_change( 0.0, tdYdx );
            }

            // check for sign switch and suppress heatload at this point
            if( tElement->segment0()->value( BELFEM_CHANNEL_DOTQ ) *
                tElement->segment1()->value( BELFEM_CHANNEL_DOTQ ) < 0 )
            {
                tElement->segment1()->set_value( BELFEM_CHANNEL_DOTQ, 0.0 );
                tElement->segment1()->set_value( BELFEM_CHANNEL_ALPHA1, 0.0 );
                if ( tElement->segment1()->num_walls() > 1 )
                {
                    tElement->segment1()->set_value( BELFEM_CHANNEL_ALPHA2, 0.0 );
                }
            }

        }

        //this->print();
    }

//------------------------------------------------------------------------------

    void
    Channel::run_simple()
    {
        this->push_heatloads();
    }

//------------------------------------------------------------------------------

    void
    Channel::run_inverse( const real & aTthroat, const real & aPthroat,
                          const real & aTt, const real & aPt )
    {
        // remember total conditions
        mTt = aTt;
        mPt = aPt;

        Vector< real > tA( 2 );
        Vector< real > tB( 2 );
        Vector< real > tC( 2 );
        Vector< real > tD( 2 );
        Vector< real > tF( 2 );
        Matrix< real > tJ( 2, 2 );
        Vector< int > tPivot( 2 );

        real tT = aTthroat;
        real tP = aPthroat;
        real tOmega = 0.9;
        real tErr = BELFEM_REAL_MAX;

        for ( uint k = 0; k < 10; ++k )
        {
            real tDeltaT = 0.01 * tT;
            real tDeltaP = 0.01 * tP;

            this->compute_inverse_step( tT - tDeltaT, tP, tA );
            this->compute_inverse_step( tT + tDeltaT, tP, tB );
            this->compute_inverse_step( tT, tP - tDeltaP, tC );
            this->compute_inverse_step( tT, tP + tDeltaP, tD );
            this->compute_inverse_step( tT, tP, tF );

            tErr = norm( tF );

            tJ( 0, 0 ) = ( tB( 0 ) - tA( 0 )) / ( 2.0 * tDeltaT );
            tJ( 1, 0 ) = ( tD( 0 ) - tC( 0 )) / ( 2.0 * tDeltaT );
            tJ( 0, 1 ) = ( tB( 1 ) - tA( 1 )) / ( 2.0 * tDeltaP );
            tJ( 1, 1 ) = ( tD( 1 ) - tC( 1 )) / ( 2.0 * tDeltaP );

            gesv( tJ, tF, tPivot );

            tT -= tOmega * tF( 0 );
            tP -= tOmega * tF( 1 );
            std::cout << k << " " << tT << " " << tP * 1e-5 << " " << tErr << std::endl;

        }
    }

//------------------------------------------------------------------------------

    void
    Channel::compute_inverse_step( const real & aTthroat, const real & aPthroat, Vector< real > & aF )
    {
        this->set_inflow_conditions( aTthroat, aPthroat, 0.999 );
        this->run();

        // get last segment
        channel::Segment * tSegment = mSegments( mSegments.size() - 1 );

        // get static flow properties
        real tT = tSegment->value( BELFEM_CHANNEL_TM );
        real tP = tSegment->value( BELFEM_CHANNEL_PM );
        real tU = tSegment->value( BELFEM_CHANNEL_UM );

        real tTt;
        real tPt;
        mGas.total( tT, tP, tU, tTt, tPt );

        // compute error
        aF( 0 ) = ( tTt - mTt ) / mTt;
        aF( 1 ) = ( tPt - mPt ) / mPt;
    }

//------------------------------------------------------------------------------

    void
    Channel::save_data( const string & aPath )
    {
        // collect the data
        Matrix< real > tData( mSegments( 0 )->data().length(), mSegments.size());

        index_t tCount = 0;

        for ( channel::Segment * tSegment : mSegments )
        {
            tData.set_col( tCount++, tSegment->data());
        }

        tData = trans( tData );

        HDF5 tFile( aPath, FileMode::NEW );

        tFile.save_data( "Data", tData );

        tFile.close();
    }

//------------------------------------------------------------------------------

    void
    Channel::print()
    {
        for ( channel::Segment * tSegment : mSegments )
        {
            tSegment->print();
        }
    }

//------------------------------------------------------------------------------

    void
    Channel::pull_temperatures()
    {
        for ( channel::Segment * tSegment : mSegments )
        {
            tSegment->pull_surface_temperatures();
        }
    }

//------------------------------------------------------------------------------

    void
    Channel::push_heatloads()
    {
        // ask boundary layer which temperature is used
        // for the reference of alpha

        for ( channel::Segment * tSegment : mSegments )
        {
            tSegment->push_heatloads();
        }
    }

//------------------------------------------------------------------------------

    void
    Channel::push_flowdata()
    {
        for ( channel::Segment * tSegment : mSegments )
        {
            tSegment->push_flowdata();
        }
    }

//------------------------------------------------------------------------------

    void
    Channel::set_reacting_flag()
    {
        mIsReacting = true;
    }

//------------------------------------------------------------------------------

    void
    Channel::unset_reacting_flag()
    {
        mIsReacting = false;
    }


//------------------------------------------------------------------------------

    real
    Channel::compute_total_enthalpy_change()
    {

        real aHeatflux = BELFEM_QUIET_NAN ;

        if( comm_rank() == 0 )
        {
            if( mIsReacting )
            {
                mGas.remix( mInitialMolarFractions , true, false );
            }

            // inflow state
            const real & tH0 = mSegments( 0 )->value( BELFEM_CHANNEL_HM );
            const real & tU0 = mSegments( 0 )->value( BELFEM_CHANNEL_UM );

            // compute massflow
            real tDotM = tU0 * mSegments( 0 )->value( BELFEM_CHANNEL_A ) * mGas.rho(
                    mSegments( 0 )->value( BELFEM_CHANNEL_TM ),
                    mSegments( 0 )->value( BELFEM_CHANNEL_PM ) );

            // outflow state
            const real & tH1 = mSegments( mSegments.size()-1 )->value( BELFEM_CHANNEL_HM );
            const real & tU1 = mSegments( mSegments.size()-1 )->value( BELFEM_CHANNEL_UM ) ;

            aHeatflux =  tDotM * ( tH1 + tU1 * tU1 - tH0 - tU0 * tU0 );
        }

        broadcast( 0, aHeatflux );

        return aHeatflux ;

    }

//------------------------------------------------------------------------------
}