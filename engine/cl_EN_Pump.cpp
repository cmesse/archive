//
// Created by Christian Messe on 19.01.21.
//

#include "cl_EN_Pump.hpp"
#include "cl_Matrix.hpp"
#include "fn_gesv.hpp"
#include "cl_GT_RefGas.hpp"
#include "constants.hpp"
#include "fn_polyval.hpp"
#include "fn_dot.hpp"

namespace belfem
{
    namespace engine
    {
//------------------------------------------------------------------------------

        Pump::Pump( Gas & aGas ) :
            mGas( aGas ),
            mIsHydrogen( aGas.number_of_components() == 1 &&
                         aGas.component( 0 )->label() == "H2" ),
            mEntry( aGas, "Entry" ),
            mExitIsotropic( aGas, "ExitIsentropic" ),
            mExit( aGas, "Exit" )
        {

        }

//------------------------------------------------------------------------------

        void
        Pump::set_mass_flux( const real & aMassFlux )
        {
            mDotM = aMassFlux ;

            mDotV = aMassFlux / mEntry.rho() ;

        }

//------------------------------------------------------------------------------

        void
        Pump::set_entry( const real & aTt, const real & aPt )
        {
            mEntry.value( BELFEM_ENGINE_STATE_TT ) = aTt ;
            mEntry.value( BELFEM_ENGINE_STATE_PT ) = aPt ;
            mEntry.value( BELFEM_ENGINE_STATE_HT ) = mGas.h( aTt, aPt );
            mEntry.value( BELFEM_ENGINE_STATE_S )  = mGas.s( aTt, aPt );

            // initial guesses for the static case
            mEntry.value( BELFEM_ENGINE_STATE_T )   = aTt ;
            mEntry.value( BELFEM_ENGINE_STATE_P )   = aPt ;
            mEntry.value( BELFEM_ENGINE_STATE_H )   = mEntry.value( BELFEM_ENGINE_STATE_HT ) ;
            mEntry.value( BELFEM_ENGINE_STATE_RHO ) = mGas.rho( aTt, aPt ) ;
            mEntry.value( BELFEM_ENGINE_STATE_U )   = 0.0 ;
            mEntry.value( BELFEM_ENGINE_STATE_MA )  = 0.0 ;

            mDotV = mDotM / mEntry.rho() ;
        }

//------------------------------------------------------------------------------

        void
        Pump::set_deltap( const real & aDeltaP )
        {

            BELFEM_ERROR( ! mPt2Flag, "You must set deltap or pt2, not both!" );

            mDeltaPFlag = true ;

            mDeltaP = aDeltaP ;

            this->compute_isotropic_exit() ;
        }

//------------------------------------------------------------------------------

        void
        Pump::set_pt2( const real & aPt2 )
        {
            BELFEM_ERROR( ! mDeltaPFlag, "You must set deltap or pt2, not both!" );

            mPt2Flag = true ;

            mPt2 = aPt2 ;

            this->compute_isotropic_exit() ;
        }

//------------------------------------------------------------------------------

        void
        Pump::set_nrpm( const real & aNRPM )
        {
            mNRPM = aNRPM ;
            mN = aNRPM / 60.0 ;
        }

//------------------------------------------------------------------------------

        void
        Pump::set_D2a( const real & aD2a, const bool aFixD2  )
        {
            mD2a = aD2a ;
            mD2aFlag = true ;
            mD2Flag = aFixD2 ;
        }

//------------------------------------------------------------------------------

        void
        Pump::set_psi( const real & aPsi  )
        {
            mPsi     = aPsi ;
            mPsiFlag = true ;
            mD2Flag  = true ;
        }

//------------------------------------------------------------------------------

        void
        Pump::set_DN( const real & aDN )
        {
            mDN = aDN ;
            mDNFlag = true ;
        }

//------------------------------------------------------------------------------

        void
        Pump::set_D1a( const real & aD1a )
        {
            mD1a = aD1a ;
            mD1aFlag = true ;
        }

//------------------------------------------------------------------------------

        /**
         * find a value for beta2 to fullfill the haller criterion
         */
        void
        Pump::set_haller( const real aHaller )
        {
            mHaller = aHaller ;
            mHallerFlag = true ;
        }

//------------------------------------------------------------------------------

        void
        Pump::set_z2( const uint & aZ2 )
        {
            mZ2Flag = true ;
            mZ2 = ( real ) aZ2 ;
        }

//------------------------------------------------------------------------------

        void
        Pump::set_z1( const uint & aZ1 )
        {
            mZ1Flag = true ;
            mZ1 = ( real ) aZ1 ;
        }

//------------------------------------------------------------------------------

        void
        Pump::set_s2( const real & aS2 )
        {
            mS2 = aS2 ;
        }

//------------------------------------------------------------------------------

        void
        Pump::set_beta2( const uint & aBeta2 )
        {
            mBeta2 = aBeta2 ;
            mBeta2Flag = true ;
        }

//------------------------------------------------------------------------------

        void
        Pump::compute()
        {
            this->compute_isotropic_exit() ;
            this->compute_rotation_specific_speeds() ;
            this->compute_exit_diameters();
            this->compute_entry_diameters();
            this->compute_npsh();
            this->compute_efficiencies();
            this->compute_triangles();
            this->compute_axial_lengths();
        }

//------------------------------------------------------------------------------

        void
        Pump::print()
        {
            std::fprintf( stdout , "Pump Head         H         : %8.3f m     [ %8.3f ft ] \n",
                          ( double ) mH,
                          ( double ) mHft );

            std::fprintf( stdout , "Pump Head         Y         : %8.3f kJ/kg \n",
                          ( double ) mY * 0.001 );

            std::fprintf( stdout , "Volume Flow       Q          : %8.3f l/s  [ %8.3f gal/min ]\n",
                ( double ) mDotV * 1000 ,
                ( double ) mQ) ;

            std::fprintf( stdout , "spec. rot         nq [ ns ]  : %8.3f      [ %8.3f ]\n",
                          ( double ) mNq ,
                          ( double ) mNs) ;

            std::fprintf( stdout , "head rise (US)    psi        : %6.5f \n",
                          ( double ) mPsi ) ;

            std::fprintf( stdout , "flow factor       phi2       : %6.5f \n",
                          ( double ) mPhi2 ) ;

            std::fprintf( stdout , "exit diameter     D2m        : %8.3f mm \n",
                          ( double ) mD2m * 1000 ) ;


            std::fprintf( stdout , "spec. exit diam.  Ds         : %8.3f  \n",
                          ( double ) mDs ) ;

            std::fprintf( stdout , "inner diameter    D1i        : %8.3f mm \n",
                          ( double ) mD1i * 1000 ) ;

            std::fprintf( stdout , "entry diameter    D1m        : %8.3f mm \n",
                          ( double ) mD1m * 1000 ) ;

            std::fprintf( stdout , "suction diameter  D1a        : %8.3f mm \n",
                          ( double ) mD1a * 1000 ) ;

            std::fprintf( stdout , "suction speed     cs         : %8.3f m/s \n",
                          ( double ) mCs ) ;

            std::fprintf( stdout , "exit width        b2         : %8.3f mm \n",
                          ( double ) mB2 * 1000 ) ;

            std::fprintf( stdout , "inducer length    l1         : %8.3f mm\n",
                          ( double ) mInducerLength * 1000) ;

            std::fprintf( stdout , "impeller length   l2         : %8.3f mm\n",
                          ( double ) mImpellerLength * 1000 ) ;


            std::fprintf( stdout , "speed of sound               : %8.3f m/s \n",
                          ( double ) mGas.c( mEntry.Tt(), mEntry.pt() ) ) ;

            std::fprintf( stdout , "suction cross section As     : %8.3f mm^2 \n",
                          ( double ) mAs *1e6 ) ;

            std::fprintf( stdout , "available  head   NPSA_a     : %8.3f m \n",
                          ( double ) mNPSHa ) ;

            std::fprintf( stdout , "required   head   NPSA_r     : %8.3f m \n",
                          ( double ) mNPSHr ) ;

            std::fprintf( stdout , "suction param.    epsilon    : %8.3f \n",
                          ( double ) mEpsilon ) ;

            std::fprintf( stdout , "volumetric eff.   eta_v      : %8.4f  \n",
                          ( double ) mEtaV ) ;


            std::fprintf( stdout , "hydraulic  eff.   eta_h      : %8.4f \n",
                          ( double ) mEtaH ) ;

            std::fprintf( stdout , "mechanic   eff.   eta_m      : %8.4f \n",
                          ( double ) mEtaM ) ;

            std::fprintf( stdout , "isentropic eff.   eta_s      : %8.4f \n",
                          ( double ) mEtaI ) ;

            std::fprintf( stdout , "total      eff.   eta        : %8.4f \n",
                          ( double ) mEta ) ;

            std::fprintf( stdout , "suction speed     c1m        : %8.3f m/s \n",
                          ( double ) mC1 ) ;

            std::fprintf( stdout , "suction speed     u1         : %8.3f m/s \n",
                          ( double ) mU1 ) ;

            std::fprintf( stdout , "suction speed     w1         : %8.3f m/s \n",
                          ( double ) mW1 ) ;

            std::fprintf( stdout , "suction speed     c2m        : %8.3f m/s \n",
                          ( double ) mC2m ) ;

            std::fprintf( stdout , "suction speed     u2         : %8.3f m/s \n",
                          ( double ) mU2 ) ;

            std::fprintf( stdout , "suction speed     w2         : %8.3f m/s \n",
                          ( double ) mW2 ) ;

            std::fprintf( stdout , "Haller            w2/w1      : %8.4f - \n",
                          ( double ) mW2 / mW1 ) ;


            std::fprintf( stdout , "num. bl. inducer   Z1        : %2u\n",
                          ( int ) mZ1 ) ;

            std::fprintf( stdout , "num. bl. impeller  Z2        : %2u\n",
                          ( int ) mZ2 ) ;

            std::fprintf( stdout , "                  beta1      : %8.4f\n",
                          ( double ) mBeta1 ) ;

            std::fprintf( stdout , "                  beta2      : %8.4f\n",
                          ( double ) mBeta2 ) ;

            std::fprintf( stdout , "                  Power      : %8.4f kW\n",
                          ( double ) mP * 0.001 ) ;
            std::fprintf( stdout , "Exit Temperature        Tt2  : %8.4f K\n",
                          ( double ) mExit.Tt()  ) ;
            std::fprintf( stdout , "Exit Pressure           pt2  : %8.4f bar\n",
                          ( double ) mExit.pt() * 1e-5  ) ;
        }

//------------------------------------------------------------------------------
// private
//------------------------------------------------------------------------------

        void
        Pump::compute_isotropic_exit()
        {
            // pressure at exit
            real & tT = mExitIsotropic.value( BELFEM_ENGINE_STATE_TT );
            real & tP = mExitIsotropic.value( BELFEM_ENGINE_STATE_PT ) ;
            real & tH = mExitIsotropic.value( BELFEM_ENGINE_STATE_HT ) ;

            // total pressure at exit
            if ( mDeltaPFlag )
            {
                tP = mEntry.value( BELFEM_ENGINE_STATE_PT ) + mDeltaP;
                mPt2 = tP;
            }
            else if ( mPt2Flag )
            {
                tP = mPt2 ;
                mDeltaP = mPt2 - mEntry.value( BELFEM_ENGINE_STATE_PT ) ;
            }
            else
            {
                BELFEM_ERROR( false, "Must set either deltap or pt2");
            }

            // total temperature at exit
            tT = mGas.isen_T( mEntry.Tt(), mEntry.pt(), tP );

            // enthalpy at exit
            tH = mGas.h( tT, tP );

            // std::cout << "test s " << tT << " " << tH << " " << mGas.s( tT, tP ) << std::endl ;

            // isentropic enthalpy difference
            mY = tH - mEntry.value( BELFEM_ENGINE_STATE_HT );

            // pump head
            mH = mY / constant::g0 ;

            // pump head in ft
            mHft = mH / constant::ft ;
        }

//------------------------------------------------------------------------------

        void
        Pump::compute_rotation_specific_speeds()
        {
            // volume flux im m^3/s
            mDotV = mDotM / mEntry.rho() ;

            // volume flux in gal/min
            mQ = mDotV * 60.0 / constant::gal ;

            // specific rotational speed in metric units
            mNq = mNRPM * std::sqrt( mDotV ) / std::pow( mH, 0.75 );

            // specific rotational speed in US units
            mNs = mNRPM * std::sqrt( mQ ) / std::pow( mHft, 0.75 );

        }

//------------------------------------------------------------------------------

        void
        Pump::compute_exit_diameters()
        {
            BELFEM_ERROR( ! ( mPsiFlag && mD2aFlag ),
                "Must prescribe either psi or D2a, not both!" );

            if( mPsiFlag )
            {
                BELFEM_ERROR( mD2Flag, "Must set D2Flag if psi is prescribed");
                mU2 = std::sqrt( mY / mPsi ) ;
                mD2m = mU2 / ( constant::pi * mN  ) ;
                mD2i = mD2m ;
                mD2a = mD2m ;
                return ;
            }
            else if ( mD2aFlag )
            {
                mDs = mD2a / constant::ft * std::sqrt( std::sqrt( mHft ) / mQ ) ;
            }
            else
            {
                // polynomial based on data from NASA SP-8109 Fig. 3
                // and NASA SP-8107 Fig 20

                mDs = std::exp( polyval( { 4.491201E-02, - 1.480553E+00 , + 6.118206E+00 },
                                         std::log( mNs ) ) );


                // compute outer diameter
                mD2a = std::sqrt( mQ / std::sqrt( mHft )) * mDs * constant::ft;
            }

            if( mD2Flag )
            {
                mD2m = mD2a ;
                mD2i = mD2a ;
            }
            else
            {
                // from Bohl, Fig. 1.35
                mD2m = polyval(  {  1.786920722671952e-09,
                                   -4.659665676034846e-07,
                                   2.271462953190475e-05,
                                   -2.095552800562279e-04,
                                   9.967931912652648e-01 }, mNq ) * mD2a ;

                mD2i = std::sqrt( 2.0 * mD2m * mD2m - mD2a * mD2a ) ;
            }

            mU2 = mD2m * constant::pi * mN ;

            // head rise, US defintion
            mPsi = mY / ( mU2 * mU2 );
        }

//------------------------------------------------------------------------------

        void
        Pump::compute_entry_diameters()
        {
            if( ! mDNFlag )
            {
                mDN = 1.35 * mDw ;
            }

            if( ! mD1aFlag )
            {
                // from Bohl Fig. 1.35
                //mD1a = polyval( { -6.446536E-09, 1.863524E-06, - 2.243131E-04,
                //                  + 1.717627E-02,  + 1.327552E-01}, mNq ) * mD2a ;

                // from Stephanov
                /*mD1a = mD2m * polyval( {   -2.914574810747589e-11,
                                           5.179037219463104e-08,
                                           -3.256064881295752e-05,
                                           9.582249629446555e-03,
                                           2.185152924363399e-01}, mNq ); */

                // chat from Bohl
                mEpsilon = polyval( {  -1.172924e-05, 4.394069e-03, -1.507149e-02  }, mNq );

                mCs = mEpsilon * std::sqrt( 2.0 * mY ) ;

                // cross section at entry
                mAs = mDotV / mCs ;

                mD1a = std::sqrt( 4.0 / constant::pi * mAs + mDN * mDN );

            }
            else
            {

                // cross section at entry
                mAs = 0.25 * constant::pi * ( mD1a * mD1a - mDN * mDN );

                // inlet flow velocity
                mCs = mDotV / mAs;

                mEpsilon = mCs / std::sqrt( 2.0 * mY );
            }
        }

//------------------------------------------------------------------------------

        void
        Pump::compute_npsh()
        {


            // Wesche Eq. ( 3.151 )
            real tVsp = 0.001 * ( 4.5 + 5. * mDN / mD2a ) * mNRPM
                        * std::pow( mNq, 0.8 ) * std::pow( mD2a, 2.5 ) / 3600 ;

            // Wesche ( Eq. 3.10 )
            mEtaV = mDotV / ( tVsp + mDotV );

            // Lomakin equation
            //mEtaV = 1.0 / ( 1.0 + 0.287 / std::pow( mNq, 2./3. ) );

            // vapor pressure
            real tPvap = mGas.eos()->p_vap( mEntry.T() ) ;


            if( mIsHydrogen )
            {
                mKnpsh = 114.7 ;

                // from fig. 27 NASA SP-8107
                mBetas = std::pow( 3.086791E-08 * tPvap, 1.538770E+00 );

            }
            else
            {
                mKnpsh = 6.35;

                // vapor state
                real tTv = mEntry.T() + 1e-6;
                real tHv = mGas.h( tTv, tPvap );
                real tRhov = mGas.rho( tTv, tPvap );

                // liquid state
                real tTl = mEntry.T() - 1e-6;
                real tHl = mGas.h( tTl, tPvap );
                real tCpl = mGas.cp( tTl, tPvap );
                real tRhol = mGas.rho( tTl, tPvap );


                // beta-parameter
                mBetas = tRhov * ( tHv - tHl ) / ( tRhol * tCpl );
            }

            // scale unit of K
            mKnpsh *=  1.8  * std::pow( 1. / constant::ft, 0.16 )  / std::pow( 1.0 / constant::ft, 1.15 ) ;

            // available suction head
            mNPSHa = ( mEntry.p() - tPvap ) / ( mEntry.rho() * constant::g0 ) ;

            // std::cout << " k " << mKnpsh << std::endl ;

            if( ! mD1aFlag )
            {
                // check for safety value
                real tF0 = this->compute_npshr( mD1a );

                if( tF0 < 0.0 )
                {
                    // find better diameter
                    real tX1 = mD1a ;
                    real tX0 = 0 ;

                    real tF = tF0 ;
                    while( tF < 0.0001 )
                    {
                        tF0 = tF ;
                        tX0 = tX1 ;
                        tX1 += 0.001 ;
                        tF = this->compute_npshr( tX1 );
                    }

                    uint tCount = 0 ;

                    real tX = 0 ;
                    while( std::abs( tF ) > 0.0001 )
                    {
                        tX = 0.5 * ( tX0 + tX1 );
                        tF = this->compute_npshr( tX );

                        if ( tF0 * tF  > 0 )
                        {
                            tX0 = tX;
                            tF0 = tF;
                        }
                        else
                        {
                            tX1 = tX;
                        }

                        ++tCount;
                        BELFEM_ERROR( tCount < 1000, "Too many iterations" );
                    }
                    mD1a = tX ;
                    this->compute_npshr( tX );
                }
            }
            else
            {
                // do not check for NPSH
                this->compute_npshr( mD1a );
            }

            // NASA SP-8107 Eq. ( 6 )
            //mSs = mNRPM * std::sqrt( mQ ) / std::pow( mNPSHr / constant::ft, 0.75 ) ;

            //real tNu = mDN / mD1a ;

            // NASA SP-8107 Eq. ( 7 )
            //mSscorr = mSs / std::sqrt( 1.0 - tNu * tNu );

            // Wesche Eq. 3.43
            /*real tW1 = std::sqrt( mC1 * mC1 + tU1a * tU1a );
            real tNPSHr = ( 1.15 * mC1 * mC1 + 0.14 * tW1 * tW1 ) / ( 2.0 * constant::g0 ) ;
            std::cout << "NPSHr " << tNPSHr << std::endl ; */


        }

//------------------------------------------------------------------------------

        real
        Pump::compute_npshr( const real & aD1a )
        {

            // cross section at entry
            mAs = 0.25 * constant::pi * ( aD1a * aD1a - mDN * mDN );

            // inlet flow velocity
            mCs = mDotV / mAs;

            // axial velocity C1
            mC1 = mCs / mEtaV ;

            // Flow coefficient
            real tU1a = constant::pi * aD1a * mN ;
            mPhi1a = mC1 / tU1a ;

            // remember suction parameter
            mEpsilon = mCs / std::sqrt( 2.0 * mY );

            // from  NASA SP-8107
            mNPSHr = ( 0.931 / std::pow( mPhi1a, 4./9. )
                       - mKnpsh * std::pow( aD1a  / mZ1, 0.16 ) * mBetas  /
                         ( mPhi1a * mPhi1a * std::pow( tU1a , 1.15 ) ) )
                     * mC1 * mC1 / ( 2.0 * constant::g0 ) ;

            return mNPSHa - mNPSHr - mNPSHs ;
        }

//------------------------------------------------------------------------------

        void
        Pump::compute_triangles()
        {
            // From Bohl, Fig. 1.35
            mD1i = polyval( {  4.745642E-09, - 8.884172E-07, 2.449606E-05,
                               1.763717E-03, 3.154522E-01} , mNq ) * mD2a ;

            // mean diameter
            mD1m = std::sqrt( 0.5 * ( mD1a * mD1a + mD1i * mD1i ) ) ;

            // velocity around 1
            mU1 = mD1m * constant::pi * mN ;

            mW1 = std::sqrt( mC1 * mC1 + mU1 * mU1 );

            mBeta1 = std::atan( mC1 / mU1 ) / constant::deg ;


            real tBeta2 ;

            if( ! mHallerFlag )
            {
                if ( mBeta2Flag )
                {
                    tBeta2 = mBeta2 * constant::deg;
                }
                else // guess value based on graph
                {
                    // from Bohl
                    tBeta2 = std::atan(
                            polyval( { -2.765261544046292e-10, -1.060736960995682e-07, 5.023919097310044e-05,
                                       -5.759522694014521e-03, 5.077866712104788e-01 }, mNq ));

                    mBeta2 = tBeta2 / constant::deg;
                }

                mHaller = this->compute_haller( tBeta2 );

            }
            else
            {
                real tDX = 1.0 * constant::deg ;

                real tX1 = 1.0 * constant::deg ;

                real tF = this->compute_haller( tX1 ) - mHaller ;

                real tF0 = 0 ;
                real tX0 = 0 ;

                // scan interval
                while( tX1 < 0.5 * constant::pi )
                {
                    tX0 = tX1 ;
                    tF0 = tF ;
                    tX1 += tDX ;
                    tF = this->compute_haller( tX1 ) - mHaller ;

                    // check for sign change
                    if ( tF0 * tF < 0 )
                    {
                        break;
                    }
                }

                // counter to avoid infinite loop
                uint tCount = 0 ;

                real tX = tX1 ;

                while ( std::abs( tF ) > 1e-6 )
                {
                    tX = 0.5 * ( tX0 + tX1 );
                    tF = this->compute_haller( tX ) - mHaller ;

                    if ( tF * tF0 > 0 )
                    {
                        tX0 = tX ;
                        tF0 = tF ;
                    }
                    else
                    {
                        tX1 = tX ;
                    }

                    BELFEM_ERROR( tCount++ < 100, "Too many iterations");
                }

                // final run to store value
                mHaller = this->compute_haller( tX );
            }



            if ( !mZ2Flag )
            {
                // from NASA SP-8109 Figure 16

                Vector< real > tC = {
                        -2.806934206930023e+00,
                        9.177720942575908e+00,
                        1.557445124240440e+01,
                        -7.346018008664592e+01,
                        3.288358055173669e+01,
                        -2.924640112867317e+01,
                        1.279785039686524e+02,
                        6.570150589745129e+00,
                        -4.548177282822825e+01,
                        2.673403497649506e+01 };

                Vector< real > tP( 10 );

                tP( 0 ) = 1.0;
                tP( 1 ) = mPhi2;
                tP( 2 ) = mPsi;
                tP( 3 ) = mPhi2 * mPhi2;
                tP( 4 ) = mPhi2 * mPsi;
                tP( 5 ) = mPsi * mPsi;
                tP( 6 ) = mPhi2 * mPhi2 * mPhi2;
                tP( 7 ) = mPhi2 * mPhi2 * mPsi;
                tP( 8 ) = mPhi2 * mPsi * mPsi;
                tP( 9 ) = mPsi * mPsi * mPsi;

                mZ2 = std::max( std::round( std::exp( dot( tC, tP ) )), 3.0 );

                // hack for oxygen pump
                if ( mGas.component( 0 )->label() == "O2" )
                {
                    mZ2 *= 2 ;
                }

                // check if Z1 was prescribed
                if( ! mZ1Flag )
                {
                    // correct both Z1 and Z2
                    index_t tIndex = int( std::round( mZ2 ) );

                    if( tIndex < mZ1Table.length() )
                    {
                        mZ1 = mZ1Table( tIndex );
                    }
                    if ( tIndex < mZ2Table.length() )
                    {
                        mZ2 = mZ2Table( tIndex );
                    }
                }
            }


            mYthinf = mC2u * mU2 ;

            mMu = mYth / mYthinf ;


            // Bohl ( 1.35 )
            // mC2m = std::tan( tBeta2 ) * ( mU2 - mC2u );


            real tV = mDotM / mExit.rho();

            // Bohl ( 2.2 )
            real tK2 = mD2a * constant::pi /
                       ( mD2a * constant::pi - mZ2 * mS2 / std::sin( mBeta2 ) );

            mB2 = tV * tK2 / ( mD2a * constant::pi * mC2m );
        }

//------------------------------------------------------------------------------

        real
        Pump::compute_haller( const real & aBeta2 )
        {
            // remember angle
            mBeta2 = aBeta2 / constant::deg ;

            // From Stephanof chart
            mC2m = mStepanoff.Km2( mNs ) * std::sqrt( 2.0 * mY );

            mPhi2 = mC2m / mU2 ;

            mC2u = mU2 - mC2m / std::tan( aBeta2 )  ;
            mW2 = std::sqrt( mC2m * mC2m + ( mU2 - mC2u ) * ( mU2 - mC2u ) ) ;

            return mW2 / mW1 ;
        }

//------------------------------------------------------------------------------

        void
        Pump::compute_static( State & aState, const real & aA )
        {
            // relaxation factor
            real tOmega = 0.1 ;

            const real & tTt  = aState.value( BELFEM_ENGINE_STATE_TT );
            const real & tPt  = aState.value( BELFEM_ENGINE_STATE_PT );
            real & tHt  = aState.value( BELFEM_ENGINE_STATE_HT );
            real & tSt  = aState.value( BELFEM_ENGINE_STATE_S );

            real & tT   = aState.value( BELFEM_ENGINE_STATE_T );
            real & tP   = aState.value( BELFEM_ENGINE_STATE_P );
            real & tH   = aState.value( BELFEM_ENGINE_STATE_H );
            real & tU   = aState.value( BELFEM_ENGINE_STATE_U );
            real & tMa  = aState.value( BELFEM_ENGINE_STATE_MA );
            real & tRho = aState.value( BELFEM_ENGINE_STATE_RHO );

            // entropy
            tHt = mGas.h( tTt, tPt );
            tSt = mGas.s( tTt, tPt );

            /*real tS ;

            // the Jacobian
            Matrix< real > tJ( 2, 2 );

            // the RHS
            Vector< real > tX( 2 );

            // pivot for gesv
            Vector< int > tPivot( 2 ); */

            // guess values
            tT = tTt ;
            tP = tPt ;

            real tError = BELFEM_REAL_MAX ;
            uint tCount = 0 ;

            real tDotV ;

            while( tError > 1e-9 )
            {
                // density at entry
                tRho = mGas.rho( tT, tP );

                // volume flux
                tDotV = mDotM / tRho ;

                // velocity
                tU = tDotV / aA ;

                // enthalpy
                tH = mGas.h( tT, tP );

                // temperature
                tT -= tOmega * ( tH + 0.5 * tU * tU - tHt ) / mGas.cp( tT, tP );

                tP = mGas.isen_p( tTt, tPt, tT );

                tError = std::abs( ( tH + 0.5 * tU * tU - tHt ) / tHt ) ;

                BELFEM_ERROR( tCount++ < 1000, "Too many iterations" );
            }

            tMa = tU / mGas.c( tT, tP );
        }

//------------------------------------------------------------------------------

        void
        Pump::compute_efficiencies()
        {
            // Wesche 3.56
            // Wesche, Eq. (3.56)
            real tEtah0 = mAhydro - 3.0 /
                                    std::pow( std::log10( 5e7 * mDotV * 3600 / mNRPM ) , 2 );


            // Hydrodonamic efficiency, Wesche Eq. (3.55)
            mEtaH = tEtah0 - 1.4e-6 * mNq * mNq / std::pow( tEtah0, 21. );

            // theoretical height, Wesche Eq. (3.11)
            mYth = mY / mEtaH ;

            // Hydraulic Power, Wesche Eq. (3.7 )
            mPh = mDotM * mYth / mEtaV ;

            // friciton power, Wesche Eq. ( 3.180 )
            mPr = 80.0 * mEntry.rho() * std::pow( mNRPM * 0.001, 3. )
                    * std::pow( mD2a, 5. );

            real tP = 0.0 ;

            // Wesche, Eq. ( 3.6 )
            mPs = mPh + mPr ;
            mP = mPs ;

            uint tCount = 0 ;
            while( std::abs( mP - tP ) > 1e-6 )
            {
                // shift power
                tP = mP ;

                // Wesche Eq. ( 3.57 )
                mPm = 200. * std::sqrt( mP * 0.001 );

                mP = mPm + mPs ;
                BELFEM_ERROR( ++tCount < 100, "Too many iterations" );
            }

            // Wesche, Eq. ( 3.5 )
            mEtaI = mDotM * mY / mPs ;

            // Wesche, Eq. ( 3.15 )
            mEtaM = mPs / mP ;

            // Wesche, Eq. ( 3.14 )
            mEta = mEtaI * mEtaM ;

            // compute exit enthalpy
            const real & tH1   = mEntry.ht();
                  real & tT2   = mExit.value( BELFEM_ENGINE_STATE_TT );
                  real & tP2   = mExit.value( BELFEM_ENGINE_STATE_PT );
                  real & tH2   = mExit.value( BELFEM_ENGINE_STATE_HT );
                  real & tRho2 = mExit.value( BELFEM_ENGINE_STATE_RHO );
            tP2 = mExitIsotropic.pt();

            tH2 = tH1 + mPs / mDotM ;
            tT2 = mGas.T_from_h( tH2, tP2 );
            tRho2 = mGas.rho( tT2, tP2 );
        }

//------------------------------------------------------------------------------
/*        void
        Pump::set_psi( const real & aPsi )
        {
            mPsi = aPsi ;
            mD2 = 0.45 / mN * std::sqrt( mY / ( 2.0 * aPsi ) ) ;

            std::cout << "D2 " << mD2 *1000 << std::endl ;

            mU2 = mD2 * constant::pi * mN ;
        }

//------------------------------------------------------------------------------

        void
        Pump::set_z( const uint aZ )
        {
            mZ = ( real ) aZ ;

        }

//------------------------------------------------------------------------------

        void
        Pump::set_beta2( const real aBeta2 )
        {
            mBeta2 = aBeta2 * constant::deg ;
        }

//------------------------------------------------------------------------------

        void
        Pump::set_D2( const real & aD2 )
        {
            mD2 = aD2 ;

            mU2 = mD2 * constant::pi * mN ;

            mPsi = mY / ( mU2 * mU2 );

            std::cout << "psi " << mPsi << std::endl ;
        }

//------------------------------------------------------------------------------

        void
        Pump::compute_static( State & aState, const real & aA )
        {
            // relaxation factor
            real tOmega = 0.9 ;

            const real & tTt  = aState.value( BELFEM_ENGINE_STATE_TT );
            const real & tPt  = aState.value( BELFEM_ENGINE_STATE_PT );
                  real & tHt  = aState.value( BELFEM_ENGINE_STATE_HT );
                  real & tSt  = aState.value( BELFEM_ENGINE_STATE_S );

                  real & tT   = aState.value( BELFEM_ENGINE_STATE_T );
                  real & tP   = aState.value( BELFEM_ENGINE_STATE_P );
                  real & tH   = aState.value( BELFEM_ENGINE_STATE_H );
                  real & tU   = aState.value( BELFEM_ENGINE_STATE_U );
                  real & tMa  = aState.value( BELFEM_ENGINE_STATE_MA );
                  real & tRho = aState.value( BELFEM_ENGINE_STATE_RHO );

            // entropy
            tHt = mGas.h( tTt, tPt );
            tSt = mGas.s( tTt, tPt );

            real tS ;

            // the jacobian
            Matrix< real > tJ( 2, 2 );

            // the RHS
            Vector< real > tX( 2 );

            // pivot for gesv
            Vector< int > tPivot( 2 );

            // guess values
            tT = tTt ;
            tP = tPt ;

            real tError = BELFEM_REAL_MAX ;
            uint tCount = 0 ;

            real tDotV ;

            while( tError > 1e-9 )
            {
                // density at entry
                tRho = mGas.rho( tT, tP );

                // volume flux
                tDotV = mDotM / tRho;

                // dhdT
                tJ( 0, 0 ) = mGas.cp( tT, tP );

                // dsdT
                tJ( 1, 0 ) = mGas.dsdT( tT, tP );

                // dhdP
                tJ( 0, 1 ) = mGas.dhdp( tT, tP );

                // dsdP
                tJ( 1, 1 ) = mGas.dsdp( tT, tP );

                tH = mGas.h( tT, tP );
                tU = tDotV / aA;
                tS = mGas.s( tT, tP );

                // rhs
                tX( 0 ) = tH + 0.5 * tU * tU - tHt;
                tX( 1 ) = tS - tSt;

                // compute error
                tError = std::sqrt( tX( 0 ) * tX( 0 ) + tX( 1 ) * tX( 1 ) );

                // solve system
                gesv( tJ, tX, tPivot );

                tT -= tOmega * tX( 0 );
                tP -= tOmega * tX( 1 );

                BELFEM_ERROR( tCount++ < 100, "Too many iterations" );
            }

            tMa = tU / mGas.c( tT, tP );
        }

//------------------------------------------------------------------------------

        void
        Pump::compute_entry_static( const real & aD0, const real & aD1 )
        {


            // cross section at entry
            real & tA = mEntry.value( BELFEM_ENGINE_STATE_A );
            tA = 0.25 * constant::pi * ( aD1 * aD1 - aD0 * aD0 );
            mDn   = aD0 ;
            mD1   = aD1 ;

            this->compute_static( mEntry, tA );

            // volume flux
            mDotV = mDotM / mEntry.rho() ;

            std::cout << "Ma " << mEntry.value( BELFEM_ENGINE_STATE_MA ) << std::endl ;
        }


//------------------------------------------------------------------------------



//------------------------------------------------------------------------------

        void
        Pump::compute_npsh()
        {
            // entry data
            // mC1 = mStepanoff.Km1( mNs ) * std::sqrt( 2.0 * mY ) ;
            // mC1 = mEntry.u() ;

            // std::cout << "Km1 " << mStepanoff.Km1( mNs ) << std::endl ;

            mC1 = mEntry.u() / mEtaV ;

            mU1 = mD1 * constant::pi * mN ;
            mW1 = std::sqrt( mC1 * mC1 + mU1 * mU1 );

            mPhi1 = mC1 / mU1 ;
            std::cout << "cs " << " " << mEntry.u() << std::endl ;
            std::cout << "c1 " << " " << mC1 << std::endl ;
            std::cout << "epsilon " << " " << mEntry.u() / std::sqrt( 2.0 * mY ) << std::endl ;
            std::cout << "u1 " << " " << mU1 << std::endl ;
            std::cout << "w1 " << " " << mW1 << std::endl ;
            std::cout << "phi1 " << " " << mPhi1 << std::endl ;

            real tPvap = mGas.eos()->p_vap( mEntry.T() ) ;

            const real & tT = mEntry.T() ;
            real tPv = tPvap * 0.9999 ;
            real tPl = tPvap * 1.0001 ;
            real tRhoL = mGas.rho( tT, tPl );
            real tHl   = mGas.h( tT, tPl );
            real tCpl  = mGas.cp( tT, tPl );
            real tRhoV = mGas.rho( tT, tPv );
            real tHv   = mGas.h( tT, tPv );

            bool tIsHydrogen = ( mGas.number_of_components() == 1 && mGas.component( 0 )->label() == "H2" ) ;

            BELFEM_ERROR( ! tIsHydrogen , "need to implement NPSH for hydrogen" );

            // this equation is for non hydrogen only
            real tBeta = ( tRhoV / tRhoL ) * ( tHv - tHl ) / tCpl ;

            const real tK = tIsHydrogen ? 77.01372 : 4.26366 ;

            std::cout << "beta " << tBeta << std::endl ;

            mNPSHr = ( 0.931 / std::pow( mPhi1, 4./9. )
                       - tK * std::pow( mD1 / mZ0, 0.16 ) * tBeta /
                         ( mPhi1 * mPhi1 * std::pow( mU1, 1.15 ) ) )
                     * mC1 * mC1 / ( 2.0 * constant::g0 ) ;

            std::cout << "NPSHr " << mNPSHr << std::endl ;

            real tHvap  = tPvap / ( mEntry.rho() * constant::g0 );

            //  Rocket Propulsion Elements, G. Sutton:Eq. 10.8
            mNPSHa = mEntry.p() / ( mEntry.rho() * constant::g0 ) - tHvap ;

            std::cout << "NPSHa " << mNPSHa << std::endl ;
        }

//------------------------------------------------------------------------------

        void
        Pump::compute_similarity()
        {
            mH = mY / constant::g0 ;

            real tH = mH / constant::ft ;

            std::cout << "Y " << mY << std::endl ;
            std::cout << "H " << mH << std::endl ;
            std::cout << "Hft " << tH << std::endl ;
            std::cout << "H(p) " << ( mDeltaP ) / ( mEntry.rho() * constant::g0 ) << std::endl ;
            mQ = mDotV * 60.0 / constant::gal ;
            std::cout << "dotV " << mDotV << std::endl ;
            std::cout << "Q " << mQ << std::endl ;

            mNq = mNRPM * std::sqrt( mDotV ) / std::pow( mH, 0.75 ) ;
            mNs = mNRPM * std::sqrt( mQ ) / std::pow( tH, 0.75 ) ;

            std::cout << "Nq " << mNq << std::endl ;
            std::cout << "Ns " << mNs << std::endl ;

            real tD2 = mD2 / constant::ft ;

            mDs = tD2 * std::pow( tH, 0.25 ) / std::sqrt( mQ ) ;

            std::cout << "Ds " << mDs << std::endl ;
            std::cout << "psi " << mPsi << std::endl ;



            // Wesche Eq. ( 3.151 )
            real tVsp = 1e-3 * ( 4.5 + 5. * mDn / mD2 ) * mNRPM
                   * std::pow( mNq, 0.8 ) * std::pow( mD2, 2.5 ) / 3600 ;

            // Wesche ( Eq. 3.10 )
            mEtaV = mDotV / ( tVsp + mDotV );



            // real tVsp = mDotV / mEtaV - mDotV ;
            // volumetric efficiency after LOMAKIN
            //mEtaV = 1.0 / ( 1.0 + 0.287 / std::pow( mNq, 2./3. ) );
            //real tVsp = mDotV / mEtaV - mDotV ;
            std::cout << "eta_V " << mEtaV << std::endl ;

            // Wesche 3.56
            // Wesche, Eq. (3.56)
            real tEtah0 = mAhydro - 3.0 /
                                    std::pow( std::log10( 5e7 * mDotV * 3600 / mNRPM ) , 2 );
            std::cout << "VN " <<  mDotV * 3600 / mNRPM << std::endl ;
            // Hydrodonamic efficiency, Wesche Eq. (3.55)
            mEtaH = tEtah0 - 1.4e-6 * mNq * mNq / std::pow( tEtah0, 21 );
            std::cout << "eta_H " << mEtaH << std::endl ;

            // Wesche (Eq. 3.180)
            mPr = 80.0 * mEntry.rho() * std::pow( 0.001 * mNRPM, 3. )
                    * std::pow( mD2, 5. ) ;


            // Wesche (Eq. 3.11 )
            real mYth = mY / mEtaH ;


            // Wesche Eq. (3.7)
            mPh = mYth * mEntry.rho() * ( mDotV + tVsp ) ;

            // Wesche Eq. (3.6)
            mPi = mPh + mPr ;

            real tP = 0.0 ;
            mP = mPi ;

            while( abs( mP - tP ) / mP > 1e-6 )
            {
                tP = mP ;
                // Wesche Eq. ( 3.57 )
                mPm = 200. * std::sqrt( mP * 0.001 );
                mP = mPm + mPi ;
            }
            // Wesche Eq. ( 3.15 )
            mEtaM = mPi / mP ;


            std::cout << "eta_m " << mEtaM << std::endl ;


            std::cout << "P_r " << mPr << std::endl ;

            // Wesche Eq. ( 3.8 )
            mEtaI = mEntry.rho() * mY * mDotV / ( mEntry.rho() * mYth * ( mDotV + tVsp ) + mPr ) ;
            std::cout << "eta_I " << mEtaI << std::endl ;

            // Wesche Eq. ( 3.14 )
            mEta = mEtaI * mEtaM ;
            std::cout << "eta " << mEta << std::endl ;

            // Bohl 2 ( 2.4 )
            real tPsiTh = 2.0 * mYth / ( mU2 * mU2 );
            std::cout << "beta2 " << mBeta2 / constant::deg << std::endl ;
            mMu = 1.0 / ( 1.0 + 2.0 * constant::pi * std::sin( mBeta2 ) / ( mZ * tPsiTh ) );
            std::cout << "mu " << mMu << std::endl ;

            mYthinf = mYth / mMu ;
        }

//------------------------------------------------------------------------------

        void
        Pump::compute_exit( const real & aT2, const real & aB2 )
        {
            real & tTt = mExit.value( BELFEM_ENGINE_STATE_TT );
            real & tPt = mExit.value( BELFEM_ENGINE_STATE_PT ) ;
            real & tHt = mExit.value( BELFEM_ENGINE_STATE_HT ) ;

            tTt = aT2 ;
            tPt = mExitIsotropic.pt() ;

            // surface at exit
            real tA = aB2 * mD2 * constant::pi ;


            this->compute_static( mExit, tA );

            // actual work
            mA = tHt - mEntry.ht() ;

            std::cout << "a " << mA << std::endl ;
            std::cout << "eta_is " << mY / mA << std::endl ;

            // Euler
            mC2u = mYthinf / mU2 ;
            mW2u = mU2 - mC2u ;
            mC2m =  mStepanoff.Km2( mNs ) * std::sqrt( 2. * mY ) ;

            std::cout << "Km2 " << mStepanoff.Km2( mNs ) << std::endl ;
            std::cout << "u2  " << mU2 << std::endl ;
            std::cout << "c2u " << mC2u << std::endl ;
            std::cout << "w2u " << mW2u << std::endl ;
            std::cout << "c2m " << mC2m << std::endl ;
            std::cout << "c2m* " << mExit.u() << std::endl ;

            // Bohl 1.33
            real tBeta1 = std::atan( mC1 / mU1 ) ;

            // Bohl 1.35
            real tBeta2 = std::atan( mC2m / ( mU2 - mC2u ) ) ;

            std::cout << "T2 " << mExit.T() << std::endl ;
            std::cout << "p2 " << mExit.p() *1e-5 << std::endl ;
            std::cout << "Ma2 " << mExit.Ma() << std::endl ;
            std::cout << "beta1 " << tBeta1 / constant::deg << std::endl ;
            std::cout << "beta2 " << tBeta2 / constant::deg << std::endl ;

        } */


//------------------------------------------------------------------------------

        void
        Pump::compute_axial_lengths()
        {
            // help constant
            real tX = std::sqrt( mAs );

            // LOX pump
            if ( mGas.number_of_components()
                == 1 && mGas.component( 0 )->label() == "O2" )
            {
                mInducerLength   = 0.6189 * tX ;
                mImpellerLength  = 0.3675 * tX ;
            }
            else // assume values for LCH4
            {
                mInducerLength  = 0.7181 * tX ;
                mImpellerLength = 0.4588 * tX ;
            }

        }

//------------------------------------------------------------------------------
    }
}