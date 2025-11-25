//
// Created by Christian Messe on 25.02.20.
//

#ifndef BELFEM_CL_CH_CHANNELODE_HPP
#define BELFEM_CL_CH_CHANNELODE_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

#include "cl_ODE.hpp"
#include "cl_CH_Geometry.hpp"
#include "cl_Gas.hpp"

namespace belfem
{

    namespace channel
    {
        class Element ;

        enum class ChannelMode
        {
            Channel,
            Combustor
        };

        class ChannelODE : public ode::ODE
        {
            const ChannelMode mMode ;

            // linked geometry
            Geometry * mGeometry = nullptr ;

            // linked gas
            Gas      & mGas;

            //
            bool mElementMode ;

            // wall temperature
            real mTwall = BELFEM_QUIET_NAN ;

            // Jacobian Matrix
            Matrix< real > mJacobi;

            // Pivot Vector for gesv
            Vector< int > mPivot;

            // value for dRdx/R ( combustion only )
            real mdRdxR = 0.0 ;

            // value for combustion heat
            real mdwdx = 0.0 ;

            // value for mass flux insertion
            real mdMdxM = 0.0 ;

            // valze for momemtum change
            real mdIdx = 0.0 ;

            // pointer to channel element, if used
            Element * mElement = nullptr ;

            // work vectors for element interpolation
            Vector< real > mWorkN ;
            Vector< real > mWorkV ;

            Vector< real > mdYdx ;

            // cross section
            real mA ;

            // cross section
            real mdAdx ;

            // hydraulic diameter
            real mDh ;

            bool mReverse = false ;

            // flag telling if combustion is set
            // ( automatically activated if set_combusiton is called )
            bool mCombust = false ;

//------------------------------------------------------------------------------

            // Funciton pointer for computation object
            void
            ( ChannelODE:: * mComputeFunction )
            (
                    const real          & aT,
                    const Vector <real> & aY,
                          Vector <real> & adYdT ) ;

//------------------------------------------------------------------------------

            // Function pointer for geometry object
            void
            ( ChannelODE:: * mGeometryFunction )
            ( const real & aX, real & aDh, real & aA, real & adAdX );

//------------------------------------------------------------------------------

            // Function pointer for compute Jacobi
            void
            ( ChannelODE:: * mJacobiFunction )
            ( const real & aV,
              const real & aU,
              const real & aT,
              const real & aP );

//------------------------------------------------------------------------------

            // Function pointer for friction
            void
            ( ChannelODE:: * mFrictionFunction )
                    ( const real & aX,
                      const real & aV,
                      const real & aU,
                      const real & aT,
                      const real & aP,
                      const real & aTwall,
                            real & aTauW,
                            real & aDotQ );

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            // constructor
            ChannelODE( Geometry & aGeometry, Gas & aGas,
                        const ChannelMode aMode = ChannelMode::Channel ) ;

//------------------------------------------------------------------------------

            // constructor
            ChannelODE( Gas & aGas, const ChannelMode aMode = ChannelMode::Channel ) ;

//------------------------------------------------------------------------------

            // destructor
            virtual
            ~ChannelODE() = default ;

//------------------------------------------------------------------------------

            virtual void
            compute(
                    const real          & aT,
                    const Vector <real> & aY,
                    Vector <real>       & adYdT ) ;

//------------------------------------------------------------------------------

            // set the wall temperature
            void
            set_wall_temperature( const real & aTw );

//------------------------------------------------------------------------------

            void
            set_combustion( const real & adRdxR, const real & adwdx );

//------------------------------------------------------------------------------

            void
            link_element( Element * aElement );

//------------------------------------------------------------------------------

            void
            link_geometry( Geometry * aGeometry, const bool aReverse=false );

//------------------------------------------------------------------------------

            void
            set_composition_change( const real & adRdxR, const Vector< real >  & adYdx );

//------------------------------------------------------------------------------

            /**
             * returns the operating mode of this ode
             */
             const ChannelMode &
             mode() const ;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            compute_channel_ode( const real          & aT,
                                 const Vector <real> & aY,
                                 Vector <real>       & adYdT );

//------------------------------------------------------------------------------

            void
            compute_combustor_ode( const real        & aT,
                                 const Vector <real> & aY,
                                 Vector <real>       & adYdT );

//------------------------------------------------------------------------------

            void
            compute_geometry_from_geometry_object(
                    const real & aX,
                          real & aDh,
                          real & aA,
                          real & adAdX );

//------------------------------------------------------------------------------

            void
            compute_geometry_from_geometry_object_reversex(
                    const real & aX,
                    real & aDh,
                    real & aA,
                    real & adAdX );

//------------------------------------------------------------------------------

            void
            compute_geometry_from_element(
                    const real & aX,
                          real & aDh,
                          real & aA,
                          real & adAdX );

//------------------------------------------------------------------------------

            void
            compute_jacobi_idgas(
                    const real & aV,
                    const real & aU,
                    const real & aT,
                    const real & aP );

//------------------------------------------------------------------------------

            void
            compute_jacobi_realgas(
                    const real & aV,
                    const real & aU,
                    const real & aT,
                    const real & aP ) ;

//------------------------------------------------------------------------------

            void
            compute_friction_dittus_boelter(
                    const real & aX,
                    const real & aV,
                    const real & aU,
                    const real & aT,
                    const real & aP,
                    const real & aTwall,
                    real & aTauW,
                    real & aDotQ ) ;

//------------------------------------------------------------------------------

            void
            compute_friction_element(
                    const real & aX,
                    const real & aV,
                    const real & aU,
                    const real & aT,
                    const real & aP,
                    const real & aTwall,
                    real & aTauW,
                    real & aDotQ ) ;

//------------------------------------------------------------------------------

            void
            compute_combustion( const real & aT, const real & aP, Vector<real> & aRHS );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------

        inline const ChannelMode &
        ChannelODE::mode() const
        {
            return mMode ;
        }

//------------------------------------------------------------------------------
    } /* namespace channel */
}  /* namespace ode */


#endif //BELFEM_CL_CH_CHANNELODE_HPP
