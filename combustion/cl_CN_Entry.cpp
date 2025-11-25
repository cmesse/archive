//
// Created by Christian Messe on 26.09.19.
//

#include "cl_CN_Entry.hpp"
#include "assert.h"
namespace belfem
{
    namespace combustion
    {
//------------------------------------------------------------------------------
// public
//------------------------------------------------------------------------------
        bool
        Entry::has_low() const
        {
            return mLow.length() > 0;
        }

//------------------------------------------------------------------------------

        bool
        Entry::has_troe() const
        {
            return mTroe.length() > 0;
        }

//------------------------------------------------------------------------------

        bool
        Entry::has_third_body_weights() const
        {
            return mThirdBodyWeights.length() > 0;
        }

//------------------------------------------------------------------------------

        const string &
        Entry::reaction() const
        {
            return mReaction;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Entry::coeffs() const
        {
            return mCoeffs;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Entry::low() const
        {
            return mLow;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Entry::troe() const
        {
            return mTroe;
        }

//------------------------------------------------------------------------------

        const Cell< string > &
        Entry::third_body_species() const
        {
            return mThirdBodySpecies;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Entry::third_body_weights() const
        {
            return mThirdBodyWeights;
        }

//------------------------------------------------------------------------------

        bool
        Entry::is_duplicate() const
        {
            return mDuplicateFlag;
        }

//------------------------------------------------------------------------------

        bool
        Entry::is_active() const
        {
            return mActiveFlag;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Entry::duplicate()  const
        {
            return mDuplicate;
        }

//------------------------------------------------------------------------------
// protected
//------------------------------------------------------------------------------

        string &
        Entry::reaction()
        {
            return mReaction;
        }

//------------------------------------------------------------------------------

        Vector< real > &
        Entry::coeffs()
        {
            return mCoeffs;
        }

//------------------------------------------------------------------------------

        Vector< real > &
        Entry::low()
        {
            return mLow;
        }

//------------------------------------------------------------------------------

        Vector< real > &
        Entry::troe()
        {
            return mTroe;
        }

//------------------------------------------------------------------------------

        Cell< string > &
        Entry::third_body_species()
        {
            return mThirdBodySpecies;
        }

//------------------------------------------------------------------------------

        Vector< real > &
        Entry::third_body_weights()
        {
            return mThirdBodyWeights;
        }

//------------------------------------------------------------------------------

        void
        Entry::set_duplicate_flag()
        {
            mDuplicateFlag = true;
        }

//------------------------------------------------------------------------------

        void
        Entry::set_duplicate( const Vector< real > & aDuplicate )
        {
            BELFEM_ASSERT( mDuplicate.length() == 0,
                "tried to link reaction %s with duplicate, but duplicat pointer is already filled. A duplicate reaction must exist twice and only twice.",
                          mReaction.c_str() );

            mDuplicate = aDuplicate;
        }

//------------------------------------------------------------------------------

        void
        Entry::deactivate()
        {
            mActiveFlag = false;
        }

//------------------------------------------------------------------------------
    }
}