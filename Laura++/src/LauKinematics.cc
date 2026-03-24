
/*
Copyright 2004 University of Warwick

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

/*
Laura++ package authors:
John Back
Paul Harrison
Thomas Latham
*/

/*! \file LauKinematics.cc
    \brief File containing implementation of LauKinematics class.
*/

#include "LauKinematics.hh"

#include "LauConstants.hh"
#include "LauRandom.hh"

#include "TF2.h"
#include "TMath.h"
#include "TRandom.h"

#include <iostream>

LauKinematics::LauKinematics( const Double_t m1,
                              const Double_t m2,
                              const Double_t m3,
                              const Double_t mParent,
                              const Bool_t calcSquareDPCoords,
                              const Bool_t symmetricalDP,
                              const Bool_t fullySymmetricDP ) :
    symmetricalDP_( symmetricalDP ),
    fullySymmetricDP_( fullySymmetricDP ),
    m1_( m1 ),
    m2_( m2 ),
    m3_( m3 ),
    mParent_( mParent ),
    m1Sq_( m1 * m1 ),
    m2Sq_( m2 * m2 ),
    m3Sq_( m3 * m3 ),
    mParentSq_( mParent * mParent ),
    mDTot_( m1 + m2 + m3 ),
    massInt_( mParent - ( m1 + m2 + m3 ) ),
    mSqDTot_( m1 * m1 + m2 * m2 + m3 * m3 ),
    m12_( 0.0 ),
    m23_( 0.0 ),
    m13_( 0.0 ),
    m12Sq_( 0.0 ),
    m23Sq_( 0.0 ),
    m13Sq_( 0.0 ),
    c12_( 0.0 ),
    c23_( 0.0 ),
    c13_( 0.0 ),
    mPrime_( 0.0 ),
    thetaPrime_( 0.0 ),
    qi_( 0.0 ),
    qk_( 0.0 ),
    p1_12_( 0.0 ),
    p3_12_( 0.0 ),
    p2_23_( 0.0 ),
    p1_23_( 0.0 ),
    p1_13_( 0.0 ),
    p2_13_( 0.0 ),
    p1_Parent_( 0.0 ),
    p2_Parent_( 0.0 ),
    p3_Parent_( 0.0 ),
    squareDP_( calcSquareDPCoords ),
    warnings_( kTRUE )
{
    // Constructor
    mass_.clear();
    mMin_.clear();
    mMax_.clear();
    mSqMin_.clear();
    mSqMax_.clear();
    mSq_.clear();
    mSqDiff_.clear();
    mDiff_.clear();

    mass_.push_back( m1_ );
    mass_.push_back( m2_ );
    mass_.push_back( m3_ );

    mSq_.push_back( m1Sq_ );
    mSq_.push_back( m2Sq_ );
    mSq_.push_back( m3Sq_ );

    // DP max and min kinematic boundary for circumscribed box
    // (see figure in PDG book).
    for ( Int_t i = 0; i < 3; i++ ) {
        mMin_.push_back( mDTot_ - mass_[i] );
        mMax_.push_back( mParent_ - mass_[i] );
        mSqMin_.push_back( mMin_[i] * mMin_[i] );
        mSqMax_.push_back( mMax_[i] * mMax_[i] );
        mSqDiff_.push_back( mSqMax_[i] - mSqMin_[i] );
        mDiff_.push_back( mMax_[i] - mMin_[i] );
    }

    if ( this->squareDP() ) {
        std::cout << "INFO in LauKinematics::LauKinematics : squareDP = kTRUE" << std::endl;
    } else {
        std::cout << "INFO in LauKinematics::LauKinematics : squareDP = kFALSE" << std::endl;
    }

    // add covariant factor calculation
}

LauKinematics::~LauKinematics()
{
    // Destructor
}

void LauKinematics::updateKinematics( const Double_t m13Sq, const Double_t m23Sq )
{
    // Sets the internal private data members
    // m13Sq and m23Sq, as well as m13 and m23.
    // Also calculate m12Sq and m12, given mParent defined in the constructor.
    // Lastly, calculate the helicity cosines.

    // Update the 3 mass-squares
    this->updateMassSquares( m13Sq, m23Sq );

    // Now update the helicity cosines
    this->calcHelicities();

    // Also calculate the m', theta' variables
    if ( squareDP_ ) {
        this->calcSqDPVars();
    }
}

void LauKinematics::updateSqDPKinematics( const Double_t mPrime, const Double_t thetaPrime )
{
    // From square DP co-ordinates, calculate remaining kinematic variables
    this->updateSqDPMassSquares( mPrime, thetaPrime );

    // Finally calculate the helicities and track momenta
    this->calcHelicities();
}

void LauKinematics::calcSqDPVars()
{
    // For given m_12 and cos(theta_12) values, calculate m' and theta' for the square Dalitz plot
    Double_t value = ( 2.0 * ( m12_ - mMin_[2] ) / mDiff_[2] ) - 1.0;
    mPrime_        = LauConstants::invPi * TMath::ACos( value );
    thetaPrime_    = LauConstants::invPi * TMath::ACos( c12_ );
    // Sometimes events are assigned exactly thetaPrime = 0.0 or 1.0
    // which gives problems with efficiency and other histograms
    if ( thetaPrime_ == 0.0 ) {
        thetaPrime_ += 1.0e-10;
    } else if ( thetaPrime_ == 1.0 ) {
        thetaPrime_ -= 1.0e-10;
    }
}

Double_t LauKinematics::calcSqDPJacobian() const
{
    return this->calcSqDPJacobian( mPrime_, thetaPrime_ );
}

Double_t LauKinematics::calcSqDPJacobian( const Double_t mPrime, const Double_t thPrime ) const
{
    // Calculate the Jacobian for the transformation
    // m23^2, m13^2 -> m', theta' (square DP)
    const Double_t m12 = 0.5 * mDiff_[2] * ( 1.0 + TMath::Cos( LauConstants::pi * mPrime ) ) +
                         mMin_[2];
    const Double_t m12Sq = m12 * m12;

    const Double_t e1Cms12 = ( m12Sq - m2Sq_ + m1Sq_ ) / ( 2.0 * m12 );
    const Double_t e3Cms12 = ( mParentSq_ - m12Sq - m3Sq_ ) / ( 2.0 * m12 );

    const Double_t p1Cms12 = this->pCalc( e1Cms12, m1Sq_ );
    const Double_t p3Cms12 = this->pCalc( e3Cms12, m3Sq_ );

    const Double_t deriv1 = LauConstants::piBy2 * mDiff_[2] * TMath::Sin( LauConstants::pi * mPrime );
    const Double_t deriv2 = LauConstants::pi * TMath::Sin( LauConstants::pi * thPrime );

    const Double_t jacobian = 4.0 * p1Cms12 * p3Cms12 * m12 * deriv1 * deriv2;

    return jacobian;
}

void LauKinematics::updateMassSquares( const Double_t m13Sq, const Double_t m23Sq )
{
    m13Sq_ = m13Sq;
    if ( m13Sq_ > 0.0 ) {
        m13_ = TMath::Sqrt( m13Sq_ );
    } else {
        m13_ = 0.0;
    }

    m23Sq_ = m23Sq;
    if ( m23Sq_ > 0.0 ) {
        m23_ = TMath::Sqrt( m23Sq_ );
    } else {
        m23_ = 0.0;
    }

    // Now calculate m12Sq and m12.
    this->calcm12Sq();

    // Calculate momenta of tracks in parent (B, D etc.) rest-frame
    this->calcParentFrameMomenta();
}

void LauKinematics::updateSqDPMassSquares( const Double_t mPrime, const Double_t thetaPrime )
{
    // From square DP co-ordinates, calculate only the mass-squares

    // First set the square-DP variables
    mPrime_     = mPrime;
    thetaPrime_ = thetaPrime;

    // Next calculate m12 and c12
    Double_t m12 = 0.5 * mDiff_[2] * ( 1.0 + TMath::Cos( LauConstants::pi * mPrime ) ) + mMin_[2];
    Double_t c12 = TMath::Cos( LauConstants::pi * thetaPrime );

    // From these calculate m13Sq and m23Sq
    this->updateMassSq_m12( m12, c12 );

    // Calculate momenta of tracks in parent (B, D etc.) rest-frame
    this->calcParentFrameMomenta();
}

void LauKinematics::calcm12Sq()
{
    // Calculate m12Sq from m13Sq and m23Sq.
    m12Sq_ = mParentSq_ + mSqDTot_ - m13Sq_ - m23Sq_;

    // If m12Sq is too low, return lower limit,
    // and modify m13Sq accordingly.
    if ( m12Sq_ < mSqMin_[2] ) {
        m12Sq_ = mSqMin_[2] + 1.0e-3;
        m13Sq_ = mParentSq_ + mSqDTot_ - m12Sq_ - m23Sq_;
    }

    if ( m12Sq_ > 0.0 ) {
        m12_ = TMath::Sqrt( m12Sq_ );
    } else {
        m12_ = 0.0;
    }
}

void LauKinematics::calcParentFrameMomenta()
{
    Double_t e1 = ( mParentSq_ + m1Sq_ - m23Sq_ ) / ( 2.0 * mParent_ );
    Double_t e2 = ( mParentSq_ + m2Sq_ - m13Sq_ ) / ( 2.0 * mParent_ );
    Double_t e3 = ( mParentSq_ + m3Sq_ - m12Sq_ ) / ( 2.0 * mParent_ );

    p1_Parent_ = this->pCalc( e1, m1Sq_ );
    p2_Parent_ = this->pCalc( e2, m2Sq_ );
    p3_Parent_ = this->pCalc( e3, m3Sq_ );
}

void LauKinematics::calcHelicities()
{
    // Calculate helicity angle cosines, given m12Sq, m13Sq and m23Sq.
    // cij_ is the cosine of the helicity angle in the rest frame of the
    // system of particles i and j.
    // It is (but note the caveat below) the angle between tracks i and k
    // in the ij rest frame with indices permuted cyclically.
    // However, it is important to note that it is not exactly a cyclic
    // permutation (there is a special treatment for c23 inside the cFromM
    // function) for reasons of preserving the symmetry about m13=m23 for
    // symmetric final states.
    // The precise definitions are:
    // theta12 is defined as the angle between 1&3 in the rest frame of 1&2
    // theta23 is defined as the angle between 3&1 in the rest frame of 2&3
    // theta13 is defined as the angle between 3&2 in the rest frame of 1&3
    //
    // It is a prerequisite that all mij_ and mijSq_ variables be correctly set.

    Int_t zero( 0 ), one( 1 ), two( 2 );

    c12_   = cFromM( m12Sq_, m13Sq_, m12_, zero, one, two );
    p1_12_ = qi_;
    p3_12_ = qk_;    // i, j = 12 (rest frame), k = 3
    c23_   = cFromM( m23Sq_, m12Sq_, m23_, one, two, zero );
    p2_23_ = qi_;
    p1_23_ = qk_;    // i, j = 23 (rest frame), k = 1
    c13_   = cFromM( m13Sq_, m23Sq_, m13_, two, zero, one );
    p1_13_ = qi_;
    p2_13_ = qk_;    // i, j = 31 (rest frame), k = 2
}

Double_t LauKinematics::cFromM( const Double_t mijSq,
                                const Double_t mikSq,
                                const Double_t mij,
                                const Int_t i,
                                const Int_t j,
                                const Int_t k ) const
{
    // Routine to calculate the cos(helicity) variables from the masses of the particles.
    // (See comment in LauKinematics::calcHelicities for futher information.)

    Double_t EiCmsij = ( mijSq - mSq_[j] + mSq_[i] ) / ( 2.0 * mij );
    Double_t EkCmsij = ( mParentSq_ - mijSq - mSq_[k] ) / ( 2.0 * mij );

    if ( EiCmsij < mass_[i] ) {
        if ( warnings_ ) {
            std::cerr << "WARNING in LauKinematics::cFromM : EiCmsij = " << EiCmsij
                      << " too small < mass_[" << i << "] = " << mass_[i]
                      << " in cFromM, i, j, k = " << i << ", " << j << ", " << k << std::endl;
            std::cerr << "                                 : mijSq = " << mijSq << "; mij = " << mij
                      << "; mSq_[" << j << "] = " << mSq_[j] << "; mSq_[" << i << "] = " << mSq_[i]
                      << std::endl;
        }
        return 0.0;
    }
    if ( EkCmsij < mass_[k] ) {
        if ( warnings_ ) {
            std::cerr << "WARNING in LauKinematics::cFromM : EkCmsij = " << EkCmsij
                      << " too small < mass_[" << k << "] = " << mass_[k]
                      << " in cFromM, i, j, k = " << i << ", " << j << ", " << k << std::endl;
            std::cerr << "                                 : mijSq = " << mijSq << "; mij = " << mij
                      << "; mSq_[" << j << "] = " << mSq_[j] << "; mSq_[" << i << "] = " << mSq_[i]
                      << std::endl;
        }
        return 0.0;
    }

    // Find track i and k momenta in ij rest frame
    qi_ = this->pCalc( EiCmsij, mSq_[i] );
    qk_ = this->pCalc( EkCmsij, mSq_[k] );

    // Find ij helicity angle
    Double_t cosHel = -( mikSq - mSq_[i] - mSq_[k] - 2.0 * EiCmsij * EkCmsij ) / ( 2.0 * qi_ * qk_ );

    if ( cosHel > 1.0 ) {
        cosHel = 1.0;
    } else if ( cosHel < -1.0 ) {
        cosHel = -1.0;
    }

    if ( i == 1 ) {
        cosHel *= -1.0;
    }

    return cosHel;
}

Double_t LauKinematics::mFromC( const Double_t mijSq,
                                const Double_t cij,
                                const Double_t mij,
                                const Int_t i,
                                const Int_t j,
                                const Int_t k ) const
{
    // Returns the mass-squared for a pair of particles, i,j, given their
    // invariant mass (squared) and the helicity angle.
    // cij is helicity angle for the pair which is made from tracks i and j.
    // It is the angle between tracks i and k in the ij rest frame.

    Double_t cosHel( cij );
    if ( i == 1 ) {
        cosHel *= -1.0;
    }

    Double_t EiCmsij = ( mijSq - mSq_[j] + mSq_[i] ) / ( 2.0 * mij );
    Double_t EkCmsij = ( mParentSq_ - mijSq - mSq_[k] ) / ( 2.0 * mij );

    if ( TMath::Abs( EiCmsij - mass_[i] ) > 1e-6 && EiCmsij < mass_[i] ) {
        if ( warnings_ ) {
            std::cerr << "WARNING in LauKinematics::mFromC : EiCmsij = " << EiCmsij << " < "
                      << mass_[i] << " in mFromC, i, j, k = " << i << ", " << j << ", " << k
                      << std::endl;
        }
        return 0.0;
    }
    if ( TMath::Abs( EkCmsij - mass_[k] ) > 1e-6 && EkCmsij < mass_[k] ) {
        if ( warnings_ ) {
            std::cerr << "WARNING in LauKinematics::mFromC : EkCmsij = " << EkCmsij << " < "
                      << mass_[k] << " in mFromC, i, j, k = " << i << ", " << j << ", " << k
                      << std::endl;
        }
        return 0.0;
    }

    // Find track i and k momenta in ij rest fram
    Double_t qi = this->pCalc( EiCmsij, mSq_[i] );
    Double_t qk = this->pCalc( EkCmsij, mSq_[k] );

    // Find mikSq
    Double_t massSq = mSq_[i] + mSq_[k] + 2.0 * EiCmsij * EkCmsij - 2.0 * qi * qk * cosHel;

    if ( massSq < mSqMin_[j] ) {
        if ( warnings_ ) {
            std::cerr << "WARNING in LauKinematics::mFromC : mFromC below bound: i, j, k = " << i
                      << ", " << j << ", " << k << std::endl;
        }
        massSq = mSqMin_[j];
    }

    return massSq;
}

void LauKinematics::genFlatPhaseSpace( Double_t& m13Sq, Double_t& m23Sq ) const
{
    // Routine to generate flat phase-space events.
    // DP max kinematic boundaries in circumscribed box
    // See DP figure in PDG book.
    // m13max=mbrec-mass(2)
    // m13sqmax=m13max*m13max
    // m23max=mbrec-mass(1)
    // m23sqmax=m23max*m23max

    // Generate m13Sq and m23Sq flat within DP overall rectangular box
    // Loop if the generated point is not within the DP kinematic boundary

    do {
        m13Sq = mSqMin_[1] + LauRandom::randomFun()->Rndm() * mSqDiff_[1];
        m23Sq = mSqMin_[0] + LauRandom::randomFun()->Rndm() * mSqDiff_[0];

    } while ( ! this->withinDPLimits( m13Sq, m23Sq ) );
}

void LauKinematics::genFlatSqDP( Double_t& mPrime, Double_t& thetaPrime ) const
{
    // Generate random event in the square Dalitz plot
    mPrime     = LauRandom::randomFun()->Rndm();
    thetaPrime = LauRandom::randomFun()->Rndm();
}

Bool_t LauKinematics::withinDPLimits( const Double_t m13Sq, const Double_t m23Sq ) const
{
    // Find out whether the point (m13Sq,m23Sq) is within the limits of the
    // Dalitz plot. The limits are specified by the invariant masses
    // of the parent (e.g. B) and its three daughters that were
    // defined in the constructor of this class. Here
    // m_13Sq = square of invariant mass of daughters 1 and 3
    // m_23Sq = square of invariant mass of daughters 2 and 3.

    Bool_t withinDP = kFALSE;

    // First check that m13Sq is within its absolute min and max
    if ( ! ( ( m13Sq > mSqMin_[1] ) && ( m13Sq < mSqMax_[1] ) ) ) {
        return kFALSE;
    }

    // Now for the given value of m13Sq calculate the local min and max of m23Sq
    Double_t m13 = TMath::Sqrt( m13Sq );

    Double_t e3Cms13 = ( m13Sq - m1Sq_ + m3Sq_ ) / ( 2.0 * m13 );
    Double_t p3Cms13 = this->pCalc( e3Cms13, m3Sq_ );

    Double_t e2Cms13 = ( mParentSq_ - m13Sq - m2Sq_ ) / ( 2.0 * m13 );
    Double_t p2Cms13 = this->pCalc( e2Cms13, m2Sq_ );

    Double_t term = 2.0 * e2Cms13 * e3Cms13 + m2Sq_ + m3Sq_;

    Double_t m23SqLocMin = term - 2.0 * p2Cms13 * p3Cms13;
    Double_t m23SqLocMax = term + 2.0 * p2Cms13 * p3Cms13;

    // Check whether the given value of m23Sq satisfies these bounds
    if ( m23Sq > m23SqLocMin && m23Sq < m23SqLocMax ) {
        withinDP = kTRUE;
    }

    return withinDP;
}

Bool_t LauKinematics::withinDPLimits2( const Double_t m13Sq, const Double_t m23Sq ) const
{
    // Same as withinDPLimits, but this time testing whether the m13Sq
    // variable is within the kinematic boundary for the given m23Sq value

    Bool_t withinDP = kFALSE;

    // First check that m13Sq is within its absolute min and max
    if ( ! ( ( m23Sq > mSqMin_[0] ) && ( m23Sq < mSqMax_[0] ) ) ) {
        return kFALSE;
    }

    // Now for the given value of m13Sq calculate the local min and max of m23Sq
    Double_t m23 = TMath::Sqrt( m23Sq );

    Double_t e3Cms23 = ( m23Sq - m2Sq_ + m3Sq_ ) / ( 2.0 * m23 );
    Double_t p3Cms23 = this->pCalc( e3Cms23, m3Sq_ );

    Double_t e1Cms23 = ( mParentSq_ - m23Sq - m1Sq_ ) / ( 2.0 * m23 );
    Double_t p1Cms23 = this->pCalc( e1Cms23, m1Sq_ );

    Double_t term = 2.0 * e3Cms23 * e1Cms23 + m1Sq_ + m3Sq_;

    Double_t m13SqLocMin = term - 2.0 * p3Cms23 * p1Cms23;
    Double_t m13SqLocMax = term + 2.0 * p3Cms23 * p1Cms23;

    // Check whether the given value of m23Sq satisfies these bounds
    if ( m13Sq > m13SqLocMin && m13Sq < m13SqLocMax ) {
        withinDP = kTRUE;
    }

    return withinDP;
}

Bool_t LauKinematics::withinSqDPLimits( const Double_t mPrime, const Double_t thetaPrime ) const
{
    // Check whether m' and theta' are between 0 and 1
    Bool_t withinDP( kFALSE );
    if ( mPrime > 0.0 && mPrime < 1.0 && thetaPrime > 0.0 && thetaPrime < 1.0 ) {
        withinDP = kTRUE;
    }

    return withinDP;
}

Double_t LauKinematics::calcThirdMassSq( const Double_t firstMassSq, const Double_t secondMassSq ) const
{
    // Calculate one massSq from the other two
    return mParentSq_ + mSqDTot_ - firstMassSq - secondMassSq;
}

Double_t LauKinematics::distanceFromDPCentre() const
{
    return this->distanceFromDPCentre( m13Sq_, m23Sq_ );
}

Double_t LauKinematics::distanceFromDPCentre( const Double_t m13Sq, const Double_t m23Sq ) const
{
    // DP centre is defined as the point where m12=m13=m23

    // First find the m^2_ij value for the centre
    Double_t centreMijSq = ( mParentSq_ + m1Sq_ + m2Sq_ + m3Sq_ ) / 3.0;

    // Then find the difference between this and the two provided co-ords
    Double_t diff13 = m13Sq - centreMijSq;
    Double_t diff23 = m23Sq - centreMijSq;

    // Calculate the total distance
    Double_t distance = TMath::Sqrt( diff13 * diff13 + diff23 * diff23 );
    return distance;
}

Double_t LauKinematics::pCalc( const Double_t energy, const Double_t massSq ) const
{
    // Routine to calculate the momentum of a particle, given its energy and
    // invariant mass (squared).
    Double_t arg = energy * energy - massSq;

    if ( arg < 0.0 ) {
        //if (warnings_) {
        //std::cerr<<"WARNING in LauKinematics::pCalc : arg < 0.0: "<<arg<<" for e = "<<energy<<", mSq = "<<massSq<<std::endl;
        //}
        arg = 0.0;
    }
    Double_t pCalcVal = TMath::Sqrt( arg );
    return pCalcVal;
}

void LauKinematics::flipAndUpdateKinematics()
{
    // Flips the DP variables m13^2 <-> m23^2.
    // Used in the case of symmetrical Dalitz plots (i.e. when two of the
    // daughter tracks are the same type) within the
    // LauIsobarDynamics::resAmp function.
    this->updateKinematics( m23Sq_, m13Sq_ );
}

void LauKinematics::rotateAndUpdateKinematics()
{
    // Cyclically rotates the DP variables (m12 -> m23, m23 -> m13, m13 -> m12)
    // Used in the case of fully symmetric Dalitz plots (i.e. when all
    // three of the daughter tracks are the same type) within the
    // LauIsobarDynamics::resAmp function.
    this->updateKinematics( m12Sq_, m13Sq_ );
}

void LauKinematics::updateMassSq_m23( const Double_t m23, const Double_t c23 )
{
    // Update the variables m12Sq_ and m13Sq_ given m23 and c23.
    m23_   = m23;
    m23Sq_ = m23 * m23;
    c23_   = c23;

    const Int_t zero( 0 ), one( 1 ), two( 2 );
    m12Sq_ = this->mFromC( m23Sq_, c23_, m23_, one, two, zero );
    m13Sq_ = mParentSq_ - m23Sq_ - m12Sq_ + mSqDTot_;

    m12_ = TMath::Sqrt( m12Sq_ );
    m13_ = TMath::Sqrt( m13Sq_ );
}

void LauKinematics::updateMassSq_m13( const Double_t m13, const Double_t c13 )
{
    // Update the variables m12Sq_ and m23Sq_ given m13 and c13.
    m13_   = m13;
    m13Sq_ = m13 * m13;
    c13_   = c13;

    const Int_t zero( 0 ), one( 1 ), two( 2 );
    m23Sq_ = this->mFromC( m13Sq_, c13_, m13_, two, zero, one );
    m12Sq_ = mParentSq_ - m23Sq_ - m13Sq_ + mSqDTot_;

    m23_ = TMath::Sqrt( m23Sq_ );
    m12_ = TMath::Sqrt( m12Sq_ );
}

void LauKinematics::updateMassSq_m12( const Double_t m12, const Double_t c12 )
{
    // Update the variables m23Sq_ and m13Sq_ given m12 and c12.
    m12_   = m12;
    m12Sq_ = m12 * m12;
    c12_   = c12;

    const Int_t zero( 0 ), one( 1 ), two( 2 );
    m13Sq_ = this->mFromC( m12Sq_, c12_, m12_, zero, one, two );
    m23Sq_ = mParentSq_ - m12Sq_ - m13Sq_ + mSqDTot_;
    m13_   = TMath::Sqrt( m13Sq_ );
    m23_   = TMath::Sqrt( m23Sq_ );
}

void LauKinematics::updateKinematicsFrom23( const Double_t m23, const Double_t c23 )
{
    // Calculate the other mass squares
    this->updateMassSq_m23( m23, c23 );

    // Calculate the remaining helicity angles
    this->calcHelicities();

    // Calculate momenta of tracks in parent (B, D etc.) rest-frame
    this->calcParentFrameMomenta();

    // Also calculate the m', theta' variables
    if ( squareDP_ ) {
        this->calcSqDPVars();
    }
}

void LauKinematics::updateKinematicsFrom13( const Double_t m13, const Double_t c13 )
{
    // Calculate the other mass squares
    this->updateMassSq_m13( m13, c13 );

    // Calculate the remaining helicity angles
    this->calcHelicities();

    // Calculate momenta of tracks in parent (B, D etc.) rest-frame
    this->calcParentFrameMomenta();

    // Also calculate the m', theta' variables
    if ( squareDP_ ) {
        this->calcSqDPVars();
    }
}

void LauKinematics::updateKinematicsFrom12( const Double_t m12, const Double_t c12 )
{
    // Calculate the other mass squares
    this->updateMassSq_m12( m12, c12 );

    // Calculate the remaining helicity angles
    this->calcHelicities();

    // Calculate momenta of tracks in parent (B, D etc.) rest-frame
    this->calcParentFrameMomenta();

    // Also calculate the m', theta' variables
    if ( squareDP_ ) {
        this->calcSqDPVars();
    }
}

Double_t LauKinematics::genm13Sq() const
{
    Double_t m13Sq = mSqMin_[1] + LauRandom::randomFun()->Rndm() * mSqDiff_[1];
    return m13Sq;
}

Double_t LauKinematics::genm23Sq() const
{
    Double_t m23Sq = mSqMin_[0] + LauRandom::randomFun()->Rndm() * mSqDiff_[0];
    return m23Sq;
}

Double_t LauKinematics::genm12Sq() const
{
    Double_t m12Sq = mSqMin_[2] + LauRandom::randomFun()->Rndm() * mSqDiff_[2];
    return m12Sq;
}

Double_t dal( Double_t* x, Double_t* par )
{
    Double_t mParent = par[0];
    Double_t mi      = par[1];
    Double_t mj      = par[2];
    Double_t mk      = par[3];

    Double_t mikSq = x[0];
    Double_t mjkSq = x[1];
    Double_t mik   = TMath::Sqrt( mikSq );
    Double_t mjk   = TMath::Sqrt( mjkSq );

    Double_t ejcmsik = ( mParent * mParent - mj * mj - mik * mik ) / ( 2.0 * mik );
    Double_t ekcmsik = ( mik * mik + mk * mk - mi * mi ) / ( 2.0 * mik );
    if ( ekcmsik < mk || ejcmsik < mj )
        return 2.0;

    Double_t pj       = TMath::Sqrt( ejcmsik * ejcmsik - mj * mj );
    Double_t pk       = TMath::Sqrt( ekcmsik * ekcmsik - mk * mk );
    Double_t coshelik = ( mjk * mjk - mk * mk - mj * mj - 2.0 * ejcmsik * ekcmsik ) /
                        ( 2.0 * pj * pk );

    Double_t coshelikSq = coshelik * coshelik;
    return coshelikSq;
}

void LauKinematics::drawDPContour( Int_t orientation, Int_t nbins )
{
    // orientation -
    // 1323 : x = m13, y = m23
    // etc.

    Double_t m1      = this->getm1();
    Double_t m2      = this->getm2();
    Double_t m3      = this->getm3();
    Double_t mParent = this->getmParent();

    Double_t m13SqMin = this->getm13SqMin();
    Double_t m23SqMin = this->getm23SqMin();
    Double_t m12SqMin = this->getm12SqMin();
    Double_t m13SqMax = this->getm13SqMax();
    Double_t m23SqMax = this->getm23SqMax();
    Double_t m12SqMax = this->getm12SqMax();

    Double_t xMin( 0.0 );
    Double_t xMax( mParent * mParent );
    Double_t yMin( 0.0 );
    Double_t yMax( mParent * mParent );
    if ( orientation == 1323 ) {
        xMin = m13SqMin - 1.0;
        xMax = m13SqMax + 1.0;
        yMin = m23SqMin - 1.0;
        yMax = m23SqMax + 1.0;
    } else if ( orientation == 2313 ) {
        xMin = m23SqMin - 1.0;
        xMax = m23SqMax + 1.0;
        yMin = m13SqMin - 1.0;
        yMax = m13SqMax + 1.0;
    } else if ( orientation == 1213 ) {
        xMin = m12SqMin - 1.0;
        xMax = m12SqMax + 1.0;
        yMin = m13SqMin - 1.0;
        yMax = m13SqMax + 1.0;
    } else if ( orientation == 1312 ) {
        xMin = m13SqMin - 1.0;
        xMax = m13SqMax + 1.0;
        yMin = m12SqMin - 1.0;
        yMax = m12SqMax + 1.0;
    } else if ( orientation == 1223 ) {
        xMin = m12SqMin - 1.0;
        xMax = m12SqMax + 1.0;
        yMin = m23SqMin - 1.0;
        yMax = m23SqMax + 1.0;
    } else if ( orientation == 2312 ) {
        xMin = m23SqMin - 1.0;
        xMax = m23SqMax + 1.0;
        yMin = m12SqMin - 1.0;
        yMax = m12SqMax + 1.0;
    } else {
        std::cerr << "ERROR in LauKinematics::drawDPContour : Unrecognised orientation, not drawing contour."
                  << std::endl;
        return;
    }

    Int_t npar = 4;
    TF2* f2    = new TF2( "contour", &dal, xMin, xMax, yMin, yMax, npar );

    // Set the parameters
    f2->SetParameter( 0, mParent );
    if ( orientation == 1323 ) {
        f2->SetParameter( 1, m1 );
        f2->SetParameter( 2, m2 );
        f2->SetParameter( 3, m3 );
    } else if ( orientation == 2313 ) {
        f2->SetParameter( 1, m2 );
        f2->SetParameter( 2, m1 );
        f2->SetParameter( 3, m3 );
    } else if ( orientation == 1213 ) {
        f2->SetParameter( 1, m2 );
        f2->SetParameter( 2, m3 );
        f2->SetParameter( 3, m1 );
    } else if ( orientation == 1312 ) {
        f2->SetParameter( 1, m3 );
        f2->SetParameter( 2, m2 );
        f2->SetParameter( 3, m1 );
    } else if ( orientation == 1223 ) {
        f2->SetParameter( 1, m1 );
        f2->SetParameter( 2, m3 );
        f2->SetParameter( 3, m2 );
    } else if ( orientation == 2312 ) {
        f2->SetParameter( 1, m3 );
        f2->SetParameter( 2, m1 );
        f2->SetParameter( 3, m2 );
    }

    // Set up the contour to be drawn when the value of the function == 1.0
    Double_t b[] = { 1.0 };
    f2->SetContour( 1, b );

    // Set the number of bins for the contour to be sampled over
    f2->SetNpx( nbins );
    f2->SetNpy( nbins );
    // and the line style
    f2->SetLineWidth( 3 );
    f2->SetLineStyle( kSolid );

    // Draw the contour on top of the histo that should already have been drawn
    f2->DrawCopy( "same" );

    delete f2;
}
