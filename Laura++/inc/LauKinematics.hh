
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

/*! \file LauKinematics.hh
    \brief File containing declaration of LauKinematics class.
*/

/*! \class LauKinematics
    \brief Class for calculating 3-body kinematic quantities.

    Class for defining the many routines related to the three body kinematics.
    Given the two DP co-ordinates, all other useful quantities e.g. the helicity angles are calculated.
    Optionally, the so-called ``square Dalitz plot'' quantities can also be calculated.
*/

#ifndef LAU_KINEMATICS
#define LAU_KINEMATICS

#include "Rtypes.h"
#include "TMath.h"

#include <vector>

class LauKinematics {

  public:
    //! Constructor
    /*!
        \param [in] m1 the first daughter mass
        \param [in] m2 the second daughter mass
        \param [in] m3 the third daughter mass
        \param [in] mParent the parent particle mass
        \param [in] calcSquareDPCoords boolean flag to enable/disable calculation of the square Dalitz plot co-ordinates
        \param [in] symmetricalDP boolean flag to indicate whether the DP is symmetric (i.e. two identical particle in final state)
        \param [in] fullySymmetricDP boolean flag to indicate whether the DP is fully symmetric (i.e. all three final-state particles are identical)
    */
    LauKinematics( const Double_t m1,
                   const Double_t m2,
                   const Double_t m3,
                   const Double_t mParent,
                   const Bool_t calcSquareDPCoords = kFALSE,
                   const Bool_t symmetricalDP      = kFALSE,
                   const Bool_t fullySymmetricDP   = kFALSE );

    //! Destructor
    virtual ~LauKinematics();

    //! Enable/disable the calculation of square Dalitz plot co-ordinates
    /*!
        \param [in] calcSquareDPCoords kTRUE/kFALSE to enable/disable calculation of the square DP co-ordinates
    */
    inline void squareDP( const Bool_t calcSquareDPCoords ) { squareDP_ = calcSquareDPCoords; }

    //! Are the square Dalitz plot co-ordinates being calculated?
    /*!
        \return kTRUE if the square Dalitz plot co-ordinates are being calculated, kFALSE otherwise
    */
    inline Bool_t squareDP() const { return squareDP_; }

    //! Is the DP symmetric?
    /*!
        \return kTRUE if the DP is symmetric (i.e. daughters 1 and 2 are identical), kFALSE otherwise
    */
    inline Bool_t gotSymmetricalDP() const { return ( symmetricalDP_ && ! fullySymmetricDP_ ); }

    //! Is the DP fully symmetric?
    /*!
        \return kTRUE if the DP is fully symmetric (i.e. daughters 1, 2 and 3 are identical), kFALSE otherwise
    */
    inline Bool_t gotFullySymmetricDP() const { return fullySymmetricDP_; }

    //! Enable/disable warning messages
    inline void warningMessages( const Bool_t boolean ) { warnings_ = boolean; }

    //! Update all kinematic quantities based on the DP co-ordinates m13Sq and m23Sq
    /*!
        It can be useful to first check that the point is within the kinematic boundary (using LauKinematics::withinDPLimits) before calling this method.

        \param [in] m13Sq the invariant mass squared of daughters 1 and 3
        \param [in] m23Sq the invariant mass squared of daughters 2 and 3
    */
    void updateKinematics( const Double_t m13Sq, const Double_t m23Sq );

    //! Update all kinematic quantities based on the square DP co-ordinates m' and theta'
    /*!
        It can be useful to first check that the point is within the kinematic boundary (using LauKinematics::withinSqDPLimits) before calling this method.

        \param [in] mPrime the m' co-ordinate
        \param [in] thetaPrime the theta' co-ordinate
    */
    void updateSqDPKinematics( const Double_t mPrime, const Double_t thetaPrime );

    //! Update all kinematic quantities based on m23 and cos(theta23)
    /*!
        \param [in] m23 the invariant mass of daughters 2 and 3
        \param [in] c23 the cosine of the helicity angle theta23, \see getc23
    */
    void updateKinematicsFrom23( const Double_t m23, const Double_t c23 );

    //! Update all kinematic quantities based on m13 and cos(theta13)
    /*!
        \param [in] m13 the invariant mass of daughters 1 and 3
        \param [in] c13 the cosine of the helicity angle theta13, \see getc13
    */
    void updateKinematicsFrom13( const Double_t m13, const Double_t c13 );

    //! Update all kinematic quantities based on m12 and cos(theta12)
    /*!
        \param [in] m12 the invariant mass of daughters 1 and 2
        \param [in] c12 the cosine of the helicity angle theta12, \see getc12
    */
    void updateKinematicsFrom12( const Double_t m12, const Double_t c12 );

    //! Calculate the Jacobian for the transformation m23^2, m13^2 -> m', theta' (square DP) at the given point in the square DP
    /*!
        \param [in] mPrime the m' co-ordinate
        \param [in] thPrime the theta' co-ordinate
        \return the jacobian of the transformation
    */
    Double_t calcSqDPJacobian( const Double_t mPrime, const Double_t thPrime ) const;

    //! Calculate the Jacobian for the transformation m23^2, m13^2 -> m', theta' (square DP) at the currently stored point in the square DP
    /*!
        \return the jacobian of the transformation
    */
    Double_t calcSqDPJacobian() const;

    //! Routine to generate events flat in phase-space
    /*!
        \param [out] m13Sq the invariant mass squared of daughters 1 and 3
        \param [out] m23Sq the invariant mass squared of daughters 2 and 3
    */
    void genFlatPhaseSpace( Double_t& m13Sq, Double_t& m23Sq ) const;

    //! Routine to generate events flat in the square Dalitz plot
    /*!
        \param [out] mPrime the m' variable
        \param [out] thetaPrime the theta' variable
    */
    void genFlatSqDP( Double_t& mPrime, Double_t& thetaPrime ) const;

    //! Check whether a given (m13Sq,m23Sq) point is within the kinematic limits of the Dalitz plot
    /*!
        This method first checks that m13Sq is within its absolute
        min and max and then for the given m13Sq calculates the
        local min and max of m23Sq and checks whether the given
        m23Sq satisfies these bounds.

        \param [in] m13Sq the invariant mass squared of daughters 1 and 3
        \param [in] m23Sq the invariant mass squared of daughters 2 and 3
        \return kTRUE if the event is inside the kinematic limit, kFALSE otherwise
    */
    Bool_t withinDPLimits( const Double_t m13Sq, const Double_t m23Sq ) const;

    //! Check whether a given (m13Sq,m23Sq) point is within the kinematic limits of the Dalitz plot (alternate method)
    /*!
        This method first checks that m23Sq is within its absolute
        min and max and then for the given m23Sq calculates the
        local min and max of m13Sq and checks whether the given
        m13Sq satisfies these bounds.

        \param [in] m13Sq the m13 invariant mass pair squared
        \param [in] m23Sq the m23 invariant mass pair squared
        \return kTRUE if the event is inside the kinematic limit, kFALSE otherwise
    */
    Bool_t withinDPLimits2( const Double_t m13Sq, const Double_t m23Sq ) const;

    //! Check whether a given (m',theta') point is within the kinematic limits of the Dalitz plot
    /*!
        \param [in] mPrime the m' co-ordinate
        \param [in] thetaPrime the theta' co-ordinate
        \return kTRUE if the event is inside the kinematic limit, kFALSE otherwise
    */
    Bool_t withinSqDPLimits( const Double_t mPrime, const Double_t thetaPrime ) const;

    //! Calculate the third invariant mass square from the two provided (e.g. mjkSq from mijSq and mikSq)
    /*!
        \param [in] firstMassSq the first invariant mass squared
        \param [in] secondMassSq the second invariant mass squared
        \return the third invariant mass square
    */
    Double_t calcThirdMassSq( const Double_t firstMassSq, const Double_t secondMassSq ) const;

    //! Calculate the distance from the currently set (m13Sq, m23Sq) point to the centre of the Dalitz plot (which is defined as the point where m12=m13=m23)
    /*!
        \return the distance to the DP centre
    */
    Double_t distanceFromDPCentre() const;

    //! Calculate the distance from a given (m13Sq, m23Sq) point to the centre of the Dalitz plot (which is defined as the point where m12=m13=m23)
    /*!
        \return the distance to the DP centre
    */
    Double_t distanceFromDPCentre( const Double_t m13Sq, const Double_t m23Sq ) const;

    //! Get the m12 invariant mass
    /*!
        \return the m12 invariant mass
    */
    inline Double_t getm12() const { return m12_; }

    //! Get the m23 invariant mass
    /*!
        \return the m23 invariant mass
    */
    inline Double_t getm23() const { return m23_; }

    //! Get the m13 invariant mass
    /*!
        \return the m13 invariant mass
    */
    inline Double_t getm13() const { return m13_; }

    //! Get the m12 invariant mass square
    /*!
        \return the m12 invariant mass square
    */
    inline Double_t getm12Sq() const { return m12Sq_; }

    //! Get the m23 invariant mass square
    /*!
        \return the m23 invariant mass square
    */
    inline Double_t getm23Sq() const { return m23Sq_; }

    //! Get the m13 invariant mass square
    /*!
        \return the m13 invariant mass square
    */
    inline Double_t getm13Sq() const { return m13Sq_; }

    //! Get the cosine of the helicity angle theta12
    /*!
        theta12 is defined as the angle between 1&3 in the rest frame of 1&2

        \return the cosine of the helicity angle theta12
    */
    inline Double_t getc12() const { return c12_; }

    //! Get the cosine of the helicity angle theta23
    /*!
        theta23 is defined as the angle between 3&1 in the rest frame of 2&3

        \return the cosine of the helicity angle theta23
    */
    inline Double_t getc23() const { return c23_; }

    //! Get the cosine of the helicity angle theta13
    /*!
        theta13 is defined as the angle between 3&2 in the rest frame of 1&3

        \return the cosine of the helicity angle theta13
    */
    inline Double_t getc13() const { return c13_; }

    //! Get m' value
    /*!
        \return m' value
    */
    inline Double_t getmPrime() const { return mPrime_; }

    //! Get theta' value
    /*!
        \return theta' value
    */
    inline Double_t getThetaPrime() const { return thetaPrime_; }

    //! Get parent mass
    /*!
        \return parent mass
    */
    inline Double_t getmParent() const { return mParent_; }

    //! Get parent mass squared
    /*!
        \return parent mass squared
    */
    inline Double_t getmParentSq() const { return mParentSq_; }

    //! Get the box area defined from the kinematic bounds
    /*!
        The box area is defined as:\n
        [(M-m3)^2 - (m1+m2)^2]*[(M-m2)^2 - (m1+m3)^2] .:. (m13SqMax - m13SqMin)*(m23SqMax - m23SqMin)

        \return the Dalitz plot box area
    */
    inline Double_t getDPBoxArea() const
    {
        return ( mSqMax_[1] - mSqMin_[1] ) * ( mSqMax_[0] - mSqMin_[0] );
    }

    //! Flips the DP variables m13^2 <-> m23^2 and recalculates all kinematic quantities
    /*!
        Useful in the case of symmetrical Dalitz plots, i.e. when two final state particles are identical
    */
    void flipAndUpdateKinematics();

    //! Cyclically rotates the DP variables (m12 -> m23, m23 -> m13, m13 -> m12) and recalculates all kinematic quantities
    /*!
        Useful in the case of a fully symmetric Dalitz plot, i.e. when all three final state particles are identical
    */
    void rotateAndUpdateKinematics();

    //! Get the m1 mass
    /*!
        \return the m1 mass
    */
    inline Double_t getm1() const { return m1_; }

    //! Get the m2 mass
    /*!
        \return the m2 mass
    */
    inline Double_t getm2() const { return m2_; }

    //! Get the m3 mass
    /*!
        \return the m3 mass
    */
    inline Double_t getm3() const { return m3_; }

    //! Get the m23 minimum defined as (m2 + m3)
    /*!
        \return the minimum value for m23
    */
    inline Double_t getm23Min() const { return TMath::Sqrt( mSqMin_[0] ); }

    //! Get the m13 minimum defined as (m1 + m3)
    /*!
        \return the minimum value for m13
    */
    inline Double_t getm13Min() const { return TMath::Sqrt( mSqMin_[1] ); }

    //! Get the m12 minimum defined as (m1 + m2)
    /*!
        \return the minimum value for m12
    */
    inline Double_t getm12Min() const { return TMath::Sqrt( mSqMin_[2] ); }

    //! Get the m23 maximum defined as (mParent - m1)
    /*!
        \return the maximum value for m23
    */
    inline Double_t getm23Max() const { return TMath::Sqrt( mSqMax_[0] ); }

    //! Get the m13 maximum defined as (mParent - m2)
    /*!
        \return the maximum value for m13
    */
    inline Double_t getm13Max() const { return TMath::Sqrt( mSqMax_[1] ); }

    //! Get the m12 maximum defined as (mParent - m3)
    /*!
        \return the maximum value for m12
    */
    inline Double_t getm12Max() const { return TMath::Sqrt( mSqMax_[2] ); }

    //! Get the m23Sq minimum, (m2 + m3)^2
    /*!
        \return the minimum value for m23Sq
    */
    inline Double_t getm23SqMin() const { return mSqMin_[0]; }

    //! Get the m13Sq minimum, (m1 + m3)^2
    /*!
        \return the minimum value for m13Sq
    */
    inline Double_t getm13SqMin() const { return mSqMin_[1]; }

    //! Get the m12Sq minimum, (m1 + m2)^2
    /*!
        \return the minimum value for m12Sq
    */
    inline Double_t getm12SqMin() const { return mSqMin_[2]; }

    //! Get the m23Sq maximum, (mParent - m1)^2
    /*!
        \return the maximum value for m23Sq
    */
    inline Double_t getm23SqMax() const { return mSqMax_[0]; }

    //! Get the m13Sq maximum, (mParent - m2)^2
    /*!
        \return the maximum value for m13Sq
    */
    inline Double_t getm13SqMax() const { return mSqMax_[1]; }

    //! Get the m12Sq maximum, (mParent - m3)^2
    /*!
        \return the maximum value for m12Sq
    */
    inline Double_t getm12SqMax() const { return mSqMax_[2]; }

    //! Get the momentum of the track 1 in 12 rest frame
    /*!
        \return the momentum of track 1 in 12 rest frame
    */
    inline Double_t getp1_12() const { return p1_12_; }

    //! Get the momentum of the track 3 in 12 rest frame
    /*!
        \return the momentum of track 3 in 12 rest frame
    */
    inline Double_t getp3_12() const { return p3_12_; }

    //! Get the momentum of the track 2 in 23 rest frame
    /*!
        \return the momentum of track 2 in 23 rest frame
    */
    inline Double_t getp2_23() const { return p2_23_; }

    //! Get the momentum of the track 1 in 23 rest frame
    /*!
        \return the momentum of track 1 in 23 rest frame
    */
    inline Double_t getp1_23() const { return p1_23_; }

    //! Get the momentum of the track 1 in 13 rest frame
    /*!
        \return the momentum of track 1 in 13 rest frame
    */
    inline Double_t getp1_13() const { return p1_13_; }

    //! Get the momentum of the track 2 in 13 rest frame
    /*!
        \return the momentum of track 2 in 13 rest frame
    */
    inline Double_t getp2_13() const { return p2_13_; }

    //! Get the momentum of the track 1 in parent rest frame
    /*!
        \return the momentum of track 1 in parent rest frame
    */
    inline Double_t getp1_Parent() const { return p1_Parent_; }

    //! Get the momentum of the track 2 in parent rest frame
    /*!
        \return the momentum of track 2 in parent rest frame
    */
    inline Double_t getp2_Parent() const { return p2_Parent_; }

    //! Get the momentum of the track 3 in parent rest frame
    /*!
        \return the momentum of track 3 in parent rest frame
    */
    inline Double_t getp3_Parent() const { return p3_Parent_; }

    //! Method to draw the Dalitz plot contours on the top of the histo previously drawn
    /*!
        \param [in] orientation orientation used for the draw, with default set to 1323 that corresponds x = m13, y = m23
        \param [in] nbins number of bins for the contour to be sampled over (default = 100)
    */
    void drawDPContour( const Int_t orientation = 1323, const Int_t nbins = 100 );

    //! Get covariant factor in 12 axis
    /*!
        \return covariant factor in 12 axis
    */
    inline Double_t getcov12() const
    {
        return ( mParentSq_ + m12Sq_ - m3Sq_ ) / ( 2. * mParent_ * m12_ );
    }

    //! Get covariant factor in 13 axis
    /*!
        \return covariant factor in 13 axis
    */
    inline Double_t getcov13() const
    {
        return ( mParentSq_ + m13Sq_ - m2Sq_ ) / ( 2. * mParent_ * m13_ );
    }

    //! Get covariant factor in 23 axis
    /*!
        \return covariant factor in 23 axis
    */
    inline Double_t getcov23() const
    {
        return ( mParentSq_ + m23Sq_ - m1Sq_ ) / ( 2. * mParent_ * m23_ );
    }

  protected:
    //! Update the variables m23Sq_ and m13Sq_ given the invariant mass m12 and the cosine of the helicity angle c12
    /*!
        \param [in] m12 the invariant mass m12
        \param [in] c12 the cosine of the helicity angle c12
    */
    void updateMassSq_m12( const Double_t m12, const Double_t c12 );

    //! Update the variables m12Sq_ and m13Sq_ given the invariant mass m23 and the cosine of the helicity angle c23
    /*!
        \param [in] m23 the invariant mass m12
        \param [in] c23 the cosine of the helicity angle c23
    */
    void updateMassSq_m23( const Double_t m23, const Double_t c23 );

    //! Update the variables m12Sq_ and m23Sq_ given the invariant mass m13 and the cosine of the helicity angle c13
    /*!
        \param [in] m13 the invariant mass m13
        \param [in] c13 the cosine of the helicity angle c13
    */
    void updateMassSq_m13( const Double_t m13, const Double_t c13 );

    //! Update some kinematic quantities based on the DP co-ordinates m13Sq and m23Sq
    /*!
        Only the three invariant masses and their squares, plus the parent rest-frame momenta are updated.

        \param [in] m13Sq the invariant mass squared of daughters 1 and 3
        \param [in] m23Sq the invariant mass squared of daughters 2 and 3
    */
    void updateMassSquares( const Double_t m13Sq, const Double_t m23Sq );

    //! Update some kinematic quantities based on the square DP co-ordinates m' and theta'
    /*!
        Only m', theta', the three invariant masses and their squares, plus the parent rest-frame momenta are updated.

        \param [in] mPrime the m' co-ordinate
        \param [in] thetaPrime the theta' co-ordinate
    */
    void updateSqDPMassSquares( const Double_t mPrime, const Double_t thetaPrime );

    //! General method to calculate the cos(helicity) variables from the masses of the particles
    /*!
        \param [in] mijSq the mij invariant mass square
        \param [in] mikSq the mik invariant mass square
        \param [in] mij the mij invariant mass
        \param [in] i the index for the first track
        \param [in] j the index for the second track
        \param [in] k the index for the third track
        \return helicity angle in the ij rest frame
    */
    Double_t cFromM( const Double_t mijSq,
                     const Double_t mikSq,
                     const Double_t mij,
                     const Int_t i,
                     const Int_t j,
                     const Int_t k ) const;

    //! General method to calculate mikSq given mijSq and cosHel_ij
    /*!
        \param [in] mijSq the mij invariant mass square
        \param [in] cij the helicity angle for the pair which is made from tracks i and j
        \param [in] mij the mij invariant mass
        \param [in] i the index for the first track
        \param [in] j the index for the second track
        \param [in] k the index for the third track
        \return the invariant mass square mikSq
    */
    Double_t mFromC( const Double_t mijSq,
                     const Double_t cij,
                     const Double_t mij,
                     const Int_t i,
                     const Int_t j,
                     const Int_t k ) const;

    //! General method to calculate the momentum of a particle, given its energy and invariant mass squared.
    /*!
        \param [in] energy the energy of the particle
        \param [in] massSq the invariant mass squared of the particle
        \return the momentum of the particle
    */
    Double_t pCalc( const Double_t energy, const Double_t massSq ) const;

    //! Randomly generate the invariant mass squared m13Sq
    /*!
        \return the invariant mass squared m13Sq
    */
    Double_t genm13Sq() const;
    //! Randomly generate the invariant mass squared m23Sq
    /*!
        \return the invariant mass squared m23Sq
    */
    Double_t genm23Sq() const;
    //! Randomly generate the invariant mass squared m12Sq
    /*!
        \return the invariant mass squared m12Sq
    */
    Double_t genm12Sq() const;

    //! Calculate m12Sq from m13Sq and m23Sq
    void calcm12Sq();

    //! Calculate cosines of the helicity angles, momenta of daughters and bachelor in various ij rest frames
    void calcHelicities();

    //! Calculate the m' and theta' variables for the square Dalitz plot
    void calcSqDPVars();

    //! Calculate the momenta of each daughter in the parent rest frame
    void calcParentFrameMomenta();

  private:
    //! Copy constructor (not implemented)
    LauKinematics( const LauKinematics& rhs );

    //! Copy assignment operator (not implemented)
    LauKinematics& operator=( const LauKinematics& rhs );

    //! Symmetrical DP
    const Bool_t symmetricalDP_;
    //! Fully-symmetrical DP
    const Bool_t fullySymmetricDP_;

    //! Mass of particle 1
    const Double_t m1_;
    //! Mass of particle 2
    const Double_t m2_;
    //! Mass of particle 3
    const Double_t m3_;
    //! Mass of parent particle
    const Double_t mParent_;

    //! Mass of particle 1 squared
    const Double_t m1Sq_;
    //! Mass of particle 2 squared
    const Double_t m2Sq_;
    //! Mass of particle 3 squared
    const Double_t m3Sq_;
    //! Mass of parent particle squared
    const Double_t mParentSq_;

    //! Vector of daughter particles masses
    std::vector<Double_t> mass_;
    //! Vector of the minimum mij values
    std::vector<Double_t> mMin_;
    //! Vector of the maximum mij values
    std::vector<Double_t> mMax_;
    //! Vector of the difference between the mMax and mMin
    std::vector<Double_t> mDiff_;

    //! Vector of daughter particles masses squared
    std::vector<Double_t> mSq_;
    //! Vector of the minimum mijSq values
    std::vector<Double_t> mSqMin_;
    //! Vector of the maximum mijSq values
    std::vector<Double_t> mSqMax_;
    //! Vector of the difference between the mSqMax and mSqMin
    std::vector<Double_t> mSqDiff_;

    //! Sum of the daughter masses
    const Double_t mDTot_;
    //! Mass difference between the parent particle and the sum of the daughter particles
    const Double_t massInt_;
    //! Sum of the squares of the daughter masses
    const Double_t mSqDTot_;

    //! Invariant mass m12
    Double_t m12_;
    //! Invariant mass m23
    Double_t m23_;
    //! Invariant mass m13
    Double_t m13_;

    //! Invariant mass m12 squared
    Double_t m12Sq_;
    //! Invariant mass m23 squared
    Double_t m23Sq_;
    //! Invariant mass m13 squared
    Double_t m13Sq_;

    //! Cosine of the helicity angle theta12, which is defined as the angle between 1&3 in the rest frame of 1&2
    Double_t c12_;
    //! Cosine of the helicity angle theta23, which is defined as the angle between 1&2 in the rest frame of 2&3
    Double_t c23_;
    //! Cosine of the helicity angle theta13, which is defined as the angle between 1&2 in the rest frame of 1&3
    Double_t c13_;

    //! m' co-ordinate
    Double_t mPrime_;
    //! theta' co-ordinate
    Double_t thetaPrime_;

    //! Momentum q of particle i
    mutable Double_t qi_;
    //! Momentum q of particle k
    mutable Double_t qk_;

    //! Momentum of track 1 in 1-2 rest frame
    Double_t p1_12_;
    //! Momentum of track 3 in 1-2 rest frame
    Double_t p3_12_;
    //! Momentum of track 2 in 2-3 rest frame
    Double_t p2_23_;
    //! Momentum of track 1 in 2-3 rest frame
    Double_t p1_23_;
    //! Momentum of track 1 in 1-3 rest frame
    Double_t p1_13_;
    //! Momentum of track 2 in 1-3 rest frame
    Double_t p2_13_;

    //! Momentum of track 1 in parent rest frame
    Double_t p1_Parent_;
    //! Momentum of track 2 in parent rest frame
    Double_t p2_Parent_;
    //! Momentum of track 3 in parent rest frame
    Double_t p3_Parent_;

    //! Should we calculate the square DP co-ordinates or not?
    Bool_t squareDP_;
    //! Enable/disable warning messages
    Bool_t warnings_;

    ClassDef( LauKinematics, 0 )
};

#endif
