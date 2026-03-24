
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

/*! \file LauVetoes.cc
    \brief File containing implementation of LauVetoes class.
*/

#include "LauVetoes.hh"

#include "LauKinematics.hh"

#include "TMath.h"

#include <iostream>

LauVetoes::LauVetoes() :
    nVetoes_( 0 )
{
}

LauVetoes::~LauVetoes()
{
}

LauVetoes::LauVetoes( const LauVetoes& other ) :
    nVetoes_( other.nVetoes_ ),
    vetoPair_( other.vetoPair_ ),
    vetoMinMass_( other.vetoMinMass_ ),
    vetoMaxMass_( other.vetoMaxMass_ )
{
}

LauVetoes& LauVetoes::operator=( const LauVetoes& other )
{
    if ( &other != this ) {
        nVetoes_     = other.nVetoes_;
        vetoPair_    = other.vetoPair_;
        vetoMinMass_ = other.vetoMinMass_;
        vetoMaxMass_ = other.vetoMaxMass_;
    }
    return *this;
}

void LauVetoes::addMassVeto( const Int_t resPairAmpInt, const Double_t minMass, const Double_t maxMass )
{

    const Double_t minMassSq = minMass * minMass;
    const Double_t maxMassSq = maxMass * maxMass;

    this->addMassSqVeto( resPairAmpInt, minMassSq, maxMassSq );
}

void LauVetoes::addMassSqVeto( const Int_t resPairAmpInt,
                               const Double_t minMassSq,
                               const Double_t maxMassSq )
{
    // Routine to add a veto in the Dalitz plot. The function takes as input the
    // bachelor track number (1, 2 or 3) and the mass-squared range of the veto.

    if ( resPairAmpInt == 1 ) {

        // The bachelor track is the first track
        std::cout << "INFO in LauVetoes::addMassSqVeto : Adding the veto for resPairAmpInt = 1, with "
                  << minMassSq << " < m^2_23 < " << maxMassSq << std::endl;

    } else if ( resPairAmpInt == 2 ) {

        // The bachelor track is the second track
        std::cout << "INFO in LauVetoes::addMassSqVeto : Adding the veto for resPairAmpInt = 2, with "
                  << minMassSq << " < m^2_13 < " << maxMassSq << std::endl;

    } else if ( resPairAmpInt == 3 ) {

        // The bachelor track is the third track
        std::cout << "INFO in LauVetoes::addMassSqVeto : Adding the veto for resPairAmpInt = 3, with "
                  << minMassSq << " < m^2_12 < " << maxMassSq << std::endl;

    } else if ( resPairAmpInt == 4 ) {

        // Special case for symmetric DPs - the veto will be applied on the minimum of m13Sq and m23Sq
        std::cout << "INFO in LauVetoes::addMassSqVeto : Adding the veto for resPairAmpInt = 4, with "
                  << minMassSq << " < m^2_min < " << maxMassSq << std::endl;

    } else if ( resPairAmpInt == 5 ) {

        // Special case for symmetric DPs - the veto will be applied on the maximum of m13Sq and m23Sq
        std::cout << "INFO in LauVetoes::addMassSqVeto : Adding the veto for resPairAmpInt = 5, with "
                  << minMassSq << " < m^2_max < " << maxMassSq << std::endl;

    } else {
        std::cerr << "ERROR in LauVetoes::addMassSqVeto : Invalid resPairAmpInt. Please use 1, 2 or 3 to specify bachelor daughter track (or 4 or 5 to specify a veto on mMinSq or mMaxSq in a symmetric DP). Veto is not added."
                  << std::endl;
        return;
    }

    // Set the veto limits
    vetoPair_.push_back( resPairAmpInt );
    vetoMinMass_.push_back( minMassSq );
    vetoMaxMass_.push_back( maxMassSq );

    // Keep track of how many vetoes we have
    ++nVetoes_;
}

Bool_t LauVetoes::passVeto( const LauKinematics* kinematics ) const
{
    // Routine to ask whether the given Dalitz plot point passes any specified vetoes.
    if ( kinematics == 0 ) {
        std::cerr << "ERROR in LauVetoes::passVeto : LauKinematics object is null." << std::endl;
        return kFALSE;
    }

    const Double_t m12Sq          = kinematics->getm12Sq();
    const Double_t m23Sq          = kinematics->getm23Sq();
    const Double_t m13Sq          = kinematics->getm13Sq();
    const Bool_t symmetricDP      = kinematics->gotSymmetricalDP();
    const Bool_t fullySymmetricDP = kinematics->gotFullySymmetricDP();

    return this->passVeto( m12Sq, m23Sq, m13Sq, symmetricDP, fullySymmetricDP );
}

Bool_t LauVetoes::passVeto( const Double_t m12Sq,
                            const Double_t m23Sq,
                            const Double_t m13Sq,
                            const Bool_t symmetricDP,
                            const Bool_t fullySymmetricDP ) const
{
    // Routine to ask whether the given Dalitz plot point passes any specified vetoes.

    // Loop over the number of possible vetoes
    for ( UInt_t i( 0 ); i < nVetoes_; ++i ) {

        if ( vetoPair_[i] == 1 ) {
            // Veto m23 combination
            if ( m23Sq > vetoMinMass_[i] && m23Sq < vetoMaxMass_[i] ) {
                return kFALSE;
            }
            // If the DP is symmetric we need to test m13 combination as well
            if ( symmetricDP || fullySymmetricDP ) {
                if ( m13Sq > vetoMinMass_[i] && m13Sq < vetoMaxMass_[i] ) {
                    return kFALSE;
                }
            }
            // If it's fully symmetric we need to check all 3 combinations
            if ( fullySymmetricDP ) {
                if ( m12Sq > vetoMinMass_[i] && m12Sq < vetoMaxMass_[i] ) {
                    return kFALSE;
                }
            }
        } else if ( vetoPair_[i] == 2 ) {
            // Veto m13 combination
            if ( m13Sq > vetoMinMass_[i] && m13Sq < vetoMaxMass_[i] ) {
                return kFALSE;
            }
            // If the DP is symmetric we need to test m23 combination as well
            if ( symmetricDP || fullySymmetricDP ) {
                if ( m23Sq > vetoMinMass_[i] && m23Sq < vetoMaxMass_[i] ) {
                    return kFALSE;
                }
            }
            // If it's fully symmetric we need to check all 3 combinations
            if ( fullySymmetricDP ) {
                if ( m12Sq > vetoMinMass_[i] && m12Sq < vetoMaxMass_[i] ) {
                    return kFALSE;
                }
            }
        } else if ( vetoPair_[i] == 3 ) {
            // Veto m12 combination
            if ( m12Sq > vetoMinMass_[i] && m12Sq < vetoMaxMass_[i] ) {
                return kFALSE;
            }
            // If it's fully symmetric we need to check all 3 combinations
            if ( fullySymmetricDP ) {
                if ( m13Sq > vetoMinMass_[i] && m13Sq < vetoMaxMass_[i] ) {
                    return kFALSE;
                }
                if ( m23Sq > vetoMinMass_[i] && m23Sq < vetoMaxMass_[i] ) {
                    return kFALSE;
                }
            }
        } else if ( vetoPair_[i] == 4 ) {
            if ( ! symmetricDP ) {
                std::cerr << "WARNING in LauVetoes::passVeto : resPairAmpInt of 4 is only valid for symmetric DPs, will ignore this veto"
                          << std::endl;
                continue;
            }
            // Veto mMin combination
            const Double_t mMinSq = TMath::Min( m13Sq, m23Sq );
            if ( mMinSq > vetoMinMass_[i] && mMinSq < vetoMaxMass_[i] ) {
                return kFALSE;
            }
        } else if ( vetoPair_[i] == 5 ) {
            if ( ! symmetricDP ) {
                std::cerr << "WARNING in LauVetoes::passVeto : resPairAmpInt of 5 is only valid for symmetric DPs, will ignore this veto"
                          << std::endl;
                continue;
            }
            // Veto mMax combination
            const Double_t mMaxSq = TMath::Max( m13Sq, m23Sq );
            if ( mMaxSq > vetoMinMass_[i] && mMaxSq < vetoMaxMass_[i] ) {
                return kFALSE;
            }
        }
    }

    return kTRUE;
}
