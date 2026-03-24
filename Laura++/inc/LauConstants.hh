
/*
Copyright 2005 University of Warwick

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

/*! \file LauConstants.hh
    \brief File containing LauConstants namespace.
*/

/*! \namespace LauConstants
    \brief Namespace to contain various constants.
*/

#ifndef LAU_CONSTANTS
#define LAU_CONSTANTS

#include "Rtypes.h"
#include "TMath.h"

namespace LauConstants {

    // Masses updated to PDG 2020 values in March 2021
    //! Mass of charged D (GeV/c^2)
    const Double_t mD = 1.86965;
    //! Mass of neutral D (GeV/c^2)
    const Double_t mD0 = 1.86483;
    //! Mass of Ds (GeV/c^2)
    const Double_t mDs = 1.96834;
    //! Mass of charged B (GeV/c^2)
    const Double_t mB = 5.27934;
    //! Mass of neutral B_d (GeV/c^2)
    const Double_t mB0 = 5.27965;
    //! Mass of neutral B_s (GeV/c^2)
    const Double_t mBs0 = 5.36688;
    //! Mass of pi+- (GeV/c^2)
    const Double_t mPi = 0.13957039;
    //! Mass of pi0 (GeV/c^2)
    const Double_t mPi0 = 0.1349768;
    //! Mass of K+- (GeV/c^2)
    const Double_t mK = 0.493677;
    //! Mass of K0 (GeV/c^2)
    const Double_t mK0 = 0.497611;
    //! Mass of eta (GeV/c^2)
    const Double_t mEta = 0.547862;
    //! Mass of eta' (GeV/c^2)
    const Double_t mEtaPrime = 0.95778;

    //! Square of charged D mass
    const Double_t mDSq = mD * mD;
    //! Square of neutral D mass
    const Double_t mD0Sq = mD0 * mD0;
    //! Square of Ds mass
    const Double_t mDsSq = mDs * mDs;
    //! Square of charged B mass
    const Double_t mBSq = mB * mB;
    //! Square of neutral B_d mass
    const Double_t mB0Sq = mB0 * mB0;
    //! Square of neutral B_s mass
    const Double_t mBs0Sq = mBs0 * mBs0;
    //! Square of pi+- mass
    const Double_t mPiSq = mPi * mPi;
    //! Square of pi0 mass
    const Double_t mPi0Sq = mPi0 * mPi0;
    //! Square of K+- mass
    const Double_t mKSq = mK * mK;
    //! Square of K0 mass
    const Double_t mK0Sq = mK0 * mK0;
    //! Square of eta mass
    const Double_t mEtaSq = mEta * mEta;
    //! Square of eta' mass
    const Double_t mEtaPrimeSq = mEtaPrime * mEtaPrime;

    //! Lifetime of the B+ in ps
    const Double_t tauB = 1.638;
    //! Lifetime of the B0 in ps
    const Double_t tauB0 = 1.519;
    //! Mass difference of the B_H and B_L - source PDG via HFLAV 2019
    const Double_t deltaMd = 0.5065;
    //! Angle beta of the unitarity triangle - source PDG via HFLAV 2019
    const Double_t beta = 22.2 * TMath::DegToRad();
    //! Sine of twice the angle beta of the unitarity triangle
    const Double_t sin2beta = TMath::Sin( 2.0 * beta );

    //! Pi
    const Double_t pi = TMath::Pi();
    //! Square root of Pi
    const Double_t rootPi = TMath::Sqrt( TMath::Pi() );
    //! Two times Pi
    const Double_t twoPi = 2.0 * TMath::Pi();
    //! Three times Pi
    const Double_t threePi = 3.0 * TMath::Pi();
    //! Six times Pi
    const Double_t sixPi = 6.0 * TMath::Pi();
    //! One over Pi to the one-fourth
    const Double_t pim4 = 1.0 / TMath::Power( TMath::Pi(), 0.25 );
    //! Pi divided by two
    const Double_t piBy2 = 0.5 * TMath::Pi();
    //! Square root of Pi divided by two
    const Double_t rootPiBy2 = TMath::Sqrt( 0.5 * TMath::Pi() );
    //! One over Pi
    const Double_t invPi = 1.0 / TMath::Pi();
    //! Square root of two
    const Double_t root2 = TMath::Sqrt( 2.0 );
    //! Logarithm of four
    const Double_t log4 = TMath::Log( 4.0 );

};    // namespace LauConstants

#endif
