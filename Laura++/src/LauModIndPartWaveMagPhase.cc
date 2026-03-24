
/*
Copyright 2015 University of Warwick

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

/*! \file LauModIndPartWaveMagPhase.cc
    \brief File containing implementation of LauModIndPartWaveMagPhase class.
*/

#include "LauModIndPartWaveMagPhase.hh"

#include "Lau1DCubicSpline.hh"
#include "LauConstants.hh"
#include "LauKinematics.hh"
#include "LauResonanceInfo.hh"

#include "TMath.h"

#include <cstdlib>
#include <iostream>

LauModIndPartWaveMagPhase::LauModIndPartWaveMagPhase( LauResonanceInfo* resInfo,
                                                      Int_t resPairAmpInt,
                                                      const LauDaughters* daughters ) :
    LauAbsModIndPartWave( resInfo, resPairAmpInt, daughters )
{
}

LauModIndPartWaveMagPhase::~LauModIndPartWaveMagPhase()
{
}

void LauModIndPartWaveMagPhase::createAmpParameters( const UInt_t iKnot )
{
    const TString& parNameBase = this->getSanitisedName();
    const Bool_t secondStage   = this->floatKnotsSecondStage();

    std::vector<LauParameter*>& magnitudePars = this->getAmp1Pars();
    std::vector<LauParameter*>& phasePars     = this->getAmp2Pars();

    TString magName( parNameBase );
    magName += "_A";
    magName += iKnot;

    magnitudePars.push_back( this->getResInfo()->getExtraParameter( magName ) );
    if ( magnitudePars[iKnot] == 0 ) {
        magnitudePars[iKnot] = new LauParameter( magName, 1.0, 0.0, 10.0, kFALSE );
        magnitudePars[iKnot]->secondStage( secondStage );
        this->getResInfo()->addExtraParameter( magnitudePars[iKnot] );
    }

    TString phaseName( parNameBase );
    phaseName += "_d";
    phaseName += iKnot;

    phasePars.push_back( this->getResInfo()->getExtraParameter( phaseName ) );
    if ( phasePars[iKnot] == 0 ) {
        phasePars[iKnot] =
            new LauParameter( phaseName, 1.0, -6.0 * LauConstants::pi, 6.0 * LauConstants::pi, kFALSE );
        phasePars[iKnot]->secondStage( secondStage );
        this->getResInfo()->addExtraParameter( phasePars[iKnot] );
    }
}

void LauModIndPartWaveMagPhase::setKnotAmp( const UInt_t knot,
                                            const Double_t magVal,
                                            const Double_t phaseVal,
                                            const Bool_t fixMagnitude,
                                            const Bool_t fixPhase )
{
    const UInt_t nknots = this->nKnots();

    if ( knot >= nknots ) {
        std::cerr << "WARNING in LauModIndPartWaveMagPhase::setKnotAmp : Index " << knot
                  << " does not correspond to an existing knot in resonance "
                  << this->getResonanceName() << std::endl;
        std::cerr << "                                                 : Index must be in range 0 to "
                  << nknots - 1 << std::endl;
        return;
    }

    const std::vector<Double_t>& masses = this->getMasses();

    std::vector<Double_t>& magnitudes = this->getAmp1Vals();
    std::vector<Double_t>& phases     = this->getAmp2Vals();

    std::vector<LauParameter*>& magnitudePars = this->getAmp1Pars();
    std::vector<LauParameter*>& phasePars     = this->getAmp2Pars();

    magnitudes[knot] = magVal;
    magnitudePars[knot]->value( magVal );
    magnitudePars[knot]->genValue( magVal );
    magnitudePars[knot]->initValue( magVal );
    magnitudePars[knot]->fixed( fixMagnitude );
    phases[knot] = phaseVal;
    phasePars[knot]->value( phaseVal );
    phasePars[knot]->genValue( phaseVal );
    phasePars[knot]->initValue( phaseVal );
    phasePars[knot]->fixed( fixPhase );

    if ( knot == 0 ) {
        std::cout << "INFO in LauModIndPartWaveMagPhase::setKnotAmp : Knot updated in resonance "
                  << this->getResonanceName() << " at lower kinematic limit (" << masses[knot]
                  << ")" << std::endl;
    } else if ( knot == nknots - 1 ) {
        std::cout << "INFO in LauModIndPartWaveMagPhase::setKnotAmp : Knot updated in resonance "
                  << this->getResonanceName() << " at upper kinematic limit (" << masses[knot]
                  << ")" << std::endl;
    } else {
        std::cout << "INFO in LauModIndPartWaveMagPhase::setKnotAmp : Knot updated in resonance "
                  << this->getResonanceName() << " at mass " << masses[knot] << std::endl;
    }

    if ( fixMagnitude ) {
        std::cout << "                                              : Magnitude fixed to " << magVal
                  << std::endl;
    } else {
        std::cout << "                                              : Magnitude set to " << magVal
                  << std::endl;
    }

    if ( fixPhase ) {
        std::cout << "                                              : Phase fixed to " << phaseVal
                  << std::endl;
    } else {
        std::cout << "                                              : Phase set to " << phaseVal
                  << std::endl;
    }
}

void LauModIndPartWaveMagPhase::evaluateAmplitude( const Double_t mass )
{
    const Lau1DCubicSpline* splineMag   = this->getSpline1();
    const Lau1DCubicSpline* splinePhase = this->getSpline2();

    const Double_t mag   = splineMag->evaluate( mass );
    const Double_t phase = splinePhase->evaluate( mass );

    this->setAmp( mag * TMath::Cos( phase ), mag * TMath::Sin( phase ) );
}
