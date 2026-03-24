
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

/*! \file LauModIndPartWaveRealImag.cc
    \brief File containing implementation of LauModIndPartWaveRealImag class.
*/

#include "LauModIndPartWaveRealImag.hh"

#include "Lau1DCubicSpline.hh"
#include "LauConstants.hh"
#include "LauKinematics.hh"
#include "LauResonanceInfo.hh"

#include <cstdlib>
#include <iostream>

LauModIndPartWaveRealImag::LauModIndPartWaveRealImag( LauResonanceInfo* resInfo,
                                                      Int_t resPairAmpInt,
                                                      const LauDaughters* daughters ) :
    LauAbsModIndPartWave( resInfo, resPairAmpInt, daughters )
{
}

LauModIndPartWaveRealImag::~LauModIndPartWaveRealImag()
{
}

void LauModIndPartWaveRealImag::createAmpParameters( const UInt_t iKnot )
{
    const TString& parNameBase = this->getSanitisedName();
    const Bool_t secondStage   = this->floatKnotsSecondStage();

    std::vector<LauParameter*>& rePars = this->getAmp1Pars();
    std::vector<LauParameter*>& imPars = this->getAmp2Pars();

    TString reName( parNameBase );
    reName += "_Re";
    reName += iKnot;

    rePars.push_back( this->getResInfo()->getExtraParameter( reName ) );
    if ( rePars[iKnot] == 0 ) {
        rePars[iKnot] = new LauParameter( reName, 1.0, -10.0, 10.0, kFALSE );
        rePars[iKnot]->secondStage( secondStage );
        this->getResInfo()->addExtraParameter( rePars[iKnot] );
    }

    TString imName( parNameBase );
    imName += "_Im";
    imName += iKnot;

    imPars.push_back( this->getResInfo()->getExtraParameter( imName ) );
    if ( imPars[iKnot] == 0 ) {
        imPars[iKnot] = new LauParameter( imName, 1.0, -10.0, 10.0, kFALSE );
        imPars[iKnot]->secondStage( secondStage );
        this->getResInfo()->addExtraParameter( imPars[iKnot] );
    }
}

void LauModIndPartWaveRealImag::setKnotAmp( const UInt_t knot,
                                            const Double_t realVal,
                                            const Double_t imagVal,
                                            const Bool_t fixRealPart,
                                            const Bool_t fixImagPart )
{
    const UInt_t nknots = this->nKnots();

    if ( knot >= nknots ) {
        std::cerr << "WARNING in LauModIndPartWaveRealImag::setKnotAmp : Index " << knot
                  << " does not correspond to an existing knot in resonance "
                  << this->getResonanceName() << std::endl;
        std::cerr << "                                                 : Index must be in range 0 to "
                  << nknots - 1 << std::endl;
        return;
    }

    const std::vector<Double_t>& masses = this->getMasses();

    std::vector<Double_t>& realVals = this->getAmp1Vals();
    std::vector<Double_t>& imagVals = this->getAmp2Vals();

    std::vector<LauParameter*>& realPars = this->getAmp1Pars();
    std::vector<LauParameter*>& imagPars = this->getAmp2Pars();

    realVals[knot] = realVal;
    realPars[knot]->value( realVal );
    realPars[knot]->genValue( realVal );
    realPars[knot]->initValue( realVal );
    realPars[knot]->fixed( fixRealPart );
    imagVals[knot] = imagVal;
    imagPars[knot]->value( imagVal );
    imagPars[knot]->genValue( imagVal );
    imagPars[knot]->initValue( imagVal );
    imagPars[knot]->fixed( fixImagPart );

    if ( knot == 0 ) {
        std::cout << "INFO in LauModIndPartWaveRealImag::setKnotAmp : Knot updated in resonance "
                  << this->getResonanceName() << " at lower kinematic limit (" << masses[knot]
                  << ")" << std::endl;
    } else if ( knot == nknots - 1 ) {
        std::cout << "INFO in LauModIndPartWaveRealImag::setKnotAmp : Knot updated in resonance "
                  << this->getResonanceName() << " at upper kinematic limit (" << masses[knot]
                  << ")" << std::endl;
    } else {
        std::cout << "INFO in LauModIndPartWaveRealImag::setKnotAmp : Knot updated in resonance "
                  << this->getResonanceName() << " at mass " << masses[knot] << std::endl;
    }

    if ( fixRealPart ) {
        std::cout << "                                              : Real part fixed to "
                  << realVal << std::endl;
    } else {
        std::cout << "                                              : Real part set to " << realVal
                  << std::endl;
    }

    if ( fixImagPart ) {
        std::cout << "                                              : Imaginary part fixed to "
                  << imagVal << std::endl;
    } else {
        std::cout << "                                              : Imaginary part set to "
                  << imagVal << std::endl;
    }
}

void LauModIndPartWaveRealImag::evaluateAmplitude( const Double_t mass )
{
    const Lau1DCubicSpline* splineReal = this->getSpline1();
    const Lau1DCubicSpline* splineImag = this->getSpline2();

    const Double_t re = splineReal->evaluate( mass );
    const Double_t im = splineImag->evaluate( mass );

    this->setAmp( re, im );
}
