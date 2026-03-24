
/*
Copyright 2010 University of Warwick

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

/*! \file LauDabbaRes.cc
    \brief File containing implementation of LauDabbaRes class.

    Formulae and data values from arXiv:0901.2217 - author D.V.Bugg
*/

#include "LauDabbaRes.hh"

#include "LauConstants.hh"
#include "LauResonanceInfo.hh"

#include "TMath.h"

#include <iostream>

LauDabbaRes::LauDabbaRes( LauResonanceInfo* resInfo,
                          const Int_t resPairAmpInt,
                          const LauDaughters* daughters ) :
    LauAbsResonance( resInfo, resPairAmpInt, daughters ),
    mSumSq_( 0.0 ),
    sAdler_( 0.0 ),
    b_( 0 ),
    alpha_( 0 ),
    beta_( 0 )
{
    // Default constant factors
    const Double_t bVal     = 24.49;
    const Double_t alphaVal = 0.1;
    const Double_t betaVal  = 0.1;

    const TString& parNameBase = this->getSanitisedName();

    TString bName( parNameBase );
    bName += "_b";
    b_     = resInfo->getExtraParameter( bName );
    if ( b_ == 0 ) {
        b_ = new LauParameter( bName, bVal, 0.0, 100.0, kTRUE );
        b_->secondStage( kTRUE );
        resInfo->addExtraParameter( b_ );
    }

    TString alphaName( parNameBase );
    alphaName += "_alpha";
    alpha_     = resInfo->getExtraParameter( alphaName );
    if ( alpha_ == 0 ) {
        alpha_ = new LauParameter( alphaName, alphaVal, 0.0, 10.0, kTRUE );
        alpha_->secondStage( kTRUE );
        resInfo->addExtraParameter( alpha_ );
    }

    TString betaName( parNameBase );
    betaName += "_beta";
    beta_     = resInfo->getExtraParameter( betaName );
    if ( beta_ == 0 ) {
        beta_ = new LauParameter( betaName, betaVal, 0.0, 10.0, kTRUE );
        beta_->secondStage( kTRUE );
        resInfo->addExtraParameter( beta_ );
    }
}

LauDabbaRes::~LauDabbaRes()
{
}

void LauDabbaRes::initialise()
{
    // check that we have a D and a pi
    this->checkDaughterTypes();

    // Initialise various constants
    Double_t massDaug1 = this->getMassDaug1();
    Double_t massDaug2 = this->getMassDaug2();

    Double_t mSum = massDaug1 + massDaug2;
    mSumSq_       = mSum * mSum;

    Double_t massDaug1Sq = massDaug1 * massDaug1;
    Double_t massDaug2Sq = massDaug2 * massDaug2;
    sAdler_              = TMath::Max( massDaug1Sq, massDaug2Sq ) -
              0.5 * TMath::Min( massDaug1Sq, massDaug2Sq );    // Adler zero at (mD)^2 - 0.5*(mpi)^2

    Int_t resSpin = this->getSpin();
    if ( resSpin != 0 ) {
        std::cerr << "WARNING in LauDabbaRes::initialise : Spin = " << resSpin
                  << " is not zero!  It will be ignored anyway!" << std::endl;
    }
}

void LauDabbaRes::checkDaughterTypes() const
{
    // Check that the daughter tracks are D and pi. Otherwise issue a warning.
    Int_t resPairAmpInt = this->getPairInt();
    if ( resPairAmpInt < 1 || resPairAmpInt > 3 ) {
        std::cerr << "WARNING in LauDabbaRes::checkDaughterTypes : resPairAmpInt = " << resPairAmpInt
                  << " is out of the range [1,2,3]." << std::endl;
        return;
    }

    // Check that daughter types agree
    const TString& nameDaug1 = this->getNameDaug1();
    const TString& nameDaug2 = this->getNameDaug2();
    if ( ! ( nameDaug1.Contains( "pi", TString::kIgnoreCase ) &&
             nameDaug2.Contains( "d", TString::kIgnoreCase ) ) ) {
        if ( ! ( nameDaug2.Contains( "pi", TString::kIgnoreCase ) &&
                 nameDaug1.Contains( "d", TString::kIgnoreCase ) ) ) {
            std::cerr << "ERROR in LauDabbaRes::checkDaughterTypes : Dabba model is using daughters \""
                      << nameDaug1 << "\" and \"" << nameDaug2
                      << "\" that are not a D meson and a pion." << std::endl;
        }
    }
}

LauComplex LauDabbaRes::resAmp( Double_t mass, Double_t spinTerm )
{
    // This function returns the complex dynamical amplitude for a Dabba distribution
    // given the invariant mass and cos(helicity) values.

    // Invariant mass squared combination for the system
    Double_t s = mass * mass;

    // Dabba is spin zero - so there are no helicity factors.
    // Just set it to 1.0 in case anyone decides to use it at a later date.
    spinTerm = 1.0;

    // Phase-space factor
    Double_t rho( 0.0 );
    Double_t sDiff = s - mSumSq_;
    if ( sDiff > 0.0 ) {
        rho = TMath::Sqrt( 1.0 - mSumSq_ / s );
    }

    const Double_t bVal     = this->getBValue();
    const Double_t alphaVal = this->getAlphaValue();
    const Double_t betaVal  = this->getBetaValue();

    Double_t realPart = 1.0 - betaVal * sDiff;

    Double_t imagPart = bVal * TMath::Exp( -alphaVal * sDiff ) * ( s - sAdler_ ) * rho;

    LauComplex resAmplitude( realPart, imagPart );

    Double_t denomFactor = realPart * realPart + imagPart * imagPart;

    Double_t invDenomFactor = 0.0;
    if ( denomFactor > 1e-10 ) {
        invDenomFactor = 1.0 / denomFactor;
    }

    resAmplitude.rescale( spinTerm * invDenomFactor );

    return resAmplitude;
}

const std::vector<LauParameter*>& LauDabbaRes::getFloatingParameters()
{
    this->clearFloatingParameters();

    if ( ! this->fixBValue() ) {
        this->addFloatingParameter( b_ );
    }

    if ( ! this->fixAlphaValue() ) {
        this->addFloatingParameter( alpha_ );
    }

    if ( ! this->fixBetaValue() ) {
        this->addFloatingParameter( beta_ );
    }

    return this->getParameters();
}

void LauDabbaRes::setResonanceParameter( const TString& name, const Double_t value )
{
    // Set various parameters for the lineshape
    if ( name == "b" ) {
        this->setBValue( value );
        std::cout << "INFO in LauDabbaRes::setResonanceParameter : Setting parameter b = "
                  << this->getBValue() << std::endl;
    } else if ( name == "alpha" ) {
        this->setAlphaValue( value );
        std::cout << "INFO in LauDabbaRes::setResonanceParameter : Setting parameter alpha = "
                  << this->getAlphaValue() << std::endl;
    } else if ( name == "beta" ) {
        this->setBetaValue( value );
        std::cout << "INFO in LauDabbaRes::setResonanceParameter : Setting parameter beta = "
                  << this->getBetaValue() << std::endl;
    } else {
        std::cerr << "WARNING in LauDabbaRes::setResonanceParameter: Parameter name not reconised.  No parameter changes made."
                  << std::endl;
    }
}

void LauDabbaRes::floatResonanceParameter( const TString& name )
{
    if ( name == "b" ) {
        if ( b_->fixed() ) {
            b_->fixed( kFALSE );
            this->addFloatingParameter( b_ );
        } else {
            std::cerr << "WARNING in LauDabbaRes::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "alpha" ) {
        if ( alpha_->fixed() ) {
            alpha_->fixed( kFALSE );
            this->addFloatingParameter( alpha_ );
        } else {
            std::cerr << "WARNING in LauDabbaRes::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "beta" ) {
        if ( beta_->fixed() ) {
            beta_->fixed( kFALSE );
            this->addFloatingParameter( beta_ );
        } else {
            std::cerr << "WARNING in LauDabbaRes::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else {
        std::cerr << "WARNING in LauDabbaRes::fixResonanceParameter: Parameter name not reconised.  No parameter changes made."
                  << std::endl;
    }
}

LauParameter* LauDabbaRes::getResonanceParameter( const TString& name )
{
    if ( name == "b" ) {
        return b_;
    } else if ( name == "alpha" ) {
        return alpha_;
    } else if ( name == "beta" ) {
        return beta_;
    } else {
        std::cerr << "WARNING in LauDabbaRes::getResonanceParameter: Parameter name not reconised."
                  << std::endl;
        return 0;
    }
}

void LauDabbaRes::setBValue( const Double_t b )
{
    b_->value( b );
    b_->genValue( b );
    b_->initValue( b );
}

void LauDabbaRes::setAlphaValue( const Double_t alpha )
{
    alpha_->value( alpha );
    alpha_->genValue( alpha );
    alpha_->initValue( alpha );
}

void LauDabbaRes::setBetaValue( const Double_t beta )
{
    beta_->value( beta );
    beta_->genValue( beta );
    beta_->initValue( beta );
}
