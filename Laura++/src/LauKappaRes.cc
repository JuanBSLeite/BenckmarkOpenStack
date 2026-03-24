
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

/*! \file LauKappaRes.cc
    \brief File containing implementation of LauKappaRes class.
*/

#include "LauKappaRes.hh"

#include "LauConstants.hh"
#include "LauResonanceInfo.hh"

#include "TMath.h"

#include <iostream>

LauKappaRes::LauKappaRes( LauResonanceInfo* resInfo,
                          const Int_t resPairAmpInt,
                          const LauDaughters* daughters ) :
    LauAbsResonance( resInfo, resPairAmpInt, daughters ),
    mSumSq_( ( LauConstants::mPi + LauConstants::mK ) * ( LauConstants::mPi + LauConstants::mK ) ),
    sAdler_( LauConstants::mKSq - 0.5 * LauConstants::mPiSq ),
    b1_( 0 ),
    b2_( 0 ),
    a_( 0 ),
    m0_( 0 )
{
    // Default constant factors from BES data
    const Double_t b1Val = 24.49;
    const Double_t b2Val = 0.0;
    const Double_t aVal  = 2.5;
    const Double_t m0Val = 3.3;

    const TString& parNameBase = this->getSanitisedName();

    TString b1Name( parNameBase );
    b1Name += "_b1";
    b1_     = resInfo->getExtraParameter( b1Name );
    if ( b1_ == 0 ) {
        b1_ = new LauParameter( b1Name, b1Val, 0.0, 100.0, kTRUE );
        b1_->secondStage( kTRUE );
        resInfo->addExtraParameter( b1_ );
    }

    TString b2Name( parNameBase );
    b2Name += "_b2";
    b2_     = resInfo->getExtraParameter( b2Name );
    if ( b2_ == 0 ) {
        b2_ = new LauParameter( b2Name, b2Val, 0.0, 100.0, kTRUE );
        b2_->secondStage( kTRUE );
        resInfo->addExtraParameter( b2_ );
    }

    TString aName( parNameBase );
    aName += "_A";
    a_     = resInfo->getExtraParameter( aName );
    if ( a_ == 0 ) {
        a_ = new LauParameter( aName, aVal, 0.0, 10.0, kTRUE );
        a_->secondStage( kTRUE );
        resInfo->addExtraParameter( a_ );
    }

    TString m0Name( parNameBase );
    m0Name += "_m0";
    m0_     = resInfo->getExtraParameter( m0Name );
    if ( m0_ == 0 ) {
        m0_ = new LauParameter( m0Name, m0Val, 0.0, 10.0, kTRUE );
        m0_->secondStage( kTRUE );
        resInfo->addExtraParameter( m0_ );
    }
}

LauKappaRes::~LauKappaRes()
{
}

void LauKappaRes::initialise()
{
    this->checkDaughterTypes();

    Double_t resSpin = this->getSpin();
    if ( resSpin != 0 ) {
        std::cerr << "WARNING in LauKappaRes::initialise : Resonance spin is " << resSpin << "."
                  << std::endl;
        std::cerr << "                                   : Kappa amplitude is only for scalers, resetting spin to 0."
                  << std::endl;
        this->changeResonance( -1.0, -1.0, 0 );
    }
}

void LauKappaRes::checkDaughterTypes() const
{
    // Check that the daughter tracks are K and pi. Otherwise issue a warning.

    Int_t resPairAmpInt = this->getPairInt();
    if ( resPairAmpInt < 1 || resPairAmpInt > 3 ) {
        std::cerr << "WARNING in LauKappaRes::checkDaughterTypes : resPairAmpInt = " << resPairAmpInt
                  << " is out of the range [1,2,3]." << std::endl;
        return;
    }

    const TString& nameDaug1 = this->getNameDaug1();
    const TString& nameDaug2 = this->getNameDaug2();
    if ( ! ( nameDaug1.Contains( "pi", TString::kIgnoreCase ) &&
             nameDaug2.Contains( "k", TString::kIgnoreCase ) ) ) {
        if ( ! ( nameDaug2.Contains( "pi", TString::kIgnoreCase ) &&
                 nameDaug1.Contains( "k", TString::kIgnoreCase ) ) ) {
            std::cerr << "ERROR in LauKappaRes::checkDaughterTypes : Kappa model is using daughters \""
                      << nameDaug1 << "\" and \"" << nameDaug2
                      << "\" that are not a kaon and a pion." << std::endl;
        }
    }
}

LauComplex LauKappaRes::resAmp( Double_t mass, Double_t spinTerm )
{
    // This function returns the complex dynamical amplitude for a Kappa distribution
    // given the invariant mass and cos(helicity) values.

    // First check that the appropriate daughters are either pi+pi- or K+K-
    // Check that the daughter tracks are the same type. Otherwise issue a warning
    // and set the type to be pion for the Kappa distribution. Returns the
    // integer defined by the enum LauKappaRes::KappaPartType.

    Double_t s = mass * mass;    // Invariant mass squared combination for the system
    Double_t rho( 0.0 );         // Phase-space factor
    if ( s > mSumSq_ ) {
        rho = TMath::Sqrt( 1.0 - mSumSq_ / s );
    }

    const Double_t m0Val = this->getM0Value();
    const Double_t m0Sq  = m0Val * m0Val;

    const Double_t aVal  = this->getAValue();
    const Double_t b1Val = this->getB1Value();
    const Double_t b2Val = this->getB2Value();

    Double_t f         = b2Val * s + b1Val;    // f(s) function
    Double_t numerator = s - sAdler_;
    Double_t denom     = m0Sq - sAdler_;
    Double_t gamma( 0.0 );
    if ( TMath::Abs( denom ) > 1e-10 && TMath::Abs( aVal ) > 1e-10 ) {
        // Decay width of the system
        gamma = rho * ( numerator / denom ) * f * TMath::Exp( -( s - m0Sq ) / aVal );
    }

    // Now form the complex amplitude - use relativistic BW form (without barrier factors)
    // Note that the M factor in the denominator is not the "pole" at ~500 MeV, but is
    // m0_ = 0.9264, the mass when the phase shift goes through 90 degrees.

    Double_t dMSq      = m0Sq - s;
    Double_t widthTerm = gamma * m0Val;
    LauComplex resAmplitude( dMSq, widthTerm );

    Double_t denomFactor = dMSq * dMSq + widthTerm * widthTerm;

    Double_t invDenomFactor = 0.0;
    if ( denomFactor > 1e-10 ) {
        invDenomFactor = 1.0 / denomFactor;
    }

    resAmplitude.rescale( spinTerm * invDenomFactor );

    return resAmplitude;
}

const std::vector<LauParameter*>& LauKappaRes::getFloatingParameters()
{
    this->clearFloatingParameters();

    if ( ! this->fixB1Value() ) {
        this->addFloatingParameter( b1_ );
    }

    if ( ! this->fixB2Value() ) {
        this->addFloatingParameter( b2_ );
    }

    if ( ! this->fixAValue() ) {
        this->addFloatingParameter( a_ );
    }

    if ( ! this->fixM0Value() ) {
        this->addFloatingParameter( m0_ );
    }

    return this->getParameters();
}

void LauKappaRes::setResonanceParameter( const TString& name, const Double_t value )
{
    // Set various parameters for the lineshape
    if ( name == "b1" ) {
        this->setB1Value( value );
        std::cout << "INFO in LauKappaRes::setResonanceParameter : Setting parameter b1 = "
                  << this->getB1Value() << std::endl;
    } else if ( name == "b2" ) {
        this->setB2Value( value );
        std::cout << "INFO in LauKappaRes::setResonanceParameter : Setting parameter b2 = "
                  << this->getB2Value() << std::endl;
    } else if ( name == "A" ) {
        this->setAValue( value );
        std::cout << "INFO in LauKappaRes::setResonanceParameter : Setting parameter A = "
                  << this->getAValue() << std::endl;
    } else if ( name == "m0" ) {
        this->setM0Value( value );
        std::cout << "INFO in LauKappaRes::setResonanceParameter : Setting parameter m0 = "
                  << this->getM0Value() << std::endl;
    } else {
        std::cerr << "WARNING in LauKappaRes::setResonanceParameter: Parameter name not reconised.  No parameter changes made."
                  << std::endl;
    }
}

void LauKappaRes::floatResonanceParameter( const TString& name )
{
    if ( name == "b1" ) {
        if ( b1_->fixed() ) {
            b1_->fixed( kFALSE );
            this->addFloatingParameter( b1_ );
        } else {
            std::cerr << "WARNING in LauKappaRes::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "b2" ) {
        if ( b2_->fixed() ) {
            b2_->fixed( kFALSE );
            this->addFloatingParameter( b2_ );
        } else {
            std::cerr << "WARNING in LauKappaRes::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "A" ) {
        if ( a_->fixed() ) {
            a_->fixed( kFALSE );
            this->addFloatingParameter( a_ );
        } else {
            std::cerr << "WARNING in LauKappaRes::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "m0" ) {
        if ( m0_->fixed() ) {
            m0_->fixed( kFALSE );
            this->addFloatingParameter( m0_ );
        } else {
            std::cerr << "WARNING in LauKappaRes::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else {
        std::cerr << "WARNING in LauKappaRes::fixResonanceParameter: Parameter name not reconised.  No parameter changes made."
                  << std::endl;
    }
}

LauParameter* LauKappaRes::getResonanceParameter( const TString& name )
{
    if ( name == "b1" ) {
        return b1_;
    } else if ( name == "b2" ) {
        return b2_;
    } else if ( name == "A" ) {
        return a_;
    } else if ( name == "m0" ) {
        return m0_;
    } else {
        std::cerr << "WARNING in LauKappaRes::getResonanceParameter: Parameter name not reconised."
                  << std::endl;
        return 0;
    }
}

void LauKappaRes::setB1Value( const Double_t b1 )
{
    b1_->value( b1 );
    b1_->genValue( b1 );
    b1_->initValue( b1 );
}

void LauKappaRes::setB2Value( const Double_t b2 )
{
    b2_->value( b2 );
    b2_->genValue( b2 );
    b2_->initValue( b2 );
}

void LauKappaRes::setAValue( const Double_t A )
{
    a_->value( A );
    a_->genValue( A );
    a_->initValue( A );
}

void LauKappaRes::setM0Value( const Double_t m0 )
{
    m0_->value( m0 );
    m0_->genValue( m0 );
    m0_->initValue( m0 );
}
