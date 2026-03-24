
/*
Copyright 2008 University of Warwick

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

/*! \file LauLASSNRRes.cc
    \brief File containing implementation of LauLASSNRRes class.
*/

#include "LauLASSNRRes.hh"

#include "LauConstants.hh"
#include "LauResonanceInfo.hh"

#include <iostream>

LauLASSNRRes::LauLASSNRRes( LauResonanceInfo* resInfo,
                            const Int_t resPairAmpInt,
                            const LauDaughters* daughters ) :
    LauAbsResonance( resInfo, resPairAmpInt, daughters ),
    r_( 0 ),
    a_( 0 ),
    cutOff_( 0.0 )
{
    // Default values for LASS parameters
    cutOff_             = 1.8;
    const Double_t rVal = 3.32;
    const Double_t aVal = 2.07;

    const TString& parNameBase = this->getSanitisedName();

    TString rName( parNameBase );
    rName += "_r";
    r_     = resInfo->getExtraParameter( rName );
    if ( r_ == 0 ) {
        r_ = new LauParameter( rName, rVal, 0.0, 10.0, kTRUE );
        r_->secondStage( kTRUE );
        resInfo->addExtraParameter( r_ );
    }

    TString aName( parNameBase );
    aName += "_a";
    a_     = resInfo->getExtraParameter( aName );
    if ( a_ == 0 ) {
        a_ = new LauParameter( aName, aVal, 0.0, 10.0, kTRUE );
        a_->secondStage( kTRUE );
        resInfo->addExtraParameter( a_ );
    }
}

LauLASSNRRes::~LauLASSNRRes()
{
}

void LauLASSNRRes::initialise()
{
    Int_t resSpin = this->getSpin();
    if ( resSpin != 0 ) {
        std::cerr << "WARNING in LauLASSNRRes::amplitude : Resonance spin is " << resSpin << "."
                  << std::endl;
        std::cerr << "                                   : LASS amplitude is only for scalers, resetting spin to 0."
                  << std::endl;
        this->changeResonance( -1.0, -1.0, 0 );
    }
}

LauComplex LauLASSNRRes::resAmp( Double_t mass, Double_t /*spinTerm*/ )
{
    LauComplex bkgAmplitude( 0.0, 0.0 );

    if ( mass < 1e-10 ) {
        std::cerr << "WARNING in LauLASSNRRes::amplitude : mass < 1e-10." << std::endl;
        return LauComplex( 0.0, 0.0 );
    }

    if ( mass < cutOff_ ) {
        // q is the momentum of either daughter in the resonance rest-frame
        const Double_t q = this->getQ();

        // Calculate the phase shift term
        const Double_t rVal       = this->getEffectiveRange();
        const Double_t aVal       = this->getScatteringLength();
        const Double_t qcotdeltaB = 1.0 / aVal + ( rVal * q * q ) / 2.0;

        // Compute the complex amplitude
        bkgAmplitude = LauComplex( qcotdeltaB, q );

        // Scale by the numerator and denominator factors
        bkgAmplitude.rescale( mass / ( qcotdeltaB * qcotdeltaB + q * q ) );
    }

    return bkgAmplitude;
}

const std::vector<LauParameter*>& LauLASSNRRes::getFloatingParameters()
{
    this->clearFloatingParameters();

    if ( ! this->fixMass() ) {
        this->addFloatingParameter( this->getMassPar() );
    }

    if ( ! this->fixWidth() ) {
        this->addFloatingParameter( this->getWidthPar() );
    }

    if ( ! this->fixEffectiveRange() ) {
        this->addFloatingParameter( r_ );
    }

    if ( ! this->fixScatteringLength() ) {
        this->addFloatingParameter( a_ );
    }

    return this->getParameters();
}

void LauLASSNRRes::setResonanceParameter( const TString& name, const Double_t value )
{
    // Set various parameters for the LASS lineshape dynamics
    if ( name == "a" ) {
        this->setScatteringLength( value );
        std::cout << "INFO in LauLASSNRRes::setResonanceParameter : Setting LASS Scattering Length = "
                  << this->getScatteringLength() << std::endl;
    } else if ( name == "r" ) {
        this->setEffectiveRange( value );
        std::cout << "INFO in LauLASSNRRes::setResonanceParameter : Setting LASS Effective Range = "
                  << this->getEffectiveRange() << std::endl;
    } else {
        std::cerr << "WARNING in LauLASSNRRes::setResonanceParameter: Parameter name not reconised.  No parameter changes made."
                  << std::endl;
    }
}

void LauLASSNRRes::floatResonanceParameter( const TString& name )
{
    if ( name == "a" ) {
        if ( a_->fixed() ) {
            a_->fixed( kFALSE );
            this->addFloatingParameter( a_ );
        } else {
            std::cerr << "WARNING in LauLASSNRRes::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "r" ) {
        if ( r_->fixed() ) {
            r_->fixed( kFALSE );
            this->addFloatingParameter( r_ );
        } else {
            std::cerr << "WARNING in LauLASSNRRes::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else {
        std::cerr << "WARNING in LauLASSNRRes::fixResonanceParameter: Parameter name not reconised.  No parameter changes made."
                  << std::endl;
    }
}

LauParameter* LauLASSNRRes::getResonanceParameter( const TString& name )
{
    if ( name == "a" ) {
        return a_;
    } else if ( name == "r" ) {
        return r_;
    } else {
        std::cerr << "WARNING in LauLASSNRRes::getResonanceParameter: Parameter name not reconised."
                  << std::endl;
        return 0;
    }
}

void LauLASSNRRes::setEffectiveRange( const Double_t r )
{
    r_->value( r );
    r_->genValue( r );
    r_->initValue( r );
}

void LauLASSNRRes::setScatteringLength( const Double_t a )
{
    a_->value( a );
    a_->genValue( a );
    a_->initValue( a );
}
