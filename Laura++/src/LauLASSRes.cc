
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

/*! \file LauLASSRes.cc
    \brief File containing implementation of LauLASSRes class.
*/

#include "LauLASSRes.hh"

#include "LauConstants.hh"
#include "LauResonanceInfo.hh"

#include "TMath.h"

#include <iostream>

LauLASSRes::LauLASSRes( LauResonanceInfo* resInfo,
                        const Int_t resPairAmpInt,
                        const LauDaughters* daughters ) :
    LauAbsResonance( resInfo, resPairAmpInt, daughters ),
    q0_( 0.0 ),
    mDaugSum_( 0.0 ),
    mDaugSumSq_( 0.0 ),
    mDaugDiff_( 0.0 ),
    mDaugDiffSq_( 0.0 ),
    resMassSq_( 0.0 ),
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

LauLASSRes::~LauLASSRes()
{
}

void LauLASSRes::initialise()
{
    // Create the mass sums and differences
    Double_t massDaug1 = this->getMassDaug1();
    Double_t massDaug2 = this->getMassDaug2();

    mDaugSum_   = massDaug1 + massDaug2;
    mDaugSumSq_ = mDaugSum_ * mDaugSum_;

    mDaugDiff_   = massDaug1 - massDaug2;
    mDaugDiffSq_ = mDaugDiff_ * mDaugDiff_;

    Int_t resSpin = this->getSpin();
    if ( resSpin != 0 ) {
        std::cerr << "WARNING in LauLASSRes::amplitude : Resonance spin is " << resSpin << "."
                  << std::endl;
        std::cerr << "                                 : LASS amplitude is only for scalers, resetting spin to 0."
                  << std::endl;
        this->changeResonance( -1.0, -1.0, 0 );
    }

    this->calcQ0();
}

void LauLASSRes::calcQ0()
{
    // Decay momentum of either daughter in the resonance rest frame
    // when resonance mass = rest-mass value, m_0 (PDG value)

    resMass_   = this->getMass();
    resMassSq_ = resMass_ * resMass_;

    q0_ = TMath::Sqrt( ( resMassSq_ - mDaugSumSq_ ) * ( resMassSq_ - mDaugDiffSq_ ) ) /
          ( 2.0 * resMass_ );
}

LauComplex LauLASSRes::resAmp( Double_t mass, Double_t /*spinTerm*/ )
{
    LauComplex resAmplitude( 0.0, 0.0 );
    LauComplex bkgAmplitude( 0.0, 0.0 );
    LauComplex totAmplitude( 0.0, 0.0 );

    if ( mass < 1e-10 ) {
        std::cerr << "WARNING in LauLASSRes::amplitude : mass < 1e-10." << std::endl;
        return LauComplex( 0.0, 0.0 );
    }

    //---------------------------
    // First do the resonant part
    //---------------------------

    // Calculate the width of the resonance (as a function of mass)
    // q is the momentum of either daughter in the resonance rest-frame
    const Double_t q        = this->getQ();
    const Double_t resMass  = this->getMass();
    const Double_t resWidth = this->getWidth();

    // If the mass is floating and their value have changed
    // we need to recalculate everything that assumes this value
    if ( ( ! this->fixMass() ) && resMass != resMass_ ) {
        this->calcQ0();
    }

    Double_t qRatio   = q / q0_;
    Double_t totWidth = resWidth * qRatio * ( resMass / mass );

    Double_t massSq     = mass * mass;
    Double_t massSqTerm = resMassSq_ - massSq;

    // Compute the complex amplitude
    resAmplitude = LauComplex( massSqTerm, resMass * totWidth );

    // Scale by the denominator factor
    resAmplitude.rescale( ( resMassSq_ * resWidth / q0_ ) /
                          ( massSqTerm * massSqTerm + resMassSq_ * totWidth * totWidth ) );

    // Calculate the phase shift term
    const Double_t rVal           = this->getEffectiveRange();
    const Double_t aVal           = this->getScatteringLength();
    const Double_t tandeltaB      = ( 2.0 * aVal * q ) / ( 2.0 + aVal * rVal * q * q );
    const Double_t tanSq          = tandeltaB * tandeltaB;
    const Double_t cos2PhaseShift = ( 1.0 - tanSq ) / ( 1.0 + tanSq );
    const Double_t sin2PhaseShift = 2.0 * tandeltaB / ( 1.0 + tanSq );
    LauComplex phaseShift( cos2PhaseShift, sin2PhaseShift );

    // Multiply by the phase shift
    resAmplitude = resAmplitude * phaseShift;

    //--------------------------------
    // Now do the effective range part
    //--------------------------------

    // Form the real and imaginary parts
    const Double_t qcotdeltaB = 1.0 / aVal + ( rVal * q * q ) / 2.0;

    // Compute the complex amplitude
    bkgAmplitude = LauComplex( qcotdeltaB, q );

    // Scale by the numerator and denominator factors
    bkgAmplitude.rescale( mass / ( qcotdeltaB * qcotdeltaB + q * q ) );

    //------------------
    // Add them together
    //------------------

    if ( mass > cutOff_ ) {
        totAmplitude = resAmplitude;
    } else {
        totAmplitude = bkgAmplitude + resAmplitude;
    }

    return totAmplitude;
}

const std::vector<LauParameter*>& LauLASSRes::getFloatingParameters()
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

void LauLASSRes::setResonanceParameter( const TString& name, const Double_t value )
{
    // Set various parameters for the LASS lineshape dynamics
    if ( name == "a" ) {
        this->setScatteringLength( value );
        std::cout << "INFO in LauLASSRes::setResonanceParameter : Setting LASS Scattering Length = "
                  << this->getScatteringLength() << std::endl;
    } else if ( name == "r" ) {
        this->setEffectiveRange( value );
        std::cout << "INFO in LauLASSRes::setResonanceParameter : Setting LASS Effective Range = "
                  << this->getEffectiveRange() << std::endl;
    } else {
        std::cerr << "WARNING in LauLASSRes::setResonanceParameter: Parameter name not reconised.  No parameter changes made."
                  << std::endl;
    }
}

void LauLASSRes::floatResonanceParameter( const TString& name )
{
    if ( name == "a" ) {
        if ( a_->fixed() ) {
            a_->fixed( kFALSE );
            this->addFloatingParameter( a_ );
        } else {
            std::cerr << "WARNING in LauLASSRes::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "r" ) {
        if ( r_->fixed() ) {
            r_->fixed( kFALSE );
            this->addFloatingParameter( r_ );
        } else {
            std::cerr << "WARNING in LauLASSRes::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else {
        std::cerr << "WARNING in LauLASSRes::fixResonanceParameter: Parameter name not reconised.  No parameter changes made."
                  << std::endl;
    }
}

LauParameter* LauLASSRes::getResonanceParameter( const TString& name )
{
    if ( name == "a" ) {
        return a_;
    } else if ( name == "r" ) {
        return r_;
    } else {
        std::cerr << "WARNING in LauLASSRes::getResonanceParameter: Parameter name not reconised."
                  << std::endl;
        return 0;
    }
}

void LauLASSRes::setEffectiveRange( const Double_t r )
{
    r_->value( r );
    r_->genValue( r );
    r_->initValue( r );
}

void LauLASSRes::setScatteringLength( const Double_t a )
{
    a_->value( a );
    a_->genValue( a );
    a_->initValue( a );
}
