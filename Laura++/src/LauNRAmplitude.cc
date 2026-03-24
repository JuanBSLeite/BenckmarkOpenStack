
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

/*! \file LauNRAmplitude.cc
    \brief File containing implementation of LauNRAmplitude class.
*/

#include "LauNRAmplitude.hh"

#include "LauKinematics.hh"
#include "LauResonanceInfo.hh"

#include "TMath.h"

#include <iostream>

LauNRAmplitude::LauNRAmplitude( LauResonanceInfo* resInfo,
                                const Int_t resPairAmpInt,
                                const LauDaughters* daughters ) :
    LauAbsResonance( resInfo, resPairAmpInt, daughters ),
    d_( 0 ),
    c1_( 0 ),
    c2_( 0 ),
    p1_( 0 ),
    p2_( 0 )
{
    // Default values for parameters, taken from arXiv:0709.0075v1 [hep-ph]
    const Double_t dVal( 1.3232e-3 );
    const Double_t c1Val( 0.65 );
    const Double_t c2Val( 0.55 );
    const Double_t p1Val( 18.0 );
    const Double_t p2Val( 15.0 );

    const TString& parNameBase = this->getSanitisedName();

    TString dName( parNameBase );
    dName += "_d";
    d_     = resInfo->getExtraParameter( dName );
    if ( d_ == 0 ) {
        d_ = new LauParameter( dName, dVal, 0.0, 1.0, kTRUE );
        d_->secondStage( kTRUE );
        resInfo->addExtraParameter( d_ );
    }

    TString c1Name( parNameBase );
    c1Name += "_c1";
    c1_     = resInfo->getExtraParameter( c1Name );
    if ( c1_ == 0 ) {
        c1_ = new LauParameter( c1Name, c1Val, 0.0, 2.0, kTRUE );
        c1_->secondStage( kTRUE );
        resInfo->addExtraParameter( c1_ );
    }

    TString c2Name( parNameBase );
    c2Name += "_c2";
    c2_     = resInfo->getExtraParameter( c2Name );
    if ( c2_ == 0 ) {
        c2_ = new LauParameter( c2Name, c2Val, 0.0, 2.0, kTRUE );
        c2_->secondStage( kTRUE );
        resInfo->addExtraParameter( c2_ );
    }

    TString p1Name( parNameBase );
    p1Name += "_p1";
    p1_     = resInfo->getExtraParameter( p1Name );
    if ( p1_ == 0 ) {
        p1_ = new LauParameter( p1Name, p1Val, 0.0, 50.0, kTRUE );
        p1_->secondStage( kTRUE );
        resInfo->addExtraParameter( p1_ );
    }

    TString p2Name( parNameBase );
    p2Name += "_p2";
    p2_     = resInfo->getExtraParameter( p2Name );
    if ( p2_ == 0 ) {
        p2_ = new LauParameter( p2Name, p2Val, 0.0, 50.0, kTRUE );
        p2_->secondStage( kTRUE );
        resInfo->addExtraParameter( p2_ );
    }
}

LauNRAmplitude::~LauNRAmplitude()
{
}

void LauNRAmplitude::initialise()
{
}

LauComplex LauNRAmplitude::amplitude( const LauKinematics* kinematics )
{
    // Get the information from the kinematics object
    const Double_t m13Sq = kinematics->getm13Sq();
    const Double_t m23Sq = kinematics->getm23Sq();
    const Double_t m13   = kinematics->getm13();
    const Double_t m23   = kinematics->getm23();

    // Calculate the magnitude
    Double_t magnitude = TMath::Sqrt(
        m13 * m23 * this->f( m23Sq, c1_->unblindValue(), p1_->unblindValue() ) *
        this->f( m13Sq, c2_->unblindValue(), p2_->unblindValue() ) *
        TMath::Exp( -1.0 * d_->unblindValue() * m13Sq * m13Sq * m23Sq * m23Sq ) );

    // return the amplitude
    LauComplex resAmplitude( magnitude, 0.0 );
    return resAmplitude;
}

LauComplex LauNRAmplitude::resAmp( Double_t mass, Double_t spinTerm )
{
    std::cerr << "ERROR in LauNRAmplitude::resAmp : This method shouldn't get called." << std::endl;
    std::cerr << "                                  Returning zero amplitude for mass = " << mass
              << " and spinTerm = " << spinTerm << "." << std::endl;
    return LauComplex( 0.0, 0.0 );
}

Double_t LauNRAmplitude::f( const Double_t s, const Double_t c, const Double_t p ) const
{
    return 1.0 / ( 1.0 + TMath::Exp( c * ( s - p ) ) );
}

const std::vector<LauParameter*>& LauNRAmplitude::getFloatingParameters()
{
    this->clearFloatingParameters();

    if ( ! this->fixdParameter() ) {
        this->addFloatingParameter( d_ );
    }

    if ( ! this->fixc1Parameter() ) {
        this->addFloatingParameter( c1_ );
    }

    if ( ! this->fixc2Parameter() ) {
        this->addFloatingParameter( c2_ );
    }

    if ( ! this->fixp1Parameter() ) {
        this->addFloatingParameter( p1_ );
    }

    if ( ! this->fixp2Parameter() ) {
        this->addFloatingParameter( p2_ );
    }

    return this->getParameters();
}

void LauNRAmplitude::setResonanceParameter( const TString& name, const Double_t value )
{
    // Set various parameters for the NRAmplitude lineshape dynamics
    if ( name == "d" ) {
        this->setdParameter( value );
        std::cout << "INFO in LauNRAmplitude::setResonanceParameter : Setting NRAmplitude d = "
                  << this->getdParameter() << std::endl;
    } else if ( name == "c1" ) {
        this->setc1Parameter( value );
        std::cout << "INFO in LauNRAmplitude::setResonanceParameter : Setting NRAmplitude c1 = "
                  << this->getc1Parameter() << std::endl;
    } else if ( name == "c2" ) {
        this->setc2Parameter( value );
        std::cout << "INFO in LauNRAmplitude::setResonanceParameter : Setting NRAmplitude c2 = "
                  << this->getc2Parameter() << std::endl;
    } else if ( name == "p1" ) {
        this->setp1Parameter( value );
        std::cout << "INFO in LauNRAmplitude::setResonanceParameter : Setting NRAmplitude p1 = "
                  << this->getp1Parameter() << std::endl;
    } else if ( name == "p2" ) {
        this->setp2Parameter( value );
        std::cout << "INFO in LauNRAmplitude::setResonanceParameter : Setting NRAmplitude p2 = "
                  << this->getp2Parameter() << std::endl;
    } else {
        std::cerr << "WARNING in LauNRAmplitude::setResonanceParameter: Parameter name not reconised.  No parameter changes made."
                  << std::endl;
    }
}

void LauNRAmplitude::floatResonanceParameter( const TString& name )
{
    if ( name == "d" ) {
        if ( d_->fixed() ) {
            d_->fixed( kFALSE );
            this->addFloatingParameter( d_ );
        } else {
            std::cerr << "WARNING in LauNRAmplitude::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "c1" ) {
        if ( c1_->fixed() ) {
            c1_->fixed( kFALSE );
            this->addFloatingParameter( c1_ );
        } else {
            std::cerr << "WARNING in LauNRAmplitude::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "c2" ) {
        if ( c2_->fixed() ) {
            c2_->fixed( kFALSE );
            this->addFloatingParameter( c2_ );
        } else {
            std::cerr << "WARNING in LauNRAmplitude::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "p1" ) {
        if ( p1_->fixed() ) {
            p1_->fixed( kFALSE );
            this->addFloatingParameter( p1_ );
        } else {
            std::cerr << "WARNING in LauNRAmplitude::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "p2" ) {
        if ( p2_->fixed() ) {
            p2_->fixed( kFALSE );
            this->addFloatingParameter( p2_ );
        } else {
            std::cerr << "WARNING in LauNRAmplitude::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else {
        std::cerr << "WARNING in LauNRAmplitude::fixResonanceParameter: Parameter name not reconised.  No parameter changes made."
                  << std::endl;
    }
}

LauParameter* LauNRAmplitude::getResonanceParameter( const TString& name )
{
    if ( name == "d" ) {
        return d_;
    } else if ( name == "c1" ) {
        return c1_;
    } else if ( name == "c2" ) {
        return c2_;
    } else if ( name == "p1" ) {
        return p1_;
    } else if ( name == "p2" ) {
        return p2_;
    } else {
        std::cerr << "WARNING in LauNRAmplitude::getResonanceParameter: Parameter name not reconised."
                  << std::endl;
        return 0;
    }
}

void LauNRAmplitude::setdParameter( const Double_t d )
{
    d_->value( d );
    d_->genValue( d );
    d_->initValue( d );
}

void LauNRAmplitude::setc1Parameter( const Double_t c1 )
{
    c1_->value( c1 );
    c1_->genValue( c1 );
    c1_->initValue( c1 );
}

void LauNRAmplitude::setc2Parameter( const Double_t c2 )
{
    c2_->value( c2 );
    c2_->genValue( c2 );
    c2_->initValue( c2 );
}

void LauNRAmplitude::setp1Parameter( const Double_t p1 )
{
    p1_->value( p1 );
    p1_->genValue( p1 );
    p1_->initValue( p1 );
}

void LauNRAmplitude::setp2Parameter( const Double_t p2 )
{
    p2_->value( p2 );
    p2_->genValue( p2 );
    p2_->initValue( p2 );
}
