
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

/*! \file LauBelleNR.cc
    \brief File containing implementation of LauBelleNR class.
*/

#include "LauBelleNR.hh"

#include "LauDaughters.hh"
#include "LauParameter.hh"
#include "LauResonanceInfo.hh"

#include "TMath.h"

#include <iostream>

LauBelleNR::LauBelleNR( LauResonanceInfo* resInfo,
                        const LauAbsResonance::LauResonanceModel resType,
                        const Int_t resPairAmpInt,
                        const LauDaughters* daughters ) :
    LauAbsResonance( resInfo, resPairAmpInt, daughters ),
    alpha_( 0 ),
    model_( resType ),
    forceLegendre_( kTRUE )
{
    TString parName  = this->getSanitisedName();
    parName         += "_alpha";
    alpha_           = resInfo->getExtraParameter( parName );
    if ( alpha_ == 0 ) {
        alpha_ = new LauParameter( parName, 0.0, -2.0, 10.0, kTRUE );
        alpha_->secondStage( kTRUE );
        resInfo->addExtraParameter( alpha_ );
    }
}

LauBelleNR::~LauBelleNR()
{
}

void LauBelleNR::initialise()
{
    const LauDaughters* daughters = this->getDaughters();
    Int_t resPairAmpInt           = this->getPairInt();
    if ( daughters->gotSymmetricalDP() && resPairAmpInt != 3 ) {
        std::cerr << "WARNING in LauBelleNR::initialise : Dalitz plot is symmetric - this lineshape is not appropriate."
                  << std::endl;
    }

    if ( model_ != LauAbsResonance::BelleNR && model_ != LauAbsResonance::PowerLawNR ) {
        std::cerr << "WARNING in LauBelleNR::initialise : Unknown model requested, defaulting to exponential."
                  << std::endl;
        model_ = LauAbsResonance::BelleNR;
    }

    // Make the spin term purely the Legendre polynomial of the cos(helicity angle)
    if ( forceLegendre_ ) {
        this->setSpinType( LauAbsResonance::Legendre );
    }
}

LauComplex LauBelleNR::resAmp( Double_t mass, Double_t spinTerm )
{
    Double_t magnitude( 1.0 );

    Double_t alpha = this->getAlpha();

    if ( model_ == LauAbsResonance::BelleNR ) {
        magnitude = spinTerm * TMath::Exp( -alpha * mass * mass );
    } else if ( model_ == LauAbsResonance::PowerLawNR ) {
        magnitude = spinTerm * TMath::Power( mass * mass, -alpha );
    }

    LauComplex resAmplitude( magnitude, 0.0 );

    return resAmplitude;
}

const std::vector<LauParameter*>& LauBelleNR::getFloatingParameters()
{
    this->clearFloatingParameters();

    if ( ! this->fixAlpha() ) {
        this->addFloatingParameter( alpha_ );
    }

    return this->getParameters();
}

void LauBelleNR::setResonanceParameter( const TString& name, const Double_t value )
{
    // Set various parameters for the lineshape
    if ( name == "alpha" ) {
        this->setAlpha( value );
        std::cout << "INFO in LauBelleNR::setResonanceParameter : Setting parameter alpha = "
                  << this->getAlpha() << std::endl;
    } else {
        std::cerr << "WARNING in LauBelleNR::setResonanceParameter: Parameter name not reconised.  No parameter changes made."
                  << std::endl;
    }
}

void LauBelleNR::floatResonanceParameter( const TString& name )
{
    if ( name == "alpha" ) {
        if ( alpha_->fixed() ) {
            alpha_->fixed( kFALSE );
            this->addFloatingParameter( alpha_ );
        } else {
            std::cerr << "WARNING in LauBelleNR::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else {
        std::cerr << "WARNING in LauBelleNR::fixResonanceParameter: Parameter name not reconised.  No parameter changes made."
                  << std::endl;
    }
}

LauParameter* LauBelleNR::getResonanceParameter( const TString& name )
{
    if ( name == "alpha" ) {
        return alpha_;
    } else {
        std::cerr << "WARNING in LauBelleNR::getResonanceParameter: Parameter name not reconised."
                  << std::endl;
        return 0;
    }
}

void LauBelleNR::setAlpha( const Double_t alpha )
{
    alpha_->value( alpha );
    alpha_->genValue( alpha );
    alpha_->initValue( alpha );
}
