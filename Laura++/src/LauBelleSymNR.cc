
/*
Copyright 2013 University of Warwick

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

/*! \file LauBelleSymNR.cc
    \brief File containing implementation of LauBelleSymNR class.
*/

#include "LauBelleSymNR.hh"

#include "LauDaughters.hh"
#include "LauKinematics.hh"
#include "LauResonanceInfo.hh"

#include "TMath.h"

#include <iostream>

LauBelleSymNR::LauBelleSymNR( LauResonanceInfo* resInfo,
                              const LauAbsResonance::LauResonanceModel resType,
                              const Int_t resPairAmpInt,
                              const LauDaughters* daughters ) :
    LauAbsResonance( resInfo, resPairAmpInt, daughters ),
    alpha_( 0 ),
    model_( resType )
{
    TString parName  = this->getSanitisedName();
    parName         += "_alpha";
    alpha_           = resInfo->getExtraParameter( parName );
    if ( alpha_ == 0 ) {
        alpha_ = new LauParameter( parName, 0.0, 0.0, 10.0, kTRUE );
        alpha_->secondStage( kTRUE );
        resInfo->addExtraParameter( alpha_ );
    }
}

LauBelleSymNR::~LauBelleSymNR()
{
}

void LauBelleSymNR::initialise()
{
    const LauDaughters* daughters = this->getDaughters();
    if ( ! daughters->gotSymmetricalDP() ) {
        std::cerr << "WARNING in LauBelleSymNR::initialise : Dalitz plot is symmetric - this lineshape is not appropriate."
                  << std::endl;
    }

    Int_t resPairAmpInt = this->getPairInt();
    if ( resPairAmpInt == 3 ) {
        std::cerr << "WARNING in LauBelleSymNR::initialise : This lineshape is intended to be on the symmetrised axes of the DP."
                  << std::endl;
    }

    if ( ( model_ != LauAbsResonance::BelleSymNR ) &&
         ( model_ != LauAbsResonance::BelleSymNRNoInter ) &&
         ( model_ != LauAbsResonance::TaylorNR ) ) {
        std::cerr << "WARNING in LauBelleSymNR::initialise : Unknown model requested, defaulting to exponential."
                  << std::endl;
        model_ = LauAbsResonance::BelleSymNR;
    }

    if ( ( model_ != LauAbsResonance::BelleSymNRNoInter ) && ( this->getSpin() != 0 ) ) {
        std::cerr << "WARNING in LauBelleSymNR::initialise : Non-zero spin will be ignored for this model - perhaps you should use LauAbsResonance::BelleSymNRNoInter instead"
                  << std::endl;
    }

    // NB we do not need to call setSpinType(LauAbsResonance::Legendre) here (as is done in LauBelleNR) since override the amplitude method and explicitly use calcLegendrePoly
}

LauComplex LauBelleSymNR::amplitude( const LauKinematics* kinematics )
{
    // This function returns the complex dynamical amplitude for a Belle Non-Resonant distribution

    // Calculate for symmetric DPs, e.g. 3pi or 3K, by using shapeNo = 1 or 2
    // Have s<->t symmetry already done in Dynamics flip function.
    // For Kpipi or similar plots, one can use the separate exponentials
    // and consider them as two separate components with their own mag and phase.
    // For this shapeNo = 3 and shapeNo = 4 need to be used to create the two
    // individual amplitudes (with the same value of alpha).

    // Calculate Mandelstam variables.
    // s = m_13^2, t = m_23^2
    const Double_t s = kinematics->getm13Sq();
    const Double_t t = kinematics->getm23Sq();

    Double_t magnitude( 1.0 );

    const Double_t alpha = this->getAlpha();

    if ( model_ == LauAbsResonance::BelleSymNR ) {

        magnitude = TMath::Exp( -alpha * s ) + TMath::Exp( -alpha * t );

    } else if ( model_ == LauAbsResonance::BelleSymNRNoInter ) {

        magnitude = ( s <= t ) ? TMath::Exp( -alpha * s ) : TMath::Exp( -alpha * t );

        const Double_t cosHel  = ( s <= t ) ? kinematics->getc13() : kinematics->getc23();
        magnitude             *= this->calcLegendrePoly( cosHel );

    } else if ( model_ == LauAbsResonance::TaylorNR ) {

        const Double_t mParentSq = kinematics->getmParentSq();
        magnitude                = alpha * ( s + t ) / mParentSq + 1.0;
    }

    LauComplex resAmplitude( magnitude, 0.0 );

    return resAmplitude;
}

LauComplex LauBelleSymNR::resAmp( Double_t mass, Double_t spinTerm )
{
    std::cerr << "ERROR in LauBelleSymNR : This method should never be called." << std::endl;
    std::cerr << "                       : Returning zero amplitude for mass = " << mass
              << " and spinTerm = " << spinTerm << "." << std::endl;
    return LauComplex( 0.0, 0.0 );
}

const std::vector<LauParameter*>& LauBelleSymNR::getFloatingParameters()
{
    this->clearFloatingParameters();

    if ( ! this->fixAlpha() ) {
        this->addFloatingParameter( alpha_ );
    }

    return this->getParameters();
}

void LauBelleSymNR::setResonanceParameter( const TString& name, const Double_t value )
{
    // Set various parameters for the lineshape
    if ( name == "alpha" ) {
        this->setAlpha( value );
        std::cout << "INFO in LauBelleSymNR::setResonanceParameter : Setting parameter alpha = "
                  << this->getAlpha() << std::endl;
    } else {
        std::cerr << "WARNING in LauBelleSymNR::setResonanceParameter : Parameter name not reconised.  No parameter changes made."
                  << std::endl;
    }
}

void LauBelleSymNR::floatResonanceParameter( const TString& name )
{
    if ( name == "alpha" ) {
        if ( alpha_->fixed() ) {
            alpha_->fixed( kFALSE );
            this->addFloatingParameter( alpha_ );
        } else {
            std::cerr << "WARNING in LauBelleSymNR::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else {
        std::cerr << "WARNING in LauBelleSymNR::fixResonanceParameter: Parameter name not reconised.  No parameter changes made."
                  << std::endl;
    }
}

LauParameter* LauBelleSymNR::getResonanceParameter( const TString& name )
{
    if ( name == "alpha" ) {
        return alpha_;
    } else {
        std::cerr << "WARNING in LauBelleSymNR::getResonanceParameter: Parameter name not reconised."
                  << std::endl;
        return 0;
    }
}

void LauBelleSymNR::setAlpha( const Double_t alpha )
{
    alpha_->value( alpha );
    alpha_->genValue( alpha );
    alpha_->initValue( alpha );
}
