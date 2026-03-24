
/*
Copyright 2018 University of Warwick

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

/*! \file LauPolarFormFactorSymNR.cc
    \brief File containing implementation of LauPolarFormFactorSymNR class.
*/

#include "LauPolarFormFactorSymNR.hh"

#include "LauConstants.hh"
#include "LauDaughters.hh"
#include "LauParameter.hh"
#include "LauResonanceInfo.hh"

#include "TMath.h"

#include <iostream>

LauPolarFormFactorSymNR::LauPolarFormFactorSymNR( LauResonanceInfo* resInfo,
                                                  const LauAbsResonance::LauResonanceModel resType,
                                                  const Int_t resPairAmpInt,
                                                  const LauDaughters* daughters ) :
    LauAbsResonance( resInfo, resPairAmpInt, daughters ),
    lambda_( 0 ),
    model_( resType )

{
    TString parName  = this->getSanitisedName();
    parName         += "_lambda";
    lambda_          = resInfo->getExtraParameter( parName );
    if ( lambda_ == 0 ) {
        lambda_ = new LauParameter( parName, 1.0, 0.0, 10.0, kTRUE );
        lambda_->secondStage( kTRUE );
        resInfo->addExtraParameter( lambda_ );
    }
}

LauPolarFormFactorSymNR::~LauPolarFormFactorSymNR()
{
}

void LauPolarFormFactorSymNR::initialise()
{
    const LauDaughters* daughters = this->getDaughters();
    if ( ! daughters->gotSymmetricalDP() ) {
        std::cerr << "WARNING in LauPolarFormFactorSymNR::initialise : Dalitz plot is symmetric - this lineshape is not appropriate."
                  << std::endl;
    }

    Int_t resPairAmpInt = this->getPairInt();
    if ( resPairAmpInt == 3 ) {
        std::cerr << "WARNING in LauPolarFormFactorSymNR::initialise : This lineshape is intended to be on the symmetrised axes of the DP."
                  << std::endl;
    }

    if ( ( model_ != LauAbsResonance::PolarFFSymNR ) &&
         ( model_ != LauAbsResonance::PolarFFSymNRNoInter ) ) {
        std::cerr << "WARNING in LauPolarFormFactorSymNR::initialise : Unknown model requested, defaulting to Polar Form Factor."
                  << std::endl;
        model_ = LauAbsResonance::PolarFFSymNR;
    }

    if ( ( model_ != LauAbsResonance::PolarFFSymNR ) && ( this->getSpin() != 0 ) ) {
        std::cerr << "WARNING in LauPolarFormFactorSymNR::initialise : Non-zero spin will be ignored for this model - perhaps you should use LauAbsResonance::PolarFFSymNRNoInter instead"
                  << std::endl;
    }

    // NB we do not need to call setSpinType(LauAbsResonance::Legendre) here (as is done in LauPolarFormFactorNR) since override the amplitude method and explicitly use calcLegendrePoly
}

LauComplex LauPolarFormFactorSymNR::resAmp( Double_t mass, Double_t spinTerm )
{
    std::cerr << "ERROR in LauPolarFormFactorSymNR : This method should never be called."
              << std::endl;
    std::cerr << "                       : Returning zero amplitude for mass = " << mass
              << " and spinTerm = " << spinTerm << "." << std::endl;
    return LauComplex( 0.0, 0.0 );
}

LauComplex LauPolarFormFactorSymNR::amplitude( const LauKinematics* kinematics )
{
    // This function returns the complex dynamical amplitude for a Polar Form Factor Non-Resonant distribution

    // Calculate for symmetric DPs, e.g. 3pi or 3K, by using shapeNo = 1 or 2
    // Have s<->t symmetry already done in Dynamics flip function.
    // For Kpipi or similar plots, one can use the separate terms
    // and consider them as two separate components with their own mag and phase.
    // For this shapeNo = 3 and shapeNo = 4 need to be used to create the two
    // individual amplitudes (with the same value of lambda).

    // Calculate Mandelstam variables.
    // s = m_13^2, t = m_23^2
    const Double_t s = kinematics->getm13Sq();
    const Double_t t = kinematics->getm23Sq();

    Double_t magnitude( 1.0 );

    const Double_t lambda = this->getLambda();

    if ( model_ == LauAbsResonance::PolarFFSymNR ) {

        magnitude = 1.0 / ( 1.0 + s / ( lambda * lambda ) ) + 1.0 / ( 1.0 + t / ( lambda * lambda ) );

    } else if ( model_ == LauAbsResonance::PolarFFSymNRNoInter ) {

        magnitude = ( s <= t ) ? 1.0 / ( 1.0 + s / ( lambda * lambda ) )
                               : 1.0 / ( 1.0 + t / ( lambda * lambda ) );
    }

    LauComplex resAmplitude( magnitude, 0.0 );

    return resAmplitude;
}

const std::vector<LauParameter*>& LauPolarFormFactorSymNR::getFloatingParameters()
{
    this->clearFloatingParameters();

    if ( ! this->fixLambda() ) {
        this->addFloatingParameter( lambda_ );
    }

    return this->getParameters();
}

void LauPolarFormFactorSymNR::setResonanceParameter( const TString& name, const Double_t value )
{
    // Set various parameters for the lineshape
    if ( name == "lambda" ) {
        this->setLambda( value );
        std::cout << "INFO in LauPolarFormFactorSymNR::setResonanceParameter : Setting parameter lambda = "
                  << this->getLambda() << std::endl;
    } else {
        std::cerr << "WARNING in LauPolarFormFactorSymNR::setResonanceParameter: Parameter name not reconised.  No parameter changes made."
                  << std::endl;
    }
}

void LauPolarFormFactorSymNR::floatResonanceParameter( const TString& name )
{
    if ( name == "lambda" ) {
        if ( lambda_->fixed() ) {
            lambda_->fixed( kFALSE );
            this->addFloatingParameter( lambda_ );
        } else {
            std::cerr << "WARNING in LauPolarFormFactorSymNR::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else {
        std::cerr << "WARNING in LauPolarFormFactorSymNR::fixResonanceParameter: Parameter name not reconised.  No parameter changes made."
                  << std::endl;
    }
}

LauParameter* LauPolarFormFactorSymNR::getResonanceParameter( const TString& name )
{
    if ( name == "lambda" ) {
        return lambda_;
    } else {
        std::cerr << "WARNING in LauPolarFormFactorSymNR::getResonanceParameter: Parameter name not reconised."
                  << std::endl;
        return 0;
    }
}

void LauPolarFormFactorSymNR::setLambda( const Double_t lambda )
{
    lambda_->value( lambda );
    lambda_->genValue( lambda );
    lambda_->initValue( lambda );
}
