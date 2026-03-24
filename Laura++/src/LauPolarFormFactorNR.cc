
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

/*! \file LauPolarFormFactorNR.cc
    \brief File containing implementation of LauPolarFormFactorNR class.
*/

#include "LauPolarFormFactorNR.hh"

#include "LauConstants.hh"
#include "LauDaughters.hh"
#include "LauParameter.hh"
#include "LauResonanceInfo.hh"

#include "TMath.h"

#include <iostream>

LauPolarFormFactorNR::LauPolarFormFactorNR( LauResonanceInfo* resInfo,
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

LauPolarFormFactorNR::~LauPolarFormFactorNR()
{
}

void LauPolarFormFactorNR::initialise()
{
    const LauDaughters* daughters = this->getDaughters();
    Int_t resPairAmpInt           = this->getPairInt();
    if ( daughters->gotSymmetricalDP() && resPairAmpInt != 3 ) {
        std::cerr << "WARNING in LauPolarFormFactorNR::initialise : Dalitz plot is symmetric - this lineshape is not appropriate."
                  << std::endl;
    }

    if ( model_ != LauAbsResonance::PolarFFNR ) {
        std::cerr << "WARNING in LauPolarFormFactorNR::initialise : Unknown model requested, defaulting to Polar Form Factor."
                  << std::endl;
        model_ = LauAbsResonance::PolarFFNR;
    }
}

LauComplex LauPolarFormFactorNR::resAmp( Double_t mass, Double_t )
{
    Double_t magnitude( 1.0 );

    Double_t lambda = this->getLambda();

    magnitude = 1.0 / ( 1.0 + mass * mass / ( lambda * lambda ) );

    LauComplex resAmplitude( magnitude, 0.0 );

    return resAmplitude;
}

const std::vector<LauParameter*>& LauPolarFormFactorNR::getFloatingParameters()
{
    this->clearFloatingParameters();

    if ( ! this->fixLambda() ) {
        this->addFloatingParameter( lambda_ );
    }

    return this->getParameters();
}

void LauPolarFormFactorNR::setResonanceParameter( const TString& name, const Double_t value )
{
    // Set various parameters for the lineshape
    if ( name == "lambda" ) {
        this->setLambda( value );
        std::cout << "INFO in LauPolarFormFactorNR::setResonanceParameter : Setting parameter lambda = "
                  << this->getLambda() << std::endl;
    } else {
        std::cerr << "WARNING in LauPolarFormFactorNR::setResonanceParameter: Parameter name not reconised.  No parameter changes made."
                  << std::endl;
    }
}

void LauPolarFormFactorNR::floatResonanceParameter( const TString& name )
{
    if ( name == "lambda" ) {
        if ( lambda_->fixed() ) {
            lambda_->fixed( kFALSE );
            this->addFloatingParameter( lambda_ );
        } else {
            std::cerr << "WARNING in LauPolarFormFactorNR::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else {
        std::cerr << "WARNING in LauPolarFormFactorNR::fixResonanceParameter: Parameter name not reconised.  No parameter changes made."
                  << std::endl;
    }
}

LauParameter* LauPolarFormFactorNR::getResonanceParameter( const TString& name )
{
    if ( name == "lambda" ) {
        return lambda_;
    } else {
        std::cerr << "WARNING in LauPolarFormFactorNR::getResonanceParameter: Parameter name not reconised."
                  << std::endl;
        return 0;
    }
}

void LauPolarFormFactorNR::setLambda( const Double_t lambda )
{
    lambda_->value( lambda );
    lambda_->genValue( lambda );
    lambda_->initValue( lambda );
}
