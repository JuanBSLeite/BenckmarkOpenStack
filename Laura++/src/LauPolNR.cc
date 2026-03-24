
/*
Copyright 2014 University of Warwick

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

/*! \file LauPolNR.cc
    \brief File containing implementation of LauPolNR class.
*/

#include "LauPolNR.hh"

#include "LauDaughters.hh"

#include "TMath.h"

#include <iostream>

LauPolNR::LauPolNR( LauResonanceInfo* resInfo,
                    const Int_t resPairAmpInt,
                    const LauDaughters* daughters ) :
    LauAbsResonance( resInfo, resPairAmpInt, daughters ),
    omega_( 0.0 ),
    order_( 0 )
{
}

LauPolNR::~LauPolNR()
{
}

void LauPolNR::initialise()
{
    const LauDaughters* daughters = this->getDaughters();
    Int_t resPairAmpInt           = this->getPairInt();
    if ( daughters->gotSymmetricalDP() && resPairAmpInt != 3 ) {
        std::cerr << "WARNING in LauPolNR::initialise : Dalitz plot is symmetric - this lineshape is not appropriate."
                  << std::endl;
    }

    // Make the spin term purely the Legendre polynomial of the cos(helicity angle)
    this->setSpinType( LauAbsResonance::Legendre );

    Double_t omega = 0.5 * ( daughters->getMassParent() +
                             ( 1.0 / 3.0 ) * ( daughters->getMassDaug1() + daughters->getMassDaug2() +
                                               daughters->getMassDaug3() ) );
    this->setOmega( omega );

    TString name = this->getResonanceName();
    if ( name.EndsWith( "0", TString::kExact ) ) {
        order_ = 0;
    } else if ( name.EndsWith( "1", TString::kExact ) ) {
        order_ = 1;
    } else if ( name.EndsWith( "2", TString::kExact ) ) {
        order_ = 2;
    } else if ( name.EndsWith( "3", TString::kExact ) ) {
        order_ = 3;
    } else {
        std::cerr << "ERROR in LauPolNR::initialise : Resonance order unknown (known orders : 0, 1, 2 and 3)"
                  << std::endl;
        std::cerr << "                              : Defaulting to 0." << std::endl;
    }
}

LauComplex LauPolNR::resAmp( Double_t mass, Double_t spinTerm )
{
    Double_t magnitude( spinTerm );

    Double_t massDiff = mass - omega_;
    for ( UInt_t i( 0 ); i < order_; ++i ) {
        magnitude *= massDiff;
    }

    LauComplex resAmplitude( magnitude, 0.0 );

    return resAmplitude;
}
