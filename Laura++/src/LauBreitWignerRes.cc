
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

/*! \file LauBreitWignerRes.cc
    \brief File containing implementation of LauBreitWignerRes class.
*/

#include "LauBreitWignerRes.hh"

LauBreitWignerRes::LauBreitWignerRes( LauResonanceInfo* resInfo,
                                      const Int_t resPairAmpInt,
                                      const LauDaughters* daughters ) :
    LauAbsResonance( resInfo, resPairAmpInt, daughters )
{
}

LauBreitWignerRes::~LauBreitWignerRes()
{
}

void LauBreitWignerRes::initialise()
{
}

LauComplex LauBreitWignerRes::resAmp( Double_t mass, Double_t spinTerm )
{
    // This function returns the complex dynamical amplitude for a Breit-Wigner resonance,
    // given the invariant mass and cos(helicity) values.

    Double_t resMass = this->getMass();
    Double_t deltaM  = mass - resMass;

    Double_t resWidthO2   = 0.5 * this->getWidth();
    Double_t resWidthO2Sq = resWidthO2 * resWidthO2;

    LauComplex resAmplitude( deltaM, resWidthO2 );

    resAmplitude.rescale( spinTerm / ( deltaM * deltaM + resWidthO2Sq ) );

    return resAmplitude;
}

const std::vector<LauParameter*>& LauBreitWignerRes::getFloatingParameters()
{
    this->clearFloatingParameters();

    if ( ! this->fixMass() ) {
        this->addFloatingParameter( this->getMassPar() );
    }
    if ( ! this->fixWidth() ) {
        this->addFloatingParameter( this->getWidthPar() );
    }

    return this->getParameters();
}
