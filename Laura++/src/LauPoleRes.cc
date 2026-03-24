
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

/*! \file LauPoleRes.cc
    \brief File containing implementation of LauPoleRes class.
*/

#include "LauPoleRes.hh"

LauPoleRes::LauPoleRes( LauResonanceInfo* resInfo,
                        const Int_t resPairAmpInt,
                        const LauDaughters* daughters ) :
    LauAbsResonance( resInfo, resPairAmpInt, daughters )
{
}

LauPoleRes::~LauPoleRes()
{
}

void LauPoleRes::initialise()
{
}

LauComplex LauPoleRes::resAmp( Double_t mass, Double_t spinTerm )
{
    // This function returns the complex dynamical amplitude for a pole. 1/((mpole*mpole)-m*m)
    // mpole==(m0-iw0)==(rePole-i*imPole), to use the already defined mass and width of the resonance

    Double_t rePole = this->getMass();
    Double_t imPole = this->getWidth();
    Double_t reTerm = rePole * rePole - imPole * imPole - mass * mass;
    Double_t imTerm = 2.0 * rePole * imPole;

    LauComplex resAmplitude( reTerm, imTerm );

    resAmplitude.rescale( 1. / ( reTerm * reTerm + imTerm * imTerm ) );

    return resAmplitude.scale( spinTerm );
}

const std::vector<LauParameter*>& LauPoleRes::getFloatingParameters()
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
