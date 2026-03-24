
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

/*! \file LauGaussIncohRes.cc
    \brief File containing implementation of LauGaussIncohRes class.
*/

#include "LauGaussIncohRes.hh"

#include "LauConstants.hh"
#include "LauKinematics.hh"

#include "TMath.h"

#include <iostream>

LauGaussIncohRes::LauGaussIncohRes( LauResonanceInfo* resInfo,
                                    const Int_t resPairAmpInt,
                                    const LauDaughters* daughters ) :
    LauAbsIncohRes( resInfo, resPairAmpInt, daughters )
{
}

LauGaussIncohRes::~LauGaussIncohRes()
{
}

void LauGaussIncohRes::initialise()
{
}

Double_t LauGaussIncohRes::intensityFactor( const LauKinematics* kinematics )
{
    Int_t pair = this->getPairInt();
    Double_t mass( 0. );

    switch ( pair ) {
        case 1 :
            mass = kinematics->getm23();
            break;
        case 2 :
            mass = kinematics->getm13();
            break;
        case 3 :
            mass = kinematics->getm12();
            break;
        default :
            std::cerr << "ERROR in LauGaussIncohRes::intensityFactor : Incorrect daughter pair defined."
                      << std::endl;
            std::cerr << "                             Returning zero intensity." << std::endl;
            return 0.;
    }

    Double_t m0 = this->getMass();
    Double_t G0 = this->getWidth();

    if ( G0 == 0 ) {
        std::cerr << "ERROR in LauGaussIncohRes::intensityFactor : Gaussian width cannot be zero."
                  << std::endl;
        std::cerr << "                             Returning zero intensity." << std::endl;
        return 0.;
    }

    Double_t exponent  = -0.5 * ( mass - m0 ) * ( mass - m0 ) / ( G0 * G0 );
    Double_t intensity = TMath::Exp( exponent );

    return intensity;
}

const std::vector<LauParameter*>& LauGaussIncohRes::getFloatingParameters()
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
