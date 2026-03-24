
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

/*! \file LauFlatNR.cc
    \brief File containing implementation of LauFlatNR class.
*/

#include "LauFlatNR.hh"

#include "LauKinematics.hh"

#include <iostream>

LauFlatNR::LauFlatNR( LauResonanceInfo* resInfo,
                      const Int_t resPairAmpInt,
                      const LauDaughters* daughters ) :
    LauAbsResonance( resInfo, resPairAmpInt, daughters )
{
}

LauFlatNR::~LauFlatNR()
{
}

void LauFlatNR::initialise()
{
}

LauComplex LauFlatNR::amplitude( const LauKinematics* /*kinematics*/ )
{
    // return the uniform amplitude
    return LauComplex( 1.0, 0.0 );
}

LauComplex LauFlatNR::resAmp( Double_t mass, Double_t spinTerm )
{
    std::cerr << "ERROR in LauFlatNR::resAmp : This method shouldn't get called." << std::endl;
    std::cerr << "                             Returning zero amplitude for mass = " << mass
              << " and spinTerm = " << spinTerm << "." << std::endl;
    return LauComplex( 0.0, 0.0 );
}
