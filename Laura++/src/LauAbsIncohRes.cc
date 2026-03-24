
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

/*! \file LauAbsIncohRes.cc
    \brief File containing implementation of LauAbsIncohRes class.
*/

#include "LauAbsIncohRes.hh"

#include "LauDaughters.hh"
#include "LauResonanceInfo.hh"

// Constructor
LauAbsIncohRes::LauAbsIncohRes( LauResonanceInfo* resInfo,
                                const Int_t resPairAmpInt,
                                const LauDaughters* daughters ) :
    LauAbsResonance( resInfo, resPairAmpInt, daughters )
{
}

// Destructor
LauAbsIncohRes::~LauAbsIncohRes()
{
}

LauComplex LauAbsIncohRes::resAmp( Double_t /*mass*/, Double_t spinTerm )
{
    return LauComplex( spinTerm, 0.0 );
}
