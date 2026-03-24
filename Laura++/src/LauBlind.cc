
/*
Copyright 2015 University of Warwick

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

/*! \file LauBlind.cc
    \brief File containing implementation of LauBlind class.
*/

#include "LauBlind.hh"

#include "TMath.h"
#include "TRandom3.h"

LauBlind::LauBlind() :
    blindingString_( "" ),
    blindingWidth_( 0.0 ),
    offset_( 0.0 )
{
}

LauBlind::LauBlind( const TString& blindingStr, const Double_t width ) :
    blindingString_( blindingStr ),
    blindingWidth_( width ),
    offset_( 0.0 )
{
    this->calcOffset();
}

LauBlind::LauBlind( const LauBlind& rhs ) :
    blindingString_( rhs.blindingString_ ),
    blindingWidth_( rhs.blindingWidth_ ),
    offset_( rhs.offset_ )
{
}

LauBlind::~LauBlind()
{
}

void LauBlind::calcOffset()
{
    //hash the blinding string to obtain a seed
    TRandom3 r( TMath::Hash( blindingString_ ) );

    //offsets are Gaussian distributed with defined width
    offset_ = blindingWidth_ * r.Gaus();
}
