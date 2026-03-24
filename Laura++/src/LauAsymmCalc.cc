
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

/*! \file LauAsymmCalc.cc
    \brief File containing implementation of LauAsymmCalc class.
*/

#include "LauAsymmCalc.hh"

#include "TMath.h"

LauAsymmCalc::LauAsymmCalc( Double_t negValue, Double_t posValue ) :
    negValue_( negValue ),
    posValue_( posValue ),
    asymm_( 0.0 )
{
    asymm_ = calcAsymmetry();
}

LauAsymmCalc::LauAsymmCalc( const LauAsymmCalc& rhs ) :
    negValue_( rhs.negValue_ ),
    posValue_( rhs.posValue_ ),
    asymm_( rhs.asymm_ )
{
}

LauAsymmCalc& LauAsymmCalc::operator=( const LauAsymmCalc& rhs )
{
    if ( &rhs != this ) {
        negValue_ = rhs.negValue_;
        posValue_ = rhs.posValue_;
        asymm_    = rhs.asymm_;
    }
    return *this;
}

LauAsymmCalc::~LauAsymmCalc()
{
}

Double_t LauAsymmCalc::calcAsymmetry()
{
    Double_t num   = negValue_ - posValue_;
    Double_t denom = negValue_ + posValue_;
    Double_t asymm( 0.0 );
    if ( TMath::Abs( denom ) > 1e-10 ) {
        asymm = num / denom;
    }

    return asymm;
}
