
/*
Copyright 2008 University of Warwick

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

/*! \file LauParticlePDG.cc
    \brief File containing implementation of LauParticlePDG class.
*/

#include <iostream>
using std::cerr;
using std::endl;

#include "LauConstants.hh"
#include "LauDatabasePDG.hh"
#include "LauParticlePDG.hh"

#include "TMath.h"
#include "TParticlePDG.h"

LauParticlePDG::LauParticlePDG( const TParticlePDG* particle ) :
    theParticle_( particle ),
    antiParticle_( 0 ),
    mass_( 0.0 ),
    kaon_( kFALSE ),
    pion_( kFALSE ),
    ks_( kFALSE ),
    pi0_( kFALSE ),
    D_( kFALSE ),
    B_( kFALSE )
{
    // store the mass and kaon? and pion?
    if ( TMath::Abs( this->code() ) == 211 ) {    // charged pion
        mass_ = LauConstants::mPi;
        pion_ = kTRUE;
    } else if ( this->code() == 111 ) {    // neutral pion
        mass_ = LauConstants::mPi0;
        pion_ = kTRUE;
        pi0_  = kTRUE;
    } else if ( TMath::Abs( this->code() ) == 321 ) {    // charged kaon
        mass_ = LauConstants::mK;
        kaon_ = kTRUE;
    } else if ( this->code() == 130 || this->code() == 310 ||
                TMath::Abs( this->code() ) == 311 ) {    // neutral kaon
        mass_ = LauConstants::mK0;
        kaon_ = kTRUE;
        if ( this->code() == 310 ) {
            ks_ = kTRUE;
        }
    } else if ( this->code() == 221 ) {    // eta
        mass_ = LauConstants::mEta;
    } else if ( this->code() == 331 ) {    // eta'
        mass_ = LauConstants::mEtaPrime;
    } else if ( TMath::Abs( this->code() ) == 411 ) {    // charged D
        mass_ = LauConstants::mD;
        D_    = kTRUE;
    } else if ( TMath::Abs( this->code() ) == 421 ) {    // neutral D
        mass_ = LauConstants::mD0;
        D_    = kTRUE;
    } else if ( TMath::Abs( this->code() ) == 431 ) {    // charged D_s
        mass_ = LauConstants::mDs;
        D_    = kTRUE;
    } else if ( TMath::Abs( this->code() ) == 511 ) {    // neutral B
        mass_ = LauConstants::mB0;
        B_    = kTRUE;
    } else if ( TMath::Abs( this->code() ) == 521 ) {    // charged B
        mass_ = LauConstants::mB;
        B_    = kTRUE;
    } else if ( TMath::Abs( this->code() ) == 531 ) {    // neutral B_s
        mass_ = LauConstants::mBs0;
        B_    = kTRUE;
    }
}

LauParticlePDG::~LauParticlePDG()
{
}

Int_t LauParticlePDG::code() const
{
    if ( theParticle_ == 0 ) {
        cerr << "ERROR in LauParticlePDG::code : Null pointer." << endl;
        return -999;
    }
    return theParticle_->PdgCode();
}

Int_t LauParticlePDG::conjugateCode() const
{
    if ( antiParticle_ == 0 ) {
        cerr << "ERROR in LauParticlePDG::conjugateCode : Null pointer." << endl;
        return -999;
    }
    return antiParticle_->code();
}

Int_t LauParticlePDG::charge() const
{
    if ( theParticle_ == 0 ) {
        cerr << "ERROR in LauParticlePDG::charge : Null pointer." << endl;
        return -999;
    }
    Double_t chargeTimesThree = theParticle_->Charge();
    Int_t theCharge           = static_cast<Int_t>( chargeTimesThree / 3.0 );
    return theCharge;
}

TString LauParticlePDG::string() const
{
    if ( theParticle_ == 0 ) {
        cerr << "ERROR in LauParticlePDG::string : Null pointer." << endl;
        return "";
    }
    return theParticle_->GetName();
}

TString LauParticlePDG::stringAlphaNum() const
{
    TString name = this->string();
    name.ReplaceAll( "*", "st" );
    name.ReplaceAll( "_", "" );
    name.ReplaceAll( "(", "" );
    name.ReplaceAll( ")", "" );
    name.ReplaceAll( "+", "p" );
    name.ReplaceAll( "-", "m" );
    return name;
}
