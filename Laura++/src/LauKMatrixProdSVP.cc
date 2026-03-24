
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

/*! \file LauKMatrixProdSVP.cc
    \brief File containing implementation of LauKMatrixProdSVP class.
*/

#include "LauKMatrixProdSVP.hh"

#include "LauKMatrixPropagator.hh"

#include "TMath.h"

#include <iostream>

LauKMatrixProdSVP::LauKMatrixProdSVP( const TString& SVPName,
                                      Int_t channelIndex,
                                      Int_t resPairAmpInt,
                                      LauKMatrixPropagator* propagator,
                                      const LauDaughters* daughters,
                                      Bool_t useProdAdler ) :
    LauAbsResonance( SVPName, resPairAmpInt, daughters, propagator->getL( propagator->getIndex() ) ),
    thePropagator_( propagator ),
    channelIndex_( channelIndex - 1 ),    // channelIndex goes from 1 to nChannels.
    useProdAdler_( useProdAdler )
{
    // Constructor
    if ( useProdAdler_ ) {
        std::cout << "Creating K matrix production SVP " << SVPName
                  << " with channelIndex = " << channelIndex
                  << " with s-dependent production Adler zero term" << std::endl;
    } else {
        std::cout << "Creating K matrix production SVP " << SVPName
                  << " with channelIndex = " << channelIndex
                  << " with production Adler zero factor = 1" << std::endl;
    }
}

LauKMatrixProdSVP::~LauKMatrixProdSVP()
{
    // Destructor
}

LauComplex LauKMatrixProdSVP::resAmp( const Double_t mass, const Double_t spinTerm )
{

    // Calculate the amplitude for the K-matrix production pole.
    // First, make sure the K-matrix propagator is up-to-date for
    // the given centre-of-mass squared value ("s") from the kinematics.
    LauComplex amp( 0.0, 0.0 );

    if ( thePropagator_ == 0 ) {
        std::cerr << "ERROR in LauKMatrixProdSVP::amplitude : The propagator is null" << std::endl;
        return amp;
    }

    // Get barrier factors ('resonance' factor is already accounted for internally via propagator 'Gamma' matrix)
    Double_t fFactorB( 1.0 );

    const Int_t resSpin  = this->getSpin();
    const Double_t pstar = this->getPstar();

    if ( resSpin > 0 ) {
        const LauBlattWeisskopfFactor* parBWFactor = this->getParBWFactor();
        if ( parBWFactor != nullptr ) {
            switch ( parBWFactor->getRestFrame() ) {
                case LauBlattWeisskopfFactor::ResonanceFrame :
                    fFactorB = parBWFactor->calcFormFactor( this->getP() );
                    break;
                case LauBlattWeisskopfFactor::ParentFrame :
                    fFactorB = parBWFactor->calcFormFactor( pstar );
                    break;
                case LauBlattWeisskopfFactor::Covariant : {
                    Double_t covFactor = this->getCovFactor();
                    if ( resSpin > 2 ) {
                        covFactor = TMath::Power( covFactor, 1.0 / resSpin );
                    } else if ( resSpin == 2 ) {
                        covFactor = TMath::Sqrt( covFactor );
                    }
                    fFactorB = parBWFactor->calcFormFactor( pstar * covFactor );
                    break;
                }
            }
        }
    }

    thePropagator_->updatePropagator( mass * mass );

    Double_t SVPTerm = thePropagator_->getProdSVPTerm();

    amp = thePropagator_->getPropTerm( channelIndex_ );

    // Include Adler zero factor if requested
    Double_t adlerZero( 1.0 );
    if ( useProdAdler_ ) {
        adlerZero = thePropagator_->getAdlerZero();
    }

    amp.rescale( SVPTerm * adlerZero );

    // Scale by the spin term
    Double_t scale = spinTerm;

    // Include Blatt-Weisskopf barrier factor for parent
    scale *= fFactorB;

    amp.rescale( scale );

    return amp;
}

const std::vector<LauParameter*>& LauKMatrixProdSVP::getFloatingParameters()
{

    this->clearFloatingParameters();

    Int_t nChannels = thePropagator_->getNChannels();

    for ( int jChannel = 0; jChannel < nChannels; jChannel++ ) {
        LauParameter& par_f_ = thePropagator_->getScatteringParameter( channelIndex_, jChannel );
        if ( ! par_f_.fixed() ) {
            this->addFloatingParameter( &par_f_ );
        }
    }

    LauParameter& par_mSq0_ = thePropagator_->getmSq0();
    if ( ! par_mSq0_.fixed() ) {
        this->addFloatingParameter( &par_mSq0_ );
    }

    LauParameter& par_s0Scatt_ = thePropagator_->gets0Scatt();
    if ( ! par_s0Scatt_.fixed() ) {
        this->addFloatingParameter( &par_s0Scatt_ );
    }

    LauParameter& par_s0Prod_ = thePropagator_->gets0Prod();
    if ( ! par_s0Prod_.fixed() ) {
        this->addFloatingParameter( &par_s0Prod_ );
    }

    LauParameter& par_sA_ = thePropagator_->getsA();
    if ( ! par_sA_.fixed() ) {
        this->addFloatingParameter( &par_sA_ );
    }

    LauParameter& par_sA0_ = thePropagator_->getsA0();
    if ( ! par_sA0_.fixed() ) {
        this->addFloatingParameter( &par_sA0_ );
    }

    return this->getParameters();
}
