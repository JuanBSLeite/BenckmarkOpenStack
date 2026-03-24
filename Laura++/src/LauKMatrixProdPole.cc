
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

/*! \file LauKMatrixProdPole.cc
    \brief File containing implementation of LauKMatrixProdPole class.
*/

#include "LauKMatrixProdPole.hh"

#include "LauKMatrixPropagator.hh"
#include "LauResonanceMaker.hh"

#include "TMath.h"

#include <iostream>

LauKMatrixProdPole::LauKMatrixProdPole( const TString& poleName,
                                        Int_t poleIndex,
                                        Int_t resPairAmpInt,
                                        LauKMatrixPropagator* propagator,
                                        const LauDaughters* daughters,
                                        Bool_t useProdAdler ) :
    LauAbsResonance( poleName, resPairAmpInt, daughters, propagator->getL( propagator->getIndex() ) ),
    thePropagator_( propagator ),
    poleIndex_( poleIndex - 1 ),    // poleIndex goes from 1 to nPoles
    useProdAdler_( useProdAdler )
{

    if ( useProdAdler_ ) {
        std::cout << "Creating K matrix production pole " << poleName
                  << " with poleIndex = " << poleIndex
                  << " with s-dependent production Adler zero term" << std::endl;
    } else {
        std::cout << "Creating K matrix production pole " << poleName
                  << " with poleIndex = " << poleIndex << " with production Adler zero factor = 1"
                  << std::endl;
    }

    // `Resonance' Blatt-Weisskopf factor is handled internally, but parent must be set here. For other lineshapes, LauResonanceMaker handles this.
    this->setBarrierRadii(
        nullptr,
        LauResonanceMaker::get().getParentBWFactor( propagator->getL( propagator->getIndex() ),
                                                    LauBlattWeisskopfFactor::BWBarrier ) );
}

LauKMatrixProdPole::~LauKMatrixProdPole()
{
}

LauComplex LauKMatrixProdPole::resAmp( const Double_t mass, const Double_t spinTerm )
{
    // Calculate the amplitude for the K-matrix production pole.
    LauComplex amp( 0.0, 0.0 );

    if ( thePropagator_ == 0 ) {
        std::cerr << "ERROR in LauKMatrixProdPole::amplitude : The propagator is null" << std::endl;
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

    // Make sure the K-matrix propagator is up-to-date for
    // the given centre-of-mass squared value ("s")
    thePropagator_->updatePropagator( mass * mass );

    // Sum the pole denominator terms over all channels j, multiplying by
    // the propagator terms. Note that we do not sum over poles, since we
    // only want one of the production pole terms.
    Int_t nChannels = thePropagator_->getNChannels();
    Int_t jChannel;

    for ( jChannel = 0; jChannel < nChannels; jChannel++ ) {

        Double_t gj         = thePropagator_->getCouplingConstant( poleIndex_, jChannel );
        LauComplex prodTerm = thePropagator_->getPropTerm( jChannel );

        prodTerm.rescale( gj );

        amp += prodTerm;
    }

    Double_t poleDenom = thePropagator_->getPoleDenomTerm( poleIndex_ );

    // Include Adler zero factor if requested
    Double_t adlerZero( 1.0 );
    if ( useProdAdler_ ) {
        adlerZero = thePropagator_->getAdlerZero();
    }

    amp.rescale( poleDenom * adlerZero );

    // Scale by the spin term
    Double_t scale = spinTerm;

    // Include Blatt-Weisskopf barrier factor for parent
    scale *= fFactorB;

    amp.rescale( scale );

    return amp;
}

const std::vector<LauParameter*>& LauKMatrixProdPole::getFloatingParameters()
{

    this->clearFloatingParameters();

    Int_t nChannels = thePropagator_->getNChannels();

    for ( int jChannel = 0; jChannel < nChannels; jChannel++ ) {
        LauParameter& par_gj_ = thePropagator_->getCouplingParameter( poleIndex_, jChannel );
        if ( ! par_gj_.fixed() ) {
            this->addFloatingParameter( &par_gj_ );
        }
    }

    LauParameter& par_polemasssq_ = thePropagator_->getPoleMassSqParameter( poleIndex_ );
    if ( ! par_polemasssq_.fixed() ) {
        this->addFloatingParameter( &par_polemasssq_ );
    }

    return this->getParameters();
}
