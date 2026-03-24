
/*
Copyright 2006 University of Warwick

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

/*! \file LauBelleCPCoeffSet.cc
    \brief File containing implementation of LauBelleCPCoeffSet class.
*/

#include "LauBelleCPCoeffSet.hh"

#include "LauComplex.hh"
#include "LauConstants.hh"
#include "LauParameter.hh"
#include "LauPrint.hh"

#include "TMath.h"
#include "TRandom.h"

#include <fstream>
#include <iostream>
#include <vector>

LauBelleCPCoeffSet::LauBelleCPCoeffSet( const TString& compName,
                                        Double_t a,
                                        Double_t delta,
                                        Double_t b,
                                        Double_t phi,
                                        Bool_t aFixed,
                                        Bool_t deltaFixed,
                                        Bool_t bFixed,
                                        Bool_t phiFixed,
                                        Bool_t bSecondStage,
                                        Bool_t phiSecondStage ) :
    LauAbsCoeffSet( compName ),
    a_( new LauParameter( "A", a, minMagnitude_, maxMagnitude_, aFixed ) ),
    b_( new LauParameter( "B", b, minMagnitude_, maxMagnitude_, bFixed ) ),
    delta_( new LauParameter( "Delta", delta, minPhase_, maxPhase_, deltaFixed ) ),
    phi_( new LauParameter( "Phi", phi, minPhase_, maxPhase_, phiFixed ) ),
    particleCoeff_( 0.0, 0.0 ),
    antiparticleCoeff_( 0.0, 0.0 ),
    acp_( "ACP", ( -2.0 * b * TMath::Cos( phi ) ) / ( 1.0 + b * b ), -1.0, 1.0, bFixed && phiFixed )
{
    if ( bSecondStage && ! bFixed ) {
        b_->secondStage( kTRUE );
        b_->initValue( 0.0 );
    }
    if ( phiSecondStage && ! phiFixed ) {
        phi_->secondStage( kTRUE );
        phi_->initValue( 0.0 );
    }
}

LauBelleCPCoeffSet::LauBelleCPCoeffSet( const LauBelleCPCoeffSet& rhs,
                                        CloneOption cloneOption,
                                        Double_t constFactor ) :
    LauAbsCoeffSet( rhs.name() ),
    a_( 0 ),
    b_( 0 ),
    delta_( 0 ),
    phi_( 0 ),
    particleCoeff_( rhs.particleCoeff_ ),
    antiparticleCoeff_( rhs.antiparticleCoeff_ ),
    acp_( rhs.acp_ )
{
    if ( cloneOption == All || cloneOption == TieMagnitude ) {
        a_ = rhs.a_->createClone( constFactor );
    } else {
        a_ = new LauParameter( "A", rhs.a_->value(), minMagnitude_, maxMagnitude_, rhs.a_->fixed() );
        if ( rhs.a_->blind() ) {
            const LauBlind* blinder = rhs.a_->blinder();
            a_->blindParameter( blinder->blindingString(), blinder->blindingWidth() );
        }
    }

    if ( cloneOption == All || cloneOption == TieCPPars ) {
        b_ = rhs.b_->createClone( constFactor );
    } else {
        b_ = new LauParameter( "B", rhs.b_->value(), minMagnitude_, maxMagnitude_, rhs.b_->fixed() );
        if ( rhs.b_->blind() ) {
            const LauBlind* blinder = rhs.b_->blinder();
            b_->blindParameter( blinder->blindingString(), blinder->blindingWidth() );
        }
    }

    if ( cloneOption == All || cloneOption == TiePhase ) {
        delta_ = rhs.delta_->createClone( constFactor );
    } else {
        delta_ =
            new LauParameter( "Delta", rhs.delta_->value(), minPhase_, maxPhase_, rhs.delta_->fixed() );
        if ( rhs.delta_->blind() ) {
            const LauBlind* blinder = rhs.delta_->blinder();
            delta_->blindParameter( blinder->blindingString(), blinder->blindingWidth() );
        }
    }

    if ( cloneOption == All || cloneOption == TieCPPars ) {
        phi_ = rhs.phi_->createClone( constFactor );
    } else {
        phi_ = new LauParameter( "Phi", rhs.phi_->value(), minPhase_, maxPhase_, rhs.phi_->fixed() );
        if ( rhs.phi_->blind() ) {
            const LauBlind* blinder = rhs.phi_->blinder();
            phi_->blindParameter( blinder->blindingString(), blinder->blindingWidth() );
        }
    }
}

std::vector<LauParameter*> LauBelleCPCoeffSet::getParameters()
{
    std::vector<LauParameter*> pars;
    pars.push_back( a_ );
    pars.push_back( b_ );
    pars.push_back( delta_ );
    pars.push_back( phi_ );
    return pars;
}

void LauBelleCPCoeffSet::printParValues() const
{
    std::cout << "INFO in LauBelleCPCoeffSet::printParValues : Component \"" << this->name()
              << "\" has ";
    std::cout << "a-magnitude = " << a_->value() << ",\t";
    std::cout << "delta = " << delta_->value() << ",\t";
    std::cout << "b-magnitude = " << b_->value() << ",\t";
    std::cout << "phi = " << phi_->value() << "." << std::endl;
}

void LauBelleCPCoeffSet::printTableHeading( std::ostream& stream ) const
{
    stream << "\\begin{tabular}{|l|c|c|c|c|}" << std::endl;
    stream << "\\hline" << std::endl;
    stream << "Component & a-Magnitude & delta & b-Magnitude & phi \\\\" << std::endl;
    stream << "\\hline" << std::endl;
}

void LauBelleCPCoeffSet::printTableRow( std::ostream& stream ) const
{
    LauPrint print;
    TString resName = this->name();
    resName         = resName.ReplaceAll( "_", "\\_" );
    stream << resName << "  &  $";
    print.printFormat( stream, a_->value() );
    stream << " \\pm ";
    print.printFormat( stream, a_->error() );
    stream << "$  &  $";
    print.printFormat( stream, delta_->value() );
    stream << " \\pm ";
    print.printFormat( stream, delta_->error() );
    stream << "$  &  $";
    print.printFormat( stream, b_->value() );
    stream << " \\pm ";
    print.printFormat( stream, b_->error() );
    stream << "$  &  $";
    print.printFormat( stream, phi_->value() );
    stream << " \\pm ";
    print.printFormat( stream, phi_->error() );
    stream << "$ \\\\" << std::endl;
}

void LauBelleCPCoeffSet::randomiseInitValues()
{
    if ( a_->fixed() == kFALSE ) {
        // Choose an a-magnitude between 0.0 and 2.0
        Double_t mag = LauAbsCoeffSet::getRandomiser()->Rndm() * 2.0;
        a_->initValue( mag );
        a_->value( mag );
    }
    if ( b_->fixed() == kFALSE && b_->secondStage() == kFALSE ) {
        // Choose a b-magnitude between 0.0 and 0.1
        Double_t mag = LauAbsCoeffSet::getRandomiser()->Rndm() * 0.1;
        b_->initValue( mag );
        b_->value( mag );
    }
    if ( delta_->fixed() == kFALSE ) {
        // Choose a phase between +- pi
        Double_t phase = LauAbsCoeffSet::getRandomiser()->Rndm() * LauConstants::twoPi -
                         LauConstants::pi;
        delta_->initValue( phase );
        delta_->value( phase );
    }
    if ( phi_->fixed() == kFALSE && phi_->secondStage() == kFALSE ) {
        // Choose a phase between +- pi
        Double_t phase = LauAbsCoeffSet::getRandomiser()->Rndm() * LauConstants::twoPi -
                         LauConstants::pi;
        phi_->initValue( phase );
        phi_->value( phase );
    }
}

void LauBelleCPCoeffSet::finaliseValues()
{
    // retrieve the current values from the parameters
    Double_t aVal     = a_->value();
    Double_t bVal     = b_->value();
    Double_t deltaVal = delta_->value();
    Double_t phiVal   = phi_->value();
    Double_t genDelta = delta_->genValue();
    Double_t genPhi   = phi_->genValue();

    // Check whether we have a negative "a" magnitude.
    // If so make it positive and add pi to the "delta" phase.
    if ( aVal < 0.0 ) {
        aVal     *= -1.0;
        deltaVal += LauConstants::pi;
    }

    // Check whether we have a negative "b" magnitude.
    // If so make it positive and add pi to the "phi" phase.
    if ( bVal < 0.0 ) {
        bVal   *= -1.0;
        phiVal += LauConstants::pi;
    }

    // Check now whether the phases lies in the right range (-pi to pi).
    Bool_t deltaWithinRange( kFALSE );
    Bool_t phiWithinRange( kFALSE );
    while ( deltaWithinRange == kFALSE && phiWithinRange == kFALSE ) {
        if ( deltaVal > -LauConstants::pi && deltaVal < LauConstants::pi ) {
            deltaWithinRange = kTRUE;
        } else {
            // Not within the specified range
            if ( deltaVal > LauConstants::pi ) {
                deltaVal -= LauConstants::twoPi;
            } else if ( deltaVal < -LauConstants::pi ) {
                deltaVal += LauConstants::twoPi;
            }
        }

        if ( phiVal > -LauConstants::pi && phiVal < LauConstants::pi ) {
            phiWithinRange = kTRUE;
        } else {
            // Not within the specified range
            if ( phiVal > LauConstants::pi ) {
                phiVal -= LauConstants::twoPi;
            } else if ( phiVal < -LauConstants::pi ) {
                phiVal += LauConstants::twoPi;
            }
        }
    }

    // A further problem can occur when the generated phase is close to -pi or pi.
    // The phase can wrap over to the other end of the scale -
    // this leads to artificially large pulls so we wrap it back.
    Double_t diff = deltaVal - genDelta;
    if ( diff > LauConstants::pi ) {
        deltaVal -= LauConstants::twoPi;
    } else if ( diff < -LauConstants::pi ) {
        deltaVal += LauConstants::twoPi;
    }

    diff = phiVal - genPhi;
    if ( diff > LauConstants::pi ) {
        phiVal -= LauConstants::twoPi;
    } else if ( diff < -LauConstants::pi ) {
        phiVal += LauConstants::twoPi;
    }

    // finally store the new values in the parameters
    // and update the pulls
    a_->value( aVal );
    a_->updatePull();
    b_->value( bVal );
    b_->updatePull();
    delta_->value( deltaVal );
    delta_->updatePull();
    phi_->value( phiVal );
    phi_->updatePull();
}

const LauComplex& LauBelleCPCoeffSet::particleCoeff()
{
    LauComplex aTerm( a_->unblindValue() * TMath::Cos( delta_->unblindValue() ),
                      a_->unblindValue() * TMath::Sin( delta_->unblindValue() ) );
    LauComplex bTerm( b_->unblindValue() * TMath::Cos( phi_->unblindValue() ),
                      b_->unblindValue() * TMath::Sin( phi_->unblindValue() ) );
    particleCoeff_.setRealImagPart( 1.0, 0.0 );
    particleCoeff_ += bTerm;
    particleCoeff_ *= aTerm;
    return particleCoeff_;
}

const LauComplex& LauBelleCPCoeffSet::antiparticleCoeff()
{
    LauComplex aTerm( a_->unblindValue() * TMath::Cos( delta_->unblindValue() ),
                      a_->unblindValue() * TMath::Sin( delta_->unblindValue() ) );
    LauComplex bTerm( b_->unblindValue() * TMath::Cos( phi_->unblindValue() ),
                      b_->unblindValue() * TMath::Sin( phi_->unblindValue() ) );
    antiparticleCoeff_.setRealImagPart( 1.0, 0.0 );
    antiparticleCoeff_ -= bTerm;
    antiparticleCoeff_ *= aTerm;
    return antiparticleCoeff_;
}

void LauBelleCPCoeffSet::setCoeffValues( const LauComplex& coeff,
                                         const LauComplex& coeffBar,
                                         Bool_t init )
{
    LauComplex sum   = coeff + coeffBar;
    LauComplex diff  = coeff - coeffBar;
    LauComplex ratio = diff / sum;

    Double_t aVal( 0.5 * sum.abs() );
    Double_t deltaVal( sum.arg() );
    Double_t bVal( ratio.abs() );
    Double_t phiVal( ratio.arg() );

    a_->value( aVal );
    delta_->value( deltaVal );
    b_->value( bVal );
    phi_->value( phiVal );

    if ( init ) {
        a_->genValue( aVal );
        delta_->genValue( deltaVal );
        b_->genValue( bVal );
        phi_->genValue( phiVal );

        a_->initValue( aVal );
        delta_->initValue( deltaVal );
        b_->initValue( bVal );
        phi_->initValue( phiVal );
    }
}

LauParameter LauBelleCPCoeffSet::acp()
{
    // set the name
    TString parName( this->baseName() );
    parName += "_ACP";
    acp_.name( parName );

    // work out the ACP value
    Double_t value = ( -2.0 * b_->value() * TMath::Cos( phi_->value() ) ) /
                     ( 1.0 + b_->value() * b_->value() );

    // is it fixed?
    Bool_t fixed = b_->fixed() && phi_->fixed();
    acp_.fixed( fixed );

    // we can't work out the error without the covariance matrix
    Double_t error( 0.0 );

    // set the value and error
    acp_.valueAndErrors( value, error );

    return acp_;
}

LauAbsCoeffSet* LauBelleCPCoeffSet::createClone( const TString& newName,
                                                 CloneOption cloneOption,
                                                 Double_t constFactor )
{
    LauAbsCoeffSet* clone( 0 );
    if ( cloneOption == All || cloneOption == TiePhase || cloneOption == TieMagnitude ||
         cloneOption == TieCPPars ) {
        clone = new LauBelleCPCoeffSet( *this, cloneOption, constFactor );
        clone->name( newName );
    } else {
        std::cerr << "ERROR in LauBelleCPCoeffSet::createClone : Invalid clone option" << std::endl;
    }
    return clone;
}
