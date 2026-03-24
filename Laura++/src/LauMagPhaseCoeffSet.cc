
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

/*! \file LauMagPhaseCoeffSet.cc
    \brief File containing implementation of LauMagPhaseCoeffSet class.
*/

#include "LauMagPhaseCoeffSet.hh"

#include "LauComplex.hh"
#include "LauConstants.hh"
#include "LauParameter.hh"
#include "LauPrint.hh"

#include "TMath.h"
#include "TRandom.h"

#include <fstream>
#include <iostream>
#include <vector>

LauMagPhaseCoeffSet::LauMagPhaseCoeffSet( const TString& compName,
                                          Double_t magnitude,
                                          Double_t phase,
                                          Bool_t magFixed,
                                          Bool_t phaseFixed ) :
    LauAbsCoeffSet( compName ),
    magnitude_( new LauParameter( "A", magnitude, minMagnitude_, maxMagnitude_, magFixed ) ),
    phase_( new LauParameter( "Delta", phase, minPhase_, maxPhase_, phaseFixed ) ),
    coeff_( magnitude * TMath::Cos( phase ), magnitude * TMath::Sin( phase ) )
{
}

LauMagPhaseCoeffSet::LauMagPhaseCoeffSet( const LauMagPhaseCoeffSet& rhs,
                                          CloneOption cloneOption,
                                          Double_t constFactor ) :
    LauAbsCoeffSet( rhs.name() ),
    magnitude_( 0 ),
    phase_( 0 ),
    coeff_( rhs.coeff_ )
{
    if ( cloneOption == All || cloneOption == TieMagnitude ) {
        magnitude_ = rhs.magnitude_->createClone( constFactor );
    } else {
        magnitude_ = new LauParameter( "A",
                                       rhs.magnitude_->value(),
                                       minMagnitude_,
                                       maxMagnitude_,
                                       rhs.magnitude_->fixed() );
        if ( rhs.magnitude_->blind() ) {
            const LauBlind* blinder = rhs.magnitude_->blinder();
            magnitude_->blindParameter( blinder->blindingString(), blinder->blindingWidth() );
        }
    }

    if ( cloneOption == All || cloneOption == TiePhase ) {
        phase_ = rhs.phase_->createClone( constFactor );
    } else {
        phase_ =
            new LauParameter( "Delta", rhs.phase_->value(), minPhase_, maxPhase_, rhs.phase_->fixed() );
        if ( rhs.phase_->blind() ) {
            const LauBlind* blinder = rhs.phase_->blinder();
            phase_->blindParameter( blinder->blindingString(), blinder->blindingWidth() );
        }
    }
}

std::vector<LauParameter*> LauMagPhaseCoeffSet::getParameters()
{
    std::vector<LauParameter*> pars;
    pars.push_back( magnitude_ );
    pars.push_back( phase_ );
    return pars;
}

void LauMagPhaseCoeffSet::printParValues() const
{
    std::cout << "INFO in LauMagPhaseCoeffSet::printParValues : Component \"" << this->name()
              << "\" has magnitude = " << magnitude_->value() << " and phase = " << phase_->value()
              << "." << std::endl;
}

void LauMagPhaseCoeffSet::printTableHeading( std::ostream& stream ) const
{
    stream << "\\begin{tabular}{|l|c|c|}" << std::endl;
    stream << "\\hline" << std::endl;
    stream << "Component & Magnitude & Phase \\\\" << std::endl;
    stream << "\\hline" << std::endl;
}

void LauMagPhaseCoeffSet::printTableRow( std::ostream& stream ) const
{
    LauPrint print;
    TString resName = this->name();
    resName         = resName.ReplaceAll( "_", "\\_" );
    stream << resName << "  &  $";
    print.printFormat( stream, magnitude_->value() );
    stream << " \\pm ";
    print.printFormat( stream, magnitude_->error() );
    stream << "$  &  $";
    print.printFormat( stream, phase_->value() );
    stream << " \\pm ";
    print.printFormat( stream, phase_->error() );
    stream << "$ \\\\" << std::endl;
}

void LauMagPhaseCoeffSet::randomiseInitValues()
{
    if ( magnitude_->fixed() == kFALSE ) {
        // Choose a magnitude between 0.0 and 2.0
        Double_t mag = LauAbsCoeffSet::getRandomiser()->Rndm() * 2.0;
        magnitude_->initValue( mag );
        magnitude_->value( mag );
    }
    if ( phase_->fixed() == kFALSE ) {
        // Choose a phase between +- pi
        Double_t phase = LauAbsCoeffSet::getRandomiser()->Rndm() * LauConstants::twoPi -
                         LauConstants::pi;
        phase_->initValue( phase );
        phase_->value( phase );
    }
}

void LauMagPhaseCoeffSet::finaliseValues()
{
    // retrieve the current values from the parameters
    Double_t mag      = magnitude_->value();
    Double_t phase    = phase_->value();
    Double_t genPhase = phase_->genValue();

    // Check whether we have a negative magnitude.
    // If so make it positive and add pi to the phase.
    if ( mag < 0.0 ) {
        mag   *= -1.0;
        phase += LauConstants::pi;
    }

    // Check now whether the phase lies in the right range (-pi to pi).
    Bool_t withinRange( kFALSE );
    while ( withinRange == kFALSE ) {
        if ( phase > -LauConstants::pi && phase <= LauConstants::pi ) {
            withinRange = kTRUE;
        } else {
            // Not within the specified range
            if ( phase > LauConstants::pi ) {
                phase -= LauConstants::twoPi;
            } else if ( phase <= -LauConstants::pi ) {
                phase += LauConstants::twoPi;
            }
        }
    }

    // A further problem can occur when the generated phase is close to -pi or pi.
    // The phase can wrap over to the other end of the scale -
    // this leads to artificially large pulls so we wrap it back.
    Double_t diff = phase - genPhase;
    if ( diff > LauConstants::pi ) {
        phase -= LauConstants::twoPi;
    } else if ( diff < -LauConstants::pi ) {
        phase += LauConstants::twoPi;
    }

    // finally store the new values in the parameters
    // and update the pulls
    magnitude_->value( mag );
    magnitude_->updatePull();
    phase_->value( phase );
    phase_->updatePull();
}

const LauComplex& LauMagPhaseCoeffSet::particleCoeff()
{
    coeff_.setRealImagPart( magnitude_->unblindValue() * TMath::Cos( phase_->unblindValue() ),
                            magnitude_->unblindValue() * TMath::Sin( phase_->unblindValue() ) );
    return coeff_;
}

const LauComplex& LauMagPhaseCoeffSet::antiparticleCoeff()
{
    return this->particleCoeff();
}

void LauMagPhaseCoeffSet::setCoeffValues( const LauComplex& coeff,
                                          const LauComplex& coeffBar,
                                          Bool_t init )
{
    LauComplex average( coeff );
    average += coeffBar;
    average.rescale( 0.5 );

    Double_t magVal( average.abs() );
    Double_t phaseVal( average.arg() );

    magnitude_->value( magVal );
    phase_->value( phaseVal );

    if ( init ) {
        magnitude_->genValue( magVal );
        phase_->genValue( phaseVal );

        magnitude_->initValue( magVal );
        phase_->initValue( phaseVal );
    }
}

LauParameter LauMagPhaseCoeffSet::acp()
{
    TString parName( this->baseName() );
    parName += "_ACP";
    return LauParameter( parName, 0.0 );
}

LauAbsCoeffSet* LauMagPhaseCoeffSet::createClone( const TString& newName,
                                                  CloneOption cloneOption,
                                                  Double_t constFactor )
{
    LauAbsCoeffSet* clone( 0 );
    if ( cloneOption == All || cloneOption == TiePhase || cloneOption == TieMagnitude ) {
        clone = new LauMagPhaseCoeffSet( *this, cloneOption, constFactor );
        clone->name( newName );
    } else {
        std::cerr << "ERROR in LauMagPhaseCoeffSet::createClone : Invalid clone option" << std::endl;
    }
    return clone;
}
