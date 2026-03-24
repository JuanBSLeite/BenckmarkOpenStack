
/*
Copyright 2011 University of Warwick

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

/*! \file LauMagPhaseCPCoeffSet.cc
    \brief File containing implementation of LauMagPhaseCPCoeffSet class.
*/

#include "LauMagPhaseCPCoeffSet.hh"

#include "LauComplex.hh"
#include "LauConstants.hh"
#include "LauParameter.hh"
#include "LauPrint.hh"

#include "TMath.h"
#include "TRandom.h"

#include <fstream>
#include <iostream>
#include <vector>

LauMagPhaseCPCoeffSet::LauMagPhaseCPCoeffSet( const TString& compName,
                                              Double_t mag,
                                              Double_t phase,
                                              Double_t magBar,
                                              Double_t phaseBar,
                                              Bool_t magFixed,
                                              Bool_t phaseFixed,
                                              Bool_t magBarFixed,
                                              Bool_t phaseBarFixed ) :
    LauAbsCoeffSet( compName ),
    mag_( new LauParameter( "Mag", mag, minMagnitude_, maxMagnitude_, magFixed ) ),
    phase_( new LauParameter( "Phase", phase, minPhase_, maxPhase_, phaseFixed ) ),
    magBar_( new LauParameter( "MagBar", magBar, minMagnitude_, maxMagnitude_, magBarFixed ) ),
    phaseBar_( new LauParameter( "PhaseBar", phaseBar, minPhase_, maxPhase_, phaseBarFixed ) ),
    particleCoeff_( mag * TMath::Cos( phase ), mag * TMath::Sin( phase ) ),
    antiparticleCoeff_( magBar * TMath::Cos( phaseBar ), magBar * TMath::Sin( phaseBar ) ),
    acp_( "ACP", ( magBar * magBar - mag * mag ) / ( magBar * magBar + mag * mag ), -1.0, 1.0 )
{
}

LauMagPhaseCPCoeffSet::LauMagPhaseCPCoeffSet( const LauMagPhaseCPCoeffSet& rhs,
                                              CloneOption cloneOption,
                                              Double_t constFactor ) :
    LauAbsCoeffSet( rhs.name() ),
    mag_( 0 ),
    phase_( 0 ),
    magBar_( 0 ),
    phaseBar_( 0 ),
    particleCoeff_( rhs.particleCoeff_ ),
    antiparticleCoeff_( rhs.antiparticleCoeff_ ),
    acp_( rhs.acp_ )
{
    if ( cloneOption == All || cloneOption == TieMagnitude ) {
        mag_    = rhs.mag_->createClone( constFactor );
        magBar_ = rhs.magBar_->createClone( constFactor );
    } else {
        mag_ = new LauParameter( "Mag",
                                 rhs.mag_->value(),
                                 minMagnitude_,
                                 maxMagnitude_,
                                 rhs.mag_->fixed() );
        if ( rhs.mag_->blind() ) {
            const LauBlind* blinder = rhs.mag_->blinder();
            mag_->blindParameter( blinder->blindingString(), blinder->blindingWidth() );
        }
        magBar_ = new LauParameter( "MagBar",
                                    rhs.magBar_->value(),
                                    minMagnitude_,
                                    maxMagnitude_,
                                    rhs.magBar_->fixed() );
        if ( rhs.magBar_->blind() ) {
            const LauBlind* blinder = rhs.magBar_->blinder();
            magBar_->blindParameter( blinder->blindingString(), blinder->blindingWidth() );
        }
    }

    if ( cloneOption == All || cloneOption == TiePhase ) {
        phase_    = rhs.phase_->createClone( constFactor );
        phaseBar_ = rhs.phaseBar_->createClone( constFactor );
    } else {
        phase_ =
            new LauParameter( "Phase", rhs.phase_->value(), minPhase_, maxPhase_, rhs.phase_->fixed() );
        if ( rhs.phase_->blind() ) {
            const LauBlind* blinder = rhs.phase_->blinder();
            phase_->blindParameter( blinder->blindingString(), blinder->blindingWidth() );
        }
        phaseBar_ = new LauParameter( "PhaseBar",
                                      rhs.phaseBar_->value(),
                                      minPhase_,
                                      maxPhase_,
                                      rhs.phaseBar_->fixed() );
        if ( rhs.phaseBar_->blind() ) {
            const LauBlind* blinder = rhs.phaseBar_->blinder();
            phaseBar_->blindParameter( blinder->blindingString(), blinder->blindingWidth() );
        }
    }
}

std::vector<LauParameter*> LauMagPhaseCPCoeffSet::getParameters()
{
    std::vector<LauParameter*> pars;
    pars.push_back( mag_ );
    pars.push_back( phase_ );
    pars.push_back( magBar_ );
    pars.push_back( phaseBar_ );
    return pars;
}

void LauMagPhaseCPCoeffSet::printParValues() const
{
    std::cout << "INFO in LauMagPhaseCPCoeffSet::printParValues : Component \"" << this->name()
              << "\" has ";
    std::cout << "mag = " << mag_->value() << ",\t";
    std::cout << "phase = " << phase_->value() << ",\t";
    std::cout << "magBar = " << magBar_->value() << ",\t";
    std::cout << "phaseBar = " << phaseBar_->value() << "." << std::endl;
}

void LauMagPhaseCPCoeffSet::printTableHeading( std::ostream& stream ) const
{
    stream << "\\begin{tabular}{|l|c|c|c|c|}" << std::endl;
    stream << "\\hline" << std::endl;
    stream << "Component & Magnitude & Phase & Magnitude_bar & Phase_bar \\\\" << std::endl;
    stream << "\\hline" << std::endl;
}

void LauMagPhaseCPCoeffSet::printTableRow( std::ostream& stream ) const
{
    LauPrint print;
    TString resName = this->name();
    resName         = resName.ReplaceAll( "_", "\\_" );
    stream << resName << "  &  $";
    print.printFormat( stream, mag_->value() );
    stream << " \\pm ";
    print.printFormat( stream, mag_->error() );
    stream << "$  &  $";
    print.printFormat( stream, phase_->value() );
    stream << " \\pm ";
    print.printFormat( stream, phase_->error() );
    stream << "$  &  $";
    print.printFormat( stream, magBar_->value() );
    stream << " \\pm ";
    print.printFormat( stream, magBar_->error() );
    stream << "$  &  $";
    print.printFormat( stream, phaseBar_->value() );
    stream << " \\pm ";
    print.printFormat( stream, phaseBar_->error() );
    stream << "$ \\\\" << std::endl;
}

void LauMagPhaseCPCoeffSet::randomiseInitValues()
{
    if ( mag_->fixed() == kFALSE ) {
        // Choose a value for "magnitude" between 0.0 and 2.0
        Double_t value = LauAbsCoeffSet::getRandomiser()->Rndm() * 2.0;
        mag_->initValue( value );
        mag_->value( value );
    }
    if ( phase_->fixed() == kFALSE ) {
        // Choose a phase between +- pi
        Double_t value = LauAbsCoeffSet::getRandomiser()->Rndm() * LauConstants::twoPi -
                         LauConstants::pi;
        phase_->initValue( value );
        phase_->value( value );
    }
    if ( magBar_->fixed() == kFALSE ) {
        // Choose a value for "magnitude" between 0.0 and 2.0
        Double_t value = LauAbsCoeffSet::getRandomiser()->Rndm() * 2.0;
        magBar_->initValue( value );
        magBar_->value( value );
    }
    if ( phaseBar_->fixed() == kFALSE ) {
        // Choose a phase between +- pi
        Double_t value = LauAbsCoeffSet::getRandomiser()->Rndm() * LauConstants::twoPi -
                         LauConstants::pi;
        phaseBar_->initValue( value );
        phaseBar_->value( value );
    }
}

void LauMagPhaseCPCoeffSet::finaliseValues()
{
    // retrieve the current values from the parameters
    Double_t mVal        = mag_->value();
    Double_t pVal        = phase_->value();
    Double_t mBarVal     = magBar_->value();
    Double_t pBarVal     = phaseBar_->value();
    Double_t genPhase    = phase_->genValue();
    Double_t genPhaseBar = phaseBar_->genValue();

    // Check whether we have a negative magnitude.
    // If so make it positive and add pi to the  phase.
    if ( mVal < 0.0 ) {
        mVal *= -1.0;
        pVal += LauConstants::pi;
    }
    if ( mBarVal < 0.0 ) {
        mBarVal *= -1.0;
        pBarVal += LauConstants::pi;
    }

    // Check now whether the phases lies in the right range (-pi to pi).
    Bool_t pWithinRange( kFALSE );
    Bool_t pBarWithinRange( kFALSE );
    while ( pWithinRange == kFALSE && pBarWithinRange == kFALSE ) {
        if ( pVal > -LauConstants::pi && pVal < LauConstants::pi ) {
            pWithinRange = kTRUE;
        } else {
            // Not within the specified range
            if ( pVal > LauConstants::pi ) {
                pVal -= LauConstants::twoPi;
            } else if ( pVal < -LauConstants::pi ) {
                pVal += LauConstants::twoPi;
            }
        }

        if ( pBarVal > -LauConstants::pi && pBarVal < LauConstants::pi ) {
            pBarWithinRange = kTRUE;
        } else {
            // Not within the specified range
            if ( pBarVal > LauConstants::pi ) {
                pBarVal -= LauConstants::twoPi;
            } else if ( pBarVal < -LauConstants::pi ) {
                pBarVal += LauConstants::twoPi;
            }
        }
    }

    // A further problem can occur when the generated phase is close to -pi or pi.
    // The phase can wrap over to the other end of the scale -
    // this leads to artificially large pulls so we wrap it back.
    Double_t diff = pVal - genPhase;
    if ( diff > LauConstants::pi ) {
        pVal -= LauConstants::twoPi;
    } else if ( diff < -LauConstants::pi ) {
        pVal += LauConstants::twoPi;
    }

    diff = pBarVal - genPhaseBar;
    if ( diff > LauConstants::pi ) {
        pBarVal -= LauConstants::twoPi;
    } else if ( diff < -LauConstants::pi ) {
        pBarVal += LauConstants::twoPi;
    }

    // finally store the new values in the parameters
    // and update the pulls
    mag_->value( mVal );
    mag_->updatePull();
    phase_->value( pVal );
    phase_->updatePull();
    magBar_->value( mBarVal );
    magBar_->updatePull();
    phaseBar_->value( pBarVal );
    phaseBar_->updatePull();
}

const LauComplex& LauMagPhaseCPCoeffSet::particleCoeff()
{
    particleCoeff_.setRealImagPart( mag_->unblindValue() * TMath::Cos( phase_->unblindValue() ),
                                    mag_->unblindValue() * TMath::Sin( phase_->unblindValue() ) );
    return particleCoeff_;
}

const LauComplex& LauMagPhaseCPCoeffSet::antiparticleCoeff()
{
    antiparticleCoeff_.setRealImagPart(
        magBar_->unblindValue() * TMath::Cos( phaseBar_->unblindValue() ),
        magBar_->unblindValue() * TMath::Sin( phaseBar_->unblindValue() ) );
    return antiparticleCoeff_;
}

void LauMagPhaseCPCoeffSet::setCoeffValues( const LauComplex& coeff,
                                            const LauComplex& coeffBar,
                                            Bool_t init )
{
    Double_t magVal( coeff.abs() );
    Double_t phaseVal( coeff.arg() );
    Double_t magBarVal( coeffBar.abs() );
    Double_t phaseBarVal( coeffBar.arg() );

    mag_->value( magVal );
    phase_->value( phaseVal );
    magBar_->value( magBarVal );
    phaseBar_->value( phaseBarVal );

    if ( init ) {
        mag_->genValue( magVal );
        phase_->genValue( phaseVal );
        magBar_->genValue( magBarVal );
        phaseBar_->genValue( phaseBarVal );

        mag_->initValue( magVal );
        phase_->initValue( phaseVal );
        magBar_->initValue( magBarVal );
        phaseBar_->initValue( phaseBarVal );
    }
}

LauParameter LauMagPhaseCPCoeffSet::acp()
{
    // set the name
    TString parName( this->baseName() );
    parName += "_ACP";
    acp_.name( parName );

    // work out the ACP value

    Double_t value( -99.0 );
    value = ( magBar_->value() * magBar_->value() - mag_->value() * mag_->value() ) /
            ( magBar_->value() * magBar_->value() + mag_->value() * mag_->value() );

    // is it fixed?
    Bool_t fixed = magBar_->fixed() && mag_->fixed();
    acp_.fixed( fixed );

    // we can't work out the error without the covariance matrix
    Double_t error( 0.0 );

    // set the value and error
    acp_.valueAndErrors( value, error );

    return acp_;
}

LauAbsCoeffSet* LauMagPhaseCPCoeffSet::createClone( const TString& newName,
                                                    CloneOption cloneOption,
                                                    Double_t constFactor )
{
    LauAbsCoeffSet* clone( 0 );
    if ( cloneOption == All || cloneOption == TiePhase || cloneOption == TieMagnitude ) {
        clone = new LauMagPhaseCPCoeffSet( *this, cloneOption, constFactor );
        clone->name( newName );
    } else {
        std::cerr << "ERROR in LauMagPhaseCPCoeffSet::createClone : Invalid clone option"
                  << std::endl;
    }
    return clone;
}
