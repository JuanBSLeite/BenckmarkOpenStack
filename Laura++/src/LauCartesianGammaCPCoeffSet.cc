
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

/*! \file LauCartesianGammaCPCoeffSet.cc
    \brief File containing implementation of LauCartesianGammaCPCoeffSet class.
*/

#include "LauCartesianGammaCPCoeffSet.hh"

#include "LauComplex.hh"
#include "LauConstants.hh"
#include "LauParameter.hh"
#include "LauPrint.hh"

#include "TMath.h"
#include "TRandom.h"

#include <fstream>
#include <iostream>
#include <vector>

LauCartesianGammaCPCoeffSet::LauCartesianGammaCPCoeffSet( const TString& compName,
                                                          const Double_t x,
                                                          const Double_t y,
                                                          const Double_t xCP,
                                                          const Double_t yCP,
                                                          const Double_t deltaXCP,
                                                          const Double_t deltaYCP,
                                                          const Bool_t xFixed,
                                                          const Bool_t yFixed,
                                                          const Bool_t xCPFixed,
                                                          const Bool_t yCPFixed,
                                                          const Bool_t deltaXCPFixed,
                                                          const Bool_t deltaYCPFixed,
                                                          const Bool_t deltaXCPSecondStage,
                                                          const Bool_t deltaYCPSecondStage ) :
    LauAbsCoeffSet( compName ),
    x_( new LauParameter( "X", x, minRealImagPart_, maxRealImagPart_, xFixed ) ),
    y_( new LauParameter( "Y", y, minRealImagPart_, maxRealImagPart_, yFixed ) ),
    xCP_( new LauParameter( "XCP", xCP, minRealImagPart_, maxRealImagPart_, xCPFixed ) ),
    yCP_( new LauParameter( "YCP", yCP, minRealImagPart_, maxRealImagPart_, yCPFixed ) ),
    deltaXCP_( new LauParameter( "DeltaXCP", deltaXCP, minDelta_, maxDelta_, deltaXCPFixed ) ),
    deltaYCP_( new LauParameter( "DeltaYCP", deltaYCP, minDelta_, maxDelta_, deltaYCPFixed ) ),
    nonCPPart_( x, y ),
    cpPart_( 1 + xCP + deltaXCP, yCP + deltaYCP ),
    cpAntiPart_( 1 + xCP - deltaXCP, yCP - deltaYCP ),
    particleCoeff_( nonCPPart_ * cpPart_ ),
    antiparticleCoeff_( nonCPPart_ * cpAntiPart_ ),
    acp_( "ACP",
          ( antiparticleCoeff_.abs2() - particleCoeff_.abs2() ) /
              ( antiparticleCoeff_.abs2() + particleCoeff_.abs2() ),
          -1.0,
          1.0,
          deltaXCPFixed && deltaYCPFixed )
{
    if ( deltaXCPSecondStage && ! deltaXCPFixed ) {
        deltaXCP_->secondStage( kTRUE );
        deltaXCP_->initValue( 0.0 );
    }
    if ( deltaYCPSecondStage && ! deltaYCPFixed ) {
        deltaYCP_->secondStage( kTRUE );
        deltaYCP_->initValue( 0.0 );
    }
}

LauCartesianGammaCPCoeffSet::LauCartesianGammaCPCoeffSet( const LauCartesianGammaCPCoeffSet& rhs,
                                                          CloneOption cloneOption,
                                                          Double_t constFactor ) :
    LauAbsCoeffSet( rhs.name() ),
    x_( 0 ),
    y_( 0 ),
    xCP_( 0 ),
    yCP_( 0 ),
    deltaXCP_( 0 ),
    deltaYCP_( 0 ),
    nonCPPart_( rhs.nonCPPart_ ),
    cpPart_( rhs.cpPart_ ),
    cpAntiPart_( rhs.cpAntiPart_ ),
    particleCoeff_( rhs.particleCoeff_ ),
    antiparticleCoeff_( rhs.antiparticleCoeff_ ),
    acp_( rhs.acp_ )
{
    if ( cloneOption == All || cloneOption == TieRealPart ) {
        x_ = rhs.x_->createClone( constFactor );
    } else {
        x_ = new LauParameter( "X",
                               rhs.x_->value(),
                               minRealImagPart_,
                               maxRealImagPart_,
                               rhs.x_->fixed() );
        if ( rhs.x_->blind() ) {
            const LauBlind* blinder = rhs.x_->blinder();
            x_->blindParameter( blinder->blindingString(), blinder->blindingWidth() );
        }
    }

    if ( cloneOption == All || cloneOption == TieImagPart ) {
        y_ = rhs.y_->createClone( constFactor );
    } else {
        y_ = new LauParameter( "Y",
                               rhs.y_->value(),
                               minRealImagPart_,
                               maxRealImagPart_,
                               rhs.y_->fixed() );
        if ( rhs.y_->blind() ) {
            const LauBlind* blinder = rhs.y_->blinder();
            y_->blindParameter( blinder->blindingString(), blinder->blindingWidth() );
        }
    }

    if ( cloneOption == All || cloneOption == TieCPPars ) {
        xCP_      = rhs.xCP_->createClone( constFactor );
        yCP_      = rhs.yCP_->createClone( constFactor );
        deltaXCP_ = rhs.deltaXCP_->createClone( constFactor );
        deltaYCP_ = rhs.deltaYCP_->createClone( constFactor );
    } else {
        xCP_ = new LauParameter( "XCP",
                                 rhs.xCP_->value(),
                                 minRealImagPart_,
                                 maxRealImagPart_,
                                 rhs.xCP_->fixed() );
        if ( rhs.xCP_->blind() ) {
            const LauBlind* blinder = rhs.xCP_->blinder();
            xCP_->blindParameter( blinder->blindingString(), blinder->blindingWidth() );
        }
        yCP_ = new LauParameter( "YCP",
                                 rhs.yCP_->value(),
                                 minRealImagPart_,
                                 maxRealImagPart_,
                                 rhs.yCP_->fixed() );
        if ( rhs.yCP_->blind() ) {
            const LauBlind* blinder = rhs.yCP_->blinder();
            yCP_->blindParameter( blinder->blindingString(), blinder->blindingWidth() );
        }
        deltaXCP_ = new LauParameter( "DeltaXCP",
                                      rhs.deltaXCP_->value(),
                                      minDelta_,
                                      maxDelta_,
                                      rhs.deltaXCP_->fixed() );
        deltaYCP_ = new LauParameter( "DeltaYCP",
                                      rhs.deltaYCP_->value(),
                                      minDelta_,
                                      maxDelta_,
                                      rhs.deltaYCP_->fixed() );
        if ( rhs.deltaXCP_->secondStage() && ! rhs.deltaXCP_->fixed() ) {
            deltaXCP_->secondStage( kTRUE );
            deltaXCP_->initValue( 0.0 );
        }
        if ( rhs.deltaYCP_->secondStage() && ! rhs.deltaYCP_->fixed() ) {
            deltaYCP_->secondStage( kTRUE );
            deltaYCP_->initValue( 0.0 );
        }
        if ( rhs.deltaXCP_->blind() ) {
            const LauBlind* blinder = rhs.deltaXCP_->blinder();
            deltaXCP_->blindParameter( blinder->blindingString(), blinder->blindingWidth() );
        }
        if ( rhs.deltaYCP_->blind() ) {
            const LauBlind* blinder = rhs.deltaYCP_->blinder();
            deltaYCP_->blindParameter( blinder->blindingString(), blinder->blindingWidth() );
        }
    }
}

std::vector<LauParameter*> LauCartesianGammaCPCoeffSet::getParameters()
{
    std::vector<LauParameter*> pars;
    pars.push_back( x_ );
    pars.push_back( y_ );
    if ( ! xCP_->fixed() )
        pars.push_back( xCP_ );
    if ( ! yCP_->fixed() )
        pars.push_back( yCP_ );
    if ( ! deltaXCP_->fixed() )
        pars.push_back( deltaXCP_ );
    if ( ! deltaYCP_->fixed() )
        pars.push_back( deltaYCP_ );
    return pars;
}

void LauCartesianGammaCPCoeffSet::printParValues() const
{
    std::cout << "INFO in LauCartesianGammaCPCoeffSet::printParValues : Component \""
              << this->name() << "\" has ";
    std::cout << "x = " << x_->value() << ",\t";
    std::cout << "y = " << y_->value() << ",\t";
    std::cout << "xCP = " << xCP_->value() << ",\t";
    std::cout << "yCP = " << yCP_->value() << ",\t";
    std::cout << "Delta xCP = " << deltaXCP_->value() << ",\t";
    std::cout << "Delta yCP = " << deltaYCP_->value() << "." << std::endl;
}

void LauCartesianGammaCPCoeffSet::printTableHeading( std::ostream& stream ) const
{
    stream << "\\begin{tabular}{|l|c|c|c|c|c|c|}" << std::endl;
    stream << "\\hline" << std::endl;
    stream << "Component & Real Part & Imaginary Part & Real CP Part & Imaginary CP Part & $\\Delta$ Real CP Part & $\\Delta$ Imaginary CP Part \\\\"
           << std::endl;
    stream << "\\hline" << std::endl;
}

void LauCartesianGammaCPCoeffSet::printTableRow( std::ostream& stream ) const
{
    LauPrint print;
    TString resName = this->name();
    resName         = resName.ReplaceAll( "_", "\\_" );
    stream << resName << "  &  $";
    print.printFormat( stream, x_->value() );
    stream << " \\pm ";
    print.printFormat( stream, x_->error() );
    stream << "$  &  $";
    print.printFormat( stream, y_->value() );
    stream << " \\pm ";
    print.printFormat( stream, y_->error() );
    stream << "$  &  $";
    print.printFormat( stream, xCP_->value() );
    stream << " \\pm ";
    print.printFormat( stream, xCP_->error() );
    stream << "$  &  $";
    print.printFormat( stream, yCP_->value() );
    stream << " \\pm ";
    print.printFormat( stream, yCP_->error() );
    stream << "$  &  $";
    print.printFormat( stream, deltaXCP_->value() );
    stream << " \\pm ";
    print.printFormat( stream, deltaXCP_->error() );
    stream << "$  &  $";
    print.printFormat( stream, deltaYCP_->value() );
    stream << " \\pm ";
    print.printFormat( stream, deltaYCP_->error() );
    stream << "$ \\\\" << std::endl;
}

void LauCartesianGammaCPCoeffSet::randomiseInitValues()
{
    if ( x_->fixed() == kFALSE ) {
        // Choose a value for "X" between -3.0 and 3.0
        Double_t value = LauAbsCoeffSet::getRandomiser()->Rndm() * 6.0 - 3.0;
        x_->initValue( value );
        x_->value( value );
    }
    if ( y_->fixed() == kFALSE ) {
        // Choose a value for "Y" between -3.0 and 3.0
        Double_t value = LauAbsCoeffSet::getRandomiser()->Rndm() * 6.0 - 3.0;
        y_->initValue( value );
        y_->value( value );
    }
    if ( xCP_->fixed() == kFALSE ) {
        // Choose a value for "XCP" between -3.0 and 3.0
        Double_t value = LauAbsCoeffSet::getRandomiser()->Rndm() * 6.0 - 3.0;
        xCP_->initValue( value );
        xCP_->value( value );
    }
    if ( yCP_->fixed() == kFALSE ) {
        // Choose a value for "YCP" between -3.0 and 3.0
        Double_t value = LauAbsCoeffSet::getRandomiser()->Rndm() * 6.0 - 3.0;
        yCP_->initValue( value );
        yCP_->value( value );
    }
    if ( deltaXCP_->fixed() == kFALSE && deltaXCP_->secondStage() == kFALSE ) {
        // Choose a value for "Delta XCP" between -0.5 and 0.5
        Double_t value = LauAbsCoeffSet::getRandomiser()->Rndm() * 1.0 - 0.5;
        deltaXCP_->initValue( value );
        deltaXCP_->value( value );
    }
    if ( deltaYCP_->fixed() == kFALSE && deltaYCP_->secondStage() == kFALSE ) {
        // Choose a value for "Delta YCP" between -0.5 and 0.5
        Double_t value = LauAbsCoeffSet::getRandomiser()->Rndm() * 1.0 - 0.5;
        deltaYCP_->initValue( value );
        deltaYCP_->value( value );
    }
}

void LauCartesianGammaCPCoeffSet::finaliseValues()
{
    // update the pulls
    x_->updatePull();
    y_->updatePull();
    xCP_->updatePull();
    yCP_->updatePull();
    deltaXCP_->updatePull();
    deltaYCP_->updatePull();
}

const LauComplex& LauCartesianGammaCPCoeffSet::particleCoeff()
{
    nonCPPart_.setRealImagPart( x_->unblindValue(), y_->unblindValue() );
    cpPart_.setRealImagPart( 1.0 + xCP_->unblindValue() + deltaXCP_->unblindValue(),
                             yCP_->unblindValue() + deltaYCP_->unblindValue() );
    particleCoeff_ = nonCPPart_ * cpPart_;
    return particleCoeff_;
}

const LauComplex& LauCartesianGammaCPCoeffSet::antiparticleCoeff()
{
    nonCPPart_.setRealImagPart( x_->unblindValue(), y_->unblindValue() );
    cpAntiPart_.setRealImagPart( 1.0 + xCP_->unblindValue() - deltaXCP_->unblindValue(),
                                 yCP_->unblindValue() - deltaYCP_->unblindValue() );
    antiparticleCoeff_ = nonCPPart_ * cpAntiPart_;
    return antiparticleCoeff_;
}

void LauCartesianGammaCPCoeffSet::setCoeffValues( const LauComplex&, const LauComplex&, Bool_t )
{
    std::cerr << "ERROR in LauCartesianGammaCPCoeffSet::setCoeffValues : Method not supported by this class - too many parameters"
              << std::endl;
}

LauParameter LauCartesianGammaCPCoeffSet::acp()
{
    // set the name
    TString parName( this->baseName() );
    parName += "_ACP";
    acp_.name( parName );

    // work out the ACP value
    const LauComplex nonCPPart( x_->value(), y_->value() );
    const LauComplex cpPart( 1.0 + xCP_->value() + deltaXCP_->value(),
                             yCP_->value() + deltaYCP_->value() );
    const LauComplex cpAntiPart( 1.0 + xCP_->value() - deltaXCP_->value(),
                                 yCP_->value() - deltaYCP_->value() );
    const LauComplex partCoeff = nonCPPart * cpPart;
    const LauComplex antiCoeff = nonCPPart * cpAntiPart;

    const Double_t numer = antiCoeff.abs2() - partCoeff.abs2();
    const Double_t denom = antiCoeff.abs2() + partCoeff.abs2();
    const Double_t value = numer / denom;

    // is it fixed?
    const Bool_t fixed = deltaXCP_->fixed() && deltaYCP_->fixed();
    acp_.fixed( fixed );

    // we can't work out the error without the covariance matrix
    const Double_t error( 0.0 );

    // set the value and error
    acp_.valueAndErrors( value, error );

    return acp_;
}

LauAbsCoeffSet* LauCartesianGammaCPCoeffSet::createClone( const TString& newName,
                                                          CloneOption cloneOption,
                                                          Double_t constFactor )
{
    LauAbsCoeffSet* clone( 0 );
    if ( cloneOption == All || cloneOption == TieRealPart || cloneOption == TieImagPart ||
         cloneOption == TieCPPars ) {
        clone = new LauCartesianGammaCPCoeffSet( *this, cloneOption, constFactor );
        clone->name( newName );
    } else {
        std::cerr << "ERROR in LauCartesianGammaCPCoeffSet::createClone : Invalid clone option"
                  << std::endl;
    }
    return clone;
}
