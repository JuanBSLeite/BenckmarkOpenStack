
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

/*! \file LauRealImagCoeffSet.cc
    \brief File containing implementation of LauRealImagCoeffSet class.
*/

#include "LauRealImagCoeffSet.hh"

#include "LauComplex.hh"
#include "LauConstants.hh"
#include "LauParameter.hh"
#include "LauPrint.hh"

#include "TMath.h"
#include "TRandom.h"

#include <fstream>
#include <iostream>
#include <vector>

LauRealImagCoeffSet::LauRealImagCoeffSet( const TString& compName,
                                          Double_t x,
                                          Double_t y,
                                          Bool_t xFixed,
                                          Bool_t yFixed ) :
    LauAbsCoeffSet( compName ),
    x_( new LauParameter( "X", x, minRealImagPart_, maxRealImagPart_, xFixed ) ),
    y_( new LauParameter( "Y", y, minRealImagPart_, maxRealImagPart_, yFixed ) ),
    coeff_( x, y )
{
}

LauRealImagCoeffSet::LauRealImagCoeffSet( const LauRealImagCoeffSet& rhs,
                                          CloneOption cloneOption,
                                          Double_t constFactor ) :
    LauAbsCoeffSet( rhs.name() ),
    x_( 0 ),
    y_( 0 ),
    coeff_( rhs.coeff_ )
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
}

std::vector<LauParameter*> LauRealImagCoeffSet::getParameters()
{
    std::vector<LauParameter*> pars;
    pars.push_back( x_ );
    pars.push_back( y_ );
    return pars;
}

void LauRealImagCoeffSet::printParValues() const
{
    std::cout << "INFO in LauRealImagCoeffSet::printParValues : Component \"" << this->name()
              << "\" has real part = " << x_->value() << " and imaginary part = " << y_->value()
              << "." << std::endl;
}

void LauRealImagCoeffSet::printTableHeading( std::ostream& stream ) const
{
    stream << "\\begin{tabular}{|l|c|c|}" << std::endl;
    stream << "\\hline" << std::endl;
    stream << "Component & Real Part & Imaginary Part \\\\" << std::endl;
    stream << "\\hline" << std::endl;
}

void LauRealImagCoeffSet::printTableRow( std::ostream& stream ) const
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
    stream << "$ \\\\" << std::endl;
}

void LauRealImagCoeffSet::randomiseInitValues()
{
    if ( x_->fixed() == kFALSE ) {
        // Choose a value between -10.0 and 10.0
        Double_t value = LauAbsCoeffSet::getRandomiser()->Rndm() * 20.0 - 10.0;
        x_->initValue( value );
        x_->value( value );
    }
    if ( y_->fixed() == kFALSE ) {
        // Choose a value between -10.0 and 10.0
        Double_t value = LauAbsCoeffSet::getRandomiser()->Rndm() * 20.0 - 10.0;
        y_->initValue( value );
        y_->value( value );
    }
}

void LauRealImagCoeffSet::finaliseValues()
{
    x_->updatePull();
    y_->updatePull();
}

const LauComplex& LauRealImagCoeffSet::particleCoeff()
{
    coeff_.setRealImagPart( x_->unblindValue(), y_->unblindValue() );
    return coeff_;
}

const LauComplex& LauRealImagCoeffSet::antiparticleCoeff()
{
    return this->particleCoeff();
}

void LauRealImagCoeffSet::setCoeffValues( const LauComplex& coeff,
                                          const LauComplex& coeffBar,
                                          Bool_t init )
{
    LauComplex average( coeff );
    average += coeffBar;
    average.rescale( 0.5 );

    Double_t xVal( average.re() );
    Double_t yVal( average.im() );

    x_->value( xVal );
    y_->value( yVal );

    if ( init ) {
        x_->genValue( xVal );
        y_->genValue( yVal );

        x_->initValue( xVal );
        y_->initValue( yVal );
    }
}

LauParameter LauRealImagCoeffSet::acp()
{
    TString parName( this->baseName() );
    parName += "_ACP";
    return LauParameter( parName, 0.0 );
}

LauAbsCoeffSet* LauRealImagCoeffSet::createClone( const TString& newName,
                                                  CloneOption cloneOption,
                                                  Double_t constFactor )
{
    LauAbsCoeffSet* clone( 0 );
    if ( cloneOption == All || cloneOption == TieRealPart || cloneOption == TieImagPart ) {
        clone = new LauRealImagCoeffSet( *this, cloneOption, constFactor );
        clone->name( newName );
    } else {
        std::cerr << "ERROR in LauRealImagCoeffSet::createClone : Invalid clone option" << std::endl;
    }
    return clone;
}
