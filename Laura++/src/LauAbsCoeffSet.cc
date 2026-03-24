
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

/*! \file LauAbsCoeffSet.cc
    \brief File containing implementation of LauAbsCoeffSet class.
*/

#include "LauAbsCoeffSet.hh"

#include "LauConstants.hh"
#include "LauParameter.hh"
#include "LauRandom.hh"

#include "TString.h"

#include <iostream>

TRandom* LauAbsCoeffSet::randomiser_      = 0;
Double_t LauAbsCoeffSet::minMagnitude_    = -10.0;
Double_t LauAbsCoeffSet::maxMagnitude_    = 10.0;
Double_t LauAbsCoeffSet::minPhase_        = -LauConstants::threePi;
Double_t LauAbsCoeffSet::maxPhase_        = LauConstants::threePi;
Double_t LauAbsCoeffSet::minRealImagPart_ = -10.0;
Double_t LauAbsCoeffSet::maxRealImagPart_ = 10.0;
Double_t LauAbsCoeffSet::minDelta_        = -2.0;
Double_t LauAbsCoeffSet::maxDelta_        = 2.0;

LauAbsCoeffSet::LauAbsCoeffSet( const TString& theName, const TString& theBaseName ) :
    name_( theName ),
    basename_( theBaseName ),
    index_( 0 )
{
}

LauAbsCoeffSet::LauAbsCoeffSet( const LauAbsCoeffSet& rhs ) :
    name_( rhs.name_ ),
    basename_( rhs.basename_ ),
    index_( rhs.index_ )
{
}

LauAbsCoeffSet& LauAbsCoeffSet::operator=( const LauAbsCoeffSet& rhs )
{
    if ( &rhs != this ) {
        name_     = rhs.name_;
        basename_ = rhs.basename_;
        index_    = rhs.index_;
    }
    return *this;
}

TRandom* LauAbsCoeffSet::getRandomiser()
{
    if ( randomiser_ == 0 ) {
        randomiser_ = LauRandom::zeroSeedRandom();
    }
    return randomiser_;
}

void LauAbsCoeffSet::index( UInt_t newIndex )
{
    index_ = newIndex;

    const TString oldBaseName( this->baseName() );
    TString basename( oldBaseName );
    if ( basename.Length() != 1 ) {
        basename.Remove( 1 );
    }
    basename += newIndex;
    basename += "_";
    this->baseName( basename );

    std::vector<LauParameter*> pars = this->getParameters();
    for ( std::vector<LauParameter*>::iterator iter = pars.begin(); iter != pars.end(); ++iter ) {
        this->adjustName( *iter, oldBaseName );
    }
}

void LauAbsCoeffSet::adjustName( LauParameter* par, const TString& oldBaseName )
{
    TString theName( par->name() );
    if ( theName.BeginsWith( oldBaseName ) && theName != oldBaseName ) {
        theName.Remove( 0, oldBaseName.Length() );
    }
    theName.Prepend( this->baseName() );
    par->name( theName );
}

void LauAbsCoeffSet::setParameterValue( const TString& parName, Double_t value, Bool_t init )
{
    LauParameter* par = this->findParameter( parName );
    if ( par == 0 ) {
        std::cerr << "ERROR in LauAbsCoeffSet::setParameterValue : Unable to find parameter \""
                  << parName << "\"" << std::endl;
        return;
    }

    par->value( value );
    if ( init ) {
        par->genValue( value );
        par->initValue( value );
    }
}

void LauAbsCoeffSet::setParameterError( const TString& parName, Double_t error )
{
    LauParameter* par = this->findParameter( parName );
    if ( par == 0 ) {
        std::cerr << "ERROR in LauAbsCoeffSet::setParameterError : Unable to find parameter \""
                  << parName << "\"" << std::endl;
        return;
    }

    par->error( error );
}

void LauAbsCoeffSet::fixParameter( const TString& parName )
{
    LauParameter* par = this->findParameter( parName );
    if ( par == 0 ) {
        std::cerr << "ERROR in LauAbsCoeffSet::fixParameter : Unable to find parameter \""
                  << parName << "\"" << std::endl;
        return;
    }

    par->fixed( kTRUE );
}

void LauAbsCoeffSet::floatParameter( const TString& parName )
{
    LauParameter* par = this->findParameter( parName );
    if ( par == 0 ) {
        std::cerr << "ERROR in LauAbsCoeffSet::floatParameter : Unable to find parameter \""
                  << parName << "\"" << std::endl;
        return;
    }

    par->fixed( kFALSE );
}

void LauAbsCoeffSet::blindParameter( const TString& parName,
                                     const TString& blindingString,
                                     const Double_t width )
{
    LauParameter* par = this->findParameter( parName );
    if ( par == 0 ) {
        std::cerr << "ERROR in LauAbsCoeffSet::blindParameter : Unable to find parameter \""
                  << parName << "\"" << std::endl;
        return;
    }

    par->blindParameter( blindingString, width );
}

void LauAbsCoeffSet::addGaussianConstraint( const TString& parName,
                                            const Double_t mean,
                                            const Double_t width )
{
    LauParameter* par = this->findParameter( parName );
    if ( par == 0 ) {
        std::cerr << "ERROR in LauAbsCoeffSet::addGaussianConstraint : Unable to find parameter \""
                  << parName << "\"" << std::endl;
        return;
    }

    par->addGaussianConstraint( mean, width );
}

void LauAbsCoeffSet::addSuffixToParameterName( const TString& parName, const TString& suffix )
{
    LauParameter* par = this->findParameter( parName );
    if ( par == 0 ) {
        std::cerr << "ERROR in LauAbsCoeffSet::addSuffixToParameterName : Unable to find parameter \""
                  << parName << "\"" << std::endl;
        return;
    }

    TString newName( par->name() );
    if ( ! suffix.BeginsWith( '_' ) ) {
        newName += "_";
    }
    newName += suffix;
    par->name( newName );
}

LauParameter* LauAbsCoeffSet::findParameter( const TString& parName )
{
    std::vector<LauParameter*> pars = this->getParameters();
    for ( std::vector<LauParameter*>::iterator iter = pars.begin(); iter != pars.end(); ++iter ) {

        const TString& iName = ( *iter )->name();

        if ( iName.EndsWith( parName ) ) {
            return ( *iter );
        }
    }

    return 0;
}
