
/*
Copyright 2015 University of Warwick

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

/*! \file LauEFKLLMRes.cc
    \brief File containing implementation of LauEFKLLMRes class.
*/

#include "LauEFKLLMRes.hh"

#include "Lau1DCubicSpline.hh"
#include "LauKinematics.hh"
#include "LauResonanceInfo.hh"
#include "LauTextFileParser.hh"

#include "TMath.h"

#include <cstdlib>

Lau1DCubicSpline* LauEFKLLMRes::magSpline_   = 0;
Lau1DCubicSpline* LauEFKLLMRes::phaseSpline_ = 0;

LauEFKLLMRes::LauEFKLLMRes( LauResonanceInfo* resInfo,
                            const Int_t resPairAmpInt,
                            const LauDaughters* daughters ) :
    LauAbsResonance( resInfo, resPairAmpInt, daughters ),
    massFactor_( 0 )
{
    const Double_t massFactorVal = 0.;

    const TString& parNameBase = this->getSanitisedName();

    TString massFactorName( parNameBase );
    massFactorName += "_massFactor";
    massFactor_     = resInfo->getExtraParameter( massFactorName );
    if ( massFactor_ == 0 ) {
        massFactor_ = new LauParameter( massFactorName, massFactorVal, -10.0, 10.0, kTRUE );
        massFactor_->secondStage( kTRUE );
        resInfo->addExtraParameter( massFactor_ );
    }
}

LauEFKLLMRes::~LauEFKLLMRes()
{
}

void LauEFKLLMRes::initialise()
{
}

void LauEFKLLMRes::setResonanceParameter( const TString& name, const Double_t value )
{
    if ( name == "massFactor" ) {
        this->setMassFactor( value );
        std::cout << "INFO in LauEFKLLMRes::setResonanceParameter: Mass factor set to " << value
                  << std::endl;
    } else {
        std::cerr << "WARNING in LauEFKLLMRes::setResonanceParameter: Parameter name not reconised."
                  << std::endl;
    }
}

void LauEFKLLMRes::floatResonanceParameter( const TString& name )
{
    if ( name == "massFactor" ) {
        if ( massFactor_->fixed() ) {
            massFactor_->fixed( kFALSE );
            this->addFloatingParameter( massFactor_ );
        } else {
            std::cerr << "WARNING in LauEFKLLMRes::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else {
        std::cerr << "WARNING in LauEFKLLMRes::floatResonanceParameter: Parameter name not reconised."
                  << std::endl;
    }
}

LauParameter* LauEFKLLMRes::getResonanceParameter( const TString& name )
{
    if ( name == "massFactor" ) {
        return massFactor_;
    } else {
        std::cerr << "WARNING in LauEFKLLMRes::getResonanceParameter: Parameter name not reconised."
                  << std::endl;
        return 0;
    }
}

const std::vector<LauParameter*>& LauEFKLLMRes::getFloatingParameters()
{
    this->clearFloatingParameters();

    if ( ! this->fixMassFactor() ) {
        this->addFloatingParameter( massFactor_ );
    }

    return this->getParameters();
}

LauComplex LauEFKLLMRes::resAmp( Double_t mass, Double_t /*spinTerm*/ )
{
    LauComplex amp( 0.0, 0.0 );

    if ( magSpline_ == 0 || phaseSpline_ == 0 ) {
        std::cerr << "ERROR in LauEFKLLMRes::resAmp : One of the splines is null." << std::endl;
        return amp;
    }

    const Double_t massSq = mass * mass;
    const Double_t mag    = magSpline_->evaluate( massSq );
    const Double_t phase  = TMath::DegToRad() * phaseSpline_->evaluate( massSq );
    LauComplex ff( mag * TMath::Cos( phase ), mag * TMath::Sin( phase ) );

    amp = ff.scale( TMath::Power( mass, this->getMassFactor() ) );

    return amp;
}

void LauEFKLLMRes::setupFormFactor( const TString& inputFile )
{
    LauTextFileParser readFile( inputFile );
    readFile.processFile();

    std::vector<Double_t> mSqVals;
    std::vector<Double_t> magVals;
    std::vector<Double_t> phaseVals;

    std::vector<std::string> line;

    line = readFile.getNextLine();
    while ( ! line.empty() ) {
        UInt_t length = line.size();
        if ( length != 3 ) {
            std::cerr << "ERROR in LauEFKLLMRes::setupFormFactor : Unexpected number of fields in text file, aborting reading of form-factor information."
                      << std::endl;
            return;
        }
        mSqVals.push_back( atof( line[0].c_str() ) );
        magVals.push_back( atof( line[1].c_str() ) );
        phaseVals.push_back( atof( line[2].c_str() ) );
        line = readFile.getNextLine();
    }

    // Destroy any splines we already had defined but issue a warning just in case
    if ( magSpline_ != 0 || phaseSpline_ != 0 ) {
        std::cerr << "WARNING in LauEFKLLMRes::setupFormFactor : Overwriting previous form-factor splines with newly read values."
                  << std::endl;
        delete magSpline_;
        delete phaseSpline_;
    }
    magSpline_   = new Lau1DCubicSpline( mSqVals, magVals, Lau1DCubicSpline::AkimaSpline );
    phaseSpline_ = new Lau1DCubicSpline( mSqVals, phaseVals, Lau1DCubicSpline::AkimaSpline );
}

void LauEFKLLMRes::setMassFactor( const Double_t massFactor )
{
    massFactor_->value( massFactor );
    massFactor_->genValue( massFactor );
    massFactor_->initValue( massFactor );
}
