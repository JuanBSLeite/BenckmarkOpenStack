
/*
Copyright 2004 University of Warwick

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

/*! \file LauAbsPdf.cc
    \brief File containing implementation of LauAbsPdf class.
*/

#include "LauAbsPdf.hh"

#include "LauIntegrals.hh"
#include "LauKinematics.hh"
#include "LauRandom.hh"

#include "TMath.h"
#include "TSystem.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <vector>

// Constructor for the abstract PDF class.
LauAbsPdf::LauAbsPdf( const TString& theVarName,
                      const std::vector<LauAbsRValue*>& params,
                      Double_t minAbscissa,
                      Double_t maxAbscissa ) :
    varNames_(),
    param_( params ),
    norm_( 0.0 ),
    maxHeight_( 0.0 ),
    heightUpToDate_( kFALSE ),
    minAbscissas_(),
    maxAbscissas_(),
    randomFun_( LauRandom::randomFun() ),
    cachePDF_( kFALSE ),
    unNormPDFVal_( 0.0 ),
    nNormPoints_( 50 ),
    integMethod_( GaussLegendre ),
    withinNormCalc_( kFALSE ),
    withinGeneration_( kFALSE ),
    normWeightsDone_( kFALSE )
{
    // Store the variable name
    varNames_.insert( std::make_pair( 0, theVarName ) );

    // Store the min and max values
    minAbscissas_.push_back( minAbscissa );
    maxAbscissas_.push_back( maxAbscissa );
}

LauAbsPdf::LauAbsPdf( const std::vector<TString>& theVarNames,
                      const std::vector<LauAbsRValue*>& params,
                      const LauFitData& minAbscissas,
                      const LauFitData& maxAbscissas ) :
    varNames_(),
    param_( params ),
    norm_( 0.0 ),
    maxHeight_( 0.0 ),
    heightUpToDate_( kFALSE ),
    minAbscissas_(),
    maxAbscissas_(),
    randomFun_( LauRandom::randomFun() ),
    cachePDF_( kFALSE ),
    unNormPDFVal_( 0.0 ),
    nNormPoints_( 50 ),
    integMethod_( GaussLegendre ),
    withinNormCalc_( kFALSE ),
    withinGeneration_( kFALSE ),
    normWeightsDone_( kFALSE )
{
    // Check that we have at least one variable
    if ( theVarNames.empty() ) {
        std::cerr << "ERROR in LauAbsPdf::LauAbsPdf : No variables supplied." << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Store the variable names
    for ( UInt_t i( 0 ); i < theVarNames.size(); ++i ) {
        varNames_.insert( std::make_pair( i, theVarNames[i] ) );
    }

    // Store the min and max abcissas for every variable
    UInt_t nVars = varNames_.size();
    minAbscissas_.resize( nVars );
    maxAbscissas_.resize( nVars );
    for ( std::map<UInt_t, TString>::const_iterator iter = varNames_.begin(); iter != varNames_.end();
          ++iter ) {
        UInt_t varIndex                     = iter->first;
        const TString& theVarName           = iter->second;
        LauFitData::const_iterator min_iter = minAbscissas.find( theVarName );
        LauFitData::const_iterator max_iter = maxAbscissas.find( theVarName );
        if ( min_iter == minAbscissas.end() || max_iter == maxAbscissas.end() ) {
            std::cerr << "ERROR in LauAbsPdf::LauAbsPdf : Either min or max abscissa not provided for variable \""
                      << theVarName << "\"." << std::endl;
            gSystem->Exit( EXIT_FAILURE );
        }
        minAbscissas_[varIndex] = min_iter->second;
        maxAbscissas_[varIndex] = max_iter->second;
    }
}

std::vector<TString> LauAbsPdf::varNames() const
{
    std::vector<TString> theVarNames;
    theVarNames.reserve( varNames_.size() );
    for ( std::map<UInt_t, TString>::const_iterator iter = varNames_.begin(); iter != varNames_.end();
          ++iter ) {
        theVarNames.push_back( iter->second );
    }
    return theVarNames;
}

UInt_t LauAbsPdf::nFixedParameters() const
{
    LauParamFixed pred;
    return std::count_if( param_.begin(), param_.end(), pred );
}

Double_t LauAbsPdf::getMinAbscissa( const TString& theVarName ) const
{
    for ( std::map<UInt_t, TString>::const_iterator iter = varNames_.begin(); iter != varNames_.end();
          ++iter ) {
        UInt_t index        = iter->first;
        const TString& name = iter->second;
        if ( name == theVarName ) {
            return minAbscissas_[index];
        }
    }

    std::cerr << "ERROR in LauAbsPdf::getMinAbscissa : Variable \"" << theVarName << "\" not found."
              << std::endl;
    return 0.0;
}

Double_t LauAbsPdf::getMaxAbscissa( const TString& theVarName ) const
{
    for ( std::map<UInt_t, TString>::const_iterator iter = varNames_.begin(); iter != varNames_.end();
          ++iter ) {
        UInt_t index        = iter->first;
        const TString& name = iter->second;
        if ( name == theVarName ) {
            return maxAbscissas_[index];
        }
    }

    std::cerr << "ERROR in LauAbsPdf::getMaxAbscissa : Variable \"" << theVarName << "\" not found."
              << std::endl;
    return 0.0;
}

LauFitData LauAbsPdf::getMinAbscissas() const
{
    LauFitData minAbscissas;
    for ( std::map<UInt_t, TString>::const_iterator iter = varNames_.begin(); iter != varNames_.end();
          ++iter ) {
        UInt_t varIndex           = iter->first;
        const TString& theVarName = iter->second;
        minAbscissas[theVarName]  = minAbscissas_[varIndex];
    }
    return minAbscissas;
}

LauFitData LauAbsPdf::getMaxAbscissas() const
{
    LauFitData maxAbscissas;
    for ( std::map<UInt_t, TString>::const_iterator iter = varNames_.begin(); iter != varNames_.end();
          ++iter ) {
        UInt_t varIndex           = iter->first;
        const TString& theVarName = iter->second;
        maxAbscissas[theVarName]  = maxAbscissas_[varIndex];
    }
    return maxAbscissas;
}

LauFitData LauAbsPdf::getRanges() const
{
    LauFitData minVals = this->getMinAbscissas();
    LauFitData maxVals = this->getMaxAbscissas();
    LauFitData range;

    for ( LauFitData::const_iterator iter = maxVals.begin(); iter != maxVals.end(); ++iter ) {
        const TString& theVarName = iter->first;
        Double_t maxVal           = iter->second;
        Double_t minVal           = minVals.find( theVarName )->second;
        range[theVarName]         = maxVal - minVal;
    }

    return range;
}

void LauAbsPdf::setMinAbscissa( const TString& theVarName, Double_t minAbscissa )
{
    for ( std::map<UInt_t, TString>::const_iterator iter = varNames_.begin(); iter != varNames_.end();
          ++iter ) {
        UInt_t index        = iter->first;
        const TString& name = iter->second;
        if ( name == theVarName ) {
            minAbscissas_[index] = minAbscissa;
            return;
        }
    }

    std::cerr << "ERROR in LauAbsPdf::setMinAbscissa : Variable \"" << theVarName << "\" not found."
              << std::endl;
}

void LauAbsPdf::setMaxAbscissa( const TString& theVarName, Double_t maxAbscissa )
{
    for ( std::map<UInt_t, TString>::const_iterator iter = varNames_.begin(); iter != varNames_.end();
          ++iter ) {
        UInt_t index        = iter->first;
        const TString& name = iter->second;
        if ( name == theVarName ) {
            maxAbscissas_[index] = maxAbscissa;
            return;
        }
    }

    std::cerr << "ERROR in LauAbsPdf::setMaxAbscissa : Variable \"" << theVarName << "\" not found."
              << std::endl;
}

Bool_t LauAbsPdf::checkRange( const LauAbscissas& abscissas ) const
{
    // This method assumes the ordering of the abscissas provided match
    // those of the min and max vectors

    // check that the size is correct
    UInt_t nVars = abscissas.size();
    if ( this->isDPDependent() ) {
        nVars = abscissas.size() -
                2;    // if we depend on the DP we'll have been provided with the DP co-ordinates as well
    }
    if ( nVars != minAbscissas_.size() || nVars != maxAbscissas_.size() ) {
        std::cerr << "ERROR in LauAbsPdf::checkRange : Unexpected number of absicssas: " << nVars
                  << std::endl;
        std::cerr << "                               : " << this->IsA()->GetName() << " expects "
                  << minAbscissas_.size() << std::endl;
        return kFALSE;
    }

    for ( UInt_t i( 0 ); i < nVars; ++i ) {
        Double_t abscissa = abscissas[i];
        Double_t minVal   = minAbscissas_[i];
        Double_t maxVal   = maxAbscissas_[i];
        if ( ( abscissa < minVal ) || ( abscissa > maxVal ) ) {
            std::cerr << "ERROR in LauAbsPdf::checkRange : " << abscissa
                      << " outside allowed range: [" << minVal << "," << maxVal << "]" << std::endl;
            return kFALSE;
        }
    }
    return kTRUE;
}

void LauAbsPdf::cacheInfo( const LauFitDataTree& inputData )
{
    Bool_t hasBranch( kTRUE );
    for ( std::map<UInt_t, TString>::const_iterator iter = varNames_.begin(); iter != varNames_.end();
          ++iter ) {
        hasBranch &= inputData.haveBranch( iter->second );
        if ( ! hasBranch ) {
            std::cerr << "ERROR in LauAbsPdf::cacheInfo : Input data does not contain variable \""
                      << iter->second << "\"." << std::endl;
            return;
        }
    }

    // determine whether we are caching our PDF value
    Bool_t doCaching( this->nFixedParameters() == this->nParameters() );
    this->cachePDF( doCaching );

    // clear the vectors and reserve enough space
    UInt_t nEvents = inputData.nEvents();
    abscissas_.clear();
    abscissas_.reserve( nEvents );
    unNormPDFValues_.clear();
    unNormPDFValues_.reserve( nEvents );

    for ( UInt_t iEvt = 0; iEvt < nEvents; ++iEvt ) {

        const LauFitData& dataValues = inputData.getData( iEvt );

        LauAbscissas myData;
        // add all our variables into the data
        for ( std::map<UInt_t, TString>::const_iterator var_iter = varNames_.begin();
              var_iter != varNames_.end();
              ++var_iter ) {
            LauFitData::const_iterator iter = dataValues.find( var_iter->second );
            myData.push_back( iter->second );
        }
        // if we're DP dependent then we'll need the DP co-ordinates as well
        if ( this->isDPDependent() ) {
            LauFitData::const_iterator iter = dataValues.find( "m13Sq" );
            myData.push_back( iter->second );
            iter = dataValues.find( "m23Sq" );
            myData.push_back( iter->second );
        }

        if ( ! this->checkRange( myData ) ) {
            gSystem->Exit( EXIT_FAILURE );
        }

        abscissas_.push_back( myData );

        if ( this->cachePDF() ) {
            this->calcLikelihoodInfo( myData );
            unNormPDFValues_.push_back( this->getUnNormLikelihood() );
        }
    }

    if ( ! this->cachePDF() ) {
        // in this case we seem to be doing a fit where the parameters are floating
        // so need to mark that the PDF height is no longer up to date
        this->heightUpToDate( kFALSE );
    }
}

LauFitData LauAbsPdf::generate( const LauKinematics* kinematics )
{
    this->withinGeneration( kTRUE );

    // Check that the PDF height is up to date
    // N.B. this must now called every time (the method will simply
    // return if there is nothing to do)
    this->calcPDFHeight( kinematics );

    Bool_t gotAbscissa( kFALSE );
    if ( randomFun_ == 0 ) {
        std::cerr << "ERROR in LauAbsPdf::generate : Please set the random number generator for this PDF by using the setRandomFun(TRandom*) function."
                  << std::endl;
        this->withinGeneration( kFALSE );
        gSystem->Exit( EXIT_FAILURE );
    }

    if ( this->isDPDependent() && ! kinematics ) {
        std::cerr << "ERROR in LauAbsPdf::generate : PDF depends on the DP and an invalid kinematics pointer has been provided."
                  << std::endl;
        this->withinGeneration( kFALSE );
        gSystem->Exit( EXIT_FAILURE );
    }

    // container for holding the generated abscissa(s)
    LauAbscissas genAbscissa( 1 );

    // Generate the value of the abscissa.
    Double_t genPDFVal( 0.0 );
    Double_t PDFheight = this->getMaxHeight() * ( 1.0 + 1e-11 );
    while ( ! gotAbscissa ) {

        if ( this->isDPDependent() ) {
            genAbscissa.resize( 3 );
            genAbscissa[1] = kinematics->getm13Sq();
            genAbscissa[2] = kinematics->getm23Sq();
        }
        genAbscissa[0] = randomFun_->Rndm() * this->getRange() + this->getMinAbscissa();

        this->calcLikelihoodInfo( genAbscissa );
        genPDFVal = this->getUnNormLikelihood();

        if ( randomFun_->Rndm() <= genPDFVal / PDFheight ) {
            gotAbscissa = kTRUE;
        }

        if ( genPDFVal > PDFheight ) {
            std::cerr << "WARNING in LauAbsPdf::generate : genPDFVal = " << genPDFVal
                      << " is larger than the maximum PDF height " << this->getMaxHeight()
                      << " for the abscissa = " << genAbscissa[0] << "." << std::endl;
            std::cerr << "                               : Need to reset height to be larger than "
                      << genPDFVal
                      << " by using the setMaxHeight(Double_t) function and re-run the Monte Carlo generation!"
                      << std::endl;
        }
    }

    LauFitData genData;
    genData[this->varName()] = genAbscissa[0];

    this->withinGeneration( kFALSE );

    return genData;
}

Double_t LauAbsPdf::getLikelihood() const
{
    if ( TMath::Abs( norm_ ) > 1e-10 ) {
        return unNormPDFVal_ / norm_;
    } else {
        return 0.0;
    }
}

Double_t LauAbsPdf::getLikelihood( const TString& theVarName ) const
{
    if ( theVarName != this->varName() ) {
        std::cerr << "ERROR in LauAbsPdf::getLikelihood : Unrecognised variable name \""
                  << theVarName << "\", cannot determine likelihood." << std::endl;
        return 0.0;
    }
    return this->getLikelihood();
}

void LauAbsPdf::calcLikelihoodInfo( UInt_t iEvt )
{
    if ( this->cachePDF() && ( unNormPDFValues_.size() == abscissas_.size() ) ) {
        unNormPDFVal_ = unNormPDFValues_[iEvt];
    } else {
        this->calcLikelihoodInfo( abscissas_[iEvt] );
    }
}

LauAbsRValue* LauAbsPdf::findParameter( const TString& parName )
{
    for ( std::vector<LauAbsRValue*>::iterator iter = param_.begin(); iter != param_.end(); ++iter ) {
        //    std::vector<LauParameter*> params = (*iter)->getPars();
        //    for (std::vector<LauParameter*>::iterator params_iter = params.begin(); params_iter != params.end(); ++params_iter ) {
        if ( ( *iter )->name().Contains( parName ) ) {
            return ( *iter );
        }
        //    }
    }

    std::cerr << "ERROR in LauAbsPdf::findParameter : Parameter \"" << parName << "\" not found."
              << std::endl;
    return 0;
}

const LauAbsRValue* LauAbsPdf::findParameter( const TString& parName ) const
{
    for ( std::vector<LauAbsRValue*>::const_iterator iter = param_.begin(); iter != param_.end();
          ++iter ) {
        //    std::vector<LauParameter*> params = (*iter)->getPars();
        //    for (std::vector<LauParameter*>::iterator params_iter = params.begin(); params_iter != params.end(); ++params_iter ) {
        if ( ( *iter )->name().Contains( parName ) ) {
            return ( *iter );
            //        }
        }
    }

    std::cerr << "ERROR in LauAbsPdf::findParameter : Parameter \"" << parName << "\" not found."
              << std::endl;
    return 0;
}

void LauAbsPdf::updatePulls()
{
    for ( std::vector<LauAbsRValue*>::iterator iter = param_.begin(); iter != param_.end(); ++iter ) {
        std::vector<LauParameter*> params = ( *iter )->getPars();
        for ( std::vector<LauParameter*>::iterator params_iter = params.begin();
              params_iter != params.end();
              ++params_iter ) {
            if ( ! ( *iter )->fixed() ) {
                ( *params_iter )->updatePull();
            }
        }
    }
}

void LauAbsPdf::addParameters( std::vector<LauAbsRValue*>& params )
{
    for ( std::vector<LauAbsRValue*>::iterator iter = params.begin(); iter != params.end(); ++iter ) {
        param_.push_back( *iter );
    }
}

void LauAbsPdf::calcNorm()
{
    this->withinNormCalc( kTRUE );

    if ( this->nInputVars() > 1 ) {
        std::cerr << "ERROR in LauAbsPdf::calcNorm : Numeric integration only works for 1D PDFs."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    IntMethod sumMethod = this->integMethod();

    Double_t normFac = ( sumMethod == GaussLegendre ) ? this->integrGaussLegendre()
                                                      : this->integTrapezoid();

    this->setNorm( normFac );

    this->withinNormCalc( kFALSE );
}

Double_t LauAbsPdf::integrGaussLegendre()
{
    if ( ! this->normWeightsDone() ) {
        this->getNormWeights();
    }

    // Now compute the integral
    Double_t norm( 0.0 );
    for ( UInt_t i = 0; i < normWeights_.size(); i++ ) {
        this->calcLikelihoodInfo( normAbscissas_[i] );
        Double_t fun        = this->getUnNormLikelihood();
        Double_t intFactor  = 0.5 * this->getRange();
        norm               += normWeights_[i] * intFactor * fun;
    }

    //std::cout<<"====================================================="<<std::endl;
    //std::cout<<"NORM = "<<norm<<std::endl;

    //std::cout<<"====================================================="<<std::endl;
    return norm;
}

void LauAbsPdf::getNormWeights()
{
    // Check whether we've already calculated the weights
    if ( this->normWeightsDone() ) {
        std::cerr << "WARNING in LauAbsPdf::getNormWeights : Already calculated weights, not doing it again."
                  << std::endl;
        return;
    }

    // Avoid integral if we have no points in x space
    if ( nNormPoints_ == 0 ) {
        std::cerr << "ERROR in LauAbsPdf::getNormWeights : Zero points specified, this is daft!"
                  << std::endl;
        return;
    }

    // Calculate the normalisation weights and abscissas
    Double_t precision( 1e-6 );

    Double_t intMean   = 0.5 * ( this->getMaxAbscissa() + this->getMinAbscissa() );
    Double_t range     = this->getMaxAbscissa() - this->getMinAbscissa();
    Double_t halfRange = 0.5 * range;

    std::vector<Double_t> abscissas;
    LauIntegrals funIntegrals( precision );
    funIntegrals.calcGaussLegendreWeights( nNormPoints_, abscissas, normWeights_ );

    //std::cout<<"====================================================="<<std::endl;
    //std::cout<<"NORM POINTS = "<<nNormPoints_<<std::endl;

    //std::cout<<"====================================================="<<std::endl;

    Int_t nWeights = static_cast<Int_t>( normWeights_.size() );
    normAbscissas_.resize( nWeights );

    // Use same number of abscissas for x and y co-ordinates
    Int_t m = ( nWeights + 1 ) / 2;
    for ( Int_t i = 0; i < m; ++i ) {

        Int_t ii = nWeights - 1 - i;    // symmetric i index

        Double_t dx     = halfRange * abscissas[i];
        Double_t tmpVal = intMean - dx;
        normAbscissas_[i].push_back( tmpVal );

        tmpVal = intMean + dx;
        normAbscissas_[ii].push_back( tmpVal );
    }

    this->normWeightsDone( kTRUE );
}

Double_t LauAbsPdf::integTrapezoid()
{
    Double_t abscVal, tnm, sum, del;
    Int_t it, j;

    static Double_t norm( 0.0 );
    Double_t range = this->getRange();

    if ( this->nNormPoints() == 1 ) {

        LauAbscissas abscissa( 1 );
        abscissa[0] = this->getMinAbscissa();
        this->calcLikelihoodInfo( abscissa );
        Double_t funAbsMin = this->getUnNormLikelihood();

        abscissa[0] = this->getMinAbscissa();
        this->calcLikelihoodInfo( abscissa );
        Double_t funAbsMax = this->getUnNormLikelihood();

        norm = 0.5 * range * ( funAbsMin + funAbsMax );
        return norm;

    } else {
        for ( it = 1, j = 1; j < this->nNormPoints() - 1; j++ ) {
            it <<= 1;
        }
        tnm     = it;
        del     = range / tnm;
        abscVal = this->getMinAbscissa() + 0.5 * del;

        for ( sum = 0.0, j = 1; j < it; j++, abscVal += del ) {

            LauAbscissas abscissa( 1 );
            abscissa[0] = abscVal;
            this->calcLikelihoodInfo( abscissa );
            Double_t funVal = this->getUnNormLikelihood();

            sum += funVal;
        }

        norm = 0.5 * ( norm + sum * range / tnm );
        return norm;
    }
}
