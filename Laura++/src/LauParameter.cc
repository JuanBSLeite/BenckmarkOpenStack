
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

/*! \file LauParameter.cc
    \brief File containing implementation of LauParameter class.
*/

#include <iostream>
#include <map>
using std::cerr;
using std::cout;
using std::endl;

#include "LauParameter.hh"
#include "LauRandom.hh"

#include "TMath.h"
#include "TRandom.h"

LauParameter::LauParameter() :
    name_( "" ),
    value_( 0.0 ),
    error_( 0.0 ),
    negError_( 0.0 ),
    posError_( 0.0 ),
    genValue_( 0.0 ),
    initValue_( 0.0 ),
    minValue_( 0.0 ),
    maxValue_( 0.0 ),
    fixed_( kTRUE ),
    secondStage_( kFALSE ),
    gaussConstraint_( kFALSE ),
    constraintTrueMean_( 0.0 ),
    constraintMean_( 0.0 ),
    constraintWidth_( 0.0 ),
    gcc_( 0.0 ),
    bias_( 0.0 ),
    pull_( 0.0 ),
    clone_( kFALSE ),
    parent_( 0 ),
    blinder_( 0 )
{
}

LauParameter::LauParameter( const TString& parName ) :
    name_( parName ),
    value_( 0.0 ),
    error_( 0.0 ),
    negError_( 0.0 ),
    posError_( 0.0 ),
    genValue_( 0.0 ),
    initValue_( 0.0 ),
    minValue_( 0.0 ),
    maxValue_( 0.0 ),
    fixed_( kTRUE ),
    secondStage_( kFALSE ),
    gaussConstraint_( kFALSE ),
    constraintTrueMean_( 0.0 ),
    constraintMean_( 0.0 ),
    constraintWidth_( 0.0 ),
    gcc_( 0.0 ),
    bias_( 0.0 ),
    pull_( 0.0 ),
    clone_( kFALSE ),
    parent_( 0 ),
    blinder_( 0 )
{
}

LauParameter::LauParameter( Double_t parValue ) :
    name_( "" ),
    value_( parValue ),
    error_( 0.0 ),
    negError_( 0.0 ),
    posError_( 0.0 ),
    genValue_( parValue ),
    initValue_( parValue ),
    minValue_( parValue - 1e-6 ),
    maxValue_( parValue + 1e-6 ),
    fixed_( kTRUE ),
    secondStage_( kFALSE ),
    gaussConstraint_( kFALSE ),
    constraintTrueMean_( 0.0 ),
    constraintMean_( 0.0 ),
    constraintWidth_( 0.0 ),
    gcc_( 0.0 ),
    bias_( 0.0 ),
    pull_( 0.0 ),
    clone_( kFALSE ),
    parent_( 0 ),
    blinder_( 0 )
{
}

LauParameter::LauParameter( const TString& parName, Double_t parValue ) :
    name_( parName ),
    value_( parValue ),
    error_( 0.0 ),
    negError_( 0.0 ),
    posError_( 0.0 ),
    genValue_( parValue ),
    initValue_( parValue ),
    minValue_( parValue - 1e-6 ),
    maxValue_( parValue + 1e-6 ),
    fixed_( kTRUE ),
    secondStage_( kFALSE ),
    gaussConstraint_( kFALSE ),
    constraintTrueMean_( 0.0 ),
    constraintMean_( 0.0 ),
    constraintWidth_( 0.0 ),
    gcc_( 0.0 ),
    bias_( 0.0 ),
    pull_( 0.0 ),
    clone_( kFALSE ),
    parent_( 0 ),
    blinder_( 0 )
{
}

LauParameter::LauParameter( Double_t parValue, Double_t min, Double_t max ) :
    name_( "" ),
    value_( parValue ),
    error_( 0.0 ),
    negError_( 0.0 ),
    posError_( 0.0 ),
    genValue_( parValue ),
    initValue_( parValue ),
    minValue_( min ),
    maxValue_( max ),
    fixed_( kTRUE ),
    secondStage_( kFALSE ),
    gaussConstraint_( kFALSE ),
    constraintTrueMean_( 0.0 ),
    constraintMean_( 0.0 ),
    constraintWidth_( 0.0 ),
    gcc_( 0.0 ),
    bias_( 0.0 ),
    pull_( 0.0 ),
    clone_( kFALSE ),
    parent_( 0 ),
    blinder_( 0 )
{
    this->checkRange();
}

LauParameter::LauParameter( Double_t parValue, Double_t parError, Double_t min, Double_t max ) :
    name_( "" ),
    value_( parValue ),
    error_( parError ),
    negError_( 0.0 ),
    posError_( 0.0 ),
    genValue_( parValue ),
    initValue_( parValue ),
    minValue_( min ),
    maxValue_( max ),
    fixed_( kTRUE ),
    secondStage_( kFALSE ),
    gaussConstraint_( kFALSE ),
    constraintTrueMean_( 0.0 ),
    constraintMean_( 0.0 ),
    constraintWidth_( 0.0 ),
    gcc_( 0.0 ),
    bias_( 0.0 ),
    pull_( 0.0 ),
    clone_( kFALSE ),
    parent_( 0 ),
    blinder_( 0 )
{
    this->checkRange();
}

LauParameter::LauParameter( Double_t parValue, Double_t min, Double_t max, Bool_t parFixed ) :
    name_( "" ),
    value_( parValue ),
    error_( 0.0 ),
    negError_( 0.0 ),
    posError_( 0.0 ),
    genValue_( parValue ),
    initValue_( parValue ),
    minValue_( min ),
    maxValue_( max ),
    fixed_( parFixed ),
    secondStage_( kFALSE ),
    gaussConstraint_( kFALSE ),
    constraintTrueMean_( 0.0 ),
    constraintMean_( 0.0 ),
    constraintWidth_( 0.0 ),
    gcc_( 0.0 ),
    bias_( 0.0 ),
    pull_( 0.0 ),
    clone_( kFALSE ),
    parent_( 0 ),
    blinder_( 0 )
{
    this->checkRange();
}

LauParameter::LauParameter( const TString& parName, Double_t parValue, Double_t min, Double_t max ) :
    name_( parName ),
    value_( parValue ),
    error_( 0.0 ),
    negError_( 0.0 ),
    posError_( 0.0 ),
    genValue_( parValue ),
    initValue_( parValue ),
    minValue_( min ),
    maxValue_( max ),
    fixed_( kTRUE ),
    secondStage_( kFALSE ),
    gaussConstraint_( kFALSE ),
    constraintTrueMean_( 0.0 ),
    constraintMean_( 0.0 ),
    constraintWidth_( 0.0 ),
    gcc_( 0.0 ),
    bias_( 0.0 ),
    pull_( 0.0 ),
    clone_( kFALSE ),
    parent_( 0 ),
    blinder_( 0 )
{
    this->checkRange();
}

LauParameter::LauParameter( const TString& parName,
                            Double_t parValue,
                            Double_t min,
                            Double_t max,
                            Bool_t parFixed ) :
    name_( parName ),
    value_( parValue ),
    error_( 0.0 ),
    negError_( 0.0 ),
    posError_( 0.0 ),
    genValue_( parValue ),
    initValue_( parValue ),
    minValue_( min ),
    maxValue_( max ),
    fixed_( parFixed ),
    secondStage_( kFALSE ),
    gaussConstraint_( kFALSE ),
    constraintTrueMean_( 0.0 ),
    constraintMean_( 0.0 ),
    constraintWidth_( 0.0 ),
    gcc_( 0.0 ),
    bias_( 0.0 ),
    pull_( 0.0 ),
    clone_( kFALSE ),
    parent_( 0 ),
    blinder_( 0 )
{
    this->checkRange();
}

LauParameter::LauParameter( const TString& parName,
                            Double_t parValue,
                            Double_t parError,
                            Double_t min,
                            Double_t max ) :
    name_( parName ),
    value_( parValue ),
    error_( parError ),
    negError_( 0.0 ),
    posError_( 0.0 ),
    genValue_( parValue ),
    initValue_( parValue ),
    minValue_( min ),
    maxValue_( max ),
    fixed_( kTRUE ),
    secondStage_( kFALSE ),
    gaussConstraint_( kFALSE ),
    constraintTrueMean_( 0.0 ),
    constraintMean_( 0.0 ),
    constraintWidth_( 0.0 ),
    gcc_( 0.0 ),
    bias_( 0.0 ),
    pull_( 0.0 ),
    clone_( kFALSE ),
    parent_( 0 ),
    blinder_( 0 )
{
    this->checkRange();
}

LauParameter::LauParameter( const LauParameter& rhs ) :
    TObject( rhs ),
    LauAbsRValue( rhs ),
    name_( rhs.name_ ),
    value_( rhs.value_ ),
    error_( rhs.error_ ),
    negError_( rhs.negError_ ),
    posError_( rhs.posError_ ),
    genValue_( rhs.genValue_ ),
    initValue_( rhs.initValue_ ),
    minValue_( rhs.minValue_ ),
    maxValue_( rhs.maxValue_ ),
    fixed_( rhs.fixed_ ),
    secondStage_( rhs.secondStage_ ),
    gaussConstraint_( rhs.gaussConstraint_ ),
    constraintTrueMean_( rhs.constraintTrueMean_ ),
    constraintMean_( rhs.constraintMean_ ),
    constraintWidth_( rhs.constraintWidth_ ),
    gcc_( rhs.gcc_ ),
    bias_( rhs.bias_ ),
    pull_( rhs.pull_ ),
    clone_( rhs.clone_ ),
    parent_( rhs.parent_ ),
    clones_( rhs.clones_ ),
    blinder_( ( rhs.blinder_ == 0 ) ? 0 : new LauBlind( *( rhs.blinder_ ) ) )
{
}

LauParameter& LauParameter::operator=( const LauParameter& rhs )
{
    if ( &rhs != this ) {
        TObject::operator=( rhs );
        LauAbsRValue::operator=( rhs );
        name_               = rhs.name_;
        value_              = rhs.value_;
        error_              = rhs.error_;
        negError_           = rhs.negError_;
        posError_           = rhs.posError_;
        genValue_           = rhs.genValue_;
        initValue_          = rhs.initValue_;
        minValue_           = rhs.minValue_;
        maxValue_           = rhs.maxValue_;
        fixed_              = rhs.fixed_;
        secondStage_        = rhs.secondStage_;
        gaussConstraint_    = rhs.gaussConstraint_;
        constraintTrueMean_ = rhs.constraintTrueMean_;
        constraintMean_     = rhs.constraintMean_;
        constraintWidth_    = rhs.constraintWidth_;
        gcc_                = rhs.gcc_;
        bias_               = rhs.bias_;
        pull_               = rhs.pull_;
        clone_              = rhs.clone_;
        parent_             = rhs.parent_;
        clones_             = rhs.clones_;
        delete blinder_;
        blinder_ = ( rhs.blinder_ == 0 ) ? 0 : new LauBlind( *( rhs.blinder_ ) );
    }
    return *this;
}

LauParameter::~LauParameter()
{
    delete blinder_;
}

std::vector<LauParameter*> LauParameter::getPars()
{
    std::vector<LauParameter*> list;
    list.push_back( this );
    return list;
}

void LauParameter::value( Double_t newValue )
{
    if ( this->clone() ) {
        parent_->value( newValue );
    } else {
        this->checkRange( newValue, this->minValue(), this->maxValue() );
        this->updateClones( kTRUE );
    }
}

void LauParameter::error( Double_t newError )
{
    if ( this->clone() ) {
        parent_->error( newError );
    } else {
        error_ = TMath::Abs( newError );
        this->updateClones( kFALSE );
    }
}

void LauParameter::negError( Double_t newNegError )
{
    if ( this->clone() ) {
        parent_->negError( newNegError );
    } else {
        negError_ = TMath::Abs( newNegError );
        this->updateClones( kFALSE );
    }
}

void LauParameter::posError( Double_t newPosError )
{
    if ( this->clone() ) {
        parent_->posError( newPosError );
    } else {
        posError_ = TMath::Abs( newPosError );
        this->updateClones( kFALSE );
    }
}

void LauParameter::errors( Double_t newError, Double_t newNegError, Double_t newPosError )
{
    if ( this->clone() ) {
        parent_->errors( newError, newNegError, newPosError );
    } else {
        error_    = TMath::Abs( newError );
        negError_ = TMath::Abs( newNegError );
        posError_ = TMath::Abs( newPosError );
        this->updateClones( kFALSE );
    }
}

void LauParameter::valueAndErrors( Double_t newValue,
                                   Double_t newError,
                                   Double_t newNegError,
                                   Double_t newPosError )
{
    if ( this->clone() ) {
        parent_->valueAndErrors( newValue, newError, newNegError, newPosError );
    } else {
        this->checkRange( newValue, this->minValue(), this->maxValue() );
        error_    = TMath::Abs( newError );
        negError_ = TMath::Abs( newNegError );
        posError_ = TMath::Abs( newPosError );
        this->updateClones( kFALSE );
    }
}

void LauParameter::globalCorrelationCoeff( Double_t newGCCValue )
{
    if ( this->clone() ) {
        parent_->globalCorrelationCoeff( newGCCValue );
    } else {
        gcc_ = newGCCValue;
        this->updateClones( kFALSE );
    }
}

void LauParameter::genValue( Double_t newGenValue )
{
    if ( this->clone() ) {
        parent_->genValue( newGenValue );
    } else {
        genValue_ = newGenValue;
        this->updateClones( kFALSE );
    }
}

void LauParameter::initValue( Double_t newInitValue )
{
    if ( this->clone() ) {
        parent_->initValue( newInitValue );
    } else {
        initValue_ = newInitValue;
        this->updateClones( kFALSE );
    }
}

void LauParameter::minValue( Double_t newMinValue )
{
    if ( this->clone() ) {
        parent_->minValue( newMinValue );
    } else {
        this->checkRange( this->value(), newMinValue, this->maxValue() );
        this->updateClones( kFALSE );
    }
}

void LauParameter::maxValue( Double_t newMaxValue )
{
    if ( this->clone() ) {
        parent_->maxValue( newMaxValue );
    } else {
        this->checkRange( this->value(), this->minValue(), newMaxValue );
        this->updateClones( kFALSE );
    }
}

void LauParameter::range( Double_t newMinValue, Double_t newMaxValue )
{
    if ( this->clone() ) {
        parent_->range( newMinValue, newMaxValue );
    } else {
        this->checkRange( this->value(), newMinValue, newMaxValue );
        this->updateClones( kFALSE );
    }
}

void LauParameter::valueAndRange( Double_t newValue, Double_t newMinValue, Double_t newMaxValue )
{
    if ( this->clone() ) {
        parent_->valueAndRange( newValue, newMinValue, newMaxValue );
    } else {
        this->checkRange( newValue, newMinValue, newMaxValue );
        this->updateClones( kFALSE );
    }
}

void LauParameter::name( const TString& newName )
{
    // no need to update clones here
    // clones are allowed to have different names
    name_ = newName;
}

void LauParameter::fixed( Bool_t parFixed )
{
    if ( this->clone() ) {
        parent_->fixed( parFixed );
    } else {
        fixed_ = parFixed;
        this->updateClones( kFALSE );
    }
}

void LauParameter::secondStage( Bool_t secondStagePar )
{
    if ( this->clone() ) {
        parent_->secondStage( secondStagePar );
    } else {
        secondStage_ = secondStagePar;
        this->updateClones( kFALSE );
    }
}

void LauParameter::addGaussianConstraint( Double_t newGaussMean, Double_t newGaussWidth )
{
    if ( this->clone() ) {
        parent_->addGaussianConstraint( newGaussMean, newGaussWidth );
    } else {
        gaussConstraint_    = kTRUE;
        constraintTrueMean_ = newGaussMean;
        constraintMean_     = newGaussMean;
        constraintWidth_    = newGaussWidth;
        this->updateClones( kFALSE );
    }
}

void LauParameter::removeGaussianConstraint()
{
    if ( this->clone() ) {
        parent_->removeGaussianConstraint();
    } else {
        gaussConstraint_ = kFALSE;
        this->updateClones( kFALSE );
    }
}

void LauParameter::generateConstraintMean()
{
    if ( this->clone() ) {
        parent_->generateConstraintMean();
    } else {
        constraintMean_ = LauRandom::randomFun()->Gaus( constraintTrueMean_, constraintWidth_ );
        this->updateClones( kFALSE );
    }
}

Double_t LauParameter::constraintPenalty() const
{
    const Double_t val { this->unblindValue() };
    const Double_t diff { val - constraintMean_ };
    const Double_t term { diff * diff };

    return term / ( 2.0 * constraintWidth_ * constraintWidth_ );
}

void LauParameter::blindParameter( const TString& blindingString, const Double_t width )
{
    if ( this->clone() ) {
        parent_->blindParameter( blindingString, width );
        return;
    }

    if ( blinder_ != 0 ) {
        std::cerr << "WARNING in LauParameter::blindParameter : blinding has already been set up for this parameter"
                  << std::endl;
        return;
    }

    blinder_ = new LauBlind( blindingString, width );

    for ( std::map<LauParameter*, Double_t>::iterator iter = clones_.begin(); iter != clones_.end();
          ++iter ) {
        LauParameter* clonePar = iter->first;
        if ( clonePar->blinder_ != 0 ) {
            std::cerr << "WARNING in LauParameter::blindParameter : blinding has already been set up for a clone of this parameter - it will be replaced!"
                      << std::endl;
            delete clonePar->blinder_;
            clonePar->blinder_ = 0;
        }
        clonePar->blinder_ = new LauBlind( *blinder_ );
    }
}

void LauParameter::updatePull()
{
    if ( this->clone() ) {
        parent_->updatePull();
        return;
    }

    // calculate the bias
    bias_ = value_ - genValue_;

    // if we have errors calculated then calculate
    // the pull using the best error available
    if ( ( bias_ > 0.0 ) && ( negError_ > 1e-10 ) ) {
        pull_ = bias_ / negError_;
    } else if ( ( bias_ < 0.0 ) && ( posError_ > 1e-10 ) ) {
        pull_ = bias_ / posError_;
    } else if ( error_ > 1e-10 ) {
        pull_ = bias_ / error_;
    } else {
        pull_ = 0.0;
    }

    this->updateClones( kFALSE );
}

void LauParameter::checkRange( Double_t val, Double_t minVal, Double_t maxVal )
{
    // first check that min is less than max (or they are the same - this is allowed)
    if ( minVal > maxVal ) {
        cerr << "ERROR in LauParameter::checkRange : minValue: " << minVal
             << " greater than maxValue: " << maxVal << endl;
        if ( minValue_ > maxValue_ ) {
            minValue_ = maxValue_;
            cerr << "                                  : Setting both to " << maxValue_ << "."
                 << endl;
        } else {
            cerr << "                                  : Not changing anything." << endl;
        }
        return;
    }

    minValue_ = minVal;
    maxValue_ = maxVal;

    // now check that the value is still within the range
    if ( ( val < minVal ) || ( val > maxVal ) ) {
        if ( name_ != "" ) {
            cerr << "ERROR in LauParameter::checkRange : value: " << val << " not in allowed range ["
                 << minVal << " , " << maxVal << "] for parameter \"" << name_ << "\"." << endl;
        } else {
            cerr << "ERROR in LauParameter::checkRange : value: " << val
                 << " not in allowed range [" << minVal << " , " << maxVal << "]." << endl;
        }
        if ( val < minVal ) {
            cerr << "                                  : Setting value to minValue: " << minVal
                 << endl;
            val = minVal;
        } else {
            cerr << "                                  : Setting value to maxValue: " << maxVal
                 << endl;
            val = maxVal;
        }
    }
    value_ = val;
}

LauParameter* LauParameter::createClone( Double_t constFactor )
{
    // if we're a clone we mustn't be cloned ourselves
    // but instead return another clone of our parent
    // (this is so that the parent knows of all its clones)
    if ( this->clone() ) {
        LauParameter* clonePar = parent_->createClone( constFactor );
        clonePar->name( this->name() );
        return clonePar;
    }

    // clone ourselves using the copy-constructor
    LauParameter* clonePar = new LauParameter( *this );
    Double_t newValue      = clonePar->value() * constFactor;
    clonePar->value( newValue );
    clonePar->wipeClones();
    clonePar->clone( this );
    clones_.insert( std::make_pair( clonePar, constFactor ) );
    return clonePar;
}

LauParameter* LauParameter::createClone( const TString& newName, Double_t constFactor )
{
    // self message to create the clone
    LauParameter* clonePar = this->createClone( constFactor );
    // set the new name
    clonePar->name( newName );
    // and return
    return clonePar;
}

void LauParameter::updateClones( Bool_t justValue )
{
    // if we don't have any clones then there's nothing to do
    if ( clones_.empty() ) {
        return;
    }

    // we have to set the values directly rather than using member functions because otherwise we'd get into an infinite loop
    if ( justValue ) {
        for ( auto& [clonePar, constFactor] : clones_ ) {
            clonePar->value_ = constFactor * value_;
        }
    } else {
        for ( auto& [clonePar, constFactor] : clones_ ) {
            clonePar->value_              = constFactor * value_;
            clonePar->error_              = constFactor * error_;
            clonePar->negError_           = constFactor * negError_;
            clonePar->posError_           = constFactor * posError_;
            clonePar->genValue_           = constFactor * genValue_;
            clonePar->initValue_          = constFactor * initValue_;
            clonePar->minValue_           = constFactor * minValue_;
            clonePar->maxValue_           = constFactor * maxValue_;
            clonePar->fixed_              = fixed_;
            clonePar->secondStage_        = secondStage_;
            clonePar->gaussConstraint_    = gaussConstraint_;
            clonePar->constraintTrueMean_ = constraintTrueMean_;
            clonePar->constraintMean_     = constraintMean_;
            clonePar->constraintWidth_    = constraintWidth_;
            clonePar->gcc_                = gcc_;
            clonePar->bias_               = bias_;
            clonePar->pull_               = pull_;
        }
    }
}

void LauParameter::randomiseValue()
{
    this->randomiseValue( this->minValue(), this->maxValue() );
}

void LauParameter::randomiseValue( Double_t minVal, Double_t maxVal )
{
    // if we're fixed then do nothing
    if ( this->fixed() ) {
        return;
    }

    // check supplied values are sensible
    if ( maxVal < minVal ) {
        cerr << "ERROR in LauParameter::randomiseValue : Supplied maximum value smaller than minimum value."
             << endl;
        return;
    }
    if ( maxVal > this->maxValue() ) {
        maxVal = this->maxValue();
    }
    if ( minVal < this->minValue() ) {
        minVal = this->minValue();
    }

    // use the zero-seed random number generator to get values that are
    // uniformly distributed over the given range
    Double_t randNo = LauRandom::zeroSeedRandom()->Rndm();
    Double_t val    = randNo * ( maxVal - minVal ) + minVal;
    this->initValue( val );
}

// ostream operator
std::ostream& operator<<( std::ostream& stream, const LauParameter& par )
{
    stream << par.value();
    return stream;
}
