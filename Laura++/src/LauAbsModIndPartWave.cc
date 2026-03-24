
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

/*! \file LauAbsModIndPartWave.cc
    \brief File containing implementation of LauAbsModIndPartWave class.
*/

#include "LauAbsModIndPartWave.hh"

#include "LauConstants.hh"
#include "LauKinematics.hh"
#include "LauResonanceInfo.hh"

#include <cstdlib>
#include <iostream>

LauAbsModIndPartWave::LauAbsModIndPartWave( LauResonanceInfo* resInfo,
                                            Int_t resPairAmpInt,
                                            const LauDaughters* daughters ) :
    LauAbsResonance( resInfo, resPairAmpInt, daughters ),
    nKnots_( 0 ),
    spline1_( 0 ),
    spline2_( 0 ),
    type1_( Lau1DCubicSpline::StandardSpline ),
    type2_( Lau1DCubicSpline::StandardSpline ),
    leftBound1_( Lau1DCubicSpline::NotAKnot ),
    rightBound1_( Lau1DCubicSpline::NotAKnot ),
    leftBound2_( Lau1DCubicSpline::NotAKnot ),
    rightBound2_( Lau1DCubicSpline::NotAKnot ),
    leftGrad1_( 0. ),
    rightGrad1_( 0. ),
    leftGrad2_( 0. ),
    rightGrad2_( 0. ),
    secondStage_( kFALSE )
{
}

LauAbsModIndPartWave::~LauAbsModIndPartWave()
{
    delete spline1_;
    delete spline2_;
}

void LauAbsModIndPartWave::floatKnotsSecondStage( const Bool_t secondStage )
{
    secondStage_ = secondStage;

    // if the parameters have not yet been created we can now just return
    if ( amp1Pars_.size() != nKnots_ ) {
        return;
    }

    // otherwise we need to toggle their second-stage parameter
    for ( UInt_t i( 0 ); i < nKnots_; ++i ) {
        amp1Pars_[i]->secondStage( secondStage_ );
        amp2Pars_[i]->secondStage( secondStage_ );
    }
}

std::set<Double_t> LauAbsModIndPartWave::checkKnots( const std::set<Double_t>& masses )
{
    std::set<Double_t> knots = masses;

    const std::set<Double_t>::const_iterator first        = knots.begin();
    const std::set<Double_t>::const_reverse_iterator last = knots.rbegin();

    const Double_t lower_limit = this->getMassDaug1() + this->getMassDaug2();
    const Double_t upper_limit = this->getMassParent() - this->getMassBachelor();

    // check whether we have been given knots at unphysical masses
    if ( *first < lower_limit ) {
        std::cerr << "WARNING in LauAbsModIndPartWave::checkKnots : Knot found at mass " << *first
                  << " is below the lower kinematic limit." << std::endl;
        std::cerr << "                                            : Lower kinematic limit is at mass "
                  << lower_limit << std::endl;
        std::cerr << "                                            : Aborting definition of knot positions."
                  << std::endl;
        knots.clear();
        return knots;
    }
    if ( *last > upper_limit ) {
        std::cerr << "WARNING in LauAbsModIndPartWave::checkKnots : Knot found at mass " << *last
                  << " is above the upper kinematic limit." << std::endl;
        std::cerr << "                                            : Upper kinematic limit is at mass "
                  << upper_limit << std::endl;
        std::cerr << "                                            : Aborting definition of knot positions."
                  << std::endl;
        knots.clear();
        return knots;
    }

    // check if we have knots at each extreme - if not, add them in
    if ( ( *first ) != lower_limit ) {
        knots.insert( lower_limit );
    }
    if ( ( *last ) != upper_limit ) {
        knots.insert( upper_limit );
    }

    return knots;
}

void LauAbsModIndPartWave::defineKnots( const std::set<Double_t>& masses )
{
    if ( ! masses_.empty() ) {
        std::cerr << "WARNING in LauAbsModIndPartWave::defineKnots : Knot positions have already been defined, not making any changes."
                  << std::endl;
        return;
    }

    const std::set<Double_t> knots = this->checkKnots( masses );

    nKnots_ = knots.size();
    if ( nKnots_ == 0 ) {
        return;
    }

    masses_.reserve( nKnots_ );
    amp1Vals_.reserve( nKnots_ );
    amp2Vals_.reserve( nKnots_ );
    amp1Pars_.reserve( nKnots_ );
    amp2Pars_.reserve( nKnots_ );

    UInt_t counter( 0 );
    for ( std::set<Double_t>::const_iterator iter = knots.begin(); iter != knots.end(); ++iter ) {
        masses_.push_back( *iter );
        amp1Vals_.push_back( 1.0 );
        amp2Vals_.push_back( 1.0 );
        this->createAmpParameters( counter );
        ++counter;
    }

    for ( std::vector<Double_t>::const_iterator iter = masses_.begin(); iter != masses_.end();
          ++iter ) {
        std::cout << "INFO in LauAbsModIndPartWave::defineKnots : Knot added to resonance "
                  << this->getResonanceName() << " at mass " << *iter << std::endl;
    }
}

void LauAbsModIndPartWave::initialise()
{
    if ( spline1_ != 0 ) {
        delete spline1_;
        spline1_ = 0;
    }
    if ( spline2_ != 0 ) {
        delete spline2_;
        spline2_ = 0;
    }

    for ( UInt_t i( 0 ); i < nKnots_; ++i ) {
        amp1Vals_[i] = amp1Pars_[i]->unblindValue();
        amp2Vals_[i] = amp2Pars_[i]->unblindValue();
    }

    spline1_ = new Lau1DCubicSpline( masses_,
                                     amp1Vals_,
                                     type1_,
                                     leftBound1_,
                                     rightBound1_,
                                     leftGrad1_,
                                     rightGrad1_ );
    spline2_ = new Lau1DCubicSpline( masses_,
                                     amp2Vals_,
                                     type2_,
                                     leftBound2_,
                                     rightBound2_,
                                     leftGrad2_,
                                     rightGrad2_ );
}

LauComplex LauAbsModIndPartWave::resAmp( Double_t mass, Double_t spinTerm )
{
    amp_.zero();

    Bool_t paramChanged1( kFALSE ), paramChanged2( kFALSE );

    for ( UInt_t i( 0 ); i < nKnots_; ++i ) {
        if ( ! amp1Pars_[i]->fixed() && amp1Pars_[i]->unblindValue() != amp1Vals_[i] ) {
            paramChanged1 = kTRUE;
            amp1Vals_[i]  = amp1Pars_[i]->unblindValue();
        }
        if ( ! amp2Pars_[i]->fixed() && amp2Pars_[i]->unblindValue() != amp2Vals_[i] ) {
            paramChanged2 = kTRUE;
            amp2Vals_[i]  = amp2Pars_[i]->unblindValue();
        }
    }

    if ( spline1_ == 0 || spline2_ == 0 ) {
        std::cerr << "ERROR in LauAbsModIndPartWave::resAmp : One of the splines is null"
                  << std::endl;
        return amp_;
    }

    if ( paramChanged1 ) {
        spline1_->updateYValues( amp1Vals_ );
    }

    if ( paramChanged2 ) {
        spline2_->updateYValues( amp2Vals_ );
    }

    this->evaluateAmplitude( mass );

    amp_.rescale( spinTerm );

    return amp_;
}

void LauAbsModIndPartWave::setSplineType( Lau1DCubicSpline::LauSplineType type1,
                                          Lau1DCubicSpline::LauSplineType type2 )
{
    type1_ = type1;
    type2_ = type2;
}

void LauAbsModIndPartWave::setSplineBoundaryConditions(
    Lau1DCubicSpline::LauSplineBoundaryType leftBound1,
    Lau1DCubicSpline::LauSplineBoundaryType rightBound1,
    Lau1DCubicSpline::LauSplineBoundaryType leftBound2,
    Lau1DCubicSpline::LauSplineBoundaryType rightBound2,
    Double_t leftGrad1,
    Double_t rightGrad1,
    Double_t leftGrad2,
    Double_t rightGrad2 )
{
    leftBound1_  = leftBound1;
    rightBound1_ = rightBound1;
    leftBound2_  = leftBound2;
    rightBound2_ = rightBound2;
    leftGrad1_   = leftGrad1;
    rightGrad1_  = rightGrad1;
    leftGrad2_   = leftGrad2;
    rightGrad2_  = rightGrad2;
}

const std::vector<LauParameter*>& LauAbsModIndPartWave::getFloatingParameters()
{
    this->clearFloatingParameters();

    for ( UInt_t i( 0 ); i < nKnots_; ++i ) {
        if ( ! amp1Pars_[i]->fixed() ) {
            this->addFloatingParameter( amp1Pars_[i] );
        }
        if ( ! amp2Pars_[i]->fixed() ) {
            this->addFloatingParameter( amp2Pars_[i] );
        }
    }

    return this->getParameters();
}
