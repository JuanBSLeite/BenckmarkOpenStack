
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

/*! \file LauEffModel.cc
    \brief File containing implementation of LauEffModel class.
*/

#include "LauEffModel.hh"

#include "Lau2DHistDP.hh"
#include "Lau2DSplineDP.hh"
#include "LauDaughters.hh"
#include "LauKinematics.hh"
#include "LauVetoes.hh"

#include "TSystem.h"

#include <cstdlib>
#include <iostream>

LauEffModel::LauEffModel( const LauDaughters* daughters, const LauVetoes* vetoes ) :
    daughters_( daughters ),
    vetoes_( vetoes ),
    effHisto_( 0 ),
    fluctuateEffHisto_( kFALSE ),
    lowBinWarningIssued_( kFALSE ),
    highBinWarningIssued_( kFALSE )
{
    if ( daughters_ == 0 ) {
        std::cerr << "ERROR in LauEffModel Constructor : invalid pointer to daughters object supplied."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }
}

LauEffModel::~LauEffModel()
{
    std::vector<Lau2DAbsDP*>::iterator it  = effHisto_.begin();
    std::vector<Lau2DAbsDP*>::iterator end = effHisto_.end();
    for ( ; it != end; ++it ) {
        delete *it;
        *it = 0;
    }
    effHisto_.clear();
}

void LauEffModel::setEffHisto( const TH2* effHisto,
                               Bool_t useInterpolation,
                               Bool_t fluctuateBins,
                               Double_t avEff,
                               Double_t absError,
                               Bool_t useUpperHalfOnly,
                               Bool_t squareDP )
{
    // Set the efficiency across the Dalitz plot using a predetermined 2D histogram
    // with x = m_13^2, y = m_23^2.
    Bool_t upperHalf( kFALSE );
    if ( daughters_->gotSymmetricalDP() && useUpperHalfOnly == kTRUE ) {
        upperHalf = kTRUE;
    }
    std::cout << "INFO in LauEffModel::setEffHisto : Efficiency histogram has upperHalf = "
              << static_cast<Int_t>( upperHalf ) << std::endl;

    std::vector<Lau2DAbsDP*>::iterator it  = effHisto_.begin();
    std::vector<Lau2DAbsDP*>::iterator end = effHisto_.end();
    for ( ; it != end; ++it ) {
        delete *it;
        *it = 0;
    }
    effHisto_.clear();

    // Copy the histogram.
    effHisto_.push_back( new Lau2DHistDP( effHisto,
                                          daughters_,
                                          useInterpolation,
                                          fluctuateBins,
                                          avEff,
                                          absError,
                                          upperHalf,
                                          squareDP ) );
    fluctuateEffHisto_ = fluctuateBins;

    if ( avEff > 0.0 && absError > 0.0 ) {
        fluctuateEffHisto_ = kTRUE;
    }
}

void LauEffModel::setEffHisto( const TH2* effHisto,
                               const TH2* errorHi,
                               const TH2* errorLo,
                               Bool_t useInterpolation,
                               Bool_t fluctuateBins,
                               Double_t avEff,
                               Double_t absError,
                               Bool_t useUpperHalfOnly,
                               Bool_t squareDP )
{
    // Set the efficiency across the Dalitz plot using a predetermined 2D histogram
    // with x = m_13^2, y = m_23^2.
    Bool_t upperHalf( kFALSE );
    if ( daughters_->gotSymmetricalDP() && useUpperHalfOnly == kTRUE ) {
        upperHalf = kTRUE;
    }
    std::cout << "INFO in LauEffModel::setEffHisto : Efficiency histogram has upperHalf = "
              << static_cast<Int_t>( upperHalf ) << std::endl;

    std::vector<Lau2DAbsDP*>::iterator it  = effHisto_.begin();
    std::vector<Lau2DAbsDP*>::iterator end = effHisto_.end();
    for ( ; it != end; ++it ) {
        delete *it;
        *it = 0;
    }
    effHisto_.clear();

    // Copy the histogram.
    effHisto_.push_back( new Lau2DHistDP( effHisto,
                                          errorHi,
                                          errorLo,
                                          daughters_,
                                          useInterpolation,
                                          fluctuateBins,
                                          avEff,
                                          absError,
                                          upperHalf,
                                          squareDP ) );
    fluctuateEffHisto_ = fluctuateBins;

    if ( avEff > 0.0 && absError > 0.0 ) {
        fluctuateEffHisto_ = kTRUE;
    }
}

void LauEffModel::setEffSpline( const TH2* effHisto,
                                Bool_t fluctuateBins,
                                Double_t avEff,
                                Double_t absError,
                                Bool_t useUpperHalfOnly,
                                Bool_t squareDP )
{
    // Set the efficiency across the Dalitz plot using a predetermined 2D histogram
    // with x = m_13^2, y = m_23^2.
    Bool_t upperHalf( kFALSE );
    if ( daughters_->gotSymmetricalDP() && useUpperHalfOnly == kTRUE ) {
        upperHalf = kTRUE;
    }
    std::cout << "INFO in LauEffModel::setEffSpline : Efficiency histogram has upperHalf = "
              << static_cast<Int_t>( upperHalf ) << std::endl;

    std::vector<Lau2DAbsDP*>::iterator it  = effHisto_.begin();
    std::vector<Lau2DAbsDP*>::iterator end = effHisto_.end();
    for ( ; it != end; ++it ) {
        delete *it;
        *it = 0;
    }
    effHisto_.clear();

    // Copy the histogram.
    effHisto_.push_back(
        new Lau2DSplineDP( effHisto, daughters_, fluctuateBins, avEff, absError, upperHalf, squareDP ) );
    fluctuateEffHisto_ = fluctuateBins;

    if ( avEff > 0.0 && absError > 0.0 ) {
        fluctuateEffHisto_ = kTRUE;
    }
}

void LauEffModel::setEffSpline( const TH2* effHisto,
                                const TH2* errorHi,
                                const TH2* errorLo,
                                Bool_t fluctuateBins,
                                Double_t avEff,
                                Double_t absError,
                                Bool_t useUpperHalfOnly,
                                Bool_t squareDP )
{
    // Set the efficiency across the Dalitz plot using a predetermined 2D histogram
    // with x = m_13^2, y = m_23^2.
    Bool_t upperHalf( kFALSE );
    if ( daughters_->gotSymmetricalDP() && useUpperHalfOnly == kTRUE ) {
        upperHalf = kTRUE;
    }
    std::cout << "INFO in LauEffModel::setEffSpline : Efficiency histogram has upperHalf = "
              << static_cast<Int_t>( upperHalf ) << std::endl;

    std::vector<Lau2DAbsDP*>::iterator it  = effHisto_.begin();
    std::vector<Lau2DAbsDP*>::iterator end = effHisto_.end();
    for ( ; it != end; ++it ) {
        delete *it;
        *it = 0;
    }
    effHisto_.clear();

    // Copy the histogram.
    effHisto_.push_back( new Lau2DSplineDP( effHisto,
                                            errorHi,
                                            errorLo,
                                            daughters_,
                                            fluctuateBins,
                                            avEff,
                                            absError,
                                            upperHalf,
                                            squareDP ) );
    fluctuateEffHisto_ = fluctuateBins;

    if ( avEff > 0.0 && absError > 0.0 ) {
        fluctuateEffHisto_ = kTRUE;
    }
}

void LauEffModel::addEffHisto( const TH2* effHisto,
                               Bool_t useInterpolation,
                               Double_t avEff,
                               Double_t absError,
                               Bool_t useUpperHalfOnly,
                               Bool_t squareDP )
{
    // Set the efficiency across the Dalitz plot using a predetermined 2D histogram
    // with x = m_13^2, y = m_23^2.
    Bool_t upperHalf( kFALSE );
    if ( daughters_->gotSymmetricalDP() && useUpperHalfOnly == kTRUE ) {
        upperHalf = kTRUE;
    }
    std::cout << "INFO in LauEffModel::addEffHisto : Efficiency histogram has upperHalf = "
              << static_cast<Int_t>( upperHalf ) << std::endl;

    // Copy the histogram.
    effHisto_.push_back( new Lau2DHistDP( effHisto,
                                          daughters_,
                                          useInterpolation,
                                          fluctuateEffHisto_,
                                          avEff,
                                          absError,
                                          upperHalf,
                                          squareDP ) );
}

void LauEffModel::addEffHisto( const TH2* effHisto,
                               const TH2* errorHi,
                               const TH2* errorLo,
                               Bool_t useInterpolation,
                               Double_t avEff,
                               Double_t absError,
                               Bool_t useUpperHalfOnly,
                               Bool_t squareDP )
{
    // Set the efficiency across the Dalitz plot using a predetermined 2D histogram
    // with x = m_13^2, y = m_23^2.
    Bool_t upperHalf( kFALSE );
    if ( daughters_->gotSymmetricalDP() && useUpperHalfOnly == kTRUE ) {
        upperHalf = kTRUE;
    }
    std::cout << "INFO in LauEffModel::addEffHisto : Efficiency histogram has upperHalf = "
              << static_cast<Int_t>( upperHalf ) << std::endl;

    // Copy the histogram.
    effHisto_.push_back( new Lau2DHistDP( effHisto,
                                          errorHi,
                                          errorLo,
                                          daughters_,
                                          useInterpolation,
                                          fluctuateEffHisto_,
                                          avEff,
                                          absError,
                                          upperHalf,
                                          squareDP ) );
}

void LauEffModel::addEffSpline( const TH2* effHisto,
                                Double_t avEff,
                                Double_t absError,
                                Bool_t useUpperHalfOnly,
                                Bool_t squareDP )
{
    // Set the efficiency across the Dalitz plot using a predetermined 2D histogram
    // with x = m_13^2, y = m_23^2.
    Bool_t upperHalf( kFALSE );
    if ( daughters_->gotSymmetricalDP() && useUpperHalfOnly == kTRUE ) {
        upperHalf = kTRUE;
    }
    std::cout << "INFO in LauEffModel::addEffSpline : Efficiency histogram has upperHalf = "
              << static_cast<Int_t>( upperHalf ) << std::endl;

    // Copy the histogram.
    effHisto_.push_back( new Lau2DSplineDP( effHisto,
                                            daughters_,
                                            fluctuateEffHisto_,
                                            avEff,
                                            absError,
                                            upperHalf,
                                            squareDP ) );
}

void LauEffModel::addEffSpline( const TH2* effHisto,
                                const TH2* errorHi,
                                const TH2* errorLo,
                                Double_t avEff,
                                Double_t absError,
                                Bool_t useUpperHalfOnly,
                                Bool_t squareDP )
{
    // Set the efficiency across the Dalitz plot using a predetermined 2D histogram
    // with x = m_13^2, y = m_23^2.
    Bool_t upperHalf( kFALSE );
    if ( daughters_->gotSymmetricalDP() && useUpperHalfOnly == kTRUE ) {
        upperHalf = kTRUE;
    }
    std::cout << "INFO in LauEffModel::addEffSpline : Efficiency histogram has upperHalf = "
              << static_cast<Int_t>( upperHalf ) << std::endl;

    // Copy the histogram.
    effHisto_.push_back( new Lau2DSplineDP( effHisto,
                                            errorHi,
                                            errorLo,
                                            daughters_,
                                            fluctuateEffHisto_,
                                            avEff,
                                            absError,
                                            upperHalf,
                                            squareDP ) );
}

Double_t LauEffModel::getEffHistValue( const LauKinematics* kinematics ) const
{
    // Get the efficiency from the 2D histo.
    Double_t eff( 1.0 );

    Double_t xVal( 0.0 );
    Double_t yVal( 0.0 );

    std::vector<Lau2DAbsDP*>::const_iterator it        = effHisto_.begin();
    const std::vector<Lau2DAbsDP*>::const_iterator end = effHisto_.end();
    for ( ; it != end; ++it ) {
        if ( ( *it )->usingSquareDP() ) {
            xVal = kinematics->getmPrime();
            yVal = kinematics->getThetaPrime();
        } else {
            xVal = kinematics->getm13Sq();
            yVal = kinematics->getm23Sq();
        }
        eff *= ( *it )->interpolateXY( xVal, yVal );
    }

    return eff;
}

Double_t LauEffModel::calcEfficiency( const LauKinematics* kinematics ) const
{
    // Routine to calculate the efficiency for the given event/point in
    // the Dalitz plot. This routine uses the 2D histogram set by the
    // setEffHisto() function.
    Double_t eff( 1.0 );

    // Check for vetoes
    Bool_t vetoOK( kTRUE );
    if ( vetoes_ != 0 ) {
        vetoOK = vetoes_->passVeto( kinematics );
    }

    if ( vetoOK == kFALSE ) {
        // We failed the veto.
        eff = 0.0;
    } else {
        // We are using a 2D histogram representation of the efficiency across the Dalitz plot.
        eff = this->getEffHistValue( kinematics );
    }

    // Check that the efficiency is in the allowed range (0-1)
    // If we're using a spline then out-of-range efficiencies can be caused by adjacent bins that all contain a value of either zero or one.
    // The spline requires the efficiency, its first derivatives and the mixed second derivative to be continuous and to match the input histogram
    // at the bin centres. Derivatives are calculated using a finite difference approximation taking the difference between the neighbouring bins.
    // If two bins are zero but the third is not then the second bin will have a positive first derivative causing the spline to dip below zero
    // between the two zero bins to remain smooth. The analogous case with adjacent maximised bins will cause peaks above one. Such dips are
    // unavoidable but are correctly removed here.
    if ( eff < 0.0 ) {
        if ( ! lowBinWarningIssued_ ) {
            std::cerr << "WARNING in LauEffModel::calcEfficiency : Efficiency " << eff
                      << " is less than 0 - setting to 0.  You may want to check your histogram!"
                      << std::endl
                      << "                                       : If you are using a spline then this could be caused by adjacent empty bins. Further warnings will be suppressed."
                      << std::endl;
            lowBinWarningIssued_ = kTRUE;
        }
        eff = 0.0;
    } else if ( eff > 1.0 ) {
        if ( ! highBinWarningIssued_ ) {
            std::cerr << "WARNING in LauEffModel::calcEfficiency : Efficiency " << eff
                      << " is greater than 1 - setting to 1.  You may want to check your histogram!"
                      << std::endl
                      << "                                       : If you are using a spline then this could be caused by adjacent full bins. Further warnings will be suppressed."
                      << std::endl;
            highBinWarningIssued_ = kTRUE;
        }
        eff = 1.0;
    }

    return eff;
}

Bool_t LauEffModel::passVeto( const LauKinematics* kinematics ) const
{
    Bool_t pass = kTRUE;
    if ( vetoes_ != 0 ) {
        pass = vetoes_->passVeto( kinematics );
    }
    return pass;
}
