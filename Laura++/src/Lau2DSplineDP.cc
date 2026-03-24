
/*
Copyright 2013 University of Warwick

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

/*! \file Lau2DSplineDP.cc
    \brief File containing implementation of Lau2DSplineDP class.
*/

#include "Lau2DSplineDP.hh"

#include "Lau2DCubicSpline.hh"
#include "LauDaughters.hh"
#include "LauKinematics.hh"
#include "LauRandom.hh"

#include "TAxis.h"
#include "TH2.h"
#include "TRandom.h"
#include "TSystem.h"

#include <iostream>

Lau2DSplineDP::Lau2DSplineDP( const TH2* hist,
                              const LauDaughters* daughters,
                              Bool_t fluctuateBins,
                              Double_t avEff,
                              Double_t avEffError,
                              Bool_t useUpperHalfOnly,
                              Bool_t squareDP ) :
    Lau2DAbsHistDP( daughters, useUpperHalfOnly, squareDP ),
    spline_( 0 )
{
    //We may need to modify the histogram so clone it
    TH2* tempHist( hist ? dynamic_cast<TH2*>( hist->Clone() ) : 0 );

    if ( ! tempHist ) {
        std::cerr << "ERROR in Lau2DSplineDP constructor : the histogram pointer is null."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    if ( fluctuateBins ) {
        this->doBinFluctuation( tempHist );
    }
    if ( avEff > 0.0 && avEffError > 0.0 ) {
        this->raiseOrLowerBins( tempHist, avEff, avEffError );
    }

    spline_ = new Lau2DCubicSpline( *tempHist );

    delete tempHist;
}

Lau2DSplineDP::Lau2DSplineDP( const TH2* hist,
                              const TH2* errorHi,
                              const TH2* errorLo,
                              const LauDaughters* daughters,
                              Bool_t fluctuateBins,
                              Double_t avEff,
                              Double_t avEffError,
                              Bool_t useUpperHalfOnly,
                              Bool_t squareDP ) :
    Lau2DAbsHistDP( daughters, useUpperHalfOnly, squareDP ),
    spline_( 0 )
{
    //We may need to modify the histogram so clone it
    TH2* tempHist( hist ? dynamic_cast<TH2*>( hist->Clone() ) : 0 );
    TH2* tempErrorHi( errorHi ? dynamic_cast<TH2*>( errorHi->Clone() ) : 0 );
    TH2* tempErrorLo( errorLo ? dynamic_cast<TH2*>( errorLo->Clone() ) : 0 );

    if ( ! tempHist ) {
        std::cerr << "ERROR in Lau2DSplineDP constructor : the histogram pointer is null."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }
    if ( ! tempErrorHi ) {
        std::cerr << "ERROR in Lau2DHistDP constructor : the upper error histogram pointer is null."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }
    if ( ! tempErrorLo ) {
        std::cerr << "ERROR in Lau2DHistDP constructor : the lower error histogram pointer is null."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    TAxis* xAxis  = tempHist->GetXaxis();
    Double_t minX = static_cast<Double_t>( xAxis->GetXmin() );
    Double_t maxX = static_cast<Double_t>( xAxis->GetXmax() );

    TAxis* yAxis  = tempHist->GetYaxis();
    Double_t minY = static_cast<Double_t>( yAxis->GetXmin() );
    Double_t maxY = static_cast<Double_t>( yAxis->GetXmax() );

    Int_t nBinsX = static_cast<Int_t>( tempHist->GetNbinsX() );
    Int_t nBinsY = static_cast<Int_t>( tempHist->GetNbinsY() );

    if ( static_cast<Int_t>( tempErrorLo->GetNbinsX() ) != nBinsX ||
         static_cast<Int_t>( tempErrorLo->GetNbinsY() ) != nBinsY ) {
        std::cerr << "ERROR in Lau2DHistDP constructor : the lower error histogram has a different number of bins to the main histogram."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    if ( static_cast<Int_t>( tempErrorHi->GetNbinsX() ) != nBinsX ||
         static_cast<Int_t>( tempErrorHi->GetNbinsY() ) != nBinsY ) {
        std::cerr << "ERROR in Lau2DHistDP constructor : the upper error histogram has a different number of bins to the main histogram."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    xAxis = tempErrorLo->GetXaxis();
    yAxis = tempErrorLo->GetYaxis();

    if ( static_cast<Double_t>( xAxis->GetXmin() ) != minX ||
         static_cast<Double_t>( xAxis->GetXmax() ) != maxX ) {
        std::cerr << "ERROR in Lau2DHistDP constructor : the lower error histogram has a different x range to the main histogram."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    if ( static_cast<Double_t>( yAxis->GetXmin() ) != minY ||
         static_cast<Double_t>( yAxis->GetXmax() ) != maxY ) {
        std::cerr << "ERROR in Lau2DHistDP constructor : the lower error histogram has a different y range to the main histogram."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    xAxis = tempErrorHi->GetXaxis();
    yAxis = tempErrorHi->GetYaxis();

    if ( static_cast<Double_t>( xAxis->GetXmin() ) != minX ||
         static_cast<Double_t>( xAxis->GetXmax() ) != maxX ) {
        std::cerr << "ERROR in Lau2DHistDP constructor : the upper error histogram has a different x range to the main histogram."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    if ( static_cast<Double_t>( yAxis->GetXmin() ) != minY ||
         static_cast<Double_t>( yAxis->GetXmax() ) != maxY ) {
        std::cerr << "ERROR in Lau2DHistDP constructor : the upper error histogram has a different y range to the main histogram."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    if ( fluctuateBins ) {
        this->doBinFluctuation( tempHist, tempErrorHi, tempErrorLo );
    }
    if ( avEff > 0.0 && avEffError > 0.0 ) {
        this->raiseOrLowerBins( tempHist, avEff, avEffError );
    }

    spline_ = new Lau2DCubicSpline( *tempHist );

    delete tempHist;
    delete tempErrorHi;
    delete tempErrorLo;
}

Lau2DSplineDP::~Lau2DSplineDP()
{
    delete spline_;
    spline_ = 0;
}

Double_t Lau2DSplineDP::interpolateXY( Double_t x, Double_t y ) const
{
    // This function returns the interpolated value of the histogram function
    // for the given values of x and y by finding the adjacent bins and extrapolating
    // using weights based on the inverse distance of the point from the adajcent
    // bin centres.
    // Here, x = m13^2, y = m23^2, or m', theta' for square DP co-ordinates

    // If we're only using one half then flip co-ordinates
    // appropriately for conventional or square DP
    this->getUpperHalf( x, y );

    // First ask whether the point is inside the kinematic region.
    if ( this->withinDPBoundaries( x, y ) == kFALSE ) {
        std::cerr << "WARNING in Lau2DSplineDP::interpolateXY : Given position is outside the DP boundary, returning 0.0."
                  << std::endl;
        return 0.0;
    }

    return spline_->evaluate( x, y );
}
