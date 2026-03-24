
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

/*! \file Lau2DHistDP.cc
    \brief File containing implementation of Lau2DHistDP class.
*/

#include "Lau2DHistDP.hh"

#include "LauDaughters.hh"
#include "LauKinematics.hh"
#include "LauRandom.hh"

#include "TAxis.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom.h"
#include "TSystem.h"

#include <iostream>

Lau2DHistDP::Lau2DHistDP( const TH2* hist,
                          const LauDaughters* daughters,
                          Bool_t useInterpolation,
                          Bool_t fluctuateBins,
                          Double_t avEff,
                          Double_t avEffError,
                          Bool_t useUpperHalfOnly,
                          Bool_t squareDP ) :
    Lau2DAbsHistDP( daughters, useUpperHalfOnly, squareDP ),
    hist_( hist ? dynamic_cast<TH2*>( hist->Clone() ) : 0 ),
    errorHi_( 0 ),
    errorLo_( 0 ),
    minX_( 0.0 ),
    maxX_( 0.0 ),
    minY_( 0.0 ),
    maxY_( 0.0 ),
    rangeX_( 0.0 ),
    rangeY_( 0.0 ),
    binXWidth_( 0.0 ),
    binYWidth_( 0.0 ),
    nBinsX_( 0 ),
    nBinsY_( 0 ),
    useInterpolation_( useInterpolation )
{
    if ( ! hist_ ) {
        std::cerr << "ERROR in Lau2DHistDP constructor : the histogram pointer is null." << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Save various attributes of the histogram
    // (axis ranges, number of bins, areas)
    TAxis* xAxis = hist_->GetXaxis();
    minX_        = static_cast<Double_t>( xAxis->GetXmin() );
    maxX_        = static_cast<Double_t>( xAxis->GetXmax() );
    rangeX_      = maxX_ - minX_;

    TAxis* yAxis = hist_->GetYaxis();
    minY_        = static_cast<Double_t>( yAxis->GetXmin() );
    maxY_        = static_cast<Double_t>( yAxis->GetXmax() );
    rangeY_      = maxY_ - minY_;

    nBinsX_ = static_cast<Int_t>( hist_->GetNbinsX() );
    nBinsY_ = static_cast<Int_t>( hist_->GetNbinsY() );

    binXWidth_ = static_cast<Double_t>( TMath::Abs( rangeX_ ) / ( nBinsX_ * 1.0 ) );
    binYWidth_ = static_cast<Double_t>( TMath::Abs( rangeY_ ) / ( nBinsY_ * 1.0 ) );

    if ( fluctuateBins ) {
        this->doBinFluctuation( hist_ );
    }
    if ( avEff > 0.0 && avEffError > 0.0 ) {
        this->raiseOrLowerBins( hist_, avEff, avEffError );
    }
}

Lau2DHistDP::Lau2DHistDP( const TH2* hist,
                          const TH2* errorHi,
                          const TH2* errorLo,
                          const LauDaughters* daughters,
                          Bool_t useInterpolation,
                          Bool_t fluctuateBins,
                          Double_t avEff,
                          Double_t avEffError,
                          Bool_t useUpperHalfOnly,
                          Bool_t squareDP ) :
    Lau2DAbsHistDP( daughters, useUpperHalfOnly, squareDP ),
    hist_( hist ? dynamic_cast<TH2*>( hist->Clone() ) : 0 ),
    errorHi_( errorHi ? dynamic_cast<TH2*>( errorHi->Clone() ) : 0 ),
    errorLo_( errorLo ? dynamic_cast<TH2*>( errorLo->Clone() ) : 0 ),
    minX_( 0.0 ),
    maxX_( 0.0 ),
    minY_( 0.0 ),
    maxY_( 0.0 ),
    rangeX_( 0.0 ),
    rangeY_( 0.0 ),
    binXWidth_( 0.0 ),
    binYWidth_( 0.0 ),
    nBinsX_( 0 ),
    nBinsY_( 0 ),
    useInterpolation_( useInterpolation )
{
    if ( ! hist_ ) {
        std::cerr << "ERROR in Lau2DHistDP constructor : the histogram pointer is null." << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }
    if ( ! errorHi_ ) {
        std::cerr << "ERROR in Lau2DHistDP constructor : the upper error histogram pointer is null."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }
    if ( ! errorLo_ ) {
        std::cerr << "ERROR in Lau2DHistDP constructor : the lower error histogram pointer is null."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Save various attributes of the histogram
    // (axis ranges, number of bins, areas)
    TAxis* xAxis = hist_->GetXaxis();
    minX_        = static_cast<Double_t>( xAxis->GetXmin() );
    maxX_        = static_cast<Double_t>( xAxis->GetXmax() );
    rangeX_      = maxX_ - minX_;

    TAxis* yAxis = hist_->GetYaxis();
    minY_        = static_cast<Double_t>( yAxis->GetXmin() );
    maxY_        = static_cast<Double_t>( yAxis->GetXmax() );
    rangeY_      = maxY_ - minY_;

    nBinsX_ = static_cast<Int_t>( hist_->GetNbinsX() );
    nBinsY_ = static_cast<Int_t>( hist_->GetNbinsY() );

    binXWidth_ = static_cast<Double_t>( TMath::Abs( rangeX_ ) / ( nBinsX_ * 1.0 ) );
    binYWidth_ = static_cast<Double_t>( TMath::Abs( rangeY_ ) / ( nBinsY_ * 1.0 ) );

    if ( static_cast<Int_t>( errorLo_->GetNbinsX() ) != nBinsX_ ||
         static_cast<Int_t>( errorLo_->GetNbinsY() ) != nBinsY_ ) {
        std::cerr << "ERROR in Lau2DHistDP constructor : the lower error histogram has a different number of bins to the main histogram."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    if ( static_cast<Int_t>( errorHi_->GetNbinsX() ) != nBinsX_ ||
         static_cast<Int_t>( errorHi_->GetNbinsY() ) != nBinsY_ ) {
        std::cerr << "ERROR in Lau2DHistDP constructor : the upper error histogram has a different number of bins to the main histogram."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    xAxis = errorLo_->GetXaxis();
    yAxis = errorLo_->GetYaxis();

    if ( static_cast<Double_t>( xAxis->GetXmin() ) != minX_ ||
         static_cast<Double_t>( xAxis->GetXmax() ) != maxX_ ) {
        std::cerr << "ERROR in Lau2DHistDP constructor : the lower error histogram has a different x range to the main histogram."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    if ( static_cast<Double_t>( yAxis->GetXmin() ) != minY_ ||
         static_cast<Double_t>( yAxis->GetXmax() ) != maxY_ ) {
        std::cerr << "ERROR in Lau2DHistDP constructor : the lower error histogram has a different y range to the main histogram."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    xAxis = errorHi_->GetXaxis();
    yAxis = errorHi_->GetYaxis();

    if ( static_cast<Double_t>( xAxis->GetXmin() ) != minX_ ||
         static_cast<Double_t>( xAxis->GetXmax() ) != maxX_ ) {
        std::cerr << "ERROR in Lau2DHistDP constructor : the upper error histogram has a different x range to the main histogram."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    if ( static_cast<Double_t>( yAxis->GetXmin() ) != minY_ ||
         static_cast<Double_t>( yAxis->GetXmax() ) != maxY_ ) {
        std::cerr << "ERROR in Lau2DHistDP constructor : the upper error histogram has a different y range to the main histogram."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    if ( fluctuateBins ) {
        this->doBinFluctuation( hist_, errorHi_, errorLo_ );
    }
    if ( avEff > 0.0 && avEffError > 0.0 ) {
        this->raiseOrLowerBins( hist_, avEff, avEffError );
    }
}

Lau2DHistDP::~Lau2DHistDP()
{
    delete hist_;
    hist_ = 0;
    if ( errorHi_ ) {
        delete errorHi_;
        errorHi_ = 0;
    }
    if ( errorLo_ ) {
        delete errorLo_;
        errorLo_ = 0;
    }
}

Double_t Lau2DHistDP::getBinHistValue( Int_t xBinNo, Int_t yBinNo ) const
{
    if ( xBinNo < 0 ) {
        xBinNo = 0;
    } else if ( xBinNo >= nBinsX_ ) {
        return 0.0;
    }

    if ( yBinNo < 0 ) {
        yBinNo = 0;
    } else if ( yBinNo >= nBinsY_ ) {
        return 0.0;
    }

    Double_t value = hist_->GetBinContent( xBinNo + 1, yBinNo + 1 );
    return value;
}

Double_t Lau2DHistDP::interpolateXY( Double_t x, Double_t y ) const
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
        std::cerr << "WARNING in Lau2DHistDP::interpolateXY : Given position is outside the DP boundary, returning 0.0."
                  << std::endl;
        return 0.0;
    }

    // Find the 2D histogram bin for x and y
    Int_t i = Int_t( ( x - minX_ ) * nBinsX_ / rangeX_ );
    Int_t j = Int_t( ( y - minY_ ) * nBinsY_ / rangeY_ );

    if ( i >= nBinsX_ ) {
        i = nBinsX_ - 1;
    }
    if ( j >= nBinsY_ ) {
        j = nBinsY_ - 1;
    }

    // Ask whether we want to do the interpolation, since this method is
    // not reliable for low statistics histograms.
    if ( useInterpolation_ == kFALSE ) {
        return this->getBinHistValue( i, j );
    }

    // Find the bin centres (actual co-ordinate positions, not histogram indices)
    Double_t cbinx = Double_t( i + 0.5 ) * rangeX_ / nBinsX_ + minX_;
    Double_t cbiny = Double_t( j + 0.5 ) * rangeY_ / nBinsY_ + minY_;

    // If bin centres are outside kinematic region, do not extrapolate
    if ( this->withinDPBoundaries( cbinx, cbiny ) == kFALSE ) {
        return this->getBinHistValue( i, j );
    }

    // Find the adjacent bins
    Double_t deltax = x - cbinx;
    Double_t deltay = y - cbiny;

    Int_t i_adj( 0 ), j_adj( 0 );
    if ( deltax > 0.0 ) {
        i_adj = i + 1;
    } else {
        i_adj = i - 1;
    }
    if ( deltay > 0.0 ) {
        j_adj = j + 1;
    } else {
        j_adj = j - 1;
    }

    Bool_t isXBoundary( kFALSE ), isYBoundary( kFALSE );

    Double_t value( 0.0 );

    if ( i_adj >= nBinsX_ || i_adj < 0 ) {
        isYBoundary = kTRUE;
    }
    if ( j_adj >= nBinsY_ || j_adj < 0 ) {
        isXBoundary = kTRUE;
    }

    // In the corners, do no interpolation. Use entry in bin (i,j)
    if ( isXBoundary == kTRUE && isYBoundary == kTRUE ) {

        value = this->getBinHistValue( i, j );

    } else if ( isXBoundary == kTRUE && isYBoundary == kFALSE ) {

        // Find the adjacent x bin centre
        Double_t cbinx_adj = Double_t( i_adj + 0.5 ) * rangeX_ / nBinsX_ + minX_;

        if ( this->withinDPBoundaries( cbinx_adj, y ) == kFALSE ) {

            // The adjacent bin is outside the DP range. Don't extrapolate.
            value = this->getBinHistValue( i, j );

        } else {

            Double_t dx0         = TMath::Abs( x - cbinx );
            Double_t dx1         = TMath::Abs( cbinx_adj - x );
            Double_t inter_denom = dx0 + dx1;

            Double_t value1 = this->getBinHistValue( i, j );
            Double_t value2 = this->getBinHistValue( i_adj, j );

            value = ( value1 * dx1 + value2 * dx0 ) / inter_denom;
        }

    } else if ( isYBoundary == kTRUE && isXBoundary == kFALSE ) {

        // Find the adjacent y bin centre
        Double_t cbiny_adj = Double_t( j_adj + 0.5 ) * rangeY_ / nBinsY_ + minY_;

        if ( this->withinDPBoundaries( x, cbiny_adj ) == kFALSE ) {

            // The adjacent bin is outside the DP range. Don't extrapolate.
            value = this->getBinHistValue( i, j );

        } else {

            Double_t dy0         = TMath::Abs( y - cbiny );
            Double_t dy1         = TMath::Abs( cbiny_adj - y );
            Double_t inter_denom = dy0 + dy1;

            Double_t value1 = this->getBinHistValue( i, j );
            Double_t value2 = this->getBinHistValue( i, j_adj );

            value = ( value1 * dy1 + value2 * dy0 ) / inter_denom;
        }

    } else {

        // 2D linear interpolation using inverse distance as weights.
        // Find the adjacent x and y bin centres
        Double_t cbinx_adj = Double_t( i_adj + 0.5 ) * rangeX_ / nBinsX_ + minX_;
        Double_t cbiny_adj = Double_t( j_adj + 0.5 ) * rangeY_ / nBinsY_ + minY_;

        if ( this->withinDPBoundaries( cbinx_adj, cbiny_adj ) == kFALSE ) {

            // The adjacent bin is outside the DP range. Don't extrapolate.
            value = this->getBinHistValue( i, j );

        } else {

            Double_t dx0 = TMath::Abs( x - cbinx );
            Double_t dx1 = TMath::Abs( cbinx_adj - x );
            Double_t dy0 = TMath::Abs( y - cbiny );
            Double_t dy1 = TMath::Abs( cbiny_adj - y );

            Double_t inter_denom = ( dx0 + dx1 ) * ( dy0 + dy1 );

            Double_t value1 = this->getBinHistValue( i, j );
            Double_t value2 = this->getBinHistValue( i_adj, j );
            Double_t value3 = this->getBinHistValue( i, j_adj );
            Double_t value4 = this->getBinHistValue( i_adj, j_adj );

            value = value1 * dx1 * dy1 + value2 * dx0 * dy1 + value3 * dx1 * dy0 + value4 * dx0 * dy0;
            value /= inter_denom;
        }
    }

    return value;
}
