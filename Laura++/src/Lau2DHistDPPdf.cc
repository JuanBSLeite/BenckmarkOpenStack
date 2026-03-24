
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

/*! \file Lau2DHistDPPdf.cc
    \brief File containing implementation of Lau2DHistDPPdf class.
*/

#include "Lau2DHistDPPdf.hh"

#include "LauKinematics.hh"
#include "LauRandom.hh"
#include "LauVetoes.hh"

#include "TAxis.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom.h"

#include <iostream>

Lau2DHistDPPdf::Lau2DHistDPPdf( const TH2* hist,
                                LauKinematics* kinematics,
                                const LauVetoes* vetoes,
                                Bool_t useInterpolation,
                                Bool_t fluctuateBins,
                                Bool_t useUpperHalfOnly,
                                Bool_t squareDP ) :
    Lau2DAbsHistDPPdf( kinematics, vetoes, useUpperHalfOnly, squareDP ),
    hist_( hist ? dynamic_cast<TH2*>( hist->Clone() ) : 0 ),
    minX_( 0.0 ),
    maxX_( 0.0 ),
    minY_( 0.0 ),
    maxY_( 0.0 ),
    rangeX_( 0.0 ),
    rangeY_( 0.0 ),
    binXWidth_( 0.0 ),
    binYWidth_( 0.0 ),
    invBinXWidth_( 0.0 ),
    invBinYWidth_( 0.0 ),
    nBinsX_( 0 ),
    nBinsY_( 0 ),
    norm_( 0.0 ),
    useInterpolation_( useInterpolation )
{
    // For square Dalitz plots, the co-ordinates must be m', theta'
    // The input data to the fit is still in m13^2, m23^2 co-ords, and a
    // transformation is applied to go from one co-ordinate system to the
    // other.

    // Save various attributes of the histogram
    // (axis ranges, number of bins, areas)
    if ( hist_ ) {
        TAxis* xAxis = hist_->GetXaxis();
        minX_        = xAxis->GetXmin();
        maxX_        = xAxis->GetXmax();

        TAxis* yAxis = hist_->GetYaxis();
        minY_        = yAxis->GetXmin();
        maxY_        = yAxis->GetXmax();

        nBinsX_ = hist_->GetNbinsX();
        nBinsY_ = hist_->GetNbinsY();
    } else {
        minX_ = getKinematics()->getm13SqMin();
        maxX_ = getKinematics()->getm13SqMax();

        minY_ = getKinematics()->getm23SqMin();
        maxY_ = getKinematics()->getm23SqMax();

        nBinsX_ = 1;
        nBinsY_ = 1;
    }

    rangeX_ = maxX_ - minX_;
    rangeY_ = maxY_ - minY_;
    if ( nBinsX_ > 0 ) {
        binXWidth_ = TMath::Abs( rangeX_ ) / ( nBinsX_ * 1.0 );
    }
    if ( nBinsY_ > 0 ) {
        binYWidth_ = TMath::Abs( rangeY_ ) / ( nBinsY_ * 1.0 );
    }
    if ( binXWidth_ > 1e-10 ) {
        invBinXWidth_ = 1.0 / binXWidth_;
    }
    if ( binYWidth_ > 1e-10 ) {
        invBinYWidth_ = 1.0 / binYWidth_;
    }

    // If the bins are to be fluctuated then do so now before
    // calculating anything that depends on the bin content.
    if ( fluctuateBins ) {
        this->doBinFluctuation( hist_ );
    }

    // Calculate the PDF normalisation.
    this->calcHistNorm();

    // And check it is OK.
    this->checkNormalisation();

    // Also obtain the maximum height
    this->calcMaxHeight( hist_ );
}

Lau2DHistDPPdf::~Lau2DHistDPPdf()
{
    // Destructor
    delete hist_;
    hist_ = 0;
}

void Lau2DHistDPPdf::calcHistNorm()
{
    // Calculate the histogram normalisation. We must integrate
    // over the allowed kinematic DP region.

    // Loop over the total x and y range to get the total area
    // under x and y. Just sum the contributions up using 1e-3 increments
    // of the range in x and y. Multiply the end result by dx and dy.

    const Double_t axisMin   = TMath::Min( minX_, minY_ );
    const Double_t axisMax   = TMath::Max( maxX_, maxY_ );
    const Double_t axisRange = axisMax - axisMin;

    const Double_t dx( 1e-3 * axisRange );
    const Double_t dy( dx );

    Double_t area( 0.0 );

    Double_t x( axisMin + dx / 2.0 );
    while ( x < axisMax ) {
        Double_t y( axisMin + dy / 2.0 );
        while ( y < axisMax ) {
            area += this->interpolateXY( x, y );
            y    += dy;
        }    // y while loop
        x += dx;
    }    // x while loop

    norm_ = area * dx * dy;

    std::cout << "INFO in Lau2DHistDPPdf::calcHistNorm : Norm = " << norm_
              << ", bX*bY = " << binXWidth_ << "*" << binYWidth_ << " = " << binXWidth_ * binYWidth_
              << std::endl;
}

Double_t Lau2DHistDPPdf::getBinHistValue( Int_t xBinNo, Int_t yBinNo ) const
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

    if ( hist_ == 0 ) {
        return 1.0;
    }

    Double_t value = hist_->GetBinContent( xBinNo + 1, yBinNo + 1 );
    return value;
}

Double_t Lau2DHistDPPdf::interpolateXYNorm( Double_t x, Double_t y ) const
{
    // Get the normalised interpolated value.
    Double_t value = this->interpolateXY( x, y );
    return value / norm_;
}

Double_t Lau2DHistDPPdf::interpolateXY( Double_t x, Double_t y ) const
{
    // This function returns the interpolated value of the histogram function
    // for the given values of x and y by finding the adjacent bins and extrapolating
    // using weights based on the inverse distance of the point from the adajcent
    // bin centres.
    // Here, x = m13^2, y = m23^2, or m', theta' for square DP co-ordinates

    // If we're only using one half then flip co-ordinates
    // appropriately for conventional or square DP
    getUpperHalf( x, y );

    // First ask whether the point is inside the kinematic region.
    if ( withinDPBoundaries( x, y ) == kFALSE ) {
        return 0.0;
    }

    // Update the kinematics to the position of interest.
    updateKinematics( x, y );

    // Check that we're not inside a veto
    Bool_t vetoOK( kTRUE );
    if ( getVetoes() ) {
        vetoOK = getVetoes()->passVeto( getKinematics() );
    }
    if ( vetoOK == kFALSE ) {
        return 0.0;
    }

    // Find the 2D histogram bin for x and y
    Int_t i = TMath::FloorNint( ( x - minX_ ) * invBinXWidth_ );
    Int_t j = TMath::FloorNint( ( y - minY_ ) * invBinYWidth_ );

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
    Double_t cbinx = ( i + 0.5 ) * binXWidth_ + minX_;
    Double_t cbiny = ( j + 0.5 ) * binYWidth_ + minY_;

    // If bin centres are outside kinematic region, do not extrapolate
    if ( withinDPBoundaries( cbinx, cbiny ) == kFALSE ) {
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
        Double_t cbinx_adj = Double_t( i_adj + 0.5 ) * binXWidth_ + minX_;

        if ( withinDPBoundaries( cbinx_adj, y ) == kFALSE ) {

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
        Double_t cbiny_adj = Double_t( j_adj + 0.5 ) * binYWidth_ + minY_;

        if ( withinDPBoundaries( x, cbiny_adj ) == kFALSE ) {

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
        Double_t cbinx_adj = Double_t( i_adj + 0.5 ) * binXWidth_ + minX_;
        Double_t cbiny_adj = Double_t( j_adj + 0.5 ) * binYWidth_ + minY_;

        if ( withinDPBoundaries( cbinx_adj, cbiny_adj ) == kFALSE ) {

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

void Lau2DHistDPPdf::checkNormalisation()
{
    // Loop over the total x and y range to get the total area
    // under x and y. Just sum the contributions up using 1e-3 increments
    // of the range in x and y. Multiply the end result by dx and dy.

    const Double_t axisMin   = TMath::Min( minX_, minY_ );
    const Double_t axisMax   = TMath::Max( maxX_, maxY_ );
    const Double_t axisRange = axisMax - axisMin;

    const Double_t dx( 1e-3 * axisRange );
    const Double_t dy( dx );

    Double_t area( 0.0 );
    Double_t areaNoNorm( 0.0 );

    // Preserve the value of a variable we change temporarily
    Bool_t interpolate = useInterpolation_;

    Double_t x( axisMin + dx / 2.0 );
    while ( x < axisMax ) {
        Double_t y( axisMin + dy / 2.0 );
        while ( y < axisMax ) {
            area              += this->interpolateXYNorm( x, y );
            useInterpolation_  = kFALSE;
            areaNoNorm        += this->interpolateXY( x, y );
            useInterpolation_  = kTRUE;
            y                 += dy;
        }    // y while loop
        x += dx;
    }    // x while loop
    Double_t norm = area * dx * dy;

    useInterpolation_ = interpolate;    //Restore old value

    std::cout << "INFO in Lau2DHistDPPdf::checkNormalisation : Area = " << area << ", dx = " << dx
              << ", dy = " << dy << ", dx*dy = " << dx * dy << std::endl;
    std::cout << "                                           : Area with no norm = " << areaNoNorm
              << "*dx*dy = " << areaNoNorm * dx * dy << std::endl;
    std::cout << "                                           : The total area of the normalised histogram PDF is "
              << norm << std::endl;
}
