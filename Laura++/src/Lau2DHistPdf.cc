
/*
Copyright 2008 University of Warwick

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

/*! \file Lau2DHistPdf.cc
    \brief File containing implementation of Lau2DHistPdf class.
*/

#include "Lau2DHistPdf.hh"

#include "Lau1DHistPdf.hh"
#include "LauRandom.hh"

#include "TAxis.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom.h"
#include "TSystem.h"

#include <cstdlib>
#include <iostream>
#include <vector>

class LauParameter;

Lau2DHistPdf::Lau2DHistPdf( const std::vector<TString>& theVarNames,
                            const TH2* hist,
                            const LauFitData& minVals,
                            const LauFitData& maxVals,
                            Bool_t useInterpolation,
                            Bool_t fluctuateBins ) :
    LauAbsPdf( theVarNames, std::vector<LauAbsRValue*>(), minVals, maxVals ),
    hist_( hist ? dynamic_cast<TH2*>( hist->Clone() ) : 0 ),
    xProj_( 0 ),
    yProj_( 0 ),
    xVarPdf_( 0 ),
    yVarPdf_( 0 ),
    xName_( "" ),
    yName_( "" ),
    nBinsX_( 0 ),
    nBinsY_( 0 ),
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
    useInterpolation_( useInterpolation ),
    fluctuateBins_( fluctuateBins )
{
    // Constructor

    // Check that we have two variables
    if ( this->nInputVars() != 2 ) {
        std::cerr << "ERROR in Lau2DHistPdf::Lau2DHistPdf : Have not been provided with exactly two variables."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Set the variable names from the abstract class
    xName_ = this->varNames()[0];
    yName_ = this->varNames()[1];

    // Set the variable limits from the abstract class
    minX_   = this->getMinAbscissa( xName_ );
    maxX_   = this->getMaxAbscissa( xName_ );
    minY_   = this->getMinAbscissa( yName_ );
    maxY_   = this->getMaxAbscissa( yName_ );
    rangeX_ = maxX_ - minX_;
    rangeY_ = maxY_ - minY_;

    // Have we got a valid histogram
    if ( hist_ == 0 ) {
        std::cerr << "ERROR in Lau2DHistPdf::Lau2DHistPdf : Histogram pointer is null." << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Set the directory for the histogram
    hist_->SetDirectory( 0 );

    // Save various attributes of the histogram
    nBinsX_ = hist_->GetNbinsX();
    nBinsY_ = hist_->GetNbinsY();

    TAxis* xAxis      = hist_->GetXaxis();
    Double_t xAxisMin = xAxis->GetXmin();
    Double_t xAxisMax = xAxis->GetXmax();

    TAxis* yAxis      = hist_->GetYaxis();
    Double_t yAxisMin = yAxis->GetXmin();
    Double_t yAxisMax = yAxis->GetXmax();

    // Check that axis ranges corresponds to ranges of abscissas
    if ( TMath::Abs( minX_ - xAxisMin ) > 1e-6 ) {
        std::cerr << "ERROR in Lau2DHistPdf::Lau2DHistPdf : Histogram x-axis minimum: " << xAxisMin
                  << " does not correspond to " << xName_ << " minimum: " << minX_ << "."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }
    if ( TMath::Abs( maxX_ - xAxisMax ) > 1e-6 ) {
        std::cerr << "ERROR in Lau2DHistPdf::Lau2DHistPdf : Histogram x-axis maximum: " << xAxisMax
                  << " does not correspond to " << xName_ << " maximum: " << maxX_ << "."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }
    if ( TMath::Abs( minY_ - yAxisMin ) > 1e-6 ) {
        std::cerr << "ERROR in Lau2DHistPdf::Lau2DHistPdf : Histogram y-axis minimum: " << yAxisMin
                  << " does not correspond to " << yName_ << " minimum: " << minY_ << "."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }
    if ( TMath::Abs( maxY_ - yAxisMax ) > 1e-6 ) {
        std::cerr << "ERROR in Lau2DHistPdf::Lau2DHistPdf : Histogram y-axis maximum: " << yAxisMax
                  << " does not correspond to " << yName_ << " maximum: " << maxY_ << "."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Calculate the bin widths and inverse bin widths
    binXWidth_ = static_cast<Double_t>( TMath::Abs( rangeX_ ) / ( nBinsX_ * 1.0 ) );
    binYWidth_ = static_cast<Double_t>( TMath::Abs( rangeY_ ) / ( nBinsY_ * 1.0 ) );
    if ( binXWidth_ > 1e-10 ) {
        invBinXWidth_ = 1.0 / binXWidth_;
    }
    if ( binYWidth_ > 1e-10 ) {
        invBinYWidth_ = 1.0 / binYWidth_;
    }

    // If the bins are to be fluctuated then do so now before
    // calculating anything that depends on the bin content.
    if ( fluctuateBins ) {
        this->doBinFluctuation();
    }

    // Create the projections and 1D PDFs from those
    xProj_ = hist_->ProjectionX();
    yProj_ = hist_->ProjectionY();
    if ( ! xProj_ || ! yProj_ ) {
        std::cerr << "ERROR in Lau2DHistPdf::Lau2DHistPdf : Can't get X or Y projection from 2D histogram."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }
    xVarPdf_ = new Lau1DHistPdf( xName_, xProj_, minX_, maxX_, useInterpolation_, fluctuateBins_ );
    yVarPdf_ = new Lau1DHistPdf( yName_, yProj_, minY_, maxY_, useInterpolation_, fluctuateBins_ );

    // Calculate the PDF normalisation.
    this->calcNorm();

    // And check it is OK.
    this->checkNormalisation();
}

Lau2DHistPdf::~Lau2DHistPdf()
{
    // Destructor
    delete hist_;
    hist_ = 0;
    delete xProj_;
    xProj_ = 0;
    delete yProj_;
    yProj_ = 0;
    delete xVarPdf_;
    xVarPdf_ = 0;
    delete yVarPdf_;
    yVarPdf_ = 0;
}

void Lau2DHistPdf::calcPDFHeight( const LauKinematics* /*kinematics*/ )
{
    if ( this->heightUpToDate() ) {
        return;
    }

    // Get the maximum height of the histogram
    Int_t maxBin    = hist_->GetMaximumBin();
    Double_t height = hist_->GetBinContent( maxBin );
    this->setMaxHeight( height );
}

void Lau2DHistPdf::calcNorm()
{
    // Calculate the histogram normalisation.

    // Loop over the total x and y range to get the total area under x
    // and y. Just sum the contributions up using 1e-3 increments of
    // the range in x and y. Multiply the end result by dx and dy.

    Double_t dx( 1e-3 * rangeX_ ), dy( 1e-3 * rangeY_ );
    Double_t area( 0.0 );

    Double_t x( minX_ + dx / 2.0 );
    while ( x > minX_ && x < maxX_ ) {
        Double_t y( minY_ + dy / 2.0 );
        while ( y > minY_ && y < maxY_ ) {
            area += this->interpolateXY( x, y );
            y    += dy;
        }    // y while loop
        x += dx;
    }    // x while loop

    Double_t norm = area * dx * dy;
    this->setNorm( norm );
}

void Lau2DHistPdf::checkNormalisation()
{
    Double_t dx( 1e-3 * rangeX_ ), dy( 1e-3 * rangeY_ );
    Double_t area( 0.0 );
    Double_t areaNoNorm( 0.0 );

    // Preserve the value of a variable we change temporarily
    Bool_t interpolate = useInterpolation_;

    Double_t x( minX_ + dx / 2.0 );
    while ( x > minX_ && x < maxX_ ) {
        Double_t y( minY_ + dy / 2.0 );
        while ( y > minY_ && y < maxY_ ) {
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

    std::cout << "INFO in Lau2DHistPdf::checkNormalisation : Area = " << area << ", dx = " << dx
              << ", dy = " << dy << ", dx*dy = " << dx * dy << std::endl;
    std::cout << "                                         : Area with no norm = " << areaNoNorm
              << "*dx*dy = " << areaNoNorm * dx * dy << std::endl;
    std::cout << "                                         : The total area of the normalised histogram PDF is "
              << norm << std::endl;
}

Double_t Lau2DHistPdf::getBinHistValue( Int_t i, Int_t j ) const
{
    // Check that x co-ord is in range [1 , nBinsX_]
    if ( ( i < 1 ) || ( i > nBinsX_ ) ) {
        return 0.0;
    }
    // Check that y co-ord is in range [1 , nBinsY_]
    if ( ( j < 1 ) || ( j > nBinsY_ ) ) {
        return 0.0;
    }
    Double_t value = hist_->GetBinContent( i, j );
    // protect against negative values
    if ( value < 0.0 ) {
        std::cerr << "WARNING in Lau2DHistPdf::getBinHistValue : Negative bin content set to zero!"
                  << std::endl;
        value = 0.0;
    }
    return value;
}

Double_t Lau2DHistPdf::interpolateXYNorm( Double_t x, Double_t y ) const
{
    // Get the normalised interpolated value.
    Double_t value = this->interpolateXY( x, y );
    Double_t norm  = this->getNorm();
    return value / norm;
}

Double_t Lau2DHistPdf::interpolateXY( Double_t x, Double_t y ) const
{
    // This function returns the interpolated value of the histogram
    // function for the given values of x and y by finding the adjacent
    // bins and extrapolating using weights based on the inverse
    // distance of the point from the adajcent bin centres.

    // Find the 2D histogram bin for x and y
    Int_t i = static_cast<Int_t>( ( x - minX_ ) * invBinXWidth_ ) + 1;
    Int_t j = static_cast<Int_t>( ( y - minY_ ) * invBinYWidth_ ) + 1;

    if ( i > nBinsX_ ) {
        i = nBinsX_;
    }
    if ( j > nBinsY_ ) {
        j = nBinsY_;
    }

    // Ask whether we want to do the interpolation, since this method is
    // not reliable for low statistics histograms.
    if ( useInterpolation_ == kFALSE ) {
        return this->getBinHistValue( i, j );
    }

    // Find the bin centres (actual co-ordinate positions, not histogram indices)
    Double_t cbinx = static_cast<Double_t>( i - 0.5 ) * binXWidth_ + minX_;
    Double_t cbiny = static_cast<Double_t>( j - 0.5 ) * binYWidth_ + minY_;

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

    if ( i_adj > nBinsX_ || i_adj < 1 ) {
        isYBoundary = kTRUE;
    }
    if ( j_adj > nBinsY_ || j_adj < 1 ) {
        isXBoundary = kTRUE;
    }

    // In the corners, do no interpolation. Use entry in bin (i,j)
    if ( isXBoundary == kTRUE && isYBoundary == kTRUE ) {

        value = this->getBinHistValue( i, j );

    } else if ( isXBoundary == kTRUE && isYBoundary == kFALSE ) {

        // Find the adjacent x bin centre
        Double_t cbinx_adj = static_cast<Double_t>( i_adj - 0.5 ) * binXWidth_ + minX_;

        Double_t dx0         = TMath::Abs( x - cbinx );
        Double_t dx1         = TMath::Abs( cbinx_adj - x );
        Double_t inter_denom = dx0 + dx1;

        Double_t value1 = this->getBinHistValue( i, j );
        Double_t value2 = this->getBinHistValue( i_adj, j );

        value = ( value1 * dx1 + value2 * dx0 ) / inter_denom;

    } else if ( isYBoundary == kTRUE && isXBoundary == kFALSE ) {

        // Find the adjacent y bin centre
        Double_t cbiny_adj = static_cast<Double_t>( j_adj - 0.5 ) * binYWidth_ + minY_;

        Double_t dy0         = TMath::Abs( y - cbiny );
        Double_t dy1         = TMath::Abs( cbiny_adj - y );
        Double_t inter_denom = dy0 + dy1;

        Double_t value1 = this->getBinHistValue( i, j );
        Double_t value2 = this->getBinHistValue( i, j_adj );

        value = ( value1 * dy1 + value2 * dy0 ) / inter_denom;

    } else {

        // 2D linear interpolation using inverse distance as weights.
        // Find the adjacent x and y bin centres
        Double_t cbinx_adj = static_cast<Double_t>( i_adj - 0.5 ) * binXWidth_ + minX_;
        Double_t cbiny_adj = static_cast<Double_t>( j_adj - 0.5 ) * binYWidth_ + minY_;

        Double_t dx0 = TMath::Abs( x - cbinx );
        Double_t dx1 = TMath::Abs( cbinx_adj - x );
        Double_t dy0 = TMath::Abs( y - cbiny );
        Double_t dy1 = TMath::Abs( cbiny_adj - y );

        Double_t inter_denom = ( dx0 + dx1 ) * ( dy0 + dy1 );

        Double_t value1 = this->getBinHistValue( i, j );
        Double_t value2 = this->getBinHistValue( i_adj, j );
        Double_t value3 = this->getBinHistValue( i, j_adj );
        Double_t value4 = this->getBinHistValue( i_adj, j_adj );

        value  = value1 * dx1 * dy1 + value2 * dx0 * dy1 + value3 * dx1 * dy0 + value4 * dx0 * dy0;
        value /= inter_denom;
    }

    return value;
}

void Lau2DHistPdf::calcLikelihoodInfo( const LauAbscissas& abscissas )
{
    // Check that the given abscissa is within the allowed range
    if ( ! this->checkRange( abscissas ) ) {
        gSystem->Exit( EXIT_FAILURE );
    }

    // Calculate the interpolated value
    Double_t x     = abscissas[0];
    Double_t y     = abscissas[1];
    Double_t value = this->interpolateXY( x, y );
    this->setUnNormPDFVal( value );

    // TODO - do we need this?  I think probably not.
    // Have the 1D PDFs calculate their likelihoods as well
    //xVarPdf_->calcLikelihoodInfo(x);
    //yVarPdf_->calcLikelihoodInfo(y);
}

void Lau2DHistPdf::calcLikelihoodInfo( UInt_t iEvt )
{
    // Call the base class method
    this->LauAbsPdf::calcLikelihoodInfo( iEvt );

    // Have the 1D PDFs retrieve their values as well
    xVarPdf_->calcLikelihoodInfo( iEvt );
    yVarPdf_->calcLikelihoodInfo( iEvt );
}

Double_t Lau2DHistPdf::getLikelihood( const TString& theVarName ) const
{
    if ( theVarName == xName_ ) {
        return xVarPdf_->getLikelihood();
    } else if ( theVarName == yName_ ) {
        return yVarPdf_->getLikelihood();
    } else {
        std::cerr << "ERROR in Lau2DHistPdf::getLikelihood : Unrecognised variable name \""
                  << theVarName << "\", cannot determine likelihood." << std::endl;
        return 0.0;
    }
}

void Lau2DHistPdf::cacheInfo( const LauFitDataTree& inputData )
{
    // Call the base class method
    this->LauAbsPdf::cacheInfo( inputData );

    // Have the 1D PDFs cache their info as well
    xVarPdf_->cacheInfo( inputData );
    yVarPdf_->cacheInfo( inputData );
}

void Lau2DHistPdf::doBinFluctuation()
{
    TRandom* random( LauRandom::randomFun() );
    for ( Int_t i( 0 ); i < nBinsX_; i++ ) {
        for ( Int_t j( 0 ); j < nBinsY_; j++ ) {
            Double_t currentContent = hist_->GetBinContent( i + 1, j + 1 );
            Double_t currentError   = hist_->GetBinError( i + 1, j + 1 );
            Double_t newContent     = random->Gaus( currentContent, currentError );
            if ( newContent < 0.0 ) {
                hist_->SetBinContent( i + 1, j + 1, 0.0 );
            } else {
                hist_->SetBinContent( i + 1, j + 1, newContent );
            }
        }
    }
}

LauFitData Lau2DHistPdf::generate( const LauKinematics* kinematics )
{
    this->withinGeneration( kTRUE );

    // Check that the PDF height is up to date
    this->calcPDFHeight( kinematics );

    // Generate the value of the abscissa.
    Bool_t gotAbscissa( kFALSE );
    if ( this->getRandomFun() == 0 ) {
        std::cerr << "ERROR in Lau2DHistPdf::generate : Please set the random number generator for this PDF by using the setRandomFun(TRandom*) function."
                  << std::endl;
        this->withinGeneration( kFALSE );
        gSystem->Exit( EXIT_FAILURE );
    }

    Double_t genPDFVal( 0.0 );
    LauAbscissas genAbscissas( 2 );

    Double_t PDFheight = this->getMaxHeight() * ( 1.0 + 1e-11 );
    while ( ! gotAbscissa ) {

        genAbscissas[0] = this->getRandomFun()->Rndm() * this->getRange( xName_ ) +
                          this->getMinAbscissa( xName_ );
        genAbscissas[1] = this->getRandomFun()->Rndm() * this->getRange( yName_ ) +
                          this->getMinAbscissa( yName_ );

        this->calcLikelihoodInfo( genAbscissas );
        genPDFVal = this->getUnNormLikelihood();

        if ( this->getRandomFun()->Rndm() <= genPDFVal / PDFheight ) {
            gotAbscissa = kTRUE;
        }

        if ( genPDFVal > PDFheight ) {
            std::cerr << "WARNING in LauAbsPdf::generate : genPDFVal = " << genPDFVal
                      << " is larger than the specified PDF height " << this->getMaxHeight()
                      << " for (x,y) = (" << genAbscissas[0] << "," << genAbscissas[1] << ")."
                      << std::endl;
            std::cerr << "                               : Need to reset height to be larger than "
                      << genPDFVal
                      << " by using the setMaxHeight(Double_t) function and re-run the Monte Carlo generation!"
                      << std::endl;
        }
    }

    LauFitData genData;
    genData[xName_] = genAbscissas[0];
    genData[yName_] = genAbscissas[1];

    this->withinGeneration( kFALSE );

    return genData;
}
