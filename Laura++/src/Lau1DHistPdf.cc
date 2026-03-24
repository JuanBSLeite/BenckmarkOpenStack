
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

/*! \file Lau1DHistPdf.cc
    \brief File containing implementation of Lau1DHistPdf class.
*/

#include "Lau1DHistPdf.hh"

#include "LauRandom.hh"

#include "TAxis.h"
#include "TH1.h"
#include "TMath.h"
#include "TRandom.h"
#include "TSystem.h"

#include <cstdlib>
#include <iostream>
#include <vector>

class LauParameter;

Lau1DHistPdf::Lau1DHistPdf( const TString& theVarName,
                            const TH1* hist,
                            Double_t minAbscissa,
                            Double_t maxAbscissa,
                            Bool_t useInterpolation,
                            Bool_t fluctuateBins ) :
    LauAbsPdf( theVarName, std::vector<LauAbsRValue*>(), minAbscissa, maxAbscissa ),
    hist_( hist ? dynamic_cast<TH1*>( hist->Clone() ) : 0 ),
    useInterpolation_( useInterpolation ),
    fluctuateBins_( fluctuateBins ),
    nBins_( 0 ),
    axisMin_( 0.0 ),
    axisMax_( 0.0 ),
    axisRange_( 0.0 )
{
    // Constructor

    // Set the directory for the histogram
    hist_->SetDirectory( 0 );

    // Save various attributes of the histogram
    nBins_       = hist_->GetNbinsX();
    TAxis* xAxis = hist_->GetXaxis();
    axisMin_     = xAxis->GetXmin();
    axisMax_     = xAxis->GetXmax();
    axisRange_   = axisMax_ - axisMin_;

    // Check that axis range corresponds to range of abscissa
    if ( TMath::Abs( this->getMinAbscissa() - axisMin_ ) > 1e-6 ) {
        std::cerr << "ERROR in Lau1DHistPdf::Lau1DHistPdf : Histogram axis minimum: " << axisMin_
                  << " does not correspond to abscissa minimum: " << this->getMinAbscissa() << "."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }
    if ( TMath::Abs( this->getMaxAbscissa() - axisMax_ ) > 1e-6 ) {
        std::cerr << "ERROR in Lau1DHistPdf::Lau1DHistPdf : Histogram axis maximum: " << axisMax_
                  << " does not correspond to abscissa maximum: " << this->getMaxAbscissa() << "."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // If the bins are to be fluctuated then do so now before
    // calculating anything that depends on the bin content.
    if ( fluctuateBins ) {
        this->doBinFluctuation();
    }

    // Calculate the PDF normalisation.
    this->calcNorm();

    // And check it is OK.
    this->checkNormalisation();
}

Lau1DHistPdf::~Lau1DHistPdf()
{
    // Destructor
    delete hist_;
    hist_ = 0;
}

void Lau1DHistPdf::calcPDFHeight( const LauKinematics* /*kinematics*/ )
{
    if ( this->heightUpToDate() ) {
        return;
    }

    // Get the maximum height of the histogram
    Int_t maxBin    = hist_->GetMaximumBin();
    Double_t height = hist_->GetBinContent( maxBin );
    this->setMaxHeight( height );
}

void Lau1DHistPdf::calcNorm()
{
    // Calculate the histogram normalisation.

    // Loop over the range to get the total area.
    // Just sum the contributions up using 1e-3 increments of the range.
    // Multiply the end result by dx.

    Double_t dx( 1e-3 * axisRange_ );
    Double_t area( 0.0 );

    Double_t x( axisMin_ + dx / 2.0 );
    while ( x > axisMin_ && x < axisMax_ ) {
        area += this->interpolate( x );
        x    += dx;
    }

    Double_t norm = area * dx;
    this->setNorm( norm );
}

void Lau1DHistPdf::checkNormalisation()
{
    Double_t dx( 1e-3 * axisRange_ );
    Double_t area( 0.0 );
    Double_t areaNoNorm( 0.0 );

    const Bool_t useInterpolationOrig { useInterpolation_ };

    Double_t x( axisMin_ + dx / 2.0 );
    while ( x > axisMin_ && x < axisMax_ ) {
        area              += this->interpolateNorm( x );
        useInterpolation_  = kFALSE;
        areaNoNorm        += this->interpolate( x );
        useInterpolation_  = useInterpolationOrig;
        x                 += dx;
    }
    Double_t norm = area * dx;

    std::cout << "INFO in Lau1DHistPdf::checkNormalisation : Area = " << area << ", dx = " << dx
              << std::endl;
    std::cout << "                                         : Area with no norm = " << areaNoNorm
              << "*dx = " << areaNoNorm * dx << std::endl;
    std::cout << "                                         : The total area of the normalised histogram PDF is "
              << norm << std::endl;
}

Double_t Lau1DHistPdf::getBinHistValue( Int_t bin ) const
{
    // Check that bin is in range [1 , nBins_]
    if ( ( bin < 1 ) || ( bin > nBins_ ) ) {
        return 0.0;
    }
    Double_t value = static_cast<Double_t>( hist_->GetBinContent( bin ) );
    // protect against negative values
    if ( value < 0.0 ) {
        std::cerr << "WARNING in Lau1DHistPdf::getBinHistValue : Negative bin content set to zero!"
                  << std::endl;
        value = 0.0;
    }
    return value;
}

Double_t Lau1DHistPdf::interpolateNorm( Double_t x ) const
{
    // Get the normalised interpolated value.
    Double_t value = this->interpolate( x );
    Double_t norm  = this->getNorm();
    return value / norm;
}

void Lau1DHistPdf::calcLikelihoodInfo( const LauAbscissas& abscissas )
{
    // Check that the given abscissa is within the allowed range
    if ( ! this->checkRange( abscissas ) ) {
        gSystem->Exit( EXIT_FAILURE );
    }

    // Get our abscissa
    Double_t abscissa = abscissas[0];

    // Calculate the interpolated value
    Double_t value = this->interpolate( abscissa );
    this->setUnNormPDFVal( value );
}

Double_t Lau1DHistPdf::interpolate( Double_t x ) const
{
    // This function returns the interpolated value of the histogram function
    // for the given value of x by finding the adjacent bins and extrapolating
    // using weights based on the inverse distance of the point from the adajcent
    // bin centres.

    // Find the histogram bin
    Int_t bin = hist_->FindFixBin( x );

    // Ask whether we want to do the interpolation, since this method is
    // not reliable for low statistics histograms.
    if ( useInterpolation_ == kFALSE ) {
        return this->getBinHistValue( bin );
    }

    // Find the bin centres (actual co-ordinate positions, not histogram indices)
    Double_t cbinx = hist_->GetBinCenter( bin );

    // Find the adjacent bins
    Double_t deltax = x - cbinx;

    Int_t bin_adj( 0 );
    if ( deltax > 0.0 ) {
        bin_adj = bin + 1;
    } else {
        bin_adj = bin - 1;
    }

    Bool_t isBoundary( kFALSE );
    if ( bin_adj > nBins_ || bin_adj < 1 ) {
        isBoundary = kTRUE;
    }

    // At the edges, do no interpolation, use entry in bin.
    if ( isBoundary == kTRUE ) {
        return this->getBinHistValue( bin );
    }

    // Linear interpolation using inverse distance as weights.
    // Find the adjacent bin centre
    Double_t cbinx_adj  = hist_->GetBinCenter( bin_adj );
    Double_t deltax_adj = cbinx_adj - x;

    Double_t dx0 = TMath::Abs( deltax );
    Double_t dx1 = TMath::Abs( deltax_adj );

    Double_t denom = dx0 + dx1;

    Double_t value0 = this->getBinHistValue( bin );
    Double_t value1 = this->getBinHistValue( bin_adj );

    Double_t value  = value0 * dx1 + value1 * dx0;
    value          /= denom;

    return value;
}

void Lau1DHistPdf::doBinFluctuation()
{
    TRandom* random = LauRandom::randomFun();
    for ( Int_t bin( 0 ); bin < nBins_; bin++ ) {
        Double_t currentContent = hist_->GetBinContent( bin + 1 );
        Double_t currentError   = hist_->GetBinError( bin + 1 );
        Double_t newContent     = random->Gaus( currentContent, currentError );
        if ( newContent < 0.0 ) {
            hist_->SetBinContent( bin + 1, 0.0 );
        } else {
            hist_->SetBinContent( bin + 1, newContent );
        }
    }
}
