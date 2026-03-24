
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

/*! \file LauScfMap.cc
    \brief File containing implementation of LauScfMap class.
*/

#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>
using std::cerr;
using std::cout;
using std::endl;
using std::flush;

#include "LauScfMap.hh"

#include "TAxis.h"
#include "TH2.h"
#include "TSystem.h"

LauScfMap::LauScfMap() :
    nBinsX_( -1 ),
    nBinsY_( -1 ),
    numberOfBins_( -1 )
{
}

LauScfMap::~LauScfMap()
{
    for ( std::vector<TH2*>::iterator iter = histos_.begin(); iter != histos_.end(); ++iter ) {
        delete ( *iter );
    }
}

void LauScfMap::setHistos( const std::vector<const TH2*>& histos )
{
    // Check that we have enough histograms
    Int_t nBinsX = histos.front()->GetNbinsX();
    Int_t nBinsY = histos.front()->GetNbinsY();
    Int_t nHist  = histos.size();

    if ( nHist != nBinsX * nBinsY ) {
        cerr << "ERROR in LauScfMap::setHistos : There should be as many histograms supplied as there are bins in each of them."
             << endl;
        cerr << "                              : There are " << nHist << " histograms and "
             << nBinsX * nBinsY << " bins." << endl;
        return;
    }

    // number of (tru, and reco as well) bins in a 2D histogram
    // Later we will create as many fake points to cache DP amplitudes
    nBinsX_       = nBinsX;
    nBinsY_       = nBinsY;
    numberOfBins_ = nBinsX_ * nBinsY_;

    // Clear out any old histos we might have lying around
    for ( std::vector<TH2*>::iterator iter = histos_.begin(); iter != histos_.end(); ++iter ) {
        delete ( *iter );
    }
    histos_.clear();

    // Make sure we have enough space for the new histograms
    histos_.reserve( nHist );

    // Loop through the supplied vector
    // Clone each histogram and store it
    for ( std::vector<const TH2*>::const_iterator iter = histos.begin(); iter != histos.end();
          ++iter ) {
        TH2* histo = dynamic_cast<TH2*>( ( *iter )->Clone() );
        if ( histo == 0 ) {
            cerr << "ERROR in LauScfMap::setHistos : Problem cloning one of the histograms." << endl;
            gSystem->Exit( EXIT_FAILURE );
        }
        histos_.push_back( histo );
    }

    if ( histos_.size() != histos.size() ) {
        cerr << "ERROR in LauScfMap::setHistos : Problem cloning the histograms." << endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Now we need to setup the map

    // First clear it
    contribs_.clear();

    // Fill in the map that links, to each reco bin, a vector with the numbers of the tru bins
    // that have non-zero contributions in that reco bin
    // Loop over each histogram bin
    for ( Int_t j( 0 ); j < nBinsY_; ++j ) {
        for ( Int_t i( 0 ); i < nBinsX_; ++i ) {
            Int_t binNo                 = histos_.front()->GetBin( i + 1, j + 1 );
            std::vector<Int_t>& truBins = contribs_[binNo];
            truBins.reserve( nHist );
            for ( Int_t k( 0 ); k < nHist; ++k ) {
                Double_t content = histos_[k]->GetBinContent( binNo );
                if ( content > 0.0 ) {
                    truBins.push_back( k );
                }
            }
        }
    }
}

void LauScfMap::listBinCentres( std::vector<Double_t>& xCoords, std::vector<Double_t>& yCoords ) const
{
    // Create the list of fake points, located at the centres of the
    // true bins, that we will add to the data points, so as to cache
    // the DP amplitudes, scfFractions and jacobians at those points.
    // that have non-zero contributions in that reco bin

    const TAxis* xAxis  = histos_.front()->GetXaxis();
    Double_t xMax       = xAxis->GetXmax();
    Double_t xMin       = xAxis->GetXmin();
    Double_t xIncrement = (Double_t)( ( xMax - xMin ) / nBinsX_ );

    const TAxis* yAxis  = histos_.front()->GetYaxis();
    Double_t yMax       = yAxis->GetXmax();
    Double_t yMin       = yAxis->GetXmin();
    Double_t yIncrement = (Double_t)( ( yMax - yMin ) / nBinsY_ );

    // Loop over each histogram bin
    Double_t yCoord = yMin + yIncrement / 2;
    for ( Int_t j( 1 ); j <= nBinsY_; ++j ) {
        Double_t xCoord = xMin + xIncrement / 2;
        for ( Int_t i( 1 ); i <= nBinsX_; ++i ) {
            xCoords.push_back( xCoord );
            yCoords.push_back( yCoord );
            xCoord += xIncrement;
        }
        yCoord += yIncrement;
    }
}

Int_t LauScfMap::binNumber( Double_t xCoord, Double_t yCoord ) const
{
    const TH2* histo = histos_.front();
    if ( histo != 0 ) {
        return histo->FindFixBin( xCoord, yCoord );
    } else {
        cerr << "ERROR in LauScfMap::recoBin : No valid histograms found." << endl;
        return -1;
    }
}

const std::vector<Int_t>* LauScfMap::trueBins( Int_t recoBin ) const
{
    std::map<Int_t, std::vector<Int_t>>::const_iterator iter = contribs_.find( recoBin );
    if ( iter != contribs_.end() ) {
        return &( iter->second );
    } else {
        return 0;
    }
}

Double_t LauScfMap::prob( Int_t recoBin, Int_t trueBin ) const
{
    const TH2* histo = histos_[trueBin];
    if ( histo != 0 ) {
        return histo->GetBinContent( recoBin );
    } else {
        cerr << "ERROR in LauScfMap::prob : No histogram found for true bin " << trueBin << endl;
        return 0.0;
    }
}

TH2* LauScfMap::trueHist( Int_t trueBin )
{
    // Need to turn the histogram bin number into the vector element
    Int_t x( 0 ), y( 0 ), z( 0 );
    histos_.front()->GetBinXYZ( trueBin, x, y, z );

    Int_t theBin = nBinsX_ * ( y - 1 ) + ( x - 1 );

    if ( theBin < 0 || static_cast<UInt_t>( theBin ) >= histos_.size() ) {
        cerr << "ERROR in LauScfMap::trueHist : No histogram found for true bin " << trueBin
             << ", which corresponds to x=" << x << ", y=" << y << ", or entry " << theBin
             << " in the vector." << endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    TH2* histo = histos_[theBin];
    if ( histo == 0 ) {
        cerr << "ERROR in LauScfMap::trueHist : Null histogram pointer found for true bin "
             << trueBin << ", which corresponds to x=" << x << ", y=" << y << ", or entry "
             << theBin << " in the vector." << endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    return histo;
}
