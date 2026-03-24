
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

/*! \file Lau2DAbsHistDPPdf.cc
    \brief File containing implementation of Lau2DAbsHistDPPdf class.
*/

#include "Lau2DAbsHistDPPdf.hh"

#include "LauKinematics.hh"
#include "LauRandom.hh"
#include "LauVetoes.hh"

#include "TAxis.h"
#include "TH2.h"
#include "TRandom.h"

#include <iostream>

Lau2DAbsHistDPPdf::Lau2DAbsHistDPPdf( LauKinematics* kinematics,
                                      const LauVetoes* vetoes,
                                      Bool_t useUpperHalfOnly,
                                      Bool_t squareDP ) :
    kinematics_( kinematics ),
    vetoes_( vetoes ),
    upperHalf_( useUpperHalfOnly ),
    squareDP_( squareDP )
{
    if ( squareDP && ! kinematics->squareDP() ) {
        // The histogram provided is defined in the square DP but the
        // kinematics object has calculation of the square DP
        // co-ordinates disabled, so need to enable it,
        std::cerr << "WARNING in Lau2DAbsHistDPPdf constructor : forcing kinematics to calculate the required square DP co-ordinates"
                  << std::endl;
        kinematics->squareDP( kTRUE );
    }
}

Lau2DAbsHistDPPdf::~Lau2DAbsHistDPPdf()
{
}

void Lau2DAbsHistDPPdf::calcMaxHeight( TH2* hist )
{
    // Get the maximum height of the 2D histogram
    maxHeight_ = 1.0;
    if ( hist ) {
        Int_t maxBin = hist->GetMaximumBin();
        maxHeight_   = hist->GetBinContent( maxBin );
    }

    std::cout << "INFO in Lau2DAbsHistDPPdf::calcMaxHeight : Max height = " << maxHeight_
              << std::endl;
}

void Lau2DAbsHistDPPdf::doBinFluctuation( TH2* hist )
{
    if ( ! hist ) {
        return;
    }

    TRandom* random = LauRandom::randomFun();

    Int_t nBinsX = static_cast<Int_t>( hist->GetNbinsX() );
    Int_t nBinsY = static_cast<Int_t>( hist->GetNbinsY() );

    for ( Int_t i( 0 ); i < nBinsX; i++ ) {
        for ( Int_t j( 0 ); j < nBinsY; j++ ) {
            Double_t currentContent = hist->GetBinContent( i + 1, j + 1 );
            Double_t currentError   = hist->GetBinError( i + 1, j + 1 );
            Double_t newContent     = random->Gaus( currentContent, currentError );
            if ( newContent < 0.0 ) {
                hist->SetBinContent( i + 1, j + 1, 0.0 );
            } else {
                hist->SetBinContent( i + 1, j + 1, newContent );
            }
        }
    }
}

Bool_t Lau2DAbsHistDPPdf::withinDPBoundaries( Double_t x, Double_t y ) const
{
    return squareDP_ ? kinematics_->withinSqDPLimits( x, y ) : kinematics_->withinDPLimits( x, y );
}

void Lau2DAbsHistDPPdf::getUpperHalf( Double_t& x, Double_t& y ) const
{
    if ( upperHalf_ == kTRUE ) {
        if ( squareDP_ == kFALSE && x > y ) {
            Double_t temp = y;
            y             = x;
            x             = temp;
        } else if ( squareDP_ == kTRUE && y > 0.5 ) {
            y = 1.0 - y;
        }
    }
}

void Lau2DAbsHistDPPdf::updateKinematics( Double_t x, Double_t y ) const
{
    if ( squareDP_ == kTRUE ) {
        kinematics_->updateSqDPKinematics( x, y );
    } else {
        kinematics_->updateKinematics( x, y );
    }
}
