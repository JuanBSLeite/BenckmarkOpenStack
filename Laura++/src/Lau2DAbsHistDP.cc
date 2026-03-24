
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

/*! \file Lau2DAbsHistDP.cc
    \brief File containing implementation of Lau2DAbsHistDP class.
*/

#include "Lau2DAbsHistDP.hh"

#include "LauBifurcatedGaussPdf.hh"
#include "LauDaughters.hh"
#include "LauKinematics.hh"
#include "LauRandom.hh"

#include "TAxis.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom.h"
#include "TSystem.h"

#include <iostream>

Lau2DAbsHistDP::Lau2DAbsHistDP( const LauDaughters* daughters,
                                Bool_t useUpperHalfOnly,
                                Bool_t squareDP ) :
    kinematics_( ( daughters != 0 ) ? daughters->getKinematics() : 0 ),
    upperHalf_( useUpperHalfOnly ),
    squareDP_( squareDP )
{
    if ( squareDP && ! daughters->squareDP() ) {
        // The histogram provided is defined in the square DP but the
        // kinematics object has calculation of the square DP
        // co-ordinates disabled, so need to enable it,
        // which requires a bit of a unpleasant const_cast
        std::cerr << "WARNING in Lau2DAbsHistDP constructor : forcing kinematics to calculate the required square DP co-ordinates"
                  << std::endl;
        LauKinematics* kine = const_cast<LauKinematics*>( kinematics_ );
        kine->squareDP( kTRUE );
    }
}

Lau2DAbsHistDP::~Lau2DAbsHistDP()
{
}

void Lau2DAbsHistDP::doBinFluctuation( TH2* hist )
{
    TRandom* random = LauRandom::randomFun();

    Int_t nBinsX = static_cast<Int_t>( hist->GetNbinsX() );
    Int_t nBinsY = static_cast<Int_t>( hist->GetNbinsY() );

    for ( Int_t i( 0 ); i < nBinsX; ++i ) {
        for ( Int_t j( 0 ); j < nBinsY; ++j ) {
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

void Lau2DAbsHistDP::doBinFluctuation( TH2* hist, const TH2* errorHi, const TH2* errorLo )
{
    LauParameter* mean = new LauParameter( "mean", 0.5, 0.0, 1.0, kFALSE );
    LauParameter* sigL = new LauParameter( "sigmaL", 0.5, 0.0, 1.0, kFALSE );
    LauParameter* sigR = new LauParameter( "sigmaR", 0.5, 0.0, 1.0, kFALSE );

    std::vector<LauAbsRValue*> pars( 3 );
    pars[0] = mean;
    pars[1] = sigL;
    pars[2] = sigR;

    const TString varName( "tmp" );

    Int_t nBinsX = static_cast<Int_t>( hist->GetNbinsX() );
    Int_t nBinsY = static_cast<Int_t>( hist->GetNbinsY() );

    LauFitData genData;

    for ( Int_t i( 0 ); i < nBinsX; ++i ) {
        for ( Int_t j( 0 ); j < nBinsY; ++j ) {
            const Double_t currentContent = hist->GetBinContent( i + 1, j + 1 );
            const Double_t currentErrorLo = errorLo->GetBinContent( i + 1, j + 1 );
            const Double_t currentErrorHi = errorHi->GetBinContent( i + 1, j + 1 );

            mean->value( currentContent );
            sigL->value( currentErrorLo );
            sigR->value( currentErrorHi );

            const Double_t minVal = TMath::Max( 0.0, currentContent - 5.0 * currentErrorLo );
            const Double_t maxVal = TMath::Min( 1.0, currentContent + 5.0 * currentErrorHi );

            LauBifurcatedGaussPdf bfgaus( varName, pars, minVal, maxVal );
            bfgaus.heightUpToDate( kFALSE );
            genData = bfgaus.generate( 0 );

            const Double_t newContent = genData[varName];

            hist->SetBinContent( i + 1, j + 1, newContent );
        }
    }

    delete pars[0];
    delete pars[1];
    delete pars[2];
    pars.clear();
}

void Lau2DAbsHistDP::raiseOrLowerBins( TH2* hist, const Double_t avEff, const Double_t avEffError )
{
    TRandom* random = LauRandom::randomFun();

    Double_t curAvg = this->computeAverageContents( hist );
    Double_t newAvg = random->Gaus( avEff, avEffError );

    hist->Scale( newAvg / curAvg );
}

Double_t Lau2DAbsHistDP::computeAverageContents( const TH2* hist ) const
{
    Double_t totalContent( 0.0 );
    Int_t binsWithinDPBoundary( 0 );

    Int_t nBinsX = static_cast<Int_t>( hist->GetNbinsX() );
    Int_t nBinsY = static_cast<Int_t>( hist->GetNbinsY() );

    // Loop through the bins and include any that have their centre or any
    // of the four corners within the kinematic boundary
    for ( Int_t i( 0 ); i < nBinsX; ++i ) {
        Double_t binXCentre    = hist->GetXaxis()->GetBinCenter( i + 1 );
        Double_t binXLowerEdge = hist->GetXaxis()->GetBinLowEdge( i + 1 );
        Double_t binXUpperEdge = hist->GetXaxis()->GetBinUpEdge( i + 1 );

        for ( Int_t j( 0 ); j < nBinsY; ++j ) {
            Double_t binYCentre    = hist->GetYaxis()->GetBinCenter( i + 1 );
            Double_t binYLowerEdge = hist->GetYaxis()->GetBinLowEdge( i + 1 );
            Double_t binYUpperEdge = hist->GetYaxis()->GetBinUpEdge( i + 1 );

            if ( this->withinDPBoundaries( binXCentre, binYCentre ) ||
                 this->withinDPBoundaries( binXLowerEdge, binYLowerEdge ) ||
                 this->withinDPBoundaries( binXUpperEdge, binYUpperEdge ) ||
                 this->withinDPBoundaries( binXLowerEdge, binYUpperEdge ) ||
                 this->withinDPBoundaries( binXUpperEdge, binYLowerEdge ) ) {

                totalContent += hist->GetBinContent( i + 1, j + 1 );
                ++binsWithinDPBoundary;
            }
        }
    }

    return totalContent / binsWithinDPBoundary;
}

Bool_t Lau2DAbsHistDP::withinDPBoundaries( Double_t x, Double_t y ) const
{
    return squareDP_ ? kinematics_->withinSqDPLimits( x, y ) : kinematics_->withinDPLimits( x, y );
}

void Lau2DAbsHistDP::getUpperHalf( Double_t& x, Double_t& y ) const
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
