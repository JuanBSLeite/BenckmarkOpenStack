
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

// Code to plot the transition amplitude T for the K-matrix pi pi S-wave

#include <iostream>
using std::cout;
using std::endl;

#include "LauComplex.hh"
#include "LauKMatrixPropagator.hh"

#include "TArrow.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

int main( const int /*argc*/, const char** /*argv*/ )
{

    // Define the "S-wave" K-matrix propagator
    Int_t nChannels                  = 5;
    Int_t nPoles                     = 5;
    Int_t resPairInt                 = 1;
    Int_t KMatrixIndex               = 1;    // for S-wave
    LauKMatrixPropagator* propagator = new LauKMatrixPropagator( "KMSWave",
                                                                 "DToKspipiKMatrixCoeff.dat",
                                                                 resPairInt,
                                                                 nChannels,
                                                                 nPoles,
                                                                 KMatrixIndex );

    // Find the transition amplitude T for the pi pi mode for
    // invariant mass squared s between (2*mpi)^2 and 4.0 GeV^2
    Double_t mMin( 0.28 );
    Double_t mMax( 2.0 );
    Int_t nS    = 2000;
    Double_t dm = ( mMax - mMin ) / ( nS * 1.0 );

    // Argand plot
    TGraph* TxyGraph = new TGraph( nS );
    // Intensity
    TGraph* TSqGraph = new TGraph( nS );
    // Phase
    TGraph* TDeltaGraph = new TGraph( nS );
    // Inelasticity eta
    TGraph* TEtaGraph = new TGraph( nS );

    Double_t radToDeg = 180.0 / TMath::Pi();
    // Keep track of the real value of T to know when it
    // changes sign from -ve to +ve => 180 phase shift
    Double_t oldReT( 0.0 );
    int n180Shift( 0 );

    for ( Int_t i = 0; i < nS; i++ ) {

        Double_t m = dm * i + mMin;
        Double_t s = m * m;

        LauComplex TAmp = propagator->getTransitionAmp( s, KMatrixIndex );

        Double_t realT = TAmp.re();
        Double_t imagT = TAmp.im();

        // Argand diagram showing the transition amplitude phase motion
        TxyGraph->SetPoint( i, realT, imagT );
        // Transition amplitude squared
        TSqGraph->SetPoint( i, m, TAmp.abs2() );

        // T - E = 0.5*i, where E = inelasticity vector, pointing to T from (0,i/2)
        Double_t ReE = realT;
        Double_t ImE = imagT - 0.5;

        LauComplex EVect( ReE, ImE );
        // Inelasticity
        Double_t EMag = EVect.abs();
        Double_t eta  = 2.0 * EMag;
        TEtaGraph->SetPoint( i, m, eta );

        // Find the phase shift angle, delta
        Double_t delta = 0.5 * TMath::ATan2( ImE, fabs( ReE ) ) * radToDeg + 45.0;
        if ( realT < 0.0 ) {
            delta = 180.0 - delta;
        }

        // Have we gone through 180 deg (or 2*delta through 360 deg)?
        if ( oldReT < 0.0 && realT > 0.0 ) {
            n180Shift += 1;
        }
        // Add 180 if required
        delta += 180.0 * n180Shift;

        TDeltaGraph->SetPoint( i, m, delta );

        // Keep track of the real value of T to see if it changes sign
        oldReT = realT;
    }

    TCanvas* theCanvas = new TCanvas( "theCanvas", "", 1400, 1000 );
    gROOT->SetStyle( "Plain" );
    gStyle->SetOptStat( 0 );
    theCanvas->Clear();
    theCanvas->UseCurrentStyle();

    theCanvas->Divide( 2, 2 );
    // Argand plot
    theCanvas->cd( 1 );
    TxyGraph->SetTitle( "Transition amplitude T#; S #equiv I + 2iT" );
    TxyGraph->GetXaxis()->SetTitle( "Re(T)" );
    TxyGraph->GetYaxis()->SetTitle( "Im(T)" );
    TxyGraph->GetXaxis()->CenterTitle( kTRUE );
    TxyGraph->GetYaxis()->CenterTitle( kTRUE );
    TxyGraph->SetLineWidth( 1 );
    TxyGraph->Draw( "apl" );

    // Overlay the elastic boundaries
    TLine lineA1( -0.5, 0.0, -0.5, 1.0 );
    lineA1.SetLineStyle( kDotted );
    lineA1.SetLineWidth( 1 );
    lineA1.Draw();
    TLine lineA2( 0.5, 0.0, 0.5, 1.0 );
    lineA2.SetLineStyle( kDotted );
    lineA2.SetLineWidth( 1 );
    lineA2.Draw();
    TLine lineA3( -0.5, 1.0, 0.5, 1.0 );
    lineA3.SetLineStyle( kDotted );
    lineA3.SetLineWidth( 1 );
    lineA3.Draw();

    // Intensity
    theCanvas->cd( 2 );
    TSqGraph->SetTitle( "Transition amplitude intensity" );
    TSqGraph->GetXaxis()->SetTitle( "m(#pi^{+}#pi^{-}) (GeV/c^{2})" );
    TSqGraph->GetYaxis()->SetTitle( "|T|^{2}" );
    TSqGraph->SetLineWidth( 1 );
    TSqGraph->GetXaxis()->CenterTitle( kTRUE );
    TSqGraph->GetYaxis()->CenterTitle( kTRUE );
    TSqGraph->Draw( "apl" );
    // Overlay the unitarity limit
    TAxis* xAxis = TSqGraph->GetXaxis();
    TLine lineB1( xAxis->GetXmin(), 1.0, xAxis->GetXmax(), 1.0 );
    lineB1.SetLineStyle( kDotted );
    lineB1.SetLineWidth( 1 );
    lineB1.Draw();

    // Also overlay "interesting regions"
    TLatex text;
    text.SetTextSize( 0.045 );

    // sigma
    text.DrawLatex( 0.65, 0.6, "f_{0}(500)" );
    text.DrawLatex( 0.68, 0.525, "or #sigma" );

    // f0(980)
    text.SetTextSize( 0.045 );
    text.DrawLatex( 0.9, 1.02, "f_{0}(980)" );
    TArrow arrow1( 0.98, 1.0, 0.98, 0.8, 0.01 );
    arrow1.SetLineWidth( 1 );
    arrow1.Draw();

    // f0(1370)
    text.DrawLatex( 1.25, 0.91, "f_{0}(1370)" );
    TArrow arrow2( 1.37, 0.9, 1.37, 0.7, 0.01 );
    arrow2.SetLineWidth( 1 );
    arrow2.Draw();

    // f0(1500)
    text.DrawLatex( 1.45, 0.76, "f_{0}(1500)" );
    TArrow arrow3( 1.50, 0.72, 1.50, 0.52, 0.01 );
    arrow3.SetLineWidth( 1 );
    arrow3.Draw();

    // f0(1710)
    text.DrawLatex( 1.60, 0.61, "f_{0}(1710)" );
    TArrow arrow4( 1.71, 0.6, 1.71, 0.45, 0.01 );
    arrow4.SetLineWidth( 1 );
    arrow4.Draw();

    // Phase
    theCanvas->cd( 3 );
    TDeltaGraph->SetTitle( "Phase shift" );
    TDeltaGraph->GetXaxis()->SetTitle( "m(#pi^{+}#pi^{-}) (GeV/c^{2})" );
    TDeltaGraph->GetYaxis()->SetTitle( "#delta (degrees)" );
    TDeltaGraph->GetXaxis()->CenterTitle( kTRUE );
    TDeltaGraph->GetYaxis()->CenterTitle( kTRUE );
    TDeltaGraph->SetLineWidth( 1 );
    TDeltaGraph->Draw( "apl" );

    // Inelasticity
    theCanvas->cd( 4 );
    TEtaGraph->SetTitle( "Inelasticity, #eta #equiv |2T - iI| = |S|" );
    TEtaGraph->GetXaxis()->SetTitle( "m(#pi^{+}#pi^{-}) (GeV/c^{2})" );
    TEtaGraph->GetYaxis()->SetTitle( "#eta" );
    TEtaGraph->GetXaxis()->CenterTitle( kTRUE );
    TEtaGraph->GetYaxis()->CenterTitle( kTRUE );
    TEtaGraph->SetLineWidth( 1 );
    TEtaGraph->Draw( "apl" );

    // Overlay the elastic limit
    xAxis = TEtaGraph->GetXaxis();
    TLine lineD1( xAxis->GetXmin(), 1.0, xAxis->GetXmax(), 1.0 );
    lineD1.SetLineStyle( kDotted );
    lineD1.SetLineWidth( 1 );
    lineD1.Draw();

    theCanvas->Print( "KMatrixTAmpPlots.png" );
    theCanvas->Print( "KMatrixTAmpPlots.eps" );
}
