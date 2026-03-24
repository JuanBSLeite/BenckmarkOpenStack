
/*
Copyright 2016 University of Warwick

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

// Code to plot amplitudes for the pi pi K matrix used for B -> 3pi:
// K-matrix itself, propagator terms and production amplitudes

#include <iostream>
using std::cout;
using std::endl;

#include "LauComplex.hh"
#include "LauDaughters.hh"
#include "LauKMatrixProdPole.hh"
#include "LauKMatrixProdSVP.hh"
#include "LauKMatrixPropagator.hh"

#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLine.h"
#include "TList.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMultiGraph.h"
#include "TPad.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TText.h"

#include <algorithm>
#include <string>
#include <vector>

Double_t radToDeg  = 180.0 / TMath::Pi();
Bool_t doEPS       = kFALSE;
Bool_t doProdAdler = kFALSE;

void doKScattPlots()
{

    // Create plots of K_ij (with i >= j, since its symmetric), showing
    // only the pole and SVP contributions separately, along with the full
    // K matrix without and with the Adler zero suppression term for low s.
    // 1st page: K_1j (5), 2nd page: K_2j (4), 3rd page: K_3j, K_4j and K_55 (6)

    Int_t resPairInt( 1 );
    Int_t nChannels( 5 );
    Int_t KMatrixIndex( 1 );    // pi pi S-wave
    Int_t nPoles( 5 );

    // Mass range for plots
    Double_t mMin( 0.28 );
    Double_t mMax( 2.0 );
    Int_t nSteps = 2000;
    Double_t dm  = ( mMax - mMin ) / ( nSteps * 1.0 );

    // Propagator with only the pole terms
    TString KFile1( "B3piKMPoles.dat" );
    LauKMatrixPropagator* propPoles =
        new LauKMatrixPropagator( "propPoles", KFile1, resPairInt, nChannels, nPoles, KMatrixIndex );

    // Propagator with only the SVP terms
    TString KFile2( "B3piKMSVPs.dat" );
    LauKMatrixPropagator* propSVPs =
        new LauKMatrixPropagator( "propSVPs", KFile2, resPairInt, nChannels, 0, KMatrixIndex );

    // Propagator with poles and SVP terms, but no Adler zero scaling factor
    TString KFile3( "B3piKMNoAdler.dat" );
    LauKMatrixPropagator* propNoAdler =
        new LauKMatrixPropagator( "propNoAdler", KFile3, resPairInt, nChannels, nPoles, KMatrixIndex );

    // Full K-matrix parameterisation with Adler zero
    TString KFile4( "B3piKMatrixCoeff.dat" );
    LauKMatrixPropagator* propFull =
        new LauKMatrixPropagator( "propFull", KFile4, resPairInt, nChannels, nPoles, KMatrixIndex );

    // For each K_ij, there are 4 graphs. Therefore, need 15 K_ij multigraphs
    std::vector<TMultiGraph*> multiGraphVect;

    // Legend for the graphs (same for all multigraphs)
    TLegend* theLegend = new TLegend( 0.2, 0.2, 0.8, 0.7 );
    theLegend->SetHeader( "Symmetric K matrix components" );

    Int_t iColours[4] = { kRed, kBlue, kMagenta, kBlack };
    Double_t yLimit( 5.0 );

    Bool_t initLegend( kFALSE );

    Int_t iChannel( 0 ), jChannel( 0 );
    for ( iChannel = 0; iChannel < nChannels; iChannel++ ) {

        for ( jChannel = iChannel; jChannel < nChannels; jChannel++ ) {

            TString label( "Graphs_" );
            label += iChannel;
            label += "_";
            label += jChannel;

            TString multiName( "KMatrix" );
            multiName += label;

            TMultiGraph* multiGraph = new TMultiGraph();
            multiGraph->SetName( multiName.Data() );

            // To limit the graphs abscissas getting to large (pole "singularities")
            multiGraph->SetMaximum( yLimit );
            multiGraph->SetMinimum( -yLimit );

            // Add to the storage vector
            multiGraphVect.push_back( multiGraph );

            // Graph with poles only
            TGraph* poleGraph = new TGraph( nSteps );
            TString poleName( "Pole" );
            poleName += label;
            poleGraph->SetName( poleName.Data() );
            poleGraph->SetMarkerColor( iColours[0] );
            poleGraph->SetLineColor( iColours[0] );
            poleGraph->SetLineWidth( 1 );
            multiGraph->Add( poleGraph, "lp" );

            // Graph with SVPs only
            TGraph* SVPGraph = new TGraph( nSteps );
            TString SVPName( "SVP" );
            SVPName += label;
            SVPGraph->SetName( SVPName.Data() );
            SVPGraph->SetMarkerColor( iColours[1] );
            SVPGraph->SetLineColor( iColours[1] );
            SVPGraph->SetLineWidth( 1 );
            multiGraph->Add( SVPGraph, "lp" );

            // Graph with poles and SVP, but no Adler zero factor
            TGraph* noAdlerGraph = new TGraph( nSteps );
            TString noAdlerName( "NoAdler" );
            noAdlerName += label;
            noAdlerGraph->SetName( noAdlerName.Data() );
            noAdlerGraph->SetMarkerColor( iColours[2] );
            noAdlerGraph->SetLineColor( iColours[2] );
            noAdlerGraph->SetLineWidth( 1 );
            multiGraph->Add( noAdlerGraph, "lp" );

            // Full K matrix graphs
            TGraph* fullGraph = new TGraph( nSteps );
            TString fullName( "Full" );
            fullName += label;
            fullGraph->SetName( fullName.Data() );
            fullGraph->SetMarkerColor( iColours[3] );
            fullGraph->SetLineColor( iColours[3] );
            fullGraph->SetLineWidth( 1 );
            multiGraph->Add( fullGraph, "lp" );

            // Initialise the graph legend
            if ( ! initLegend ) {

                theLegend->AddEntry( poleGraph, "Poles only", "l" );
                theLegend->AddEntry( SVPGraph, "SVPs only (f_{1j} #neq 0)", "l" );
                theLegend->AddEntry( noAdlerGraph, "Poles & SVPs, no Adler zero", "l" );
                theLegend->AddEntry( fullGraph, "Poles & SVPs with Adler zero", "l" );

                initLegend = kTRUE;
            }
        }
    }

    for ( Int_t iM = 0; iM < nSteps; iM++ ) {

        Double_t m = dm * iM + mMin;
        Double_t s = m * m;

        // Find the K matrices for each pole/SVP/Adler condition
        propPoles->updatePropagator( s );
        propSVPs->updatePropagator( s );
        propNoAdler->updatePropagator( s );
        propFull->updatePropagator( s );

        TMatrixD KMPoles   = propPoles->getKMatrix();
        TMatrixD KMSVPs    = propSVPs->getKMatrix();
        TMatrixD KMNoAdler = propNoAdler->getKMatrix();
        TMatrixD KMFull    = propFull->getKMatrix();

        Int_t vectIndex( -1 );

        for ( iChannel = 0; iChannel < nChannels; iChannel++ ) {

            for ( jChannel = iChannel; jChannel < nChannels; jChannel++ ) {

                vectIndex++;
                TMultiGraph* theGraphs = multiGraphVect[vectIndex];
                TList* graphList       = theGraphs->GetListOfGraphs();

                TGraph* poleGraph = dynamic_cast<TGraph*>( graphList->At( 0 ) );
                poleGraph->SetPoint( iM, m, KMPoles[iChannel][jChannel] );

                TGraph* SVPGraph = dynamic_cast<TGraph*>( graphList->At( 1 ) );
                SVPGraph->SetPoint( iM, m, KMSVPs[iChannel][jChannel] );

                TGraph* noAdlerGraph = dynamic_cast<TGraph*>( graphList->At( 2 ) );
                noAdlerGraph->SetPoint( iM, m, KMNoAdler[iChannel][jChannel] );

                TGraph* fullGraph = dynamic_cast<TGraph*>( graphList->At( 3 ) );
                fullGraph->SetPoint( iM, m, KMFull[iChannel][jChannel] );
            }
        }
    }

    // Create plots
    TCanvas* theCanvas = new TCanvas( "theCanvas", "", 1400, 1000 );
    gROOT->SetStyle( "Plain" );
    gStyle->SetOptStat( 0 );
    theCanvas->Clear();
    theCanvas->UseCurrentStyle();

    // K_1j graphs, j = 1 to 5
    theCanvas->Divide( 2, 3 );

    for ( Int_t iL = 0; iL < 5; iL++ ) {

        theCanvas->cd( iL + 1 );
        TMultiGraph* graphs = multiGraphVect[iL];
        graphs->Draw( "a" );

        // Update the axes, which can only be done after drawing the graphs
        TAxis* xAxis = graphs->GetXaxis();
        xAxis->SetTitle( "m(#pi^{+}#pi^{-}) (GeV/c^{2})" );
        xAxis->SetLabelSize( 0.05 );
        xAxis->SetTitleSize( 0.05 );
        xAxis->CenterTitle( kTRUE );
        xAxis->SetTitleOffset( 0.9 );

        TString yLabel( "K(1," );
        yLabel += ( iL + 1 );
        yLabel += ")";

        TAxis* yAxis = graphs->GetYaxis();
        yAxis->SetTitleOffset( 0.9 );
        yAxis->SetTitle( yLabel.Data() );
        yAxis->SetLabelSize( 0.05 );
        yAxis->SetTitleSize( 0.05 );
        yAxis->CenterTitle( kTRUE );

        // Update the canvas with the new axes settings
        theCanvas->Update();
    }

    // Print the graph legend
    theCanvas->cd( 6 );
    theLegend->Draw();

    theCanvas->Print( "K_1jPlots.png" );
    if ( doEPS ) {
        theCanvas->Print( "K_1jPlots.eps" );
    }

    // K_2j graphs, j = 2 to 5
    theCanvas->Clear();
    theCanvas->Divide( 2, 2 );

    for ( Int_t iL = 0; iL < 4; iL++ ) {

        theCanvas->cd( iL + 1 );

        TMultiGraph* graphs = multiGraphVect[iL + 5];
        graphs->Draw( "a" );

        // Update the axes, which can only be done after drawing the graphs
        TAxis* xAxis = graphs->GetXaxis();
        xAxis->SetTitle( "m(#pi^{+}#pi^{-}) (GeV/c^{2})" );
        xAxis->SetLabelSize( 0.05 );
        xAxis->SetTitleSize( 0.05 );
        xAxis->CenterTitle( kTRUE );
        xAxis->SetTitleOffset( 0.9 );

        TString yLabel( "K(2," );
        yLabel += ( iL + 2 );
        yLabel += ")";

        TAxis* yAxis = graphs->GetYaxis();
        yAxis->SetTitleOffset( 0.9 );
        yAxis->SetTitle( yLabel.Data() );
        yAxis->SetLabelSize( 0.05 );
        yAxis->SetTitleSize( 0.05 );
        yAxis->CenterTitle( kTRUE );

        // Update the canvas with the new axes settings
        theCanvas->Update();
    }

    theCanvas->Print( "K_2jPlots.png" );
    if ( doEPS ) {
        theCanvas->Print( "K_2jPlots.eps" );
    }

    // K_3j (j = 3,4,5), K_4j (j = 4,5) and K_55
    theCanvas->Clear();
    theCanvas->Divide( 2, 3 );

    for ( Int_t iL = 0; iL < 6; iL++ ) {

        theCanvas->cd( iL + 1 );
        TMultiGraph* graphs = multiGraphVect[iL + 9];
        graphs->Draw( "a" );

        // Update the axes, which can only be done after drawing the graphs
        TAxis* xAxis = graphs->GetXaxis();
        xAxis->SetTitle( "m(#pi^{+}#pi^{-}) (GeV/c^{2})" );
        xAxis->SetLabelSize( 0.05 );
        xAxis->SetTitleSize( 0.05 );
        xAxis->CenterTitle( kTRUE );
        xAxis->SetTitleOffset( 0.9 );

        TString yLabel( "K(" );
        if ( iL == 0 ) {
            yLabel += "3,3)";
        } else if ( iL == 1 ) {
            yLabel += "3,4)";
        } else if ( iL == 2 ) {
            yLabel += "3,5)";
        } else if ( iL == 3 ) {
            yLabel += "4,4)";
        } else if ( iL == 4 ) {
            yLabel += "4,5)";
        } else {
            yLabel += "5,5)";
        }

        TAxis* yAxis = graphs->GetYaxis();
        yAxis->SetTitleOffset( 0.9 );
        yAxis->SetTitle( yLabel.Data() );
        yAxis->SetLabelSize( 0.05 );
        yAxis->SetTitleSize( 0.05 );
        yAxis->CenterTitle( kTRUE );

        // Update the canvas with the new axes settings
        theCanvas->Update();
    }

    theCanvas->Print( "K_345jPlots.png" );
    if ( doEPS ) {
        theCanvas->Print( "K_345jPlots.eps" );
    }

    // Cleanup
    delete propPoles;
    delete propSVPs;
    delete propNoAdler;
    delete propFull;

    std::vector<TMultiGraph*>::iterator multiIter;
    for ( multiIter = multiGraphVect.begin(); multiIter != multiGraphVect.end(); ++multiIter ) {

        TMultiGraph* multiGraph = *multiIter;
        delete multiGraph;
    }

    delete theLegend;
    delete theCanvas;
}

void doPropagatorPlots()
{

    // Create plots of the s variation of the elements of the
    // propagator (I - i K rho)^-1

    Int_t resPairInt( 1 );
    Int_t nChannels( 5 );
    Int_t KMatrixIndex( 1 );    // pi pi S-wave
    Int_t nPoles( 5 );

    // Mass range for plots
    Double_t mMin( 0.28 );
    Double_t mMax( 2.0 );
    Int_t nSteps = 2000;
    Double_t dm  = ( mMax - mMin ) / ( nSteps * 1.0 );

    // Find the point numbers for m = mMin, 0.5, 0.75, 1, 1.25, 1.5, 1.75 and 2 GeV
    std::vector<Int_t> points;
    std::vector<Int_t> colours, styles;
    std::vector<Double_t> massPoints;
    massPoints.push_back( mMin );
    massPoints.push_back( 0.5 );
    massPoints.push_back( 0.75 );
    massPoints.push_back( 1.0 );
    massPoints.push_back( 1.25 );
    massPoints.push_back( 1.5 );
    massPoints.push_back( 1.75 );
    massPoints.push_back( 1.99999 );

    Int_t NMain = static_cast<Int_t>( massPoints.size() );
    for ( Int_t iP = 0; iP < NMain; iP++ ) {

        Int_t pointInt = static_cast<Int_t>( ( massPoints[iP] - mMin ) / dm );
        points.push_back( pointInt );

        // Red or Blue
        Int_t iColour = kRed;
        if ( iP > 3 ) {
            iColour = kBlue;
        }
        Int_t iStyle = kFullCircle;
        Int_t iPMod  = iP % 4;
        if ( iPMod == 1 ) {
            iStyle = kOpenCircle;
        } else if ( iPMod == 2 ) {
            iStyle = kFullSquare;
        } else if ( iPMod == 3 ) {
            iStyle = kOpenSquare;
        }

        colours.push_back( iColour );
        styles.push_back( iStyle );
    }

    // Full K-matrix parameterisation
    TString KFile( "B3piKMatrixCoeff.dat" );
    LauKMatrixPropagator* propFull =
        new LauKMatrixPropagator( "propFull", KFile, resPairInt, nChannels, nPoles, KMatrixIndex );

    // Maps of <channel row, vector of graph for channel column>

    // Argand diagram graphs (full plot, then others showing only selected points)
    std::map<int, std::vector<TMultiGraph*>> argMap;
    // Real, imaginary and magnitude amplitude multigraphs
    std::map<int, std::vector<TMultiGraph*>> ampMap;

    // TLegend for argand graphs
    TLegend* argLegend = new TLegend( 0.1, 0.1, 0.7, 0.995 );
    argLegend->SetHeader( "m(#pi^{+}#pi^{-}) points (GeV/c^{2})" );
    argLegend->SetTextSize( 0.075 );
    Bool_t initArgLegend( kFALSE );

    // Legend for the amplitude graphs (same for all multigraphs)
    TLegend* ampLegend = new TLegend( 0.2, 0.2, 0.8, 0.7 );
    ampLegend->SetHeader( "Propagator (I - iK#rho)^{-1} components" );
    Bool_t initAmpLegend( kFALSE );

    Int_t iChannel( 0 ), jChannel( 0 );

    for ( iChannel = 0; iChannel < nChannels; iChannel++ ) {

        std::vector<TMultiGraph*> argGraphVect;
        std::vector<TMultiGraph*> ampGraphVect;

        for ( jChannel = 0; jChannel < nChannels; jChannel++ ) {

            // Common numbering scheme
            TString numbers( "" );
            numbers += iChannel;
            numbers += "_";
            numbers += jChannel;

            // Argand plot showing phase motion
            TMultiGraph* argGraphs = new TMultiGraph();
            TString argMultiLabel( "PropMultiArg_" );
            argMultiLabel += numbers;
            argGraphs->SetName( argMultiLabel.Data() );
            argGraphVect.push_back( argGraphs );

            TGraph* argGraph = new TGraph( nSteps );
            TString argLabel( "PropArg_" );
            argLabel += numbers;
            argGraph->SetName( argLabel.Data() );
            argGraph->SetLineWidth( 1 );
            argGraphs->Add( argGraph, "lp" );

            // Argand graphs showing only one highlighted point
            for ( Int_t iG = 0; iG < NMain; iG++ ) {

                TGraph* pointGraph = new TGraph( 1 );
                TString pointLabel( "PropPoint_" );
                pointLabel += iG;
                pointLabel += "_";
                pointLabel += numbers;

                pointGraph->SetName( pointLabel.Data() );
                pointGraph->SetMarkerColor( colours[iG] );
                pointGraph->SetMarkerSize( 1.1 );
                pointGraph->SetMarkerStyle( styles[iG] );
                argGraphs->Add( pointGraph, "p" );
            }

            if ( ! initArgLegend ) {
                TList* argList = argGraphs->GetListOfGraphs();
                for ( Int_t iG = 0; iG < NMain; iG++ ) {
                    TGraph* pointGraph = dynamic_cast<TGraph*>( argList->At( iG + 1 ) );
                    argLegend->AddEntry( pointGraph, Form( "%.2f", massPoints[iG] ), "p" );
                    initArgLegend = kTRUE;
                }
            }

            // Real, imaginary and magnitude amplitude
            TMultiGraph* ampGraphs = new TMultiGraph();
            TString multiLabel( "PropMultiAmp_" );
            multiLabel += numbers;
            ampGraphs->SetName( multiLabel.Data() );
            ampGraphVect.push_back( ampGraphs );

            TGraph* realGraph = new TGraph( nSteps );
            TString realLabel( "PropReal_" );
            realLabel += numbers;
            realGraph->SetName( realLabel.Data() );
            realGraph->SetLineColor( kRed );
            realGraph->SetMarkerColor( kRed );
            realGraph->SetLineWidth( 1 );
            ampGraphs->Add( realGraph, "lp" );

            TGraph* imagGraph = new TGraph( nSteps );
            TString imagLabel( "PropImag_" );
            imagLabel += numbers;
            imagGraph->SetName( imagLabel.Data() );
            imagGraph->SetLineColor( kBlue );
            imagGraph->SetMarkerColor( kBlue );
            imagGraph->SetLineWidth( 1 );
            ampGraphs->Add( imagGraph, "lp" );

            TGraph* magGraph = new TGraph( nSteps );
            TString magLabel( "PropMag_" );
            magLabel += numbers;
            magGraph->SetName( magLabel.Data() );
            magGraph->SetLineColor( kBlack );
            magGraph->SetMarkerColor( kBlack );
            magGraph->SetLineWidth( 1 );
            ampGraphs->Add( magGraph, "lp" );

            if ( ! initAmpLegend ) {

                ampLegend->AddEntry( realGraph, "Real part", "l" );
                ampLegend->AddEntry( imagGraph, "Imag part", "l" );
                ampLegend->AddEntry( magGraph, "Magnitude", "l" );

                initAmpLegend = kTRUE;
            }
        }

        argMap[iChannel] = argGraphVect;
        ampMap[iChannel] = ampGraphVect;
    }

    // Fill the graphs
    for ( Int_t iM = 0; iM < nSteps; iM++ ) {

        Double_t m = dm * iM + mMin;
        Double_t s = m * m;

        propFull->updatePropagator( s );

        TMatrixD realPropMatrix    = propFull->getRealPropMatrix();
        TMatrixD negImagPropMatrix = propFull->getNegImagPropMatrix();

        for ( iChannel = 0; iChannel < nChannels; iChannel++ ) {

            std::vector<TMultiGraph*> argGraphVect = argMap[iChannel];
            std::vector<TMultiGraph*> ampGraphVect = ampMap[iChannel];

            for ( jChannel = 0; jChannel < nChannels; jChannel++ ) {

                // Get the real and imaginary parts of the propagator term
                Double_t realProp = realPropMatrix[iChannel][jChannel];
                Double_t imagProp = -negImagPropMatrix[iChannel][jChannel];
                Double_t magProp  = sqrt( realProp * realProp + imagProp * imagProp );

                // Argand diagram: full plot and only selected points
                TMultiGraph* argGraphs = argGraphVect[jChannel];
                TList* argList         = argGraphs->GetListOfGraphs();

                TGraph* argGraph = dynamic_cast<TGraph*>( argList->At( 0 ) );
                argGraph->SetPoint( iM, realProp, imagProp );

                // Graphs showing specific points
                for ( Int_t iP = 0; iP < NMain; iP++ ) {

                    if ( points[iP] == iM ) {
                        TGraph* pointGraph = dynamic_cast<TGraph*>( argList->At( iP + 1 ) );
                        pointGraph->SetPoint( 0, realProp, imagProp );
                    }
                }

                // Amplitude multigraphs
                TMultiGraph* ampGraphs = ampGraphVect[jChannel];
                TList* ampList         = ampGraphs->GetListOfGraphs();

                TGraph* realGraph = dynamic_cast<TGraph*>( ampList->At( 0 ) );
                realGraph->SetPoint( iM, m, realProp );

                TGraph* imagGraph = dynamic_cast<TGraph*>( ampList->At( 1 ) );
                imagGraph->SetPoint( iM, m, imagProp );

                TGraph* magGraph = dynamic_cast<TGraph*>( ampList->At( 2 ) );
                magGraph->SetPoint( iM, m, magProp );
            }
        }
    }

    // Create plots
    TCanvas* theCanvas = new TCanvas( "theCanvas", "", 1000, 1000 );
    gROOT->SetStyle( "Plain" );
    gStyle->SetOptStat( 0 );

    for ( iChannel = 0; iChannel < nChannels; iChannel++ ) {

        std::vector<TMultiGraph*> argGraphVect = argMap[iChannel];
        std::vector<TMultiGraph*> ampGraphVect = ampMap[iChannel];

        theCanvas->Clear();
        theCanvas->UseCurrentStyle();
        theCanvas->Divide( 2, 3 );

        // Argand plots
        for ( jChannel = 0; jChannel < nChannels; jChannel++ ) {

            theCanvas->cd( jChannel + 1 );
            TMultiGraph* argGraphs = argGraphVect[jChannel];
            argGraphs->Draw( "a" );

            TString title( "Propagator(" );
            title += iChannel + 1;
            title += ",";
            title += jChannel + 1;
            title += ")";
            argGraphs->SetTitle( title.Data() );

            TAxis* xAxis = argGraphs->GetXaxis();
            xAxis->SetTitle( "Real" );
            xAxis->SetLabelSize( 0.05 );
            xAxis->SetTitleSize( 0.05 );
            xAxis->CenterTitle( kTRUE );
            xAxis->SetTitleOffset( 0.9 );

            TAxis* yAxis = argGraphs->GetYaxis();
            yAxis->SetTitle( "Imag" );
            yAxis->SetLabelSize( 0.05 );
            yAxis->SetTitleSize( 0.05 );
            yAxis->CenterTitle( kTRUE );

            theCanvas->Update();
        }

        theCanvas->cd( 6 );
        argLegend->Draw();

        TString argPlotName( "PropArg_" );
        argPlotName += iChannel;
        argPlotName += ".png";
        theCanvas->Print( argPlotName.Data() );

        if ( doEPS ) {
            TString argEpsPlotName( argPlotName );
            argEpsPlotName.ReplaceAll( ".png", ".eps" );
            theCanvas->Print( argEpsPlotName.Data() );
        }

        // Amplitude plots
        theCanvas->Clear();
        theCanvas->UseCurrentStyle();
        theCanvas->Divide( 2, 3 );

        for ( jChannel = 0; jChannel < nChannels; jChannel++ ) {

            theCanvas->cd( jChannel + 1 );

            TMultiGraph* ampGraphs = ampGraphVect[jChannel];
            ampGraphs->Draw( "a" );

            TAxis* xAxis = ampGraphs->GetXaxis();
            xAxis->SetTitle( "m(#pi^{+}#pi^{-}) (GeV/c^{2})" );
            xAxis->SetLabelSize( 0.05 );
            xAxis->SetTitleSize( 0.05 );
            xAxis->CenterTitle( kTRUE );
            xAxis->SetTitleOffset( 0.9 );

            TString yLabel( "Propagator(" );
            yLabel += iChannel + 1;
            yLabel += ",";
            yLabel += jChannel + 1;
            yLabel += ")";

            TAxis* yAxis = ampGraphs->GetYaxis();
            yAxis->SetTitleOffset( 0.9 );
            yAxis->SetTitle( yLabel.Data() );
            yAxis->SetLabelSize( 0.05 );
            yAxis->SetTitleSize( 0.05 );
            yAxis->CenterTitle( kTRUE );

            theCanvas->Update();
        }

        theCanvas->cd( 6 );
        ampLegend->Draw();

        TString ampPlotName( "PropAmp_" );
        ampPlotName += iChannel;
        ampPlotName += ".png";
        theCanvas->Print( ampPlotName.Data() );

        if ( doEPS ) {
            TString ampEpsPlotName( ampPlotName );
            ampEpsPlotName.ReplaceAll( ".png", ".eps" );
            theCanvas->Print( ampEpsPlotName.Data() );
        }
    }

    // Cleanup
    delete propFull;

    std::map<int, std::vector<TMultiGraph*>>::iterator argIter;
    for ( argIter = argMap.begin(); argIter != argMap.end(); ++argIter ) {

        std::vector<TMultiGraph*> argGraphs = argIter->second;

        std::vector<TMultiGraph*>::iterator theIter;
        for ( theIter = argGraphs.begin(); theIter != argGraphs.end(); ++theIter ) {

            TMultiGraph* aGraph = *theIter;
            delete aGraph;
        }

        argGraphs.clear();
    }

    std::map<int, std::vector<TMultiGraph*>>::iterator ampIter;
    for ( ampIter = ampMap.begin(); ampIter != ampMap.end(); ++ampIter ) {

        std::vector<TMultiGraph*> ampGraphs = ampIter->second;

        std::vector<TMultiGraph*>::iterator theIter;
        for ( theIter = ampGraphs.begin(); theIter != ampGraphs.end(); ++theIter ) {

            TMultiGraph* aGraph = *theIter;
            delete aGraph;
        }

        ampGraphs.clear();
    }

    delete ampLegend;
    delete theCanvas;
}

void doProdPolePlots( Bool_t useProdAdler )
{

    // Create plots of the pole production amplitudes using the first row of
    // the propagator (for pi pi S-wave). For each pole we sum over all channels

    Int_t resPairInt( 1 );
    Int_t nChannels( 5 );
    Int_t KMatrixIndex( 1 );    // pi pi S-wave
    Int_t nPoles( 5 );

    // Mass range for plots
    Double_t mMin( 0.28 );
    Double_t mMax( 2.0 );
    Int_t nSteps = 2000;
    Double_t dm  = ( mMax - mMin ) / ( nSteps * 1.0 );

    // Find the point numbers for m = mMin, 0.5, 0.75, 1, 1.25, 1.5, 1.75 and 2 GeV
    std::vector<Int_t> points;
    std::vector<Int_t> colours, styles;
    std::vector<Double_t> massPoints;
    massPoints.push_back( mMin );
    massPoints.push_back( 0.5 );
    massPoints.push_back( 0.75 );
    massPoints.push_back( 1.0 );
    massPoints.push_back( 1.25 );
    massPoints.push_back( 1.5 );
    massPoints.push_back( 1.75 );
    massPoints.push_back( 1.99999 );

    Int_t NMain = static_cast<Int_t>( massPoints.size() );
    for ( Int_t iP = 0; iP < NMain; iP++ ) {

        Int_t pointInt = static_cast<Int_t>( ( massPoints[iP] - mMin ) / dm );
        points.push_back( pointInt );

        // Red or Blue
        Int_t iColour = kRed;
        if ( iP > 3 ) {
            iColour = kBlue;
        }
        Int_t iStyle = kFullCircle;
        Int_t iPMod  = iP % 4;
        if ( iPMod == 1 ) {
            iStyle = kOpenCircle;
        } else if ( iPMod == 2 ) {
            iStyle = kFullSquare;
        } else if ( iPMod == 3 ) {
            iStyle = kOpenSquare;
        }

        colours.push_back( iColour );
        styles.push_back( iStyle );
    }

    Double_t mPi     = 0.13957018;
    Double_t mPiSq   = mPi * mPi;
    Double_t mB      = 5.279;
    Double_t mBSq    = mB * mB;
    Double_t mSqTerm = 2.0 * ( ( mBSq / 3.0 ) + mPiSq );

    // Full K-matrix parameterisation
    TString KFile( "B3piKMatrixCoeff.dat" );
    LauKMatrixPropagator* propFull =
        new LauKMatrixPropagator( "propFull", KFile, resPairInt, nChannels, nPoles, KMatrixIndex );

    // Define kinematics pointer
    Bool_t squareDP            = kTRUE;
    LauDaughters* theDaughters = new LauDaughters( "B+", "pi+", "pi+", "pi-", squareDP );

    // Kinematics
    LauKinematics* kinematics = theDaughters->getKinematics();

    std::vector<LauKMatrixProdPole*> poleVect;

    Int_t iPole( 0 );
    for ( iPole = 0; iPole < nPoles; iPole++ ) {

        TString poleName( "Pole" );
        poleName += iPole;

        LauKMatrixProdPole* thePole = new LauKMatrixProdPole( poleName,
                                                              iPole + 1,
                                                              resPairInt,
                                                              propFull,
                                                              theDaughters,
                                                              useProdAdler );

        poleVect.push_back( thePole );
    }

    // Argand diagram graphs (full plot, then others showing only selected points)
    std::vector<TMultiGraph*> argGraphVect;
    // Real, imaginary and magnitude amplitude multigraphs
    std::vector<TMultiGraph*> ampGraphVect;

    // TLegend for argand graphs
    TLegend* argLegend = new TLegend( 0.1, 0.1, 0.7, 0.995 );
    argLegend->SetHeader( "m(#pi^{+}#pi^{-}) points (GeV/c^{2})" );
    argLegend->SetTextSize( 0.075 );
    Bool_t initArgLegend( kFALSE );

    // Legend for the graphs (same for all multigraphs)
    TLegend* ampLegend = new TLegend( 0.2, 0.2, 0.8, 0.7 );
    ampLegend->SetHeader( "Production pole components" );
    Bool_t initLegend( kFALSE );

    for ( iPole = 0; iPole < nPoles; iPole++ ) {

        // Argand plot showing phase motion
        TMultiGraph* argGraphs = new TMultiGraph();
        TString argMultiLabel( "PoleMultiArg_" );
        argMultiLabel += iPole;
        argGraphs->SetName( argMultiLabel.Data() );
        argGraphVect.push_back( argGraphs );

        TGraph* argGraph = new TGraph( nSteps );
        TString argLabel( "PoleArg_" );
        argLabel += iPole;
        argGraph->SetName( argLabel.Data() );
        argGraph->SetLineWidth( 1 );
        argGraphs->Add( argGraph, "lp" );

        // Argand graphs showing only one highlighted point
        for ( Int_t iG = 0; iG < NMain; iG++ ) {

            TGraph* pointGraph = new TGraph( 1 );
            TString pointLabel( "PolePoint_" );
            pointLabel += iG;

            pointGraph->SetName( pointLabel.Data() );
            pointGraph->SetMarkerColor( colours[iG] );
            pointGraph->SetMarkerSize( 1.1 );
            pointGraph->SetMarkerStyle( styles[iG] );
            argGraphs->Add( pointGraph, "p" );
        }

        if ( ! initArgLegend ) {
            TList* argList = argGraphs->GetListOfGraphs();
            for ( Int_t iG = 0; iG < NMain; iG++ ) {
                TGraph* pointGraph = dynamic_cast<TGraph*>( argList->At( iG + 1 ) );
                argLegend->AddEntry( pointGraph, Form( "%.2f", massPoints[iG] ), "p" );
                initArgLegend = kTRUE;
            }
        }

        // Real, imaginary and magnitude amplitude
        TMultiGraph* ampGraphs = new TMultiGraph();
        TString multiLabel( "PoleMulti_" );
        multiLabel += iPole;
        ampGraphs->SetName( multiLabel.Data() );
        ampGraphVect.push_back( ampGraphs );

        TGraph* realGraph = new TGraph( nSteps );
        TString realLabel( "PoleReal_" );
        realLabel += iPole;
        realGraph->SetName( realLabel.Data() );
        realGraph->SetLineColor( kRed );
        realGraph->SetMarkerColor( kRed );
        realGraph->SetLineWidth( 1 );
        ampGraphs->Add( realGraph, "lp" );

        TGraph* imagGraph = new TGraph( nSteps );
        TString imagLabel( "PoleImag_" );
        imagLabel += iPole;
        imagGraph->SetName( imagLabel.Data() );
        imagGraph->SetLineColor( kBlue );
        imagGraph->SetMarkerColor( kBlue );
        imagGraph->SetLineWidth( 1 );
        ampGraphs->Add( imagGraph, "lp" );

        TGraph* magGraph = new TGraph( nSteps );
        TString magLabel( "PoleMag_" );
        magLabel += iPole;
        magGraph->SetName( magLabel.Data() );
        magGraph->SetLineColor( kBlack );
        magGraph->SetMarkerColor( kBlack );
        magGraph->SetLineWidth( 1 );
        ampGraphs->Add( magGraph, "lp" );

        if ( ! initLegend ) {

            ampLegend->AddEntry( realGraph, "Real part", "l" );
            ampLegend->AddEntry( imagGraph, "Imag part", "l" );
            ampLegend->AddEntry( magGraph, "Magnitude", "l" );

            initLegend = kTRUE;
        }
    }

    for ( Int_t iM = 0; iM < nSteps; iM++ ) {

        Double_t m = dm * iM + mMin;
        Double_t s = m * m;

        // Update kinematics. The K-matrix only really needs s = m23Sq, i.e. m13Sq is "ignored".
        // Here, take the diagonal where m12Sq = 0.5*(mB^2 + 3mpi^2), the DP centre and use
        // m13Sq = mBSq + 3mPiSq - m12Sq - m23Sq
        Double_t m13Sq = mSqTerm - s;

        kinematics->updateKinematics( m13Sq, s );

        for ( iPole = 0; iPole < nPoles; iPole++ ) {

            // Get the amplitude for the ith pole
            LauComplex amp  = poleVect[iPole]->amplitude( kinematics );
            Double_t reAmp  = amp.re();
            Double_t imAmp  = amp.im();
            Double_t magAmp = sqrt( reAmp * reAmp + imAmp * imAmp );

            // Argand diagrams: full plot and only selected points
            TMultiGraph* argGraphs = argGraphVect[iPole];
            TList* argList         = argGraphs->GetListOfGraphs();

            TGraph* argGraph = dynamic_cast<TGraph*>( argList->At( 0 ) );
            argGraph->SetPoint( iM, reAmp, imAmp );

            // Graphs showing specific points
            for ( Int_t iP = 0; iP < NMain; iP++ ) {

                if ( points[iP] == iM ) {
                    TGraph* pointGraph = dynamic_cast<TGraph*>( argList->At( iP + 1 ) );
                    pointGraph->SetPoint( 0, reAmp, imAmp );
                }
            }

            // Amplitude multigraphs
            TMultiGraph* ampGraphs = ampGraphVect[iPole];
            TList* graphList       = ampGraphs->GetListOfGraphs();

            TGraph* realGraph = dynamic_cast<TGraph*>( graphList->At( 0 ) );
            realGraph->SetPoint( iM, m, reAmp );

            TGraph* imagGraph = dynamic_cast<TGraph*>( graphList->At( 1 ) );
            imagGraph->SetPoint( iM, m, imAmp );

            TGraph* magGraph = dynamic_cast<TGraph*>( graphList->At( 2 ) );
            magGraph->SetPoint( iM, m, magAmp );
        }
    }

    // Create plots
    TCanvas* theCanvas = new TCanvas( "theCanvas", "", 1000, 1000 );
    gROOT->SetStyle( "Plain" );
    gStyle->SetOptStat( 0 );

    theCanvas->Clear();
    theCanvas->UseCurrentStyle();
    theCanvas->Divide( 2, 3 );

    // Argand plots
    for ( iPole = 0; iPole < nPoles; iPole++ ) {

        theCanvas->cd( iPole + 1 );
        TMultiGraph* argGraphs = argGraphVect[iPole];
        argGraphs->Draw( "a" );

        TString title( "Production Pole " );
        title += iPole + 1;
        argGraphs->SetTitle( title.Data() );

        TAxis* xAxis = argGraphs->GetXaxis();
        xAxis->SetTitle( "Real" );
        xAxis->SetLabelSize( 0.05 );
        xAxis->SetTitleSize( 0.05 );
        xAxis->CenterTitle( kTRUE );
        xAxis->SetTitleOffset( 0.9 );

        TAxis* yAxis = argGraphs->GetYaxis();
        yAxis->SetTitle( "Imag" );
        yAxis->SetLabelSize( 0.05 );
        yAxis->SetTitleSize( 0.05 );
        yAxis->CenterTitle( kTRUE );

        theCanvas->Update();
    }

    theCanvas->cd( 6 );
    argLegend->Draw();

    theCanvas->Print( "ProdPoleArgand.png" );
    if ( doEPS ) {
        theCanvas->Print( "ProdPoleArgand.eps" );
    }

    // Amplitude plots
    theCanvas->Clear();
    theCanvas->UseCurrentStyle();
    theCanvas->Divide( 2, 3 );

    for ( iPole = 0; iPole < nPoles; iPole++ ) {

        theCanvas->cd( iPole + 1 );

        TMultiGraph* ampGraphs = ampGraphVect[iPole];
        ampGraphs->Draw( "a" );

        TAxis* xAxis = ampGraphs->GetXaxis();
        xAxis->SetTitle( "m(#pi^{+}#pi^{-}) (GeV/c^{2})" );
        xAxis->SetLabelSize( 0.05 );
        xAxis->SetTitleSize( 0.05 );
        xAxis->CenterTitle( kTRUE );
        xAxis->SetTitleOffset( 0.9 );

        TString yLabel( "Production Pole " );
        yLabel += iPole + 1;

        TAxis* yAxis = ampGraphs->GetYaxis();
        yAxis->SetTitleOffset( 0.9 );
        yAxis->SetTitle( yLabel.Data() );
        yAxis->SetLabelSize( 0.05 );
        yAxis->SetTitleSize( 0.05 );
        yAxis->CenterTitle( kTRUE );

        theCanvas->Update();
    }

    // Legend
    theCanvas->cd( 6 );
    ampLegend->Draw();

    theCanvas->Print( "ProdPoleAmp.png" );
    if ( doEPS ) {
        theCanvas->Print( "ProdPoleAmp.eps" );
    }

    // Cleanup
    std::vector<TMultiGraph*>::iterator argIter;
    for ( argIter = argGraphVect.begin(); argIter != argGraphVect.end(); ++argIter ) {

        TMultiGraph* theGraph = *argIter;
        delete theGraph;
    }

    std::vector<TMultiGraph*>::iterator ampIter;
    for ( ampIter = ampGraphVect.begin(); ampIter != ampGraphVect.end(); ++ampIter ) {

        TMultiGraph* theGraph = *ampIter;
        delete theGraph;
    }

    std::vector<LauKMatrixProdPole*>::iterator poleIter;
    for ( poleIter = poleVect.begin(); poleIter != poleVect.end(); ++poleIter ) {

        LauKMatrixProdPole* thePole = *poleIter;
        delete thePole;
    }

    delete theDaughters;
    delete propFull;

    delete ampLegend;
    delete theCanvas;
}

void doProdSVPPlots( Bool_t useProdAdler )
{

    // Create plots of the SVP production amplitudes using the first row
    // of the propagator (for pipi). Each SVP is a pipi -> channel mode

    Int_t resPairInt( 1 );
    Int_t nChannels( 5 );
    Int_t KMatrixIndex( 1 );    // pi pi S-wave
    Int_t nSVPs( 5 );

    // Mass range for plots
    Double_t mMin( 0.28 );
    Double_t mMax( 2.0 );
    Int_t nSteps = 2000;
    Double_t dm  = ( mMax - mMin ) / ( nSteps * 1.0 );

    // Find the point numbers for m = mMin, 0.5, 0.75, 1, 1.25, 1.5, 1.75 and 2 GeV
    std::vector<Int_t> points;
    std::vector<Int_t> colours, styles;
    std::vector<Double_t> massPoints;
    massPoints.push_back( mMin );
    massPoints.push_back( 0.5 );
    massPoints.push_back( 0.75 );
    massPoints.push_back( 1.0 );
    massPoints.push_back( 1.25 );
    massPoints.push_back( 1.5 );
    massPoints.push_back( 1.75 );
    massPoints.push_back( 1.99999 );

    Int_t NMain = static_cast<Int_t>( massPoints.size() );
    for ( Int_t iP = 0; iP < NMain; iP++ ) {

        Int_t pointInt = static_cast<Int_t>( ( massPoints[iP] - mMin ) / dm );
        points.push_back( pointInt );

        // Red or Blue
        Int_t iColour = kRed;
        if ( iP > 3 ) {
            iColour = kBlue;
        }
        Int_t iStyle = kFullCircle;
        Int_t iPMod  = iP % 4;
        if ( iPMod == 1 ) {
            iStyle = kOpenCircle;
        } else if ( iPMod == 2 ) {
            iStyle = kFullSquare;
        } else if ( iPMod == 3 ) {
            iStyle = kOpenSquare;
        }

        colours.push_back( iColour );
        styles.push_back( iStyle );
    }

    Double_t mPi     = 0.13957018;
    Double_t mPiSq   = mPi * mPi;
    Double_t mB      = 5.279;
    Double_t mBSq    = mB * mB;
    Double_t mSqTerm = 2.0 * ( ( mBSq / 3.0 ) + mPiSq );

    // Full K-matrix parameterisation
    TString KFile( "B3piKMatrixCoeff.dat" );
    LauKMatrixPropagator* propFull =
        new LauKMatrixPropagator( "propFull", KFile, resPairInt, nChannels, nSVPs, KMatrixIndex );

    // Define kinematics pointer
    Bool_t squareDP            = kTRUE;
    LauDaughters* theDaughters = new LauDaughters( "B+", "pi+", "pi+", "pi-", squareDP );

    // Kinematics
    LauKinematics* kinematics = theDaughters->getKinematics();

    std::vector<LauKMatrixProdSVP*> SVPVect;

    Int_t iSVP( 0 );
    for ( iSVP = 0; iSVP < nSVPs; iSVP++ ) {

        TString SVPName( "SVP" );
        SVPName += iSVP;

        LauKMatrixProdSVP* theSVP =
            new LauKMatrixProdSVP( SVPName, iSVP + 1, resPairInt, propFull, theDaughters, useProdAdler );

        SVPVect.push_back( theSVP );
    }

    // Argand diagram graphs (full plot, then others showing only selected points)
    std::vector<TMultiGraph*> argGraphVect;
    // Real, imaginary and magnitude amplitude multigraphs
    std::vector<TMultiGraph*> ampGraphVect;

    // TLegend for argand graphs
    TLegend* argLegend = new TLegend( 0.1, 0.1, 0.7, 0.995 );
    argLegend->SetHeader( "m(#pi^{+}#pi^{-}) points (GeV/c^{2})" );
    argLegend->SetTextSize( 0.075 );
    Bool_t initArgLegend( kFALSE );

    // Legend for the graphs (same for all multigraphs)
    TLegend* ampLegend = new TLegend( 0.2, 0.2, 0.8, 0.7 );
    ampLegend->SetHeader( "Production SVP components" );
    Bool_t initLegend( kFALSE );

    for ( iSVP = 0; iSVP < nSVPs; iSVP++ ) {

        // Argand plot showing phase motion
        TMultiGraph* argGraphs = new TMultiGraph();
        TString argMultiLabel( "SVPMultiArg_" );
        argMultiLabel += iSVP;
        argGraphs->SetName( argMultiLabel.Data() );
        argGraphVect.push_back( argGraphs );

        TGraph* argGraph = new TGraph( nSteps );
        TString argLabel( "SVPArg_" );
        argLabel += iSVP;
        argGraph->SetName( argLabel.Data() );
        argGraph->SetLineWidth( 1 );
        argGraphs->Add( argGraph, "lp" );

        // Argand graphs showing only one highlighted point
        for ( Int_t iG = 0; iG < NMain; iG++ ) {

            TGraph* pointGraph = new TGraph( 1 );
            TString pointLabel( "SVPPoint_" );
            pointLabel += iG;

            pointGraph->SetName( pointLabel.Data() );
            pointGraph->SetMarkerColor( colours[iG] );
            pointGraph->SetMarkerSize( 1.1 );
            pointGraph->SetMarkerStyle( styles[iG] );
            argGraphs->Add( pointGraph, "p" );
        }

        if ( ! initArgLegend ) {
            TList* argList = argGraphs->GetListOfGraphs();
            for ( Int_t iG = 0; iG < NMain; iG++ ) {
                TGraph* pointGraph = dynamic_cast<TGraph*>( argList->At( iG + 1 ) );
                argLegend->AddEntry( pointGraph, Form( "%.2f", massPoints[iG] ), "p" );
                initArgLegend = kTRUE;
            }
        }

        // Real, imaginary and magnitude amplitude
        TMultiGraph* ampGraphs = new TMultiGraph();
        TString multiLabel( "SVPMulti_" );
        multiLabel += iSVP;
        ampGraphs->SetName( multiLabel.Data() );
        ampGraphVect.push_back( ampGraphs );

        TGraph* realGraph = new TGraph( nSteps );
        TString realLabel( "SVPReal_" );
        realLabel += iSVP;
        realGraph->SetName( realLabel.Data() );
        realGraph->SetLineColor( kRed );
        realGraph->SetMarkerColor( kRed );
        realGraph->SetLineWidth( 1 );
        ampGraphs->Add( realGraph, "lp" );

        TGraph* imagGraph = new TGraph( nSteps );
        TString imagLabel( "SVPImag_" );
        imagLabel += iSVP;
        imagGraph->SetName( imagLabel.Data() );
        imagGraph->SetLineColor( kBlue );
        imagGraph->SetMarkerColor( kBlue );
        imagGraph->SetLineWidth( 1 );
        ampGraphs->Add( imagGraph, "lp" );

        TGraph* magGraph = new TGraph( nSteps );
        TString magLabel( "SVPMag_" );
        magLabel += iSVP;
        magGraph->SetName( magLabel.Data() );
        magGraph->SetLineColor( kBlack );
        magGraph->SetMarkerColor( kBlack );
        magGraph->SetLineWidth( 1 );
        ampGraphs->Add( magGraph, "lp" );

        if ( ! initLegend ) {

            ampLegend->AddEntry( realGraph, "Real part", "l" );
            ampLegend->AddEntry( imagGraph, "Imag part", "l" );
            ampLegend->AddEntry( magGraph, "Magnitude", "l" );

            initLegend = kTRUE;
        }
    }

    for ( Int_t iM = 0; iM < nSteps; iM++ ) {

        Double_t m = dm * iM + mMin;
        Double_t s = m * m;

        // Update kinematics. The K-matrix only really needs s = m23Sq, i.e. m13Sq is "ignored".
        // Here, take the diagonal where m12Sq = 0.5*(mB^2 + 3mpi^2), the DP centre and use
        // m13Sq = mBSq + 3mPiSq - m12Sq - m23Sq
        Double_t m13Sq = mSqTerm - s;

        kinematics->updateKinematics( m13Sq, s );

        for ( iSVP = 0; iSVP < nSVPs; iSVP++ ) {

            // Get the amplitude for the ith SVP
            LauComplex amp  = SVPVect[iSVP]->amplitude( kinematics );
            Double_t reAmp  = amp.re();
            Double_t imAmp  = amp.im();
            Double_t magAmp = sqrt( reAmp * reAmp + imAmp * imAmp );

            // Argand diagrams: full plot and only selected points
            TMultiGraph* argGraphs = argGraphVect[iSVP];
            TList* argList         = argGraphs->GetListOfGraphs();

            TGraph* argGraph = dynamic_cast<TGraph*>( argList->At( 0 ) );
            argGraph->SetPoint( iM, reAmp, imAmp );

            // Graphs showing specific points
            for ( Int_t iP = 0; iP < NMain; iP++ ) {

                if ( points[iP] == iM ) {
                    TGraph* pointGraph = dynamic_cast<TGraph*>( argList->At( iP + 1 ) );
                    pointGraph->SetPoint( 0, reAmp, imAmp );
                }
            }

            // Amplitude multigraphs
            TMultiGraph* ampGraphs = ampGraphVect[iSVP];
            TList* graphList       = ampGraphs->GetListOfGraphs();

            TGraph* realGraph = dynamic_cast<TGraph*>( graphList->At( 0 ) );
            realGraph->SetPoint( iM, m, reAmp );

            TGraph* imagGraph = dynamic_cast<TGraph*>( graphList->At( 1 ) );
            imagGraph->SetPoint( iM, m, imAmp );

            TGraph* magGraph = dynamic_cast<TGraph*>( graphList->At( 2 ) );
            magGraph->SetPoint( iM, m, magAmp );
        }
    }

    // Create plots
    TCanvas* theCanvas = new TCanvas( "theCanvas", "", 1000, 1000 );
    gROOT->SetStyle( "Plain" );
    gStyle->SetOptStat( 0 );

    theCanvas->Clear();
    theCanvas->UseCurrentStyle();
    theCanvas->Divide( 2, 3 );

    // Argand plots
    for ( iSVP = 0; iSVP < nSVPs; iSVP++ ) {

        theCanvas->cd( iSVP + 1 );
        TMultiGraph* argGraphs = argGraphVect[iSVP];
        argGraphs->Draw( "a" );

        TString title( "Production SVP " );
        title += iSVP + 1;
        argGraphs->SetTitle( title.Data() );

        TAxis* xAxis = argGraphs->GetXaxis();
        xAxis->SetTitle( "Real" );
        xAxis->SetLabelSize( 0.05 );
        xAxis->SetTitleSize( 0.05 );
        xAxis->CenterTitle( kTRUE );
        xAxis->SetTitleOffset( 0.9 );

        TAxis* yAxis = argGraphs->GetYaxis();
        yAxis->SetTitle( "Imag" );
        yAxis->SetLabelSize( 0.05 );
        yAxis->SetTitleSize( 0.05 );
        yAxis->CenterTitle( kTRUE );

        theCanvas->Update();
    }

    theCanvas->cd( 6 );
    argLegend->Draw();

    theCanvas->Print( "ProdSVPArgand.png" );
    if ( doEPS ) {
        theCanvas->Print( "ProdSVPArgand.eps" );
    }

    // Amplitude plots
    theCanvas->Clear();
    theCanvas->UseCurrentStyle();
    theCanvas->Divide( 2, 3 );

    for ( iSVP = 0; iSVP < nSVPs; iSVP++ ) {

        theCanvas->cd( iSVP + 1 );

        TMultiGraph* ampGraphs = ampGraphVect[iSVP];
        ampGraphs->Draw( "a" );

        TAxis* xAxis = ampGraphs->GetXaxis();
        xAxis->SetTitle( "m(#pi^{+}#pi^{-}) (GeV/c^{2})" );
        xAxis->SetLabelSize( 0.05 );
        xAxis->SetTitleSize( 0.05 );
        xAxis->CenterTitle( kTRUE );
        xAxis->SetTitleOffset( 0.9 );

        TAxis* yAxis = ampGraphs->GetYaxis();

        TString yLabel( "Production SVP " );
        yLabel += iSVP + 1;

        yAxis->SetTitleOffset( 0.9 );
        yAxis->SetTitle( yLabel.Data() );
        yAxis->SetLabelSize( 0.05 );
        yAxis->SetTitleSize( 0.05 );
        yAxis->CenterTitle( kTRUE );

        theCanvas->Update();
    }

    // Legend
    theCanvas->cd( 6 );
    ampLegend->Draw();

    theCanvas->Print( "ProdSVPAmp.png" );
    if ( doEPS ) {
        theCanvas->Print( "ProdSVPAmp.eps" );
    }

    // Cleanup
    std::vector<TMultiGraph*>::iterator argIter;
    for ( argIter = argGraphVect.begin(); argIter != argGraphVect.end(); ++argIter ) {

        TMultiGraph* theGraph = *argIter;
        delete theGraph;
    }

    std::vector<TMultiGraph*>::iterator ampIter;
    for ( ampIter = ampGraphVect.begin(); ampIter != ampGraphVect.end(); ++ampIter ) {

        TMultiGraph* theGraph = *ampIter;
        delete theGraph;
    }

    std::vector<LauKMatrixProdSVP*>::iterator SVPIter;
    for ( SVPIter = SVPVect.begin(); SVPIter != SVPVect.end(); ++SVPIter ) {

        LauKMatrixProdSVP* theSVP = *SVPIter;
        delete theSVP;
    }

    delete theDaughters;
    delete propFull;

    delete ampLegend;
    delete theCanvas;
}

int main( const int argc, const char** argv )
{

    Int_t prodAdlerInt( 0 ), epsInt( 0 );

    if ( argc > 1 ) {
        prodAdlerInt = std::atoi( argv[1] );
        if ( prodAdlerInt != 0 ) {
            doProdAdler = kTRUE;
        }
    }

    if ( argc > 2 ) {
        epsInt = std::atoi( argv[2] );
        if ( epsInt != 0 ) {
            doEPS = kTRUE;
        }
    }

    // Scattering K-matrix components
    doKScattPlots();

    // First row of the propagator (pi pi channel)
    doPropagatorPlots();

    // SVP and SVP production amplitudes (pi pi channel)
    doProdPolePlots( doProdAdler );
    doProdSVPPlots( doProdAdler );

    return 0;
}
