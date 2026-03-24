
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

// Plot the mass projections of K-matrix amplitude terms generated over the Dalitz plot so that
// we can see what shapes the poles and SVP components have. Also show the rho for comparison

#include <iostream>
using std::cout;
using std::endl;

#include "LauCartesianCPCoeffSet.hh"
#include "LauDaughters.hh"
#include "LauEffModel.hh"
#include "LauIsobarDynamics.hh"
#include "LauKMatrixPropagator.hh"
#include "LauParameter.hh"
#include "LauResonanceMaker.hh"
#include "LauSimpleFitModel.hh"
#include "LauVetoes.hh"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include <map>
#include <vector>

void createHistos( const std::string& rootFileName );
void plotHistos( const std::string& rootFileName );
std::vector<TH1*> genMassHistos( Int_t index, Int_t NSignal );
Bool_t useProdAdler = kFALSE;

int main( const int argc, const char** argv )
{

    int histFlag( 0 ), prodFlag( 0 );
    if ( argc > 1 ) {
        histFlag = std::atoi( argv[1] );
    }

    if ( argc > 2 ) {
        prodFlag = std::atoi( argv[2] );
        if ( prodFlag != 0 ) {
            useProdAdler = kTRUE;
        }
    }

    std::string rootFileName( "KMatrixMassHistos.root" );
    if ( histFlag ) {
        createHistos( rootFileName );
    }

    plotHistos( rootFileName );
}

void plotHistos( const std::string& rootFileName )
{

    cout << "Running plotHistos for " << rootFileName << endl;

    TCanvas* theCanvas = new TCanvas( "theCanvas", "", 1600, 1200 );
    gROOT->SetStyle( "Plain" );
    gStyle->SetOptStat( 0 );

    theCanvas->Clear();
    theCanvas->UseCurrentStyle();

    TFile* rootFile = TFile::Open( rootFileName.c_str(), "read" );
    rootFile->cd();

    // Plot rho with K-matrix poles
    TLegend* poleLegend = new TLegend( 0.7, 0.5, 0.90, 0.90, "" );
    poleLegend->SetFillColor( 0 );
    poleLegend->SetTextSize( 0.03 );

    THStack* poleStack = new THStack();
    TH1* rhoHist       = dynamic_cast<TH1*>( rootFile->Get( "Rho_m13" ) );
    rhoHist->SetLineColor( kBlack );
    rhoHist->SetLineWidth( 2 );
    rhoHist->Scale( 1.0 / rhoHist->Integral() );

    poleStack->Add( rhoHist );
    poleLegend->AddEntry( rhoHist, "#rho(770)", "l" );

    Int_t colours[5] = { kRed, kViolet, kBlue, kGreen + 2, kOrange + 2 };

    Int_t i( 0 ), iN( 5 );
    for ( i = 0; i < iN; i++ ) {

        Int_t j = i + 1;
        TString poleName( "KMPole" );
        poleName += j;
        poleName += "_m13";

        TH1* poleHist = dynamic_cast<TH1*>( rootFile->Get( poleName.Data() ) );
        poleHist->SetLineColor( colours[i] );
        poleHist->SetLineWidth( 2 );

        // Normalise histogram to unit area
        poleHist->Scale( 1.0 / poleHist->Integral() );
        poleStack->Add( poleHist );

        TString poleLabel( "Pole" );
        poleLabel += j;
        poleLegend->AddEntry( poleHist, poleLabel, "l" );
    }

    poleStack->Draw( "nostackc" );

    TAxis* poleXAxis = poleStack->GetXaxis();
    poleXAxis->SetTitle( "m(#pi^{+}#pi^{-}) (GeV/c^{2})" );
    poleXAxis->CenterTitle( kTRUE );
    TAxis* poleYAxis = poleStack->GetYaxis();
    poleYAxis->SetTitle( "Normalised intensity" );
    poleYAxis->CenterTitle( kTRUE );
    poleYAxis->SetTitleOffset( 1.25 );

    poleLegend->Draw();
    theCanvas->Update();

    theCanvas->Print( "B3piKMatrixProjPoles.png" );
    theCanvas->Print( "B3piKMatrixProjPoles.eps" );

    // Plot NR with K-matrix SVPs
    TLegend* SVPLegend = new TLegend( 0.7, 0.5, 0.90, 0.90, "" );
    SVPLegend->SetFillColor( 0 );
    SVPLegend->SetTextSize( 0.03 );

    THStack* SVPStack = new THStack();
    SVPStack->Add( rhoHist );
    SVPLegend->AddEntry( rhoHist, "rho(770)", "l" );

    for ( i = 0; i < iN; i++ ) {

        Int_t j = i + 1;
        TString SVPName( "KMSVP" );
        SVPName += j;
        SVPName += "_m13";

        TH1* SVPHist = dynamic_cast<TH1*>( rootFile->Get( SVPName.Data() ) );
        SVPHist->SetLineColor( colours[i] );
        SVPHist->SetLineWidth( 2 );

        SVPHist->Scale( 1.0 / SVPHist->Integral() );
        SVPStack->Add( SVPHist );

        TString SVPLabel( "SVP" );
        SVPLabel += j;
        SVPLegend->AddEntry( SVPHist, SVPLabel, "l" );
    }

    SVPStack->Draw( "nostackc" );

    TAxis* SVPXAxis = SVPStack->GetXaxis();
    SVPXAxis->SetTitle( "m(#pi^{+}#pi^{-}) (GeV/c^{2})" );
    SVPXAxis->CenterTitle( kTRUE );
    TAxis* SVPYAxis = SVPStack->GetYaxis();
    SVPYAxis->SetTitle( "Normalised intensity" );
    SVPYAxis->CenterTitle( kTRUE );
    SVPYAxis->SetTitleOffset( 1.25 );

    SVPLegend->Draw();
    theCanvas->Update();

    theCanvas->Print( "B3piKMatrixProjSVPs.png" );
    theCanvas->Print( "B3piKMatrixProjSVPs.eps" );

    delete SVPStack;
    delete SVPLegend;
    delete poleStack;
    delete poleLegend;

    rootFile->Close();
}

void createHistos( const std::string& rootFileName )
{

    cout << "rootFileName = " << rootFileName << endl;

    Int_t i( 0 ), iN( 12 );
    Int_t NSignal = 100000;

    std::map<int, std::vector<TH1*>> histMap;

    // Create mass/DP histograms for the different amplitude terms
    for ( i = 0; i < iN; i++ ) {

        std::vector<TH1*> histVect = genMassHistos( i, NSignal );
        histMap[i]                 = histVect;
    }

    TFile* rootFile = TFile::Open( rootFileName.c_str(), "recreate" );
    rootFile->cd();

    for ( i = 0; i < iN; i++ ) {

        cout << "Writing histograms for i = " << i << endl;

        std::vector<TH1*> histVect = histMap[i];
        std::vector<TH1*>::iterator iter;
        for ( iter = histVect.begin(); iter != histVect.end(); ++iter ) {

            TH1* theHist = *iter;
            theHist->Write();
        }
    }

    rootFile->Close();
}

std::vector<TH1*> genMassHistos( Int_t index, Int_t NSignal )
{

    // Vector to store output mass/DP histograms
    std::vector<TH1*> histVect;

    // Define Dalitz plot: kinematics and resonances
    Bool_t squareDP            = kTRUE;
    LauDaughters* theDaughters = new LauDaughters( "B+", "pi+", "pi+", "pi-", squareDP );

    // Apply some vetoes to the DP
    LauVetoes* vetoes = new LauVetoes();

    // Define the efficiency model
    LauEffModel* theEffModel = new LauEffModel( theDaughters, vetoes );

    // Create the isobar model
    LauResonanceMaker& resMaker = LauResonanceMaker::get();
    resMaker.setDefaultBWRadius( LauBlattWeisskopfFactor::Parent, 4.0 );
    resMaker.setDefaultBWRadius( LauBlattWeisskopfFactor::Light, 4.0 );
    resMaker.fixBWRadius( LauBlattWeisskopfFactor::Parent, kTRUE );
    resMaker.fixBWRadius( LauBlattWeisskopfFactor::Light, kTRUE );

    // Define the "S-wave" K-matrix propagator
    Int_t nChannels    = 5;
    Int_t nPoles       = 5;
    Int_t resPairInt   = 1;
    Int_t KMatrixIndex = 1;    // for S-wave

    // Add one pole or SVP amplitude for the K-matrix and plot its mass projections
    // so that we can see the amplitude/intensity shape of the different terms
    Double_t aSqMaxValue = 0.3;

    Double_t nSigEvents  = NSignal * 1.0;
    Bool_t fixNSigEvents = kTRUE;
    Int_t nExpt          = 1;
    Int_t firstExpt      = 0;
    Bool_t isToy         = kTRUE;

    LauIsobarDynamics* theSigModel = new LauIsobarDynamics( theDaughters, theEffModel );
    theSigModel->setASqMaxValue( aSqMaxValue );

    TString theName( "KMPole" );
    if ( index >= 5 ) {
        theName = "KMSVP";
    }
    TString parName( "" );

    if ( index == 11 ) {

        // Simple Non-resonant
        cout << "Non-resonant" << endl;
        theName = "NR";
        parName = "NonReson";
        theSigModel->addResonance( "NonReson", 0, LauAbsResonance::FlatNR );

    } else if ( index == 10 ) {

        // Rho resonance
        cout << "Using rho(770)" << endl;

        theName = "Rho";
        parName = "rho0(770)";
        theSigModel->addResonance( "rho0(770)", 1, LauAbsResonance::GS );

    } else {

        // K-matrix amplitude
        theSigModel->defineKMatrixPropagator( "KMSWave",
                                              "B3piKMatrixCoeff.dat",
                                              resPairInt,
                                              nChannels,
                                              nPoles,
                                              KMatrixIndex );

        cout << "Implementing K matrix term number " << index << ", name = " << theName << endl;

        Int_t j  = index % 5 + 1;
        theName += j;
        parName  = theName;

        if ( index < 5 ) {
            theSigModel->addKMatrixProdPole( theName.Data(), "KMSWave", j, useProdAdler );
        } else {
            theSigModel->addKMatrixProdSVP( theName.Data(), "KMSWave", j, useProdAdler );
        }
    }

    // Create the fit model
    LauSimpleFitModel* fitModel = new LauSimpleFitModel( theSigModel );

    std::vector<LauCartesianCPCoeffSet*> coeffSet;
    Bool_t fixPar        = kTRUE;
    Bool_t fixCP         = kTRUE;
    Bool_t doTwoStageFit = kFALSE;

    coeffSet.push_back( new LauCartesianCPCoeffSet( parName.Data(),
                                                    1.0,
                                                    0.0,
                                                    0.0,
                                                    0.0,
                                                    fixPar,
                                                    fixPar,
                                                    fixCP,
                                                    fixCP,
                                                    doTwoStageFit,
                                                    doTwoStageFit ) );

    std::vector<LauCartesianCPCoeffSet*>::iterator coeffIter;
    for ( coeffIter = coeffSet.begin(); coeffIter != coeffSet.end(); ++coeffIter ) {
        fitModel->setAmpCoeffSet( *coeffIter );
    }

    // Set the number of signal events and the number of experiments
    LauParameter* signalEvents = new LauParameter( "signalEvents",
                                                   nSigEvents,
                                                   -2.0 * nSigEvents,
                                                   2.0 * nSigEvents,
                                                   fixNSigEvents );
    fitModel->setNSigEvents( signalEvents );
    fitModel->setNExpts( nExpt, firstExpt, isToy );

    // Execute the generation/fit
    TString treeName( "" );
    TString dummyName( "" );
    TString tableFileName( "" );

    TString dataFileName( "Gen" );
    dataFileName += theName;
    dataFileName += ".root";

    fitModel->run( "gen", dataFileName, treeName, dummyName, tableFileName );

    // Create DP and mass projection plots so that we can see the K matrix
    // amplitude term shapes

    // Histogram of the invariant mass for pi+ pi- pairs: m13 & m23
    TString m13Name( theName );
    m13Name += "_m13";

    TH1D* m13Hist = new TH1D( m13Name.Data(), "", 100, 0.0, 5.5 );
    m13Hist->SetXTitle( "m(#pi^{+}#pi^{-})" );
    m13Hist->SetDirectory( 0 );
    histVect.push_back( m13Hist );

    // Histogram of the wrong sign invariant mass pi+ pi+: m12
    TString m12Name( theName );
    m12Name += "_m12";

    TH1D* m12Hist = new TH1D( m12Name.Data(), "", 100, 0.0, 5.5 );
    m12Hist->SetXTitle( "m(#pi^{+}#pi^{+})" );
    m12Hist->SetDirectory( 0 );
    histVect.push_back( m12Hist );

    // Dalitz plot 2d histogram
    TString DPName( theName );
    DPName += "_DP";

    TH2D* DPHist = new TH2D( DPName.Data(), "", 100, 0.0, 14.0, 100, 0.0, 28.0 );
    //TH2D* DPHist = new TH2D("DPHist", "", 100, 0.0, 28.0, 100, 0.0, 28.0);
    DPHist->SetXTitle( "m^{2}(#pi^{+}#pi^{-})_{min}" );
    DPHist->SetYTitle( "m^{2}(#pi^{+}#pi^{-})_{max}" );
    DPHist->SetDirectory( 0 );
    histVect.push_back( DPHist );

    // Simulation results file
    TFile* genFile = TFile::Open( dataFileName.Data(), "read" );
    genFile->cd();
    TTree* genTree = dynamic_cast<TTree*>( genFile->Get( "genResults" ) );
    genTree->SetBranchStatus( "*", 0 );
    genTree->SetBranchStatus( "m12", 1 );
    genTree->SetBranchStatus( "m23", 1 );
    genTree->SetBranchStatus( "m13", 1 );

    Double_t m12, m23, m13;
    genTree->SetBranchAddress( "m12", &m12 );
    genTree->SetBranchAddress( "m23", &m23 );
    genTree->SetBranchAddress( "m13", &m13 );

    Int_t nEntries = genTree->GetEntries();
    cout << "nEntries = " << nEntries << endl;
    Int_t k;

    // Fill mass and DP histograms
    for ( k = 0; k < nEntries; k++ ) {

        genTree->GetEntry( k );

        Double_t mMin = m13;
        Double_t mMax = m23;
        if ( m23 < m13 ) {
            mMin = m23;
            mMax = m13;
        }

        m13Hist->Fill( mMin );
        m13Hist->Fill( mMax );

        m12Hist->Fill( m12 );

        DPHist->Fill( mMin * mMin, mMax * mMax );
    }

    // Cleanup
    genFile->Close();

    delete signalEvents;

    for ( coeffIter = coeffSet.begin(); coeffIter != coeffSet.end(); ++coeffIter ) {
        delete *coeffIter;
    }

    delete fitModel;
    delete theSigModel;
    delete theEffModel;
    delete vetoes;
    delete theDaughters;

    // Return histograms
    return histVect;
}
