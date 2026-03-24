
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

#include "LauBkgndDPModel.hh"
#include "LauDaughters.hh"
#include "LauEffModel.hh"
#include "LauIsobarDynamics.hh"
#include "LauMagPhaseCoeffSet.hh"
#include "LauRandom.hh"
#include "LauSimpleFitModel.hh"
#include "LauVetoes.hh"

#include "TFile.h"
#include "TH2.h"
#include "TString.h"
#include "TTree.h"

#include <cstdlib>
#include <iostream>
#include <vector>

void usage( std::ostream& out, const TString& progName )
{
    out << "Usage:\n";
    out << progName << " gen <category = DD or LL> [nExpt = 1] [firstExpt = 0]\n";
    out << "or\n";
    out << progName
        << " fit <category = DD or LL> <iFit> <port> [hostname] [nExpt = 1] [firstExpt = 0]"
        << std::endl;
}

int main( int argc, char** argv )
{
    // Process command-line arguments
    // Usage:
    // ./SimFitTask gen <category = DD or LL> [nExpt = 1] [firstExpt = 0]
    // or
    // ./SimFitTask fit <category = DD or LL> <iFit> <port> [nExpt = 1] [firstExpt = 0]
    if ( argc < 3 ) {
        usage( std::cerr, argv[0] );
        return EXIT_FAILURE;
    }

    TString command = argv[1];
    command.ToLower();

    TString category = argv[2];
    if ( category != "DD" && category != "LL" ) {
        usage( std::cerr, argv[0] );
        return EXIT_FAILURE;
    }

    Int_t iFit( 0 );
    UInt_t port( 5000 );
    TString hostname( "localhost" );
    Int_t nExpt( 1 );
    Int_t firstExpt( 0 );
    Bool_t isToy( kTRUE );
    if ( command == "gen" ) {
        if ( argc > 3 ) {
            nExpt = atoi( argv[3] );
            if ( argc > 4 ) {
                firstExpt = atoi( argv[4] );
            }
        }
    } else if ( command == "fit" ) {
        if ( argc < 5 ) {
            usage( std::cerr, argv[0] );
            return EXIT_FAILURE;
        }
        iFit = atoi( argv[3] );
        port = atoi( argv[4] );
        if ( argc > 5 ) {
            hostname = argv[5];
            if ( argc > 6 ) {
                nExpt = atoi( argv[6] );
                if ( argc > 7 ) {
                    firstExpt = atoi( argv[7] );
                }
            }
        }
    } else {
        usage( std::cerr, argv[0] );
        return EXIT_FAILURE;
    }

    // If you want to use square DP histograms for efficiency,
    // backgrounds or you just want the square DP co-ordinates
    // stored in the toy MC ntuple then set this to kTRUE
    Bool_t squareDP = kFALSE;

    // This defines the DP => decay is B0 -> pi0 pi0 K_S0
    // Particle 1 = pi0
    // Particle 2 = pi0
    // Particle 3 = KS_0
    // The DP is defined in terms of m13Sq and m23Sq
    LauDaughters* daughters = new LauDaughters( "B0", "pi0", "pi0", "K_S0", squareDP );

    // Optionally apply some vetoes to the DP
    LauVetoes* vetoes = new LauVetoes();

    // Define the efficiency model (defaults to unity everywhere)
    // Can optionally provide a histogram to model variation over DP
    // (example syntax given in commented-out section)
    LauEffModel* effModel = new LauEffModel( daughters, vetoes );
    //TFile *effHistFile = TFile::Open("histoFiles/efficiency.root", "read");
    //TH2* effHist = dynamic_cast<TH2*>(effHistFile->Get("effHist"));
    //Bool_t useInterpolation = kTRUE;
    //Bool_t fluctuateBins = kFALSE;
    //Bool_t useUpperHalf = kTRUE;
    //effModel->setEffHisto(effHist, useInterpolation, fluctuateBins, 0.0, 0.0, useUpperHalf, squareDP);

    // Create the isobar model
    LauIsobarDynamics* sigModel = new LauIsobarDynamics( daughters, effModel );
    LauAbsResonance* res( 0 );
    // resPairAmpInt = 1 => resonance mass is m23.
    // resPairAmpInt = 3 => resonance mass is m12.
    res = sigModel->addResonance( "f_0(980)", 3, LauAbsResonance::Flatte );
    res = sigModel->addResonance( "f_2(1270)", 3, LauAbsResonance::RelBW );
    res = sigModel->addResonance( "K*0(892)", 1, LauAbsResonance::RelBW );
    res->fixMass( kFALSE );
    res = sigModel->addResonance( "K*0_0(1430)", 1, LauAbsResonance::LASS );

    // Set the maximum signal DP ASq value
    // If you do not provide a value, one will be determined automatically,
    // which should be close to the true maximum but is not guaranteed to
    // be optimal.
    // Any value, whether manually provided or automatically determined,
    // will be automatically adjusted to avoid bias or extreme inefficiency
    // but it is best to set this by hand once you've found the right value
    // through some trial and error.
    sigModel->setASqMaxValue( 1.25 );

    TString integralsFileName( "integ-" );
    integralsFileName += category;
    integralsFileName += ".dat";
    sigModel->setIntFileName( integralsFileName );

    // Create the fit model
    LauSimpleFitModel* fitModel = new LauSimpleFitModel( sigModel );

    // Create the complex coefficients for the isobar model
    // Here we're using the magnitude and phase form:
    // c_j = a_j exp(i*delta_j)
    std::vector<LauAbsCoeffSet*> coeffset;
    coeffset.push_back( new LauMagPhaseCoeffSet( "f_0(980)", 1.00, 0.00, kTRUE, kTRUE ) );
    coeffset.push_back( new LauMagPhaseCoeffSet( "f_2(1270)", 0.53, 1.39, kFALSE, kFALSE ) );
    coeffset.push_back( new LauMagPhaseCoeffSet( "K*0(892)", 0.87, 1.99, kFALSE, kFALSE ) );
    coeffset.push_back( new LauMagPhaseCoeffSet( "K*0_0(1430)", 1.17, -1.59, kFALSE, kFALSE ) );
    for ( std::vector<LauAbsCoeffSet*>::iterator iter = coeffset.begin(); iter != coeffset.end();
          ++iter ) {
        fitModel->setAmpCoeffSet( *iter );
    }

    // Set the signal yield and define whether it is fixed or floated
    TString sigEventsName = "signalEvents" + category;
    Double_t nSig( 500.0 );
    if ( category == "DD" ) {
        nSig = 750.0;
    }
    LauParameter* nSigEvents = new LauParameter( sigEventsName, nSig, -2.0 * nSig, 2.0 * nSig, kFALSE );
    fitModel->setNSigEvents( nSigEvents );

    // Set the number of experiments to generate or fit and which
    // experiment to start with
    fitModel->setNExpts( nExpt, firstExpt, isToy );

    // Optionally load in continuum background DP model histogram
    // (example syntax given in commented-out section)
    std::vector<TString> bkgndNames( 1 );
    bkgndNames[0] = "comb" + category;
    fitModel->setBkgndClassNames( bkgndNames );
    Double_t nBkgnd = 1200.0;
    if ( category == "DD" ) {
        nBkgnd = 2500.0;
    }
    LauParameter* nBkgndEvents =
        new LauParameter( bkgndNames[0], nBkgnd, -2.0 * nBkgnd, 2.0 * nBkgnd, kFALSE );
    fitModel->setNBkgndEvents( nBkgndEvents );
    //TString bkgndFileName("histoFiles/bkgndDPs.root");
    //TFile* bkgndFile = TFile::Open(bkgndFileName.Data(), "read");
    //TH2* bkgndDP = dynamic_cast<TH2*>(bkgndFile->Get("AllmTheta")); // m', theta'
    LauBkgndDPModel* bkgndModel = new LauBkgndDPModel( daughters, vetoes );
    //bkgndModel->setBkgndHisto(bkgndDP, useInterpolation, fluctuateBins, useUpperHalf, squareDP);
    fitModel->setBkgndDPModel( bkgndNames[0], bkgndModel );

    // Randomise initial fit values for the signal mode
    fitModel->useRandomInitFitPars( kTRUE );

    const Bool_t haveBkgnds = ( fitModel->nBkgndClasses() > 0 );

    // Switch on/off Poissonian smearing of total number of events
    fitModel->doPoissonSmearing( haveBkgnds );

    // Switch on/off Extended ML Fit option
    fitModel->doEMLFit( haveBkgnds );

    // Activate two-stage fit
    fitModel->twoStageFit( kTRUE );

    // Generate toy from the fitted parameters
    //TString fitToyFileName("fitToyMC_");
    //fitToyFileName += category;
    //fitToyFileName += "-Task_";
    //fitToyFileName += iFit;
    //fitToyFileName += ".root";
    //fitModel->compareFitData(100, fitToyFileName);

    // Write out per-event likelihoods and sWeights
    //TString splotFileName("splot_");
    //splotFileName += category;
    //splotFileName += "-Task_";
    //splotFileName += iFit;
    //splotFileName += ".root";
    //fitModel->writeSPlotData(splotFileName, "splot", kFALSE);

    // Set the names of the files to read/write
    TString dataFile( "gen-" );
    dataFile += category;
    dataFile += "-Task.root";
    TString treeName( "genResults" );
    TString rootFileName( "" );
    TString tableFileName( "" );
    if ( command == "fit" ) {
        rootFileName   = "fit";
        rootFileName  += category;
        rootFileName  += "-Task_";
        rootFileName  += iFit;
        rootFileName  += "_expt_";
        rootFileName  += firstExpt;
        rootFileName  += "-";
        rootFileName  += ( firstExpt + nExpt - 1 );
        rootFileName  += ".root";
        tableFileName  = "fit";
        tableFileName += category;
        tableFileName += "TaskResults_";
        tableFileName += iFit;
    } else {
        rootFileName   = "dummy.root";
        tableFileName  = "gen";
        tableFileName += category;
        tableFileName += "TaskResults_";
    }

    // Execute the generation/fit
    if ( command == "fit" ) {
        fitModel->runTask( dataFile, treeName, rootFileName, tableFileName, hostname, port );
    } else {
        fitModel->run( command, dataFile, treeName, rootFileName, tableFileName );
    }

    return EXIT_SUCCESS;
}
