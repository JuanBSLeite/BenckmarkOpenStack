
/*
Copyright 2005 University of Warwick

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
#include "LauResonanceMaker.hh"
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
    out << progName << " gen [nExpt = 1] [firstExpt = 0]\n";
    out << "or\n";
    out << progName << " fit <iFit> [nExpt = 1] [firstExpt = 0]" << std::endl;
}

int main( int argc, char** argv )
{
    // Process command-line arguments
    // Usage:
    // ./GenFit3pi gen [nExpt = 1] [firstExpt = 0]
    // or
    // ./GenFit3pi fit <iFit> [nExpt = 1] [firstExpt = 0]
    if ( argc < 2 ) {
        usage( std::cerr, argv[0] );
        return EXIT_FAILURE;
    }

    TString command = argv[1];
    command.ToLower();
    Int_t iFit( 0 );
    Int_t nExpt( 1 );
    Int_t firstExpt( 0 );
    Bool_t isToy( kTRUE );
    if ( command == "gen" ) {
        if ( argc > 2 ) {
            nExpt = atoi( argv[2] );
            if ( argc > 3 ) {
                firstExpt = atoi( argv[3] );
            }
        }
    } else if ( command == "fit" ) {
        if ( argc < 3 ) {
            usage( std::cerr, argv[0] );
            return EXIT_FAILURE;
        }
        iFit = atoi( argv[2] );
        if ( argc > 3 ) {
            nExpt = atoi( argv[3] );
            if ( argc > 4 ) {
                firstExpt = atoi( argv[4] );
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

    // This defines the DP => decay is B+ -> pi+ pi+ pi-
    // Particle 1 = pi+
    // Particle 2 = pi+
    // Particle 3 = pi-
    // The DP is defined in terms of m13Sq and m23Sq
    LauDaughters* daughters = new LauDaughters( "B+", "pi+", "pi+", "pi-", squareDP );

    // Optionally apply some vetoes to the DP
    // (example syntax given but commented-out)
    LauVetoes* vetoes = new LauVetoes();
    //Double_t DMin = 1.70;
    //Double_t DMax = 1.925;
    //Double_t JpsiMin = 3.051;
    //Double_t JpsiMax = 3.222;
    //Double_t psi2SMin = 3.676;
    //Double_t psi2SMax = 3.866;
    //vetoes->addMassVeto(1, DMin, DMax); // D0 veto, m23 (and automatically m13)
    //vetoes->addMassVeto(1, JpsiMin, JpsiMax); // J/psi veto, m23 (and automatically m13)
    //vetoes->addMassVeto(1, psi2SMin, psi2SMax); // psi(2S) veto, m23 (and automatically m13)

    // Define the efficiency model (defaults to unity everywhere)
    LauEffModel* effModel = new LauEffModel( daughters, vetoes );

    // Can optionally provide a histogram to model variation over DP
    // (example syntax given but commented-out)
    //TFile *effHistFile = TFile::Open("histoFiles/B3piNRDPEff.root", "read");
    //TH2* effHist = dynamic_cast<TH2*>(effHistFile->Get("effHist"));
    //Bool_t useInterpolation = kTRUE;
    //Bool_t fluctuateBins = kFALSE;
    //Bool_t useUpperHalf = kTRUE;
    //effModel->setEffHisto(effHist, useInterpolation, fluctuateBins, 0.0, 0.0, useUpperHalf, squareDP);

    // Create the isobar model

    // Set the values of the Blatt-Weisskopf barrier radii and whether they are fixed or floating
    LauResonanceMaker& resMaker = LauResonanceMaker::get();
    resMaker.setDefaultBWRadius( LauBlattWeisskopfFactor::Parent, 5.0 );
    resMaker.setDefaultBWRadius( LauBlattWeisskopfFactor::Light, 4.0 );
    resMaker.fixBWRadius( LauBlattWeisskopfFactor::Parent, kTRUE );
    resMaker.fixBWRadius( LauBlattWeisskopfFactor::Light, kTRUE );

    LauIsobarDynamics* sigModel = new LauIsobarDynamics( daughters, effModel );
    // Add various components to the isobar model,
    // modifying some resonance parameters
    LauAbsResonance* reson( 0 );
    // addResonance arguments: resName, resPairAmpInt, resType
    // resPairAmpInt = 1 => resonance mass is m23.
    reson = sigModel->addResonance( "rho0(770)", 1, LauAbsResonance::GS );
    reson = sigModel->addResonance( "rho0(1450)", 1, LauAbsResonance::RelBW );
    reson = sigModel->addResonance( "f_0(980)", 1, LauAbsResonance::Flatte );
    reson->setResonanceParameter( "g1", 0.2 );
    reson->setResonanceParameter( "g2", 1.0 );
    reson = sigModel->addResonance( "f_2(1270)", 1, LauAbsResonance::RelBW );
    const TString nrName { "BelleNR_Swave" };
    reson = sigModel->addResonance( nrName, 0, LauAbsResonance::BelleSymNRNoInter );
    reson->setResonanceParameter( "alpha", 0.2 );

    // Set the maximum signal DP ASq value
    // If you do not provide a value, one will be determined automatically,
    // which should be close to the true maximum but is not guaranteed to
    // be optimal.
    // Any value, whether manually provided or automatically determined,
    // will be automatically adjusted to avoid bias or extreme inefficiency
    // but it is best to set this by hand once you've found the right value
    // through some trial and error.
    sigModel->setASqMaxValue( 0.35 );

    // Create the fit model
    LauSimpleFitModel* fitModel = new LauSimpleFitModel( sigModel );

    // Create the complex coefficients for the isobar model
    // Here we're using the magnitude and phase form:
    // c_j = a_j exp(i*delta_j)
    std::vector<LauAbsCoeffSet*> coeffset;
    coeffset.push_back( new LauMagPhaseCoeffSet( "rho0(770)", 1.00, 0.00, kTRUE, kTRUE ) );
    coeffset.push_back( new LauMagPhaseCoeffSet( "rho0(1450)", 0.37, 1.99, kFALSE, kFALSE ) );
    coeffset.push_back( new LauMagPhaseCoeffSet( "f_0(980)", 0.27, -1.59, kFALSE, kFALSE ) );
    coeffset.push_back( new LauMagPhaseCoeffSet( "f_2(1270)", 0.53, 1.39, kFALSE, kFALSE ) );
    coeffset.push_back( new LauMagPhaseCoeffSet( nrName, 0.54, -0.84, kFALSE, kFALSE ) );
    for ( std::vector<LauAbsCoeffSet*>::iterator iter = coeffset.begin(); iter != coeffset.end();
          ++iter ) {
        fitModel->setAmpCoeffSet( *iter );
    }

    // Set the signal yield and define whether it is fixed or floated
    const Double_t nSigEvents = 1500.0;
    LauParameter* signalEvents =
        new LauParameter( "signalEvents", nSigEvents, -1.0 * nSigEvents, 2.0 * nSigEvents, kFALSE );
    fitModel->setNSigEvents( signalEvents );

    // Set the number of experiments to generate or fit and which
    // experiment to start with
    fitModel->setNExpts( nExpt, firstExpt, isToy );

    // Set up a background model

    // First declare the names of the background class(es)
    std::vector<TString> bkgndNames( 1 );
    bkgndNames[0] = "qqbar";
    fitModel->setBkgndClassNames( bkgndNames );

    // Define and set the yield parameter for the background
    const Double_t nBkg        = 1250.0;
    LauParameter* nBkgndEvents = new LauParameter( "qqbar", nBkg, -2.0 * nBkg, 2.0 * nBkg, kFALSE );
    fitModel->setNBkgndEvents( nBkgndEvents );

    // Create the background DP model
    LauBkgndDPModel* qqbarModel = new LauBkgndDPModel( daughters, vetoes );

    // Load in background DP model histogram
    // (example syntax given but commented-out - the background will be treated as being uniform in the DP in the absence of a histogram)
    //TString qqFileName("histoFiles/offResDP.root");
    //TFile* qqFile = TFile::Open(qqFileName.Data(), "read");
    //TH2* qqDP = dynamic_cast<TH2*>(qqFile->Get("AllmTheta")); // m', theta'
    //qqbarModel->setBkgndHisto(qqDP, useInterpolation, fluctuateBins, useUpperHalf, squareDP);

    // Add the background DP model into the fit model
    fitModel->setBkgndDPModel( "qqbar", qqbarModel );

    // Configure various fit options

    // Switch on/off calculation of asymmetric errors.
    fitModel->useAsymmFitErrors( kFALSE );

    // Randomise initial fit values for the signal mode
    fitModel->useRandomInitFitPars( kTRUE );

    const Bool_t haveBkgnds = ( fitModel->nBkgndClasses() > 0 );

    // Switch on/off Poissonian smearing of total number of events
    fitModel->doPoissonSmearing( haveBkgnds );

    // Switch on/off Extended ML Fit option
    fitModel->doEMLFit( haveBkgnds );

    // Generate toy from the fitted parameters
    //TString fitToyFileName("fitToyMC_3pi_");
    //fitToyFileName += iFit;
    //fitToyFileName += ".root";
    //fitModel->compareFitData(100, fitToyFileName);

    // Write out per-event likelihoods and sWeights
    //TString splotFileName("splot_3pi_");
    //splotFileName += iFit;
    //splotFileName += ".root";
    //fitModel->writeSPlotData(splotFileName, "splot", kFALSE);

    // Set the names of the files to read/write
    TString dataFile( "gen-3pi.root" );
    TString treeName( "genResults" );
    TString rootFileName( "" );
    TString tableFileName( "" );
    if ( command == "fit" ) {
        rootFileName   = "fit3pi_";
        rootFileName  += iFit;
        rootFileName  += "_expt_";
        rootFileName  += firstExpt;
        rootFileName  += "-";
        rootFileName  += ( firstExpt + nExpt - 1 );
        rootFileName  += ".root";
        tableFileName  = "fit3piResults_";
        tableFileName += iFit;
    } else {
        rootFileName  = "dummy.root";
        tableFileName = "gen3piResults";
    }

    // Execute the generation/fit
    fitModel->run( command, dataFile, treeName, rootFileName, tableFileName );

    return EXIT_SUCCESS;
}
