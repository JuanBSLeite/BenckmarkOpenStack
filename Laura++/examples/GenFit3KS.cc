
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

#include "LauBkgndDPModel.hh"
#include "LauDaughters.hh"
#include "LauEffModel.hh"
#include "LauIsobarDynamics.hh"
#include "LauMagPhaseCoeffSet.hh"
#include "LauSimpleFitModel.hh"
#include "LauVetoes.hh"

#include "TFile.h"
#include "TH2.h"
#include "TString.h"
#include "TTree.h"

#include <cstdlib>
#include <iostream>
#include <vector>

void usage( std::ostream& out, const char* progName )
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
    // ./GenFit3KS gen [nExpt = 1] [firstExpt = 0]
    // or
    // ./GenFit3KS fit <iFit> [nExpt = 1] [firstExpt = 0]
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

    // This defines the DP => decay is B0 -> KS KS KS
    // Particle 1 = KS
    // Particle 2 = KS
    // Particle 3 = KS
    // The DP is defined in terms of m13Sq and m23Sq
    LauDaughters* daughters = new LauDaughters( "B0", "K_S0", "K_S0", "K_S0", squareDP );

    // Optionally apply some vetoes to the DP
    LauVetoes* vetoes = new LauVetoes();

    // Define the efficiency model (defaults to unity everywhere)
    // Can optionally provide a histogram to model variation over DP
    LauEffModel* effModel = new LauEffModel( daughters, vetoes );

    // Create the isobar model
    LauIsobarDynamics* sigModel = new LauIsobarDynamics( daughters, effModel );
    //LauAbsResonance* reson(0);
    /*reson =*/sigModel->addResonance( "f_0(980)", 3, LauAbsResonance::Flatte );
    /*reson =*/sigModel->addResonance( "f_0(1710)", 3, LauAbsResonance::RelBW );
    /*reson =*/sigModel->addResonance( "f_2(2010)", 3, LauAbsResonance::RelBW );
    /*reson =*/sigModel->addResonance( "chi_c0", 3, LauAbsResonance::RelBW );

    // Set the maximum signal DP ASq value
    // If you do not provide a value, one will be determined automatically,
    // which should be close to the true maximum but is not guaranteed to
    // be optimal.
    // Any value, whether manually provided or automatically determined,
    // will be automatically adjusted to avoid bias or extreme inefficiency
    // but it is best to set this by hand once you've found the right value
    // through some trial and error.
    sigModel->setASqMaxValue( 0.285 );

    // Create the fit model
    LauSimpleFitModel* fitModel = new LauSimpleFitModel( sigModel );

    // Create the complex coefficients for the isobar model
    // Here we're using the magnitude and phase form:
    // c_j = a_j exp(i*delta_j)
    std::vector<LauAbsCoeffSet*> coeffset;
    coeffset.push_back( new LauMagPhaseCoeffSet( "f_0(980)", 1.00, 0.00, kTRUE, kTRUE ) );
    coeffset.push_back( new LauMagPhaseCoeffSet( "f_0(1710)", 0.40, 1.11, kFALSE, kFALSE ) );
    coeffset.push_back( new LauMagPhaseCoeffSet( "f_2(2010)", 0.45, 2.50, kFALSE, kFALSE ) );
    coeffset.push_back( new LauMagPhaseCoeffSet( "chi_c0", 0.40, 0.63, kFALSE, kFALSE ) );
    for ( std::vector<LauAbsCoeffSet*>::iterator iter = coeffset.begin(); iter != coeffset.end();
          ++iter ) {
        fitModel->setAmpCoeffSet( *iter );
    }

    // Set the signal yield and define whether it is fixed or floated
    LauParameter* nSigEvents = new LauParameter( "nSigEvents", 10000.0, -50000.0, 50000.0, kFALSE );
    fitModel->setNSigEvents( nSigEvents );

    // Set the number of experiments to generate or fit and which
    // experiment to start with
    fitModel->setNExpts( nExpt, firstExpt, isToy );

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
    //TString fitToyFileName("fitToyMC_3KS_");
    //fitToyFileName += iFit;
    //fitToyFileName += ".root";
    //fitModel->compareFitData(100, fitToyFileName);

    // Write out per-event likelihoods and sWeights
    //TString splotFileName("splot_3KS_");
    //splotFileName += iFit;
    //splotFileName += ".root";
    //fitModel->writeSPlotData(splotFileName, "splot", kFALSE);

    // Set the names of the files to read/write
    TString dataFile( "gen-3KS.root" );
    TString treeName( "genResults" );
    TString rootFileName( "" );
    TString tableFileName( "" );
    if ( command == "fit" ) {
        rootFileName   = "fit3KS_";
        rootFileName  += iFit;
        rootFileName  += "_expt_";
        rootFileName  += firstExpt;
        rootFileName  += "-";
        rootFileName  += ( firstExpt + nExpt - 1 );
        rootFileName  += ".root";
        tableFileName  = "fit3KSResults_";
        tableFileName += iFit;
    } else {
        rootFileName  = "dummy.root";
        tableFileName = "gen3KSResults";
    }

    // Execute the generation/fit
    fitModel->run( command, dataFile, treeName, rootFileName, tableFileName );

    return EXIT_SUCCESS;
}
