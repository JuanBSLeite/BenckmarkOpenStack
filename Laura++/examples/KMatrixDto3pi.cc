
/*
Copyright 2008 University of Warwick

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

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

#include "LauAbsCoeffSet.hh"
#include "LauCartesianCPCoeffSet.hh"
#include "LauDaughters.hh"
#include "LauEffModel.hh"
#include "LauIsobarDynamics.hh"
#include "LauSimpleFitModel.hh"
#include "LauVetoes.hh"

#include "TString.h"

#include <map>
#include <vector>

void usage( std::ostream& out, const TString& progName )
{
    out << "Usage:\n";
    out << progName << " gen [nExpt = 1] [firstExpt = 0]\n";
    out << "or\n";
    out << progName << " fit <iFit> [nExpt = 1] [firstExpt = 0]" << std::endl;
}

int main( const int argc, const char** argv )
{

    // Process command-line arguments
    // Usage:
    // ./KMatrixDto3pi gen [nExpt = 1] [firstExpt = 0]
    // or
    // ./KMatrixDto3pi fit <iFit> [nExpt = 1] [firstExpt = 0]
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

    // Generate K-matrix amplitude distribution across the DP
    LauDaughters* daughters = new LauDaughters( "D+", "pi+", "pi+", "pi-", kFALSE );
    LauVetoes* vetoes       = new LauVetoes();
    LauEffModel* effModel   = new LauEffModel( daughters, vetoes );

    LauIsobarDynamics* sigModel = new LauIsobarDynamics( daughters, effModel );

    // Define the "S-wave" K-matrix propagator
    Int_t nChannels    = 5;
    Int_t nPoles       = 5;
    Int_t resPairInt   = 1;
    Int_t KMatrixIndex = 1;    // for S-wave
    sigModel->defineKMatrixPropagator( "KMSWave",
                                       "FOCUSD3pi.dat",
                                       resPairInt,
                                       nChannels,
                                       nPoles,
                                       KMatrixIndex );
    // Now define the pole and "slowly-varying part" (SVP) amplitudes for the K-matrix.
    // Each part will have their own complex coefficients which can (should) be fitted.
    sigModel->addKMatrixProdPole( "KMPole1", "KMSWave", 1 );
    sigModel->addKMatrixProdPole( "KMPole2", "KMSWave", 2 );
    sigModel->addKMatrixProdPole( "KMPole3", "KMSWave", 3 );
    //sigModel->addKMatrixProdPole("KMPole4", "KMSWave", 4);
    //sigModel->addKMatrixProdPole("KMPole5", "KMSWave", 5);

    sigModel->addKMatrixProdSVP( "KMSVP1", "KMSWave", 1 );
    sigModel->addKMatrixProdSVP( "KMSVP2", "KMSWave", 2 );
    //sigModel->addKMatrixProdSVP("KMSVP3", "KMSWave", 3);
    //sigModel->addKMatrixProdSVP("KMSVP4", "KMSWave", 4);
    //sigModel->addKMatrixProdSVP("KMSVP5", "KMSWave", 5);

    sigModel->addResonance( "rho0(770)", 1, LauAbsResonance::RelBW );
    sigModel->addResonance( "f_2(1270)", 1, LauAbsResonance::RelBW );

    sigModel->setASqMaxValue( 2.1 );

    LauSimpleFitModel* fitModel = new LauSimpleFitModel( sigModel );

    std::vector<LauAbsCoeffSet*> coeffset;

    coeffset.push_back(
        new LauCartesianCPCoeffSet( "KMPole1", 0.85, -0.55, 0.0, 0.0, kFALSE, kFALSE, kTRUE, kTRUE ) );
    coeffset.push_back(
        new LauCartesianCPCoeffSet( "KMPole2", -1.10, -0.30, 0.0, 0.0, kFALSE, kFALSE, kTRUE, kTRUE ) );
    coeffset.push_back(
        new LauCartesianCPCoeffSet( "KMPole3", -0.84, 0.35, 0.0, 0.0, kFALSE, kFALSE, kTRUE, kTRUE ) );
    //coeffset.push_back(new LauCartesianCPCoeffSet("KMPole4",   -0.004, .46,  0.0, 0.0, kFALSE, kFALSE, kTRUE, kTRUE));
    //coeffset.push_back(new LauCartesianCPCoeffSet("KMPole5",   0.00, 0.00,  0.0, 0.0, kFALSE, kFALSE, kTRUE, kTRUE));

    coeffset.push_back(
        new LauCartesianCPCoeffSet( "KMSVP1", 0.230, 0.352, 0.0, 0.0, kFALSE, kFALSE, kTRUE, kTRUE ) );
    coeffset.push_back(
        new LauCartesianCPCoeffSet( "KMSVP2", 0.330, 0.42, 0.0, 0.0, kFALSE, kFALSE, kTRUE, kTRUE ) );
    //coeffset.push_back(new LauCartesianCPCoeffSet("KMSVP3",   -0.48, 1.25,  0.0, 0.0, kFALSE, kFALSE, kTRUE, kTRUE));
    //coeffset.push_back(new LauCartesianCPCoeffSet("KMSVP4",   0.30, 1.32,  0.0, 0.0, kFALSE, kFALSE, kTRUE, kTRUE));
    //coeffset.push_back(new LauCartesianCPCoeffSet("KMSVP5",   0.00, 0.00,  0.0, 0.0, kFALSE, kFALSE, kTRUE, kTRUE));

    coeffset.push_back(
        new LauCartesianCPCoeffSet( "rho0(770)", 1.0, 0.00, 0.0, 0.0, kTRUE, kTRUE, kTRUE, kTRUE ) );
    coeffset.push_back(
        new LauCartesianCPCoeffSet( "f_2(1270)", 0.45, 0.35, 0.0, 0.0, kFALSE, kFALSE, kTRUE, kTRUE ) );

    for ( std::vector<LauAbsCoeffSet*>::iterator iter = coeffset.begin(); iter != coeffset.end();
          ++iter ) {
        fitModel->setAmpCoeffSet( *iter );
    }

    Double_t nSigEvents = 20000;
    Bool_t fixSigEvents = kTRUE;
    LauParameter* nSig  = new LauParameter( "signalEvents",
                                           nSigEvents,
                                           -2.0 * nSigEvents,
                                           2.0 * nSigEvents,
                                           fixSigEvents );
    fitModel->setNSigEvents( nSig );
    fitModel->setNExpts( nExpt, firstExpt, isToy );

    // Do not calculate asymmetric errors
    fitModel->useAsymmFitErrors( kFALSE );

    // Randomise initial fit values for the signal mode
    fitModel->useRandomInitFitPars( kTRUE );

    // Switch off Poissonian smearing of total number of events
    fitModel->doPoissonSmearing( kFALSE );

    // Switch on Extended ML Fit option - only really use if Poisson smearing is also enabled
    fitModel->doEMLFit( kFALSE );

    // Set the names of the files to read/write
    TString dataFile( "gen-KMatrixDto3pi.root" );
    TString treeName( "genResults" );
    TString rootFileName( "" );
    TString tableFileName( "" );
    if ( command == "fit" ) {
        rootFileName   = "fitKMatrixDto3pi_";
        rootFileName  += iFit;
        rootFileName  += "_expt_";
        rootFileName  += firstExpt;
        rootFileName  += "-";
        rootFileName  += ( firstExpt + nExpt - 1 );
        rootFileName  += ".root";
        tableFileName  = "fitKMatrixDto3piResults_";
        tableFileName += iFit;
    } else {
        rootFileName  = "dummy.root";
        tableFileName = "genKMatrixDto3piResults";
    }

    // Execute the generation/fit
    fitModel->run( command, dataFile, treeName, rootFileName, tableFileName );
}
