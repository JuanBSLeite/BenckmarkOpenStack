
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
#include "LauDaughters.hh"
#include "LauEffModel.hh"
#include "LauIsobarDynamics.hh"
#include "LauRealImagCoeffSet.hh"
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

int main( int argc, char** argv )
{
    // Process command-line arguments
    // Usage:
    // ./KMatrixExample gen [nExpt = 1] [firstExpt = 0]
    // or
    // ./KMatrixExample fit <iFit> [nExpt = 1] [firstExpt = 0]
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

    // This defines the DP => decay is B0 -> pi+ pi- K_S0
    // Particle 1 = pi+
    // Particle 2 = pi-
    // Particle 3 = K_S0
    // The DP is defined in terms of m13Sq and m23Sq
    LauDaughters* daughters = new LauDaughters( "B0", "pi+", "pi-", "K_S0", squareDP );

    LauVetoes* vetoes     = new LauVetoes();
    LauEffModel* effModel = new LauEffModel( daughters, vetoes );

    LauIsobarDynamics* sigModel = new LauIsobarDynamics( daughters, effModel );

    // Add relativistic BWs for the rho and K* resonances
    sigModel->addResonance( "rho0(770)", 3, LauAbsResonance::RelBW );
    sigModel->addResonance( "K*+(892)", 2, LauAbsResonance::RelBW );
    sigModel->addResonance( "K*+_2(1430)", 2, LauAbsResonance::RelBW );

    // Define the pi+ pi- "S-wave" K-matrix propagator
    // The defineKMatrixPropagator requires four arguments
    // name of propagator, parameter file, mass pair index (1,2,3) for "s" variable, number of channels, number of
    // mass poles and the row index defining the K-matrix state (e.g. index = 1 for the first row -> S-wave).
    Int_t resPairInt   = 3;
    Int_t nChannels    = 5;
    Int_t nPoles       = 5;
    Int_t KMatrixIndex = 1;
    sigModel->defineKMatrixPropagator( "KM_pipi",
                                       "LauKMatrixCoeff.dat",
                                       resPairInt,
                                       nChannels,
                                       nPoles,
                                       KMatrixIndex );

    // Now define the pole and "slowly-varying part" (SVP) amplitudes for this K-matrix.
    // Each part will have their own complex coefficients which can (should) be fitted.
    // The addKMatrixProd functions need: name of pole/SVP, name of propagator and the index
    // number of the amplitude (from 1 to nPoles or nChannels)
    sigModel->addKMatrixProdPole( "KM_pipi_Pole1", "KM_pipi", 1 );
    sigModel->addKMatrixProdPole( "KM_pipi_Pole2", "KM_pipi", 2 );
    sigModel->addKMatrixProdPole( "KM_pipi_Pole3", "KM_pipi", 3 );
    sigModel->addKMatrixProdPole( "KM_pipi_Pole4", "KM_pipi", 4 );
    sigModel->addKMatrixProdPole( "KM_pipi_Pole5", "KM_pipi", 5 );

    sigModel->addKMatrixProdSVP( "KM_pipi_SVP1", "KM_pipi", 1 );
    sigModel->addKMatrixProdSVP( "KM_pipi_SVP2", "KM_pipi", 2 );
    sigModel->addKMatrixProdSVP( "KM_pipi_SVP3", "KM_pipi", 3 );

    sigModel->setASqMaxValue( 3.1 );

    LauSimpleFitModel* fitModel = new LauSimpleFitModel( sigModel );

    std::vector<LauAbsCoeffSet*> coeffset;

    // Define the coefficients for the P-wave rho and K* resonances
    coeffset.push_back( new LauRealImagCoeffSet( "rho0(770)", 1.00, 0.00, kTRUE, kTRUE ) );
    coeffset.push_back( new LauRealImagCoeffSet( "K*+(892)", 0.97, -1.10, kFALSE, kFALSE ) );
    coeffset.push_back( new LauRealImagCoeffSet( "K*+_2(1430)", 0.38, -0.54, kFALSE, kFALSE ) );

    // Define the coefficients for the K-matrix amplitudes defined above ("beta" and "f" production coeffs)
    coeffset.push_back( new LauRealImagCoeffSet( "KM_pipi_Pole1", -0.18, -0.91, kFALSE, kFALSE ) );
    coeffset.push_back( new LauRealImagCoeffSet( "KM_pipi_Pole2", -1.02, -0.38, kFALSE, kFALSE ) );
    coeffset.push_back( new LauRealImagCoeffSet( "KM_pipi_Pole3", -0.37, 0.50, kFALSE, kFALSE ) );
    coeffset.push_back( new LauRealImagCoeffSet( "KM_pipi_Pole4", -0.01, 0.92, kFALSE, kFALSE ) );
    coeffset.push_back( new LauRealImagCoeffSet( "KM_pipi_Pole5", 0.01, -0.05, kFALSE, kFALSE ) );

    coeffset.push_back( new LauRealImagCoeffSet( "KM_pipi_SVP1", 0.22, 0.76, kFALSE, kFALSE ) );
    coeffset.push_back( new LauRealImagCoeffSet( "KM_pipi_SVP2", 0.89, -0.31, kFALSE, kFALSE ) );
    coeffset.push_back( new LauRealImagCoeffSet( "KM_pipi_SVP3", -0.97, 0.93, kFALSE, kFALSE ) );

    for ( std::vector<LauAbsCoeffSet*>::iterator iter = coeffset.begin(); iter != coeffset.end();
          ++iter ) {
        fitModel->setAmpCoeffSet( *iter );
    }

    Double_t nSigEvents = 5000;
    Bool_t fixSigEvents = kTRUE;
    LauParameter* nSig  = new LauParameter( "signalEvents",
                                           nSigEvents,
                                           -1.0 * nSigEvents,
                                           2.0 * nSigEvents,
                                           fixSigEvents );
    fitModel->setNSigEvents( nSig );
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
    //TString fitToyFileName("fitToyMC_KMatrixExample_");
    //fitToyFileName += iFit;
    //fitToyFileName += ".root";
    //fitModel->compareFitData(100, fitToyFileName);

    // Write out per-event likelihoods and sWeights
    //TString splotFileName("splot_KMatrixExample_");
    //splotFileName += iFit;
    //splotFileName += ".root";
    //fitModel->writeSPlotData(splotFileName, "splot", kFALSE);

    // Set the names of the files to read/write
    TString dataFile( "gen-KMatrixExample.root" );
    TString treeName( "genResults" );
    TString rootFileName( "" );
    TString tableFileName( "" );
    if ( command == "fit" ) {
        rootFileName   = "fitKMatrixExample_";
        rootFileName  += iFit;
        rootFileName  += "_expt_";
        rootFileName  += firstExpt;
        rootFileName  += "-";
        rootFileName  += ( firstExpt + nExpt - 1 );
        rootFileName  += ".root";
        tableFileName  = "fitKMatrixExampleResults_";
        tableFileName += iFit;
    } else {
        rootFileName  = "dummy.root";
        tableFileName = "genKMatrixExampleResults";
    }

    // Execute the generation/fit
    fitModel->run( command, dataFile, treeName, rootFileName, tableFileName );

    return EXIT_SUCCESS;
}
