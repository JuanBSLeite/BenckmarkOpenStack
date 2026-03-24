
/*
Copyright 2010 University of Warwick

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

#include "LauArgusPdf.hh"
#include "LauDaughters.hh"
#include "LauEffModel.hh"
#include "LauGaussPdf.hh"
#include "LauIsobarDynamics.hh"
#include "LauLinearPdf.hh"
#include "LauRealImagCoeffSet.hh"
#include "LauSimpleFitModel.hh"
#include "LauSumPdf.hh"
#include "LauVetoes.hh"

#include "TFile.h"
#include "TH2.h"
#include "TString.h"

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
    // ./GenFitNoDPMultiDim gen [nExpt = 1] [firstExpt = 0]
    // or
    // ./GenFitNoDPMultiDim fit <iFit> [nExpt = 1] [firstExpt = 0]
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

    // Still need to define the DP, so just make one up!
    LauDaughters* daughters = new LauDaughters( "B+", "K+", "pi+", "pi-" );

    // Need a dummy vetoes object
    LauVetoes* vetoes = new LauVetoes();

    // Need a dummy efficiency object
    LauEffModel* effModel = new LauEffModel( daughters, vetoes );

    // Create the signal dynamics object and just give it a non-resonant component
    LauIsobarDynamics* sigModel = new LauIsobarDynamics( daughters, effModel );
    sigModel->setIntFileName( "integ.dat" );
    sigModel->addResonance( "NonReson", 0, LauAbsResonance::FlatNR );

    // Create the fit model, passing it the signal dynamics model
    LauSimpleFitModel* fitModel = new LauSimpleFitModel( sigModel );
    fitModel->useDP( kFALSE );

    // Set isobar coefficients for the NR to be 1+0i (and fixed)
    LauAbsCoeffSet* coeffset = new LauRealImagCoeffSet( "NonReson", 1.0, 0.0, kTRUE, kTRUE );
    fitModel->setAmpCoeffSet( coeffset );

    // Define the background categories
    std::vector<TString> bkgndClasses( 1 );
    bkgndClasses[0] = "qqbar";
    fitModel->setBkgndClassNames( bkgndClasses );

    // Signal and background yields
    const Double_t nSigEvents     = 5000.0;
    const Bool_t fixNSigEvents    = kFALSE;
    const Double_t nCombBgEvents  = 7000.0;
    const Bool_t fixNCombBgEvents = kFALSE;

    // Set the number of signal and background events and the number of experiments
    LauParameter* nSig = new LauParameter( "signalEvents",
                                           nSigEvents,
                                           -1.0 * nSigEvents,
                                           2.0 * nSigEvents,
                                           fixNSigEvents );
    LauParameter* nBkg = new LauParameter( "qqbar",
                                           nCombBgEvents,
                                           -1.0 * nCombBgEvents,
                                           2.0 * nCombBgEvents,
                                           fixNCombBgEvents );

    fitModel->setNSigEvents( nSig );
    fitModel->setNBkgndEvents( nBkg );
    fitModel->setNExpts( nExpt, firstExpt, isToy );

    // Signal mES PDF parameter values
    Double_t sig_mes_mean1_value  = 5.279;
    Double_t sig_mes_sigma1_value = 0.002;
    Double_t sig_mes_sigma2_value = 0.007;
    Double_t sig_mes_frac_value   = 0.90;

    // mES PDFs
    Double_t mesMin = 5.260;
    Double_t mesMax = 5.286;
    std::vector<LauAbsRValue*> mesPars;
    mesPars.reserve( 2 );

    // Signal PDF is a double Gaussian with the means constrained to be the same
    LauParameter* sig_mes_mean1 =
        new LauParameter( "sig_mes_mean1", sig_mes_mean1_value, 5.2, 5.3, kTRUE );
    LauParameter* sig_mes_sigma1 =
        new LauParameter( "sig_mes_sigma1", sig_mes_sigma1_value, 0.0, 0.1, kTRUE );
    LauParameter* sig_mes_mean2 = sig_mes_mean1->createClone();
    LauParameter* sig_mes_sigma2 =
        new LauParameter( "sig_mes_sigma2", sig_mes_sigma2_value, 0.0, 0.2, kTRUE );
    LauParameter* sig_mes_frac =
        new LauParameter( "sig_mes_frac", sig_mes_frac_value, 0.0, 1.0, kTRUE );

    mesPars.clear();
    mesPars.push_back( sig_mes_mean1 );
    mesPars.push_back( sig_mes_sigma1 );
    LauAbsPdf* sigMESPdf1 = new LauGaussPdf( "mES", mesPars, mesMin, mesMax );

    mesPars.clear();
    mesPars.push_back( sig_mes_mean2 );
    mesPars.push_back( sig_mes_sigma2 );
    LauAbsPdf* sigMESPdf2 = new LauGaussPdf( "mES", mesPars, mesMin, mesMax );

    LauAbsPdf* sigMESPdf = new LauSumPdf( sigMESPdf1, sigMESPdf2, sig_mes_frac );
    fitModel->setSignalPdf( sigMESPdf );

    // Combinatoric background PDF is an ARGUS function
    LauParameter* comb_mes_m0 = new LauParameter( "comb_mes_m0", 5.289 );
    LauParameter* comb_mes_xi = new LauParameter( "comb_mes_xi", 20.0, 0.0, 50.0, kTRUE );

    mesPars.clear();
    mesPars.push_back( comb_mes_m0 );
    mesPars.push_back( comb_mes_xi );
    LauAbsPdf* combMESPdf = new LauArgusPdf( "mES", mesPars, mesMin, mesMax );

    fitModel->setBkgndPdf( "qqbar", combMESPdf );

    // Signal DeltaE PDF parameter values
    Double_t sig_de_mean1_value  = -0.005;
    Double_t sig_de_sigma1_value = 0.020;
    Double_t sig_de_sigma2_value = 0.050;
    Double_t sig_de_frac_value   = 0.90;

    // DeltaE PDFs
    Double_t deMin = -0.2;
    Double_t deMax = 0.2;
    std::vector<LauAbsRValue*> dePars;
    dePars.reserve( 2 );

    // Signal PDF is a double Gaussian with the means constrained to be the same
    LauParameter* sig_de_mean1 =
        new LauParameter( "sig_de_mean1", sig_de_mean1_value, -0.1, 0.1, kTRUE );
    LauParameter* sig_de_sigma1 =
        new LauParameter( "sig_de_sigma1", sig_de_sigma1_value, 0.0, 0.2, kTRUE );
    LauParameter* sig_de_mean2 = sig_de_mean1->createClone();
    LauParameter* sig_de_sigma2 =
        new LauParameter( "sig_de_sigma2", sig_de_sigma2_value, 0.0, 0.2, kTRUE );
    LauParameter* sig_de_frac = new LauParameter( "sig_de_frac", sig_de_frac_value, 0.0, 1.0, kTRUE );

    dePars.clear();
    dePars.push_back( sig_de_mean1 );
    dePars.push_back( sig_de_sigma1 );
    LauAbsPdf* sigDEPdf1 = new LauGaussPdf( "DeltaE", dePars, deMin, deMax );

    dePars.clear();
    dePars.push_back( sig_de_mean2 );
    dePars.push_back( sig_de_sigma2 );
    LauAbsPdf* sigDEPdf2 = new LauGaussPdf( "DeltaE", dePars, deMin, deMax );

    LauAbsPdf* sigDEPdf = new LauSumPdf( sigDEPdf1, sigDEPdf2, sig_de_frac );
    fitModel->setSignalPdf( sigDEPdf );

    // Combinatoric background PDF is a linear function
    LauParameter* comb_de_slope = new LauParameter( "comb_de_slope", -5.0, -10.0, 10.0, kTRUE );

    dePars.clear();
    dePars.push_back( comb_de_slope );
    LauAbsPdf* combDEPdf = new LauLinearPdf( "DeltaE", dePars, deMin, deMax );

    fitModel->setBkgndPdf( "qqbar", combDEPdf );

    // Configure various fit options

    // Do not calculate asymmetric errors.
    fitModel->useAsymmFitErrors( kFALSE );

    // Randomise initial fit values for the signal isobar parameters
    fitModel->useRandomInitFitPars( kTRUE );

    // Switch on/off Poissonian smearing of total number of events
    fitModel->doPoissonSmearing( kTRUE );

    // Switch on/off Extended ML Fit option
    fitModel->doEMLFit( kTRUE );

    // Switch on/off two-stage fitting for CPV parameters
    fitModel->twoStageFit( kFALSE );

    // Generate toy from the fitted parameters
    //TString fitToyFileName("fitToyMC_NoDPMultiDim_");
    //fitToyFileName += iFit;
    //fitToyFileName += ".root";
    //fitModel->compareFitData(100, fitToyFileName);

    // Write out per-event likelihoods and sWeights
    //TString splotFileName("splot_NoDPMultiDim_");
    //splotFileName += iFit;
    //splotFileName += ".root";
    //fitModel->writeSPlotData(splotFileName, "splot", kFALSE);

    TString dataFile( "gen-NoDPMultiDim.root" );
    TString treeName( "genResults" );
    TString rootFileName( "" );
    TString tableFileName( "" );
    if ( command == "fit" ) {
        rootFileName   = "fitNoDPMultiDim_";
        rootFileName  += iFit;
        rootFileName  += "_expt_";
        rootFileName  += firstExpt;
        rootFileName  += "-";
        rootFileName  += firstExpt + nExpt - 1;
        rootFileName  += ".root";
        tableFileName  = "fitNoDPMultiDimResults_";
        tableFileName += iFit;
    } else {
        rootFileName  = "dummy.root";
        tableFileName = "genNoDPMultiDimResults";
    }

    // Execute the generation/fit
    fitModel->run( command, dataFile, treeName, rootFileName, tableFileName );

    return EXIT_SUCCESS;
}
