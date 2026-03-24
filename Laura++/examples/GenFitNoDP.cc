
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

#include "LauCrystalBallPdf.hh"
#include "LauDaughters.hh"
#include "LauEffModel.hh"
#include "LauExponentialPdf.hh"
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
    // ./GenFitNoDP gen [nExpt = 1] [firstExpt = 0]
    // or
    // ./GenFitNoDP fit <iFit> [nExpt = 1] [firstExpt = 0]
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

    // Set up the additional background configuration
    const Int_t nBkgnds( 2 );
    std::vector<TString> bkgndNames( nBkgnds );
    bkgndNames[0] = "Combinatorial";
    bkgndNames[1] = "PartialReco";
    fitModel->setBkgndClassNames( bkgndNames );

    // Signal and background yields
    Double_t nSigEvents         = 5000.0;
    Bool_t fixNSigEvents        = kFALSE;
    Double_t nCombBgEvents      = 7000.0;
    Bool_t fixNCombBgEvents     = kFALSE;
    Double_t nPartRecoBgEvents  = 3000.0;
    Bool_t fixNPartRecoBgEvents = kFALSE;
    Double_t nTotEvents         = nSigEvents + nCombBgEvents + nPartRecoBgEvents;

    // Set the number of signal and background events and the number of experiments
    LauParameter* nSig     = new LauParameter( "signalEvents",
                                           nSigEvents,
                                           -1.0 * nTotEvents,
                                           2.0 * nTotEvents,
                                           fixNSigEvents );
    LauParameter* nCombBkg = new LauParameter( "Combinatorial",
                                               nCombBgEvents,
                                               -1.0 * nTotEvents,
                                               2.0 * nTotEvents,
                                               fixNCombBgEvents );
    LauParameter* nPartBkg = new LauParameter( "PartialReco",
                                               nPartRecoBgEvents,
                                               -1.0 * nTotEvents,
                                               2.0 * nTotEvents,
                                               fixNPartRecoBgEvents );
    // Optionally blind the yield parameters
    //nSig->blindParameter("something1",2000.0);
    //nCombBkg->blindParameter("something2",2000.0);
    //nPartBkg->blindParameter("something3",2000.0);
    fitModel->setNSigEvents( nSig );
    fitModel->setNBkgndEvents( nCombBkg );
    fitModel->setNBkgndEvents( nPartBkg );
    fitModel->setNExpts( nExpt, firstExpt, isToy );

    // Signal mB PDF parameter values
    Double_t sig_mb_mean1_value  = 5.279;
    Double_t sig_mb_sigma1_value = 0.02;
    Double_t sig_mb_sigma2_value = 0.07;
    Double_t sig_mb_frac_value   = 0.90;

    // m_B PDFs
    Double_t mbMin = 5.150;
    Double_t mbMax = 5.600;
    std::vector<LauAbsRValue*> mbPars;
    mbPars.reserve( 2 );

    // Signal PDF is a double Gaussian with the means constrained to be the same
    LauParameter* sig_mb_mean1 =
        new LauParameter( "sig_mb_mean1", sig_mb_mean1_value, 5.2, 5.3, kTRUE );
    LauParameter* sig_mb_sigma1 =
        new LauParameter( "sig_mb_sigma1", sig_mb_sigma1_value, 0.0, 0.1, kTRUE );
    LauParameter* sig_mb_mean2 = sig_mb_mean1->createClone();
    LauParameter* sig_mb_sigma2 =
        new LauParameter( "sig_mb_sigma2", sig_mb_sigma2_value, 0.0, 0.2, kTRUE );
    LauParameter* sig_mb_frac = new LauParameter( "sig_mb_frac", sig_mb_frac_value, 0.0, 1.0, kTRUE );

    mbPars.clear();
    mbPars.push_back( sig_mb_mean1 );
    mbPars.push_back( sig_mb_sigma1 );
    LauAbsPdf* sigMBPdf1 = new LauGaussPdf( "mB", mbPars, mbMin, mbMax );

    mbPars.clear();
    mbPars.push_back( sig_mb_mean2 );
    mbPars.push_back( sig_mb_sigma2 );
    LauAbsPdf* sigMBPdf2 = new LauGaussPdf( "mB", mbPars, mbMin, mbMax );

    LauAbsPdf* sigMBPdf = new LauSumPdf( sigMBPdf1, sigMBPdf2, sig_mb_frac );
    fitModel->setSignalPdf( sigMBPdf );

    // Combinatoric background PDF is linear function
    LauParameter* comb_mb_slope = new LauParameter( "comb_mb_slope", -0.05, -1.0, 1.0, kTRUE );

    mbPars.clear();
    mbPars.push_back( comb_mb_slope );
    LauAbsPdf* combMBPdf = new LauLinearPdf( "mB", mbPars, mbMin, mbMax );

    fitModel->setBkgndPdf( bkgndNames[0], combMBPdf );

    // Partially reconstructed background PDF is a crystal ball function
    LauParameter* pr_mb_mean  = new LauParameter( "pr_mb_mean", 5.200, 5.1, 5.3, kTRUE );
    LauParameter* pr_mb_sigma = new LauParameter( "pr_mb_sigma", 0.050, 0.0, 0.2, kTRUE );
    LauParameter* pr_mb_alpha = new LauParameter( "pr_mb_alpha", 0.100, -5.0, 5.0, kTRUE );
    LauParameter* pr_mb_order = new LauParameter( "pr_mb_order", 4.000, 0.0, 5.0, kTRUE );

    mbPars.clear();
    mbPars.push_back( pr_mb_mean );
    mbPars.push_back( pr_mb_sigma );
    mbPars.push_back( pr_mb_alpha );
    mbPars.push_back( pr_mb_order );
    LauAbsPdf* prbgMBPdf = new LauCrystalBallPdf( "mB", mbPars, mbMin, mbMax );

    fitModel->setBkgndPdf( bkgndNames[1], prbgMBPdf );

    // Configure various fit options

    // Do not calculate asymmetric errors.
    fitModel->useAsymmFitErrors( kTRUE );

    // Randomise initial fit values for the signal isobar parameters
    fitModel->useRandomInitFitPars( kTRUE );

    // Switch on/off Poissonian smearing of total number of events
    fitModel->doPoissonSmearing( kTRUE );

    // Switch on/off Extended ML Fit option
    fitModel->doEMLFit( kTRUE );

    // Switch on/off two-stage fitting for CPV parameters
    fitModel->twoStageFit( kFALSE );

    // Generate toy from the fitted parameters
    //TString fitToyFileName("fitToyMC_NoDP_");
    //fitToyFileName += iFit;
    //fitToyFileName += ".root";
    //fitModel->compareFitData(100, fitToyFileName);

    // Write out per-event likelihoods and sWeights
    //TString splotFileName("splot_NoDP_");
    //splotFileName += iFit;
    //splotFileName += ".root";
    //fitModel->writeSPlotData(splotFileName, "splot", kFALSE);

    TString dataFile( "gen-NoDP.root" );
    TString treeName( "genResults" );
    TString rootFileName( "" );
    TString tableFileName( "" );
    if ( command == "fit" ) {
        rootFileName   = "fitNoDP_";
        rootFileName  += iFit;
        rootFileName  += "_expt_";
        rootFileName  += firstExpt;
        rootFileName  += "-";
        rootFileName  += firstExpt + nExpt - 1;
        rootFileName  += ".root";
        tableFileName  = "fitNoDPResults_";
        tableFileName += iFit;
    } else {
        rootFileName  = "dummy.root";
        tableFileName = "genNoDPResults";
    }

    // Execute the generation/fit
    fitModel->run( command, dataFile, treeName, rootFileName, tableFileName );

    return EXIT_SUCCESS;
}
