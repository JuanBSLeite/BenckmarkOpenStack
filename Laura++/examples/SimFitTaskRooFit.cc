
/*
Copyright 2017 University of Warwick

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

#include "LauRooFitTask.hh"

#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooRealVar.h"
#include "TString.h"

#include <cstdlib>
#include <iostream>
#include <vector>

void usage( std::ostream& out, const TString& progName )
{
    out << "Usage:\n";
    out << progName << " gen <category = DD or LL> [nExpt = 1] [firstExpt = 0]\n";
    out << "or\n";
    out << progName << " fit <category = DD or LL> <port> [nExpt = 1] [firstExpt = 0]" << std::endl;
}

int main( int argc, char** argv )
{
    // Process command-line arguments
    // Usage:
    // ./SimFitTaskRooFit gen <category = DD or LL> [nExpt = 1] [firstExpt = 0]
    // or
    // ./SimFitTaskRooFit fit <category = DD or LL> <port> [nExpt = 1] [firstExpt = 0]
    if ( argc < 3 ) {
        usage( std::cerr, argv[0] );
        return EXIT_FAILURE;
    }

    TString command = argv[1];
    command.ToLower();
    TString category = argv[2];
    category.ToUpper();
    if ( category != "DD" && category != "LL" ) {
        usage( std::cerr, argv[0] );
        return EXIT_FAILURE;
    }
    UInt_t port( 5000 );
    Int_t nExpt( 1 );
    Int_t firstExpt( 0 );
    if ( command == "gen" ) {
        if ( argc > 3 ) {
            nExpt = atoi( argv[3] );
            if ( argc > 4 ) {
                firstExpt = atoi( argv[4] );
            }
        }
    } else if ( command == "fit" ) {
        if ( argc < 4 ) {
            usage( std::cerr, argv[0] );
            return EXIT_FAILURE;
        }
        port = atoi( argv[3] );
        if ( argc > 4 ) {
            nExpt = atoi( argv[4] );
            if ( argc > 5 ) {
                firstExpt = atoi( argv[5] );
            }
        }
    } else {
        usage( std::cerr, argv[0] );
        return EXIT_FAILURE;
    }

    // Set the names of the files to read/write
    TString dataFile = TString::Format( "data-%s.root", category.Data() );
    TString treeName( "fitTree" );
    TString rootFileName( "" );
    TString tableFileName( "" );
    if ( command == "fit" ) {
        rootFileName   = "fit";
        rootFileName  += category;
        rootFileName  += "_expt_";
        rootFileName  += firstExpt;
        rootFileName  += "-";
        rootFileName  += ( firstExpt + nExpt - 1 );
        rootFileName  += ".root";
        tableFileName  = "fitResults";
    } else {
        rootFileName  = "dummy.root";
        tableFileName = "genResults";
    }

    // Create the RooFit model
    std::cout << "Building the model..." << std::endl;
    RooRealVar mB( "mB", "m_{B}", 5.150, 5.600, "GeV/c^{2}" );

    Double_t sig_mb_mean1_value  = 5.279;
    Double_t sig_mb_sigma1_value = 0.02;
    Double_t sig_mb_sigma2_value = 0.07;
    if ( category == "DD" ) {
        sig_mb_sigma1_value *= 1.05;
        sig_mb_sigma2_value *= 1.05;
    }
    Double_t sig_mb_frac_value = 0.90;

    RooRealVar sig_mb_mean1( "sig_mb_mean1", "", sig_mb_mean1_value, 5.200, 5.300 );
    RooRealVar sig_mb_sigma1( "sig_mb_sigma1", "", sig_mb_sigma1_value, 0.000, 0.100 );
    RooRealVar sig_mb_sigma2( "sig_mb_sigma2", "", sig_mb_sigma2_value, 0.000, 0.200 );
    RooRealVar sig_mb_frac( "sig_mb_frac", "", sig_mb_frac_value, 0.000, 1.000 );
    sig_mb_mean1.setConstant( kFALSE );
    sig_mb_sigma1.setConstant( kTRUE );
    sig_mb_sigma2.setConstant( kTRUE );
    sig_mb_frac.setConstant( kTRUE );

    RooGaussian gauss1( "gauss1", "", mB, sig_mb_mean1, sig_mb_sigma1 );
    RooGaussian gauss2( "gauss2", "", mB, sig_mb_mean1, sig_mb_sigma2 );
    RooAddPdf sigPDF( "sigPDF", "", gauss1, gauss2, sig_mb_frac );

    RooRealVar comb_mb_slope( "comb_mb_slope", "", -1.5, -10.0, 2.0 );
    comb_mb_slope.setConstant( kFALSE );

    RooExponential combPDF( "combPDF", "", mB, comb_mb_slope );

    Double_t nSigEvents = 5000.0;
    if ( category == "DD" ) {
        nSigEvents *= 1.5;
    }
    Double_t nCombBgEvents = 7000.0;
    if ( category == "DD" ) {
        nCombBgEvents *= 1.7;
    }
    Double_t nTotEvents = nSigEvents + nCombBgEvents;

    RooRealVar sigYield( "signalEvents" + category, "", nSigEvents, -2.0 * nTotEvents, 2.0 * nTotEvents );
    RooRealVar combYield( "CombinatorialEvents" + category,
                          "",
                          nCombBgEvents,
                          -2.0 * nTotEvents,
                          2.0 * nTotEvents );

    RooAddPdf model( "model", "", RooArgList( sigPDF, combPDF ), RooArgList( sigYield, combYield ) );

    LauRooFitTask fitModel( model, kTRUE, RooArgSet( mB ) );

    // Execute the generation/fit
    if ( command == "fit" ) {
        // Create the task instance
        fitModel.runTask( dataFile, treeName, rootFileName, tableFileName, "localhost", port );
    } else {
        fitModel.generate( dataFile, treeName, nExpt, firstExpt );
    }

    return EXIT_SUCCESS;
}
