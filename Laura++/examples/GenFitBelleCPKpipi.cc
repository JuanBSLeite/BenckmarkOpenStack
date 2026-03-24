
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

#include "LauBelleCPCoeffSet.hh"
#include "LauBkgndDPModel.hh"
#include "LauCPFitModel.hh"
#include "LauDaughters.hh"
#include "LauEffModel.hh"
#include "LauIsobarDynamics.hh"
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
    // ./GenFitBelleCPKpipi gen [nExpt = 1] [firstExpt = 0]
    // or
    // ./GenFitBelleCPKpipi fit <iFit> [nExpt = 1] [firstExpt = 0]
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

    // These define the DP => decay is B+ -> K+ pi+ pi- and it's charge conjugate
    // The DP is defined in terms of m13Sq and m23Sq
    LauDaughters* negDaughters = new LauDaughters( "B-", "K-", "pi-", "pi+", squareDP );
    LauDaughters* posDaughters = new LauDaughters( "B+", "K+", "pi+", "pi-", squareDP );

    // Create the isobar models
    LauAbsResonance* res( 0 );
    LauIsobarDynamics* negSigModel = new LauIsobarDynamics( negDaughters, 0 );
    // resPairAmpInt = 1 => resonance mass is m23.
    // resPairAmpInt = 2 => resonance mass is m13.
    res = negSigModel->addResonance( "K*0(892)", 2, LauAbsResonance::RelBW );
    res = negSigModel->addResonance( "K*0_0(1430)", 2, LauAbsResonance::LASS_BW );
    res->floatResonanceParameter( "a" );
    res->floatResonanceParameter( "r" );
    res = negSigModel->addResonance( "LASSNR0", 2, LauAbsResonance::LASS_NR );
    res = negSigModel->addResonance( "rho0(770)", 1, LauAbsResonance::RelBW );
    res = negSigModel->addResonance( "f_0(980)", 1, LauAbsResonance::Flatte );
    res = negSigModel->addResonance( "chi_c0", 1, LauAbsResonance::RelBW );
    res = negSigModel->addResonance( "NonReson", 1, LauAbsResonance::BelleNR );
    res->setResonanceParameter( "alpha", 0.50 );
    res->floatResonanceParameter( "alpha" );

    LauIsobarDynamics* posSigModel = new LauIsobarDynamics( posDaughters, 0 );
    res = posSigModel->addResonance( "K*0(892)", 2, LauAbsResonance::RelBW );
    res = posSigModel->addResonance( "K*0_0(1430)", 2, LauAbsResonance::LASS_BW );
    res->floatResonanceParameter( "a" );
    res->floatResonanceParameter( "r" );
    res = posSigModel->addResonance( "LASSNR0", 2, LauAbsResonance::LASS_NR );
    res = posSigModel->addResonance( "rho0(770)", 1, LauAbsResonance::RelBW );
    res = posSigModel->addResonance( "f_0(980)", 1, LauAbsResonance::Flatte );
    res = posSigModel->addResonance( "chi_c0", 1, LauAbsResonance::RelBW );
    res = posSigModel->addResonance( "NonReson", 1, LauAbsResonance::BelleNR );
    res->setResonanceParameter( "alpha", 0.50 );
    res->floatResonanceParameter( "alpha" );

    // Set the file names for the integrals information (can be useful for debugging)
    negSigModel->setIntFileName( "integ_neg.dat" );
    posSigModel->setIntFileName( "integ_pos.dat" );

    // Set the maximum signal DP ASq value
    // If you do not provide a value, one will be determined automatically,
    // which should be close to the true maximum but is not guaranteed to
    // be optimal.
    // Any value, whether manually provided or automatically determined,
    // will be automatically adjusted to avoid bias or extreme inefficiency
    // but it is best to set this by hand once you've found the right value
    // through some trial and error.
    const Double_t aSqMaxValue { 1.62 };
    negSigModel->setASqMaxValue( aSqMaxValue );
    posSigModel->setASqMaxValue( aSqMaxValue );

    // Create the fit model, giving it both isobar models
    LauCPFitModel* fitModel = new LauCPFitModel( negSigModel, posSigModel );

    // Create the complex coefficients for the isobar model
    // Here we're using the "Belle" form:
    // The amplitude has the form a * exp(i*delta) * ( 1 +/- b * exp(i*phi) ) where
    // a is a CP conserving magnitude,
    // b is a CP violating magnitude,
    // delta is the strong phase
    // and phi is the weak phase.
    // The last two arguments indicate that the CP-violating parameters
    // should only be floated in the second-stage of the two-stage fit
    LauAbsCoeffSet* a0CoeffSet = new LauBelleCPCoeffSet( "K*0(892)",
                                                         1.00,
                                                         0.00,
                                                         0.00,
                                                         0.00,
                                                         kTRUE,
                                                         kTRUE,
                                                         kFALSE,
                                                         kTRUE,
                                                         kTRUE,
                                                         kTRUE );

    // Here we show how to use the coefficient cloning mechanism to clone
    // only the CP-violating parameters.  This means that the K*0_0(1430)
    // component will have identical CPV parameters to the K*0(892).
    // But the CP-conserving parameters will be independent and as such we
    // need to set their values and, in this case because we're cloning
    // from the reference amplitude, also set them to float.
    //LauAbsCoeffSet* a1CoeffSet = new LauBelleCPCoeffSet("K*0_0(1430)", 2.00,  3.00, 0.00, 0.00, kFALSE, kFALSE, kFALSE, kFALSE, kTRUE, kTRUE);
    LauAbsCoeffSet* a1CoeffSet = a0CoeffSet->createClone( "K*0_0(1430)", LauAbsCoeffSet::TieCPPars );
    a1CoeffSet->setParameterValue( "A", 2.0, kTRUE );
    a1CoeffSet->setParameterValue( "Delta", 3.0, kTRUE );
    a1CoeffSet->floatParameter( "A" );
    a1CoeffSet->floatParameter( "Delta" );
    LauAbsCoeffSet* a2CoeffSet = a1CoeffSet->createClone( "LASSNR0", LauAbsCoeffSet::TieCPPars );

    LauAbsCoeffSet* a3CoeffSet = new LauBelleCPCoeffSet( "rho0(770)",
                                                         0.66,
                                                         1.00,
                                                         0.00,
                                                         0.00,
                                                         kFALSE,
                                                         kFALSE,
                                                         kFALSE,
                                                         kFALSE,
                                                         kTRUE,
                                                         kTRUE );
    LauAbsCoeffSet* a4CoeffSet = new LauBelleCPCoeffSet( "f_0(980)",
                                                         1.00,
                                                         -1.00,
                                                         0.00,
                                                         0.00,
                                                         kFALSE,
                                                         kFALSE,
                                                         kFALSE,
                                                         kFALSE,
                                                         kTRUE,
                                                         kTRUE );
    LauAbsCoeffSet* a5CoeffSet = new LauBelleCPCoeffSet( "chi_c0",
                                                         0.33,
                                                         0.50,
                                                         0.00,
                                                         0.00,
                                                         kFALSE,
                                                         kFALSE,
                                                         kFALSE,
                                                         kFALSE,
                                                         kTRUE,
                                                         kTRUE );
    LauAbsCoeffSet* a6CoeffSet = new LauBelleCPCoeffSet( "NonReson",
                                                         0.50,
                                                         1.50,
                                                         0.00,
                                                         0.00,
                                                         kFALSE,
                                                         kFALSE,
                                                         kFALSE,
                                                         kFALSE,
                                                         kTRUE,
                                                         kTRUE );

    std::vector<LauAbsCoeffSet*> coeffset;
    coeffset.push_back( a0CoeffSet );
    coeffset.push_back( a1CoeffSet );
    coeffset.push_back( a2CoeffSet );
    coeffset.push_back( a3CoeffSet );
    coeffset.push_back( a4CoeffSet );
    coeffset.push_back( a5CoeffSet );
    coeffset.push_back( a6CoeffSet );
    for ( std::vector<LauAbsCoeffSet*>::iterator iter = coeffset.begin(); iter != coeffset.end();
          ++iter ) {
        fitModel->setAmpCoeffSet( *iter );
    }

    // Set the signal yield and define whether it is fixed or floated
    const Double_t nSigEvents  = 5000.0;
    Bool_t fixNSigEvents       = kFALSE;
    LauParameter* signalEvents = new LauParameter( "signalEvents",
                                                   nSigEvents,
                                                   -1.0 * nSigEvents,
                                                   2.0 * nSigEvents,
                                                   fixNSigEvents );
    fitModel->setNSigEvents( signalEvents );

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

    // Switch on/off two-stage fit
    fitModel->twoStageFit( kTRUE );

    // Generate toy from the fitted parameters
    //TString fitToyFileName("fitToyMC_BelleCPKpipi_");
    //fitToyFileName += iFit;
    //fitToyFileName += ".root";
    //fitModel->compareFitData(100, fitToyFileName);

    // Write out per-event likelihoods and sWeights
    //TString splotFileName("splot_BelleCPKpipi_");
    //splotFileName += iFit;
    //splotFileName += ".root";
    //fitModel->writeSPlotData(splotFileName, "splot", kFALSE);

    // Set the names of the files to read/write
    TString dataFile( "gen-BelleCPKpipi.root" );
    TString treeName( "genResults" );
    TString rootFileName( "" );
    TString tableFileName( "" );
    if ( command == "fit" ) {
        rootFileName   = "fitBelleCPKpipi_";
        rootFileName  += iFit;
        rootFileName  += "_expt_";
        rootFileName  += firstExpt;
        rootFileName  += "-";
        rootFileName  += ( firstExpt + nExpt - 1 );
        rootFileName  += ".root";
        tableFileName  = "fitBelleCPKpipiResults_";
        tableFileName += iFit;
    } else {
        rootFileName  = "dummy.root";
        tableFileName = "genBelleCPKpipiResults";
    }

    // Execute the generation/fit
    fitModel->run( command, dataFile, treeName, rootFileName, tableFileName );

    return EXIT_SUCCESS;
}
