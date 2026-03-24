
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

#include "Lau1DHistPdf.hh"
#include "LauArgusPdf.hh"
#include "LauBkgndDPModel.hh"
#include "LauCPFitModel.hh"
#include "LauCartesianCPCoeffSet.hh"
#include "LauDaughters.hh"
#include "LauEffModel.hh"
#include "LauGaussPdf.hh"
#include "LauIsobarDynamics.hh"
#include "LauLinearPdf.hh"
#include "LauResonanceMaker.hh"
#include "LauSumPdf.hh"
#include "LauVetoes.hh"

#include "TFile.h"
#include "TH2.h"
#include "TString.h"

#include <cstdlib>
#include <iostream>
#include <vector>

/*
 * Histogram ROOT files for this example are available at the following URL:
 * https://laura.hepforge.org/example-histograms/GenFitKpipi-histos.tar.gz
 *\  Please copy the "histos" directory into the working directory.
 */

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
    // ./GenFitKpipi gen [nExpt = 1] [firstExpt = 0]
    // or
    // ./GenFitKpipi fit <iFit> [nExpt = 1] [firstExpt = 0]
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
    const Bool_t squareDP = kTRUE;

    // Set this to kFALSE if you want to remove the DP from the fit
    const Bool_t doDP = kTRUE;

    // Set this to kTRUE if you want to fix the CPV parameters
    const Bool_t fixCP         = kFALSE;
    const Bool_t doTwoStageFit = ! fixCP;

    // General histogram booleans
    const Bool_t useInterpolation = kTRUE;
    const Bool_t fluctuateBins    = kFALSE;
    const Bool_t useUpperHalfOnly = kFALSE;

    // Signal and continuum yields - from Phys.Rev.D78:012004,2008
    const Double_t nSigEvents  = 4585.0;
    const Bool_t fixNSigEvents = kFALSE;
    const Double_t nBgEvents   = 6830.0;
    const Bool_t fixNBgEvents  = kFALSE;

    // Signal and continuum asymmetries
    // NB the signal asymmetry value here is only used if the DP is NOT
    // in the fit, otherwise the asymmetry is included within the
    // isobar model.
    const Double_t sigAsym  = 0.028;
    const Bool_t fixSigAsym = kFALSE;
    const Double_t bgAsym   = -0.028;
    const Bool_t fixBgAsym  = kFALSE;

    // Define the DP (for both B+ and B- candidates)
    LauDaughters* negDaughters = new LauDaughters( "B-", "K-", "pi-", "pi+", squareDP );
    LauDaughters* posDaughters = new LauDaughters( "B+", "K+", "pi+", "pi-", squareDP );

    // Add some vetoes
    LauVetoes* vetoes        = new LauVetoes();
    const Double_t D_KpiMin  = 1.756;
    const Double_t D_KpiMax  = 1.931;
    const Double_t D_pipiMin = 1.660;
    const Double_t D_pipiMax = 1.800;
    const Double_t JpsiMin   = 3.019;
    const Double_t JpsiMax   = 3.179;
    const Double_t psi2SMin  = 3.627;
    const Double_t psi2SMax  = 3.747;
    vetoes->addMassVeto( 2, D_KpiMin, D_KpiMax );      // D0 veto, m13
    vetoes->addMassVeto( 1, D_pipiMin, D_pipiMax );    // D0 veto, m23
    vetoes->addMassVeto( 1, JpsiMin, JpsiMax );        // J/psi veto, m23
    vetoes->addMassVeto( 1, psi2SMin, psi2SMax );      // psi(2S) veto, m23

    // Create the efficiency models
    TFile* effHistoFile = TFile::Open( "histos/eff/Kpipi_50x50.root", "read" );

    LauEffModel* negEffModel = new LauEffModel( negDaughters, vetoes );
    TString negEffHistoName( "eff_all_neg" );
    TH2* negEffHisto = dynamic_cast<TH2*>( effHistoFile->Get( negEffHistoName ) );
    negEffModel->setEffHisto( negEffHisto,
                              useInterpolation,
                              fluctuateBins,
                              0.0,
                              0.0,
                              useUpperHalfOnly,
                              squareDP );

    LauEffModel* posEffModel = new LauEffModel( posDaughters, vetoes );
    TString posEffHistoName( "eff_all_pos" );
    TH2* posEffHisto = dynamic_cast<TH2*>( effHistoFile->Get( posEffHistoName ) );
    posEffModel->setEffHisto( posEffHisto,
                              useInterpolation,
                              fluctuateBins,
                              0.0,
                              0.0,
                              useUpperHalfOnly,
                              squareDP );

    // Create the signal dynamics
    LauResonanceMaker& resMaker = LauResonanceMaker::get();
    resMaker.setDefaultBWRadius( LauBlattWeisskopfFactor::Parent, 4.0 );
    resMaker.setDefaultBWRadius( LauBlattWeisskopfFactor::Kstar, 4.0 );
    resMaker.setDefaultBWRadius( LauBlattWeisskopfFactor::Light, 4.0 );
    resMaker.setDefaultBWRadius( LauBlattWeisskopfFactor::Charmonium, 4.0 );

    LauIsobarDynamics* negSigModel = new LauIsobarDynamics( negDaughters, negEffModel );
    negSigModel->setIntFileName( "integ_neg.dat" );
    negSigModel->setASqMaxValue( 0.45 );

    LauAbsResonance* res( 0 );
    res = negSigModel->addResonance( "K*0(892)", 2, LauAbsResonance::RelBW );
    res = negSigModel->addResonance( "K*0_0(1430)", 2, LauAbsResonance::LASS );
    res = negSigModel->addResonance( "rho0(770)", 1, LauAbsResonance::RelBW );
    res = negSigModel->addResonance( "f_0(980)", 1, LauAbsResonance::Flatte );
    res->changeResonance( 0.965, -1.0, -1 );
    res = negSigModel->addResonance( "chi_c0", 1, LauAbsResonance::RelBW );
    res = negSigModel->addResonance( "NonReson", 0, LauAbsResonance::FlatNR );

    res = negSigModel->addResonance( "K*0_2(1430)", 2, LauAbsResonance::RelBW );
    res = negSigModel->addResonance( "omega(782)", 1, LauAbsResonance::RelBW );
    res = negSigModel->addResonance( "f_2(1270)", 1, LauAbsResonance::RelBW );
    res = negSigModel->addResonance( "f'_0(1300)", 1, LauAbsResonance::RelBW );
    res->changeResonance( 1.479, 0.080, -1 );

    LauIsobarDynamics* posSigModel = new LauIsobarDynamics( posDaughters, posEffModel );
    posSigModel->setIntFileName( "integ_pos.dat" );
    posSigModel->setASqMaxValue( 0.45 );

    res = posSigModel->addResonance( "K*0(892)", 2, LauAbsResonance::RelBW );
    res = posSigModel->addResonance( "K*0_0(1430)", 2, LauAbsResonance::LASS );
    res = posSigModel->addResonance( "rho0(770)", 1, LauAbsResonance::RelBW );
    res = posSigModel->addResonance( "f_0(980)", 1, LauAbsResonance::Flatte );
    res->changeResonance( 0.965, -1.0, -1 );
    res = posSigModel->addResonance( "chi_c0", 1, LauAbsResonance::RelBW );
    res = posSigModel->addResonance( "NonReson", 0, LauAbsResonance::FlatNR );

    res = posSigModel->addResonance( "K*0_2(1430)", 2, LauAbsResonance::RelBW );
    res = posSigModel->addResonance( "omega(782)", 1, LauAbsResonance::RelBW );
    res = posSigModel->addResonance( "f_2(1270)", 1, LauAbsResonance::RelBW );
    res = posSigModel->addResonance( "f'_0(1300)", 1, LauAbsResonance::RelBW );
    res->changeResonance( 1.479, 0.080, -1 );

    // Create the fit model, passing it the two signal dynamics models
    LauCPFitModel* fitModel = new LauCPFitModel( negSigModel, posSigModel );
    fitModel->useDP( doDP );

    // Set isobar coefficients - taken from Phys.Rev.D78:012004,2008
    std::vector<LauAbsCoeffSet*> coeffset;
    coeffset.push_back( new LauCartesianCPCoeffSet( "K*0(892)",
                                                    1.000,
                                                    0.000,
                                                    -0.017,
                                                    -0.238,
                                                    kTRUE,
                                                    kTRUE,
                                                    fixCP,
                                                    fixCP,
                                                    doTwoStageFit,
                                                    doTwoStageFit ) );
    coeffset.push_back( new LauCartesianCPCoeffSet( "K*0_0(1430)",
                                                    1.718,
                                                    -0.727,
                                                    -0.154,
                                                    -0.285,
                                                    kFALSE,
                                                    kFALSE,
                                                    fixCP,
                                                    fixCP,
                                                    doTwoStageFit,
                                                    doTwoStageFit ) );
    coeffset.push_back( new LauCartesianCPCoeffSet( "rho0(770)",
                                                    0.683,
                                                    -0.025,
                                                    -0.160,
                                                    0.169,
                                                    kFALSE,
                                                    kFALSE,
                                                    fixCP,
                                                    fixCP,
                                                    doTwoStageFit,
                                                    doTwoStageFit ) );
    coeffset.push_back( new LauCartesianCPCoeffSet( "f_0(980)",
                                                    -0.220,
                                                    1.203,
                                                    -0.109,
                                                    0.047,
                                                    kFALSE,
                                                    kFALSE,
                                                    fixCP,
                                                    fixCP,
                                                    doTwoStageFit,
                                                    doTwoStageFit ) );
    coeffset.push_back( new LauCartesianCPCoeffSet( "chi_c0",
                                                    -0.263,
                                                    0.180,
                                                    -0.033,
                                                    -0.007,
                                                    kFALSE,
                                                    kFALSE,
                                                    fixCP,
                                                    fixCP,
                                                    doTwoStageFit,
                                                    doTwoStageFit ) );
    coeffset.push_back(
        new LauCartesianCPCoeffSet( "NonReson", -0.594, 0.068, 0.000, 0.000, kFALSE, kFALSE, kTRUE, kTRUE ) );
    coeffset.push_back( new LauCartesianCPCoeffSet( "K*0_2(1430)",
                                                    -0.301,
                                                    0.424,
                                                    0.032,
                                                    0.007,
                                                    kFALSE,
                                                    kFALSE,
                                                    fixCP,
                                                    fixCP,
                                                    doTwoStageFit,
                                                    doTwoStageFit ) );
    coeffset.push_back(
        new LauCartesianCPCoeffSet( "omega(782)", -0.058, 0.100, 0.000, 0.000, kFALSE, kFALSE, kTRUE, kTRUE ) );
    coeffset.push_back( new LauCartesianCPCoeffSet( "f_2(1270)",
                                                    -0.193,
                                                    0.110,
                                                    -0.089,
                                                    0.125,
                                                    kFALSE,
                                                    kFALSE,
                                                    fixCP,
                                                    fixCP,
                                                    doTwoStageFit,
                                                    doTwoStageFit ) );
    coeffset.push_back( new LauCartesianCPCoeffSet( "f'_0(1300)",
                                                    -0.290,
                                                    -0.136,
                                                    0.024,
                                                    0.056,
                                                    kFALSE,
                                                    kFALSE,
                                                    fixCP,
                                                    fixCP,
                                                    doTwoStageFit,
                                                    doTwoStageFit ) );
    for ( std::vector<LauAbsCoeffSet*>::iterator iter = coeffset.begin(); iter != coeffset.end();
          ++iter ) {
        fitModel->setAmpCoeffSet( *iter );
    }

    // Set the number of signal events and the number of experiments
    LauParameter* signalEvents = new LauParameter( "signalEvents",
                                                   nSigEvents,
                                                   -1.0 * nSigEvents,
                                                   2.0 * nSigEvents,
                                                   fixNSigEvents );
    LauParameter* signalAsym   = new LauParameter( "signalAsym", sigAsym, -1.0, 1.0, fixSigAsym );

    fitModel->setNSigEvents( signalEvents, signalAsym );
    fitModel->setNExpts( nExpt, firstExpt, isToy );

    // Set up the background configuration
    const Int_t nBBs( 1 );
    TString vetoCorrHistName( "" );
    TString qqDPName = "";
    std::vector<TString> BBNames( nBBs );
    std::vector<TString> BBDPNames( nBBs );
    std::vector<TString> BBMESNames( nBBs );
    std::vector<TString> BBDENames( nBBs );
    if ( squareDP ) {
        vetoCorrHistName = "correctionSqDPHisto";
        qqDPName         = "onres_sqDP_histo";
        for ( Int_t i( 0 ); i < nBBs; ++i ) {
            TString name( "BB" );
            name += i + 1;
            TString dpName( name );
            TString mesName( name );
            TString deName( name );
            dpName        += "_sqDP_histo";
            mesName       += "_mES_histo";
            deName        += "_deltaEsig_histo";
            BBNames[i]     = name;
            BBDPNames[i]   = dpName;
            BBMESNames[i]  = mesName;
            BBDENames[i]   = deName;
        }
    } else {
        vetoCorrHistName = "correctionDPHisto";
        qqDPName         = "onres_DP_histo";
        for ( Int_t i( 0 ); i < nBBs; ++i ) {
            TString name( "BB" );
            name += i + 1;
            TString dpName( name );
            TString mesName( name );
            TString deName( name );
            dpName        += "_DP_histo";
            mesName       += "_mES_histo";
            deName        += "_deltaEsig_histo";
            BBNames[i]     = name;
            BBDPNames[i]   = dpName;
            BBMESNames[i]  = mesName;
            BBDENames[i]   = deName;
        }
    }

    std::vector<TString> bkgndNames( nBBs + 1 );
    bkgndNames[0] = "qqbar";
    for ( Int_t i( 0 ); i < nBBs; ++i ) {
        bkgndNames[i + 1] = BBNames[i];
    }
    fitModel->setBkgndClassNames( bkgndNames );

    // Get the histogram for correcting the other histograms for the
    // effect of the vetoes on bins on a veto boundary
    TString vetoCorrFileName( "histos/veto-correction/vetoCorrectionHisto_50.root" );
    TFile* vetoCorrFile = TFile::Open( vetoCorrFileName, "read" );
    TH2* vetoCorrHist   = dynamic_cast<TH2*>( vetoCorrFile->Get( vetoCorrHistName ) );

    // Load the qqbar DP histograms
    TString negqqFileName( "histos/qqbar/qqbarBackgroundHistos_50_neg.root" );
    TFile* negqqFile = TFile::Open( negqqFileName, "read" );
    TH2* negqqDP     = dynamic_cast<TH2*>( negqqFile->Get( qqDPName ) );
    negqqDP->Divide( vetoCorrHist );
    LauBkgndDPModel* negqqbarModel = new LauBkgndDPModel( negDaughters, vetoes );
    negqqbarModel->setBkgndHisto( negqqDP, useInterpolation, fluctuateBins, useUpperHalfOnly, squareDP );

    TString posqqFileName( "histos/qqbar/qqbarBackgroundHistos_50_pos.root" );
    TFile* posqqFile = TFile::Open( posqqFileName, "read" );
    TH2* posqqDP     = dynamic_cast<TH2*>( posqqFile->Get( qqDPName ) );
    posqqDP->Divide( vetoCorrHist );
    LauBkgndDPModel* posqqbarModel = new LauBkgndDPModel( posDaughters, vetoes );
    posqqbarModel->setBkgndHisto( posqqDP, useInterpolation, fluctuateBins, useUpperHalfOnly, squareDP );

    LauParameter* backgroundEvents =
        new LauParameter( "qqbar", nBgEvents, -1.0 * nBgEvents, 2.0 * nBgEvents, fixNBgEvents );
    LauParameter* backgroundAsym = new LauParameter( "qqbar", bgAsym, -1.0, 1.0, fixBgAsym );

    fitModel->setNBkgndEvents( backgroundEvents, backgroundAsym );
    fitModel->setBkgndDPModels( "qqbar", negqqbarModel, posqqbarModel );

    // Load the BBbar DP histograms
    TString negBBFileName( "histos/BBbar/BbackgroundHistos_50_neg_categories.root" );
    TFile* negBBFile = TFile::Open( negBBFileName, "read" );
    TString posBBFileName( "histos/BBbar/BbackgroundHistos_50_pos_categories.root" );
    TFile* posBBFile = TFile::Open( posBBFileName, "read" );

    std::vector<TH2*> negBBDP( nBBs );
    std::vector<TH2*> posBBDP( nBBs );
    std::vector<LauParameter*> nBBEvents( nBBs );
    std::vector<LauParameter*> BBAsym( nBBs );
    std::vector<LauBkgndDPModel*> negBBModels( nBBs );
    std::vector<LauBkgndDPModel*> posBBModels( nBBs );

    for ( Int_t i( 0 ); i < nBBs; ++i ) {
        negBBFile->cd();
        negBBDP[i] = dynamic_cast<TH2*>( negBBFile->Get( BBDPNames[i] ) );
        if ( negBBDP[i] == 0 ) {
            std::cerr << "Problem finding histogram \"" << BBDPNames[i] << "\" in B- file."
                      << std::endl;
            return EXIT_FAILURE;
        }
        Double_t negBBEvents = negBBDP[i]->Integral();
        negBBDP[i]->Divide( vetoCorrHist );
        negBBModels[i] = new LauBkgndDPModel( negDaughters, vetoes );
        negBBModels[i]->setBkgndHisto( negBBDP[i],
                                       useInterpolation,
                                       fluctuateBins,
                                       useUpperHalfOnly,
                                       squareDP );

        posBBFile->cd();
        posBBDP[i] = dynamic_cast<TH2*>( posBBFile->Get( BBDPNames[i] ) );
        if ( posBBDP[i] == 0 ) {
            std::cerr << "Problem finding histogram \"" << BBDPNames[i] << "\" in B+ file."
                      << std::endl;
            return EXIT_FAILURE;
        }
        Double_t posBBEvents = posBBDP[i]->Integral();
        posBBDP[i]->Divide( vetoCorrHist );
        posBBModels[i] = new LauBkgndDPModel( posDaughters, vetoes );
        posBBModels[i]->setBkgndHisto( posBBDP[i],
                                       useInterpolation,
                                       fluctuateBins,
                                       useUpperHalfOnly,
                                       squareDP );

        Double_t totalBBEvents = negBBEvents + posBBEvents;
        Double_t asym          = ( negBBEvents - posBBEvents ) / totalBBEvents;

        nBBEvents[i] = new LauParameter( BBNames[i],
                                         totalBBEvents,
                                         -2.0 * totalBBEvents,
                                         2.0 * totalBBEvents,
                                         kTRUE );
        BBAsym[i]    = new LauParameter( BBNames[i], asym, -1.0, 1.0, kTRUE );

        fitModel->setBkgndDPModels( BBNames[i], negBBModels[i], posBBModels[i] );
        fitModel->setNBkgndEvents( nBBEvents[i], BBAsym[i] );
    }

    // Signal mES and DeltaEsig PDF parameter values, shifts and scale factors
    Double_t mes_mean1_value  = 5.279443;
    Double_t mes_mean2_value  = 5.275484;
    Double_t mes_sigma1_value = 0.0023974;
    Double_t mes_frac_value   = 0.9094;
    Double_t de_mean1_value   = -0.00955;
    Double_t de_mean2_value   = -0.632;
    Double_t de_sigma1_value  = 1.0771;
    Double_t de_sigma2_value  = 2.258;
    Double_t de_frac_value    = 0.8069;

    // mES PDFs
    // PDFs for B+ and B- are set to be identical here but could be different in general
    Double_t mesMin = 5.272;
    Double_t mesMax = 5.286;

    // Signal PDF is a double Gaussian with the widths constrained to be the same
    LauParameter* mes_mean1_neg =
        new LauParameter( "mes_mean1_neg", mes_mean1_value, 5.272, 5.286, kTRUE );
    LauParameter* mes_sigma1_neg =
        new LauParameter( "mes_sigma1_neg", mes_sigma1_value, 0.000, 0.010, kTRUE );
    LauParameter* mes_mean2_neg =
        new LauParameter( "mes_mean2_neg", mes_mean2_value, 5.272, 5.286, kTRUE );
    LauParameter* mes_sigma2_neg = mes_sigma1_neg->createClone();
    LauParameter* mes_frac_neg =
        new LauParameter( "mes_frac_neg", mes_frac_value, 0.000, 1.000, kTRUE );
    std::vector<LauAbsRValue*> mesPars;
    mesPars.reserve( 2 );
    mesPars.push_back( mes_mean1_neg );
    mesPars.push_back( mes_sigma1_neg );
    LauAbsPdf* sigNegMESPdf1 = new LauGaussPdf( "mES", mesPars, mesMin, mesMax );
    mesPars.clear();
    mesPars.push_back( mes_mean2_neg );
    mesPars.push_back( mes_sigma2_neg );
    LauAbsPdf* sigNegMESPdf2 = new LauGaussPdf( "mES", mesPars, mesMin, mesMax );
    LauAbsPdf* sigNegMESPdf  = new LauSumPdf( sigNegMESPdf1, sigNegMESPdf2, mes_frac_neg );
    mesPars.clear();
    LauParameter* mes_mean1_pos  = mes_mean1_neg->createClone();
    LauParameter* mes_sigma1_pos = mes_sigma1_neg->createClone();
    LauParameter* mes_mean2_pos  = mes_mean2_neg->createClone();
    LauParameter* mes_sigma2_pos = mes_sigma1_pos->createClone();
    LauParameter* mes_frac_pos   = mes_frac_neg->createClone();
    mesPars.push_back( mes_mean1_pos );
    mesPars.push_back( mes_sigma1_pos );
    LauAbsPdf* sigPosMESPdf1 = new LauGaussPdf( "mES", mesPars, mesMin, mesMax );
    mesPars.clear();
    mesPars.push_back( mes_mean2_pos );
    mesPars.push_back( mes_sigma2_pos );
    LauAbsPdf* sigPosMESPdf2 = new LauGaussPdf( "mES", mesPars, mesMin, mesMax );
    LauAbsPdf* sigPosMESPdf  = new LauSumPdf( sigPosMESPdf1, sigPosMESPdf2, mes_frac_pos );

    fitModel->setSignalPdfs( sigNegMESPdf, sigPosMESPdf );

    // Continuum PDF is an ARGUS function
    mesPars.clear();
    LauParameter* mes_qqbar_m0_neg = new LauParameter( "mes_qqbar_m0", 5.289 );
    LauParameter* mes_qqbar_xi_neg = new LauParameter( "mes_qqbar_xi_neg", 13.3, 0.0, 50.0, kTRUE );
    mesPars.push_back( mes_qqbar_m0_neg );
    mesPars.push_back( mes_qqbar_xi_neg );
    LauAbsPdf* qqbarNegMESPdf = new LauArgusPdf( "mES", mesPars, mesMin, mesMax );

    mesPars.clear();
    LauParameter* mes_qqbar_m0_pos = mes_qqbar_m0_neg->createClone();
    LauParameter* mes_qqbar_xi_pos = mes_qqbar_xi_neg->createClone();
    mesPars.push_back( mes_qqbar_m0_pos );
    mesPars.push_back( mes_qqbar_xi_pos );
    LauAbsPdf* qqbarPosMESPdf = new LauArgusPdf( "mES", mesPars, mesMin, mesMax );

    fitModel->setBkgndPdfs( "qqbar", qqbarNegMESPdf, qqbarPosMESPdf );

    // BB background PDFs are 1D histograms
    std::vector<TH1F*> negBBMES( nBBs );
    std::vector<TH1F*> posBBMES( nBBs );
    std::vector<LauAbsPdf*> negBBMESPdfs( nBBs );
    std::vector<LauAbsPdf*> posBBMESPdfs( nBBs );
    for ( Int_t i( 0 ); i < nBBs; ++i ) {
        negBBFile->cd();
        negBBMES[i] = dynamic_cast<TH1F*>( negBBFile->Get( BBMESNames[i] ) );
        if ( negBBMES[i] == 0 ) {
            std::cerr << "Problem finding histogram \"" << BBMESNames[i] << "\" in B- file."
                      << std::endl;
            return EXIT_FAILURE;
        }
        negBBMESPdfs[i] =
            new Lau1DHistPdf( "mES", negBBMES[i], mesMin, mesMax, useInterpolation, fluctuateBins );

        posBBFile->cd();
        posBBMES[i] = dynamic_cast<TH1F*>( posBBFile->Get( BBMESNames[i] ) );
        if ( posBBMES[i] == 0 ) {
            std::cerr << "Problem finding histogram \"" << BBMESNames[i] << "\" in B+ file."
                      << std::endl;
            return EXIT_FAILURE;
        }
        posBBMESPdfs[i] =
            new Lau1DHistPdf( "mES", posBBMES[i], mesMin, mesMax, useInterpolation, fluctuateBins );

        fitModel->setBkgndPdfs( BBNames[i], negBBMESPdfs[i], posBBMESPdfs[i] );
    }

    // deltaE significance PDFs
    // PDFs for B+ and B- are set to be identical here but could be different in general
    Double_t deMin = -4.0;
    Double_t deMax = 4.0;

    // Signal PDF is a double Gaussian
    LauParameter* de_sig_mean1_neg =
        new LauParameter( "de_sig_mean1_neg", de_mean1_value, -4.0, 4.0, kTRUE );
    LauParameter* de_sig_sigma1_neg =
        new LauParameter( "de_sig_sigma1_neg", de_sigma1_value, 0.0, 4.0, kTRUE );
    LauParameter* de_sig_mean2_neg =
        new LauParameter( "de_sig_mean2_neg", de_mean2_value, -4.0, 4.0, kTRUE );
    LauParameter* de_sig_sigma2_neg =
        new LauParameter( "de_sig_sigma2_neg", de_sigma2_value, 0.0, 4.0, kTRUE );
    LauParameter* de_sig_frac_neg =
        new LauParameter( "de_sig_frac_neg", de_frac_value, 0.0, 1.0, kTRUE );
    std::vector<LauAbsRValue*> dePars;
    dePars.reserve( 2 );
    dePars.push_back( de_sig_mean1_neg );
    dePars.push_back( de_sig_sigma1_neg );
    LauAbsPdf* sigNegDE1Pdf = new LauGaussPdf( "deltaEsig", dePars, deMin, deMax );
    dePars.clear();
    dePars.push_back( de_sig_mean2_neg );
    dePars.push_back( de_sig_sigma2_neg );
    LauAbsPdf* sigNegDE2Pdf = new LauGaussPdf( "deltaEsig", dePars, deMin, deMax );
    LauAbsPdf* sigNegDEPdf  = new LauSumPdf( sigNegDE1Pdf, sigNegDE2Pdf, de_sig_frac_neg );

    dePars.clear();
    LauParameter* de_sig_mean1_pos  = de_sig_mean1_neg->createClone();
    LauParameter* de_sig_sigma1_pos = de_sig_sigma1_neg->createClone();
    LauParameter* de_sig_mean2_pos  = de_sig_mean2_neg->createClone();
    LauParameter* de_sig_sigma2_pos = de_sig_sigma2_neg->createClone();
    LauParameter* de_sig_frac_pos   = de_sig_frac_neg->createClone();
    dePars.push_back( de_sig_mean1_pos );
    dePars.push_back( de_sig_sigma1_pos );
    LauAbsPdf* sigPosDE1Pdf = new LauGaussPdf( "deltaEsig", dePars, deMin, deMax );
    dePars.clear();
    dePars.push_back( de_sig_mean2_pos );
    dePars.push_back( de_sig_sigma2_pos );
    LauAbsPdf* sigPosDE2Pdf = new LauGaussPdf( "deltaEsig", dePars, deMin, deMax );
    LauAbsPdf* sigPosDEPdf  = new LauSumPdf( sigPosDE1Pdf, sigPosDE2Pdf, de_sig_frac_pos );

    fitModel->setSignalPdfs( sigNegDEPdf, sigPosDEPdf );

    // Continuum PDF is a 1st order polynomial
    dePars.clear();
    LauParameter* de_qqbar_slope_neg =
        new LauParameter( "de_qqbar_slope_neg", -0.00241, -0.01, 0.01, kFALSE );
    dePars.push_back( de_qqbar_slope_neg );
    LauAbsPdf* qqbarNegDEPdf = new LauLinearPdf( "deltaEsig", dePars, deMin, deMax );

    dePars.clear();
    LauParameter* de_qqbar_slope_pos = de_qqbar_slope_neg->createClone();
    dePars.push_back( de_qqbar_slope_pos );
    LauAbsPdf* qqbarPosDEPdf = new LauLinearPdf( "deltaEsig", dePars, deMin, deMax );

    fitModel->setBkgndPdfs( "qqbar", qqbarNegDEPdf, qqbarPosDEPdf );

    // BB background PDFs are 1D histograms
    std::vector<TH1F*> negBBDE( nBBs );
    std::vector<TH1F*> posBBDE( nBBs );
    std::vector<LauAbsPdf*> negBBDEPdfs( nBBs );
    std::vector<LauAbsPdf*> posBBDEPdfs( nBBs );
    for ( Int_t i( 0 ); i < nBBs; ++i ) {
        negBBFile->cd();
        negBBDE[i] = dynamic_cast<TH1F*>( negBBFile->Get( BBDENames[i] ) );
        if ( negBBDE[i] == 0 ) {
            std::cerr << "Problem finding histogram \"" << BBDENames[i] << "\" in B- file."
                      << std::endl;
            return EXIT_FAILURE;
        }
        negBBDEPdfs[i] =
            new Lau1DHistPdf( "deltaEsig", negBBDE[i], deMin, deMax, useInterpolation, fluctuateBins );

        posBBFile->cd();
        posBBDE[i] = dynamic_cast<TH1F*>( posBBFile->Get( BBDENames[i] ) );
        if ( posBBDE[i] == 0 ) {
            std::cerr << "Problem finding histogram \"" << BBDENames[i] << "\" in B+ file."
                      << std::endl;
            return EXIT_FAILURE;
        }
        posBBDEPdfs[i] =
            new Lau1DHistPdf( "deltaEsig", posBBDE[i], deMin, deMax, useInterpolation, fluctuateBins );

        fitModel->setBkgndPdfs( BBNames[i], negBBDEPdfs[i], posBBDEPdfs[i] );
    }

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
    fitModel->twoStageFit( doTwoStageFit );

    // Generate toy from the fitted parameters
    //TString fitToyFileName("fitToyMC_Kpipi_");
    //fitToyFileName += iFit;
    //fitToyFileName += ".root";
    //fitModel->compareFitData(100, fitToyFileName);

    // Write out per-event likelihoods and sWeights
    //TString splotFileName("splot_Kpipi_");
    //splotFileName += iFit;
    //splotFileName += ".root";
    //fitModel->writeSPlotData(splotFileName, "splot", kFALSE);

    // Set the names of the files to read/write
    TString dataFile( "gen-Kpipi.root" );
    TString treeName( "genResults" );
    TString rootFileName( "" );
    TString tableFileName( "" );
    if ( command == "fit" ) {
        rootFileName   = "fitKpipi_";
        rootFileName  += iFit;
        rootFileName  += "_expts";
        rootFileName  += firstExpt;
        rootFileName  += "-";
        rootFileName  += firstExpt + nExpt - 1;
        rootFileName  += ".root";
        tableFileName  = "fitKpipiResults_";
        tableFileName += iFit;
    } else {
        rootFileName  = "dummy.root";
        tableFileName = "genKpipiResults";
    }

    // Execute the generation/fit
    fitModel->run( command, dataFile, treeName, rootFileName, tableFileName );

    return EXIT_SUCCESS;
}
