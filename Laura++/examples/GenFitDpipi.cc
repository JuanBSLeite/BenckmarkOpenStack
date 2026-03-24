
#include "LauAbsModIndPartWave.hh"
#include "LauBkgndDPModel.hh"
#include "LauDaughters.hh"
#include "LauEffModel.hh"
#include "LauIsobarDynamics.hh"
#include "LauMagPhaseCoeffSet.hh"
#include "LauModIndPartWaveMagPhase.hh"
#include "LauModIndPartWaveRealImag.hh"
#include "LauRandom.hh"
#include "LauRealImagCoeffSet.hh"
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
    out << "or\n";
    out << progName << " plot" << std::endl;
}

int main( int argc, char** argv )
{
    // Process command-line arguments
    // Usage:
    // ./GenFitDpipi gen [nExpt = 1] [firstExpt = 0]
    // or
    // ./GenFitDpipi fit <iFit> [nExpt = 1] [firstExpt = 0]
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
    Bool_t squareDP = kTRUE;

    // This defines the DP => decay is B+ -> pi+ pi+ D-
    // Particle 1 = pi+
    // Particle 2 = pi+
    // Particle 3 = D-
    // The DP is defined in terms of m13Sq and m23Sq
    LauDaughters* daughters = new LauDaughters( "B+", "pi+", "pi+", "D-", squareDP );

    // Optionally apply some vetoes to the DP
    LauVetoes* vetoes = new LauVetoes();

    LauEffModel* effModel = new LauEffModel( daughters, vetoes );

    // Set the values of the Blatt-Weisskopf barrier radii and whether they are fixed or floating
    LauResonanceMaker& resMaker = LauResonanceMaker::get();
    resMaker.setDefaultBWRadius( LauBlattWeisskopfFactor::Parent, 4.0 );
    resMaker.setDefaultBWRadius( LauBlattWeisskopfFactor::Charm, 4.0 );
    resMaker.setDefaultBWRadius( LauBlattWeisskopfFactor::Light, 4.0 );
    resMaker.setDefaultBWRadius( LauBlattWeisskopfFactor::Beauty, 4.0 );
    resMaker.fixBWRadius( LauBlattWeisskopfFactor::Parent, kTRUE );
    resMaker.fixBWRadius( LauBlattWeisskopfFactor::Charm, kTRUE );
    resMaker.fixBWRadius( LauBlattWeisskopfFactor::Light, kTRUE );
    resMaker.fixBWRadius( LauBlattWeisskopfFactor::Beauty, kTRUE );
    // Create the isobar model

    std::vector<LauAbsCoeffSet*> coeffset;
    LauIsobarDynamics* sigModel = new LauIsobarDynamics( daughters, effModel );
    LauAbsResonance* reson( 0 );

    reson = sigModel->addResonance( "D*0", 2, LauAbsResonance::RelBW );
    reson = sigModel->addResonance( "D*0_0", 2, LauAbsResonance::MIPW_MagPhase );
    LauModIndPartWaveMagPhase* mipw = dynamic_cast<LauModIndPartWaveMagPhase*>( reson );
    if ( mipw == nullptr ) {
        std::cout << "MIPW pointer is null" << std::endl;
        return 0;
    }

    //vector of knot masses - ignore ends as they are dealt with internally
    std::set<Double_t> knot_mass {
        2.10,
        2.20,
        2.30,
        2.40,
        2.50,
        2.60,
        2.70,
        2.80,
        2.90,
        3.10,
        4.10,
    };

    mipw->defineKnots( knot_mass );
    mipw->floatKnotsSecondStage( kFALSE );

    //Set magnitude and phase for the knots (including the end points here)
    mipw->setKnotAmp( 0, 0.12, -2.82, kFALSE, kFALSE );
    mipw->setKnotAmp( 1, 0.58, -1.56, kFALSE, kFALSE );
    mipw->setKnotAmp( 2, 0.73, -1.00, kFALSE, kFALSE );
    mipw->setKnotAmp( 3, 0.68, -0.42, kFALSE, kFALSE );
    mipw->setKnotAmp( 4, 0.5, 0.0, kTRUE, kTRUE );
    mipw->setKnotAmp( 5, 0.23, -0.00, kFALSE, kFALSE );
    mipw->setKnotAmp( 6, 0.23, -0.42, kFALSE, kFALSE );
    mipw->setKnotAmp( 7, 0.15, -0.31, kFALSE, kFALSE );
    mipw->setKnotAmp( 8, 0.17, -0.63, kFALSE, kFALSE );
    mipw->setKnotAmp( 9, 0.20, -0.87, kFALSE, kFALSE );
    mipw->setKnotAmp( 10, 0.14, -1.16, kFALSE, kFALSE );
    mipw->setKnotAmp( 11, 0.08, 1.02, kFALSE, kFALSE );
    mipw->setKnotAmp( 12, 0.0, 0.0, kTRUE, kTRUE );

    reson = sigModel->addResonance( "D*0_2", 2, LauAbsResonance::RelBW );
    reson = sigModel->addResonance( "D*0_1(2680)", 2, LauAbsResonance::RelBW );
    reson = sigModel->addResonance( "B*0", 2, LauAbsResonance::RelBW );
    reson = sigModel->addResonance( "D*0_3(2760)", 2, LauAbsResonance::RelBW );
    reson = sigModel->addResonance( "D0(3000)", 2, LauAbsResonance::RelBW );

    sigModel->setASqMaxValue( 1.0 );

    // Create the fit model
    LauSimpleFitModel* fitModel = new LauSimpleFitModel( sigModel );

    coeffset.push_back( new LauMagPhaseCoeffSet( "D*0", 0.55, -0.38, kFALSE, kFALSE ) );
    // Still need a total mag/phase for the MIPW shape
    coeffset.push_back( new LauMagPhaseCoeffSet( "D*0_0", 1.26, -0.28, kFALSE, kFALSE ) );
    coeffset.push_back( new LauMagPhaseCoeffSet( "D*0_2", 1.00, 0.00, kTRUE, kTRUE ) );
    coeffset.push_back( new LauMagPhaseCoeffSet( "D*0_1(2680)", 0.48, 2.47, kFALSE, kFALSE ) );
    coeffset.push_back( new LauMagPhaseCoeffSet( "B*0", 0.27, 0.14, kFALSE, kFALSE ) );
    coeffset.push_back( new LauMagPhaseCoeffSet( "D*0_3(2760)", 0.17, 0.01, kFALSE, kFALSE ) );
    coeffset.push_back( new LauMagPhaseCoeffSet( "D0(3000)", 0.08, -0.84, kFALSE, kFALSE ) );

    for ( std::vector<LauAbsCoeffSet*>::iterator iter = coeffset.begin(); iter != coeffset.end();
          ++iter ) {
        fitModel->setAmpCoeffSet( *iter );
    }

    Double_t nSig( 50000 );

    TString sigEventsName = "signalEvents";
    LauParameter* nSigEvents = new LauParameter( sigEventsName, nSig, -2.0 * nSig, 2.0 * nSig, kTRUE );
    fitModel->setNSigEvents( nSigEvents );

    // Set the number of experiments to generate or fit and which
    // experiment to start with
    fitModel->setNExpts( nExpt, firstExpt, isToy );

    // Switch on/off calculation of asymmetric errors.
    fitModel->useAsymmFitErrors( kFALSE );

    // Randomise initial fit values for the signal mode
    fitModel->useRandomInitFitPars( kTRUE );

    // Switch on/off Poissonian smearing of total number of events
    fitModel->doPoissonSmearing( kTRUE );

    // Switch on/off Extended ML Fit option
    fitModel->doEMLFit( kFALSE );

    // Switch on the two-stage fit (for the resonance parameters)
    fitModel->twoStageFit( kFALSE );

    // Generate toy from the fitted parameters
    //TString fitToyFileName("fitToyMC_Dpipi_");
    //fitToyFileName += iFit;
    //fitToyFileName += ".root";
    //fitModel->compareFitData(100, fitToyFileName);

    // Set the names of the files to read/write
    TString dataFile( "gen-Dpipi.root" );
    TString treeName( "genResults" );
    TString rootFileName( "" );
    TString tableFileName( "" );
    if ( command == "fit" ) {
        rootFileName   = "fitDpipi_";
        rootFileName  += iFit;
        rootFileName  += "_expt_";
        rootFileName  += firstExpt;
        rootFileName  += "-";
        rootFileName  += ( firstExpt + nExpt - 1 );
        rootFileName  += ".root";
        tableFileName  = "fitDpipiResults_";
        tableFileName += iFit;
    } else {
        rootFileName  = "dummy.root";
        tableFileName = "genDpipiResults";
    }

    // Execute the generation/fit
    fitModel->run( command, dataFile, treeName, rootFileName, tableFileName );

    return EXIT_SUCCESS;
}
