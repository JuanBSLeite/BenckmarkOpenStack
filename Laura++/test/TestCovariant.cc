
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

#include "LauBelleNR.hh"
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
    out << progName << " <spin>" << std::endl;
}

int main( int argc, char** argv )
{
    // Process command-line arguments
    // Usage:
    // ./TestCovariant <spin>
    if ( argc < 2 ) {
        usage( std::cerr, argv[0] );
        return EXIT_FAILURE;
    }

    const Int_t spin = atoi( argv[1] );
    if ( spin < 0 || spin > 3 ) {
        std::cerr << "Only valid for spin 0, 1, 2 or 3" << std::endl;
        return EXIT_FAILURE;
    }

    // If you want to use square DP histograms for efficiency,
    // backgrounds or you just want the square DP co-ordinates
    // stored in the toy MC ntuple then set this to kTRUE
    Bool_t squareDP = kFALSE;

    // This defines the DP => decay is B+ -> K+ pi+ pi-
    // Particle 1 = K+
    // Particle 2 = pi+
    // Particle 3 = pi-
    // The DP is defined in terms of m13Sq and m23Sq
    LauDaughters* daughters = new LauDaughters( "B+", "K+", "pi+", "pi-", squareDP );

    // Optionally apply some vetoes to the DP
    LauVetoes* vetoes = new LauVetoes();

    // Define the efficiency model (defaults to unity everywhere)
    // Can optionally provide a histogram to model variation over DP
    // (example syntax given in commented-out section)
    LauEffModel* effModel = new LauEffModel( daughters, vetoes );

    // Create the isobar model

    // Set the values of the Blatt-Weisskopf barrier radii and whether they are fixed or floating
    LauResonanceMaker& resMaker = LauResonanceMaker::get();
    resMaker.setSpinFormalism( LauAbsResonance::Covariant );
    resMaker.setDefaultBWRadius( LauBlattWeisskopfFactor::Parent, 5.0 );
    resMaker.setDefaultBWRadius( LauBlattWeisskopfFactor::Light, 4.0 );
    resMaker.fixBWRadius( LauBlattWeisskopfFactor::Parent, kTRUE );
    resMaker.fixBWRadius( LauBlattWeisskopfFactor::Light, kFALSE );

    LauIsobarDynamics* sigModel = new LauIsobarDynamics( daughters, effModel );
    LauAbsResonance* reson( 0 );
    std::map<Int_t, TString> waves;
    waves[0]             = "S";
    waves[1]             = "P";
    waves[2]             = "D";
    waves[3]             = "F";
    const TString nrName = TString::Format( "BelleNR_%swave", waves[spin].Data() );
    reson                = sigModel->addResonance( nrName, 2, LauAbsResonance::BelleNR );
    reson->setResonanceParameter( "alpha", 0.0 );
    LauBelleNR* belleNR = dynamic_cast<LauBelleNR*>( reson );
    belleNR->enforceLegendreSpinFactors( kFALSE );

    // Reset the maximum signal DP ASq value
    // This will be automatically adjusted to avoid bias or extreme
    // inefficiency if you get the value wrong but best to set this by
    // hand once you've found the right value through some trial and
    // error.
    sigModel->setASqMaxValue( 0.32 );

    // Create the fit model
    LauSimpleFitModel* fitModel = new LauSimpleFitModel( sigModel );

    // Create the complex coefficients for the isobar model
    // Here we're using the magnitude and phase form:
    // c_j = a_j exp(i*delta_j)
    std::vector<LauAbsCoeffSet*> coeffset;
    coeffset.push_back( new LauMagPhaseCoeffSet( nrName, 1.00, 0.00, kTRUE, kTRUE ) );
    for ( std::vector<LauAbsCoeffSet*>::iterator iter = coeffset.begin(); iter != coeffset.end();
          ++iter ) {
        fitModel->setAmpCoeffSet( *iter );
    }

    // Set the signal yield and define whether it is fixed or floated
    const Double_t nSig = 100000.0;
    LauParameter* nSigEvents = new LauParameter( "nSigEvents", nSig, -2.0 * nSig, 2.0 * nSig, kFALSE );
    fitModel->setNSigEvents( nSigEvents );

    // Set the number of experiments to generate or fit and which
    // experiment to start with
    Int_t nExpt( 1 );
    Int_t firstExpt( 0 );
    fitModel->setNExpts( nExpt, firstExpt, kTRUE );

    // Switch on/off calculation of asymmetric errors.
    fitModel->useAsymmFitErrors( kFALSE );

    // Randomise initial fit values for the signal mode
    fitModel->useRandomInitFitPars( kTRUE );

    const Bool_t haveBkgd = ( fitModel->nBkgndClasses() > 0 );

    // Switch on/off Poissonian smearing of total number of events
    fitModel->doPoissonSmearing( haveBkgd );

    // Switch on/off Extended ML Fit option
    fitModel->doEMLFit( haveBkgd );

    // Set the names of the files to read/write
    TString dataFile = TString::Format( "gen-covtest-spin%d.root", spin );
    TString treeName( "genResults" );
    TString rootFileName( "dummy.root" );
    TString tableFileName( "genResults" );

    // Execute the generation/fit
    fitModel->run( "gen", dataFile, treeName, rootFileName, tableFileName );

    return EXIT_SUCCESS;
}
