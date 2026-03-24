
/*
Copyright 2022 University of Warwick

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

#include "LauDaughters.hh"
#include "LauKinematics.hh"
#include "LauRandom.hh"

#include "TFile.h"
#include "TRandom.h"
#include "TTree.h"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

/*
 * Generates toy MC uniformly in the Dalitz plot
 *
 * Usage: FlatPhaseSpace <parent> <child_1> <child_2> <child_3> <nevents>
 */

int main( int argc, char** argv )
{
    // Convert the command-line arguments into a useful format
    const std::vector<TString> cmdLineArgs { argv, argv + argc };

    if ( cmdLineArgs.size() != 6 ) {
        std::cout << "Usage: " << cmdLineArgs[0]
                  << " <parent> <child_1> <child_2> <child_3> <nevents>" << std::endl;
        return 0;
    }

    const TString& parent { cmdLineArgs[1] };
    const TString& child_1 { cmdLineArgs[2] };
    const TString& child_2 { cmdLineArgs[3] };
    const TString& child_3 { cmdLineArgs[4] };
    const ULong_t nEvents { std::stoul( cmdLineArgs[5].Data() ) };

    const Bool_t squareDP { kTRUE };

    LauDaughters* daughters { new LauDaughters( parent, child_1, child_2, child_3, squareDP ) };
    LauKinematics* kinematics { daughters->getKinematics() };

    Double_t m12Sq { 0.0 }, m13Sq { 0.0 }, m23Sq { 0.0 }, mPrime { 0.0 }, thPrime { 0.0 };

    const TString outputFileName { TString::Format( "%s_%s_%s_%s.root",
                                                    daughters->getSanitisedNameParent().Data(),
                                                    daughters->getSanitisedNameDaug1().Data(),
                                                    daughters->getSanitisedNameDaug2().Data(),
                                                    daughters->getSanitisedNameDaug3().Data() ) };
    TFile* file { TFile::Open( outputFileName, "recreate" ) };
    TTree* tree { new TTree( "genResults", "genResults" ) };

    tree->Branch( "m12Sq", &m12Sq, "m12Sq/D" );
    tree->Branch( "m13Sq", &m13Sq, "m13Sq/D" );
    tree->Branch( "m23Sq", &m23Sq, "m23Sq/D" );
    if ( squareDP ) {
        tree->Branch( "mPrime", &mPrime, "mPrime/D" );
        tree->Branch( "thPrime", &thPrime, "thPrime/D" );
    }

    for ( ULong_t i { 0 }; i < nEvents; ++i ) {

        if ( nEvents > 1000 && i % ( nEvents / 1000 ) == 0 ) {
            std::cout << "Generating event " << i << std::endl;
        }

        kinematics->genFlatPhaseSpace( m13Sq, m23Sq );
        kinematics->updateKinematics( m13Sq, m23Sq );
        m12Sq = kinematics->getm12Sq();
        if ( squareDP ) {
            mPrime  = kinematics->getmPrime();
            thPrime = kinematics->getThetaPrime();
        }

        tree->Fill();
    }

    file->Write();
    file->Close();
    delete file;

    return EXIT_SUCCESS;
}
