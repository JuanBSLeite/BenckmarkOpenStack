
/*
Copyright 2013 University of Warwick

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

#include "LauResultsExtractor.hh"

#include "TString.h"

#include <cstdlib>

int main( const int argc, const char** argv )
{
    Int_t nExpt( 1 );
    TString inputFileName( "input-list.txt" );
    TString outputFileName( "output.root" );
    TString treeName( "fitResults" );

    if ( argc > 1 ) {
        nExpt = atoi( argv[1] );
    }
    if ( argc > 2 ) {
        inputFileName = argv[2];
    }
    if ( argc > 3 ) {
        outputFileName = argv[3];
    }
    if ( argc > 4 ) {
        treeName = argv[4];
    }

    LauResultsExtractor myExtractor( inputFileName, outputFileName, treeName );
    myExtractor.process( nExpt );

    return EXIT_SUCCESS;
}
