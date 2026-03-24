
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

#include "LauCalcChiSq.hh"

#include "TString.h"

int main( const int argc, const char** argv )
{

    TString inputFile( "chiSqInput.txt" );

    if ( argc > 1 ) {
        inputFile = TString( argv[1] );
    }

    LauCalcChiSq a( inputFile );

    //a.setVerbose(kTRUE);
    a.run();

    return 0;
}
