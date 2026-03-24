
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

/*! \file LauRandom.cc
    \brief File containing implementation of LauRandom methods.
*/

#include "LauRandom.hh"

#include "TRandom3.h"

TRandom* LauRandom::randomFun()
{
    // Returns a pointer to a singleton random-number generator implementation.
    // Creates the object the first time it is called.

    static TRandom* theGenerator = 0;
    if ( theGenerator == 0 ) {
        theGenerator = new TRandom3( 65539 );
    }
    return theGenerator;
}

TRandom* LauRandom::zeroSeedRandom()
{
    // Returns a pointer to a singleton random-number generator implementation.
    // Creates the object the first time it is called.

    static TRandom* theGenerator = 0;
    if ( theGenerator == 0 ) {
        theGenerator = new TRandom3( 0 );
    }
    // Use TRandom3(0) to set seed from machine clock time (within +-1 sec).
    return theGenerator;
}

void LauRandom::setSeed( UInt_t seed )
{
    TRandom* theGenerator = randomFun();
    theGenerator->SetSeed( seed );
}
