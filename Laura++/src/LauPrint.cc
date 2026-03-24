
/*
Copyright 2004 University of Warwick

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

/*! \file LauPrint.cc
    \brief File containing implementation of LauPrint class.
*/

#include "LauPrint.hh"

#include "TMath.h"

#include <iostream>

LauPrint::LauPrint()
{
}

LauPrint::~LauPrint()
{
}

void LauPrint::printFormat( std::ostream& stream, Double_t value ) const
{
    // Print double value to a specified level of precision
    char* string = new char[100];

    if ( TMath::Abs( value ) < 1e-30 ) {
        sprintf( string, "%5.3f", 0.0 );
    } else if ( TMath::Abs( value ) < 0.01 ) {
        sprintf( string, "%5.3e", value );
    } else if ( TMath::Abs( value ) < 0.1 ) {
        sprintf( string, "%5.3f", value );
    } else {
        sprintf( string, "%5.3f", value );
    }

    stream << string;

    delete[] string;
    string = 0;
}
