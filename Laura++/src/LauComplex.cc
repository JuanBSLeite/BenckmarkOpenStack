
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

/*! \file LauComplex.cc
    \brief File containing implementation of LauComplex class.
*/

/*****************************************************************************
 * Class based on RooFit/RooComplex.                                         *
 * Original copyright given below.                                           *
 *****************************************************************************
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#include "LauComplex.hh"

#include <iostream>

void LauComplex::print() const
{
    std::cout << ( *this ) << std::endl;
}

std::ostream& operator<<( std::ostream& os, const LauComplex& z )
{
    return os << "(" << z.re() << "," << z.im() << ")";
}

std::istream& operator>>( std::istream& ios, LauComplex& z )
{
    Double_t a( 0.0 ), b( 0.0 );
    char temp;
    ios >> temp >> a >> temp >> b >> temp;
    z.setRealImagPart( a, b );
    return ios;
}
