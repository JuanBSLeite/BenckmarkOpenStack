
/*
Copyright 2009 University of Warwick

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

/*! \file LauChebychevPdf.cc
    \brief File containing implementation of LauChebychevPdf class.
*/

/*****************************************************************************
 * Class based on RooFit/RooChebychev.                                       *
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

#include "LauChebychevPdf.hh"

#include "TMath.h"
#include "TSystem.h"

#include <iostream>
#include <vector>

LauChebychevPdf::LauChebychevPdf( const TString& theVarName,
                                  const std::vector<LauAbsRValue*>& params,
                                  Double_t minAbscissa,
                                  Double_t maxAbscissa ) :
    LauAbsPdf( theVarName, params, minAbscissa, maxAbscissa ),
    coeffs_( params.size(), 0 )
{
    // Constructor for the Chebychev Polynomial PDF.
    //
    // The parameters in params are the coefficients of the polynomial
    // (polynomial can be anything from 0 to 7 orders).
    // The last two arguments specify the range in which the PDF is defined, and the PDF
    // will be normalised w.r.t. these limits.

    if ( ( this->nParameters() > 7 ) ) {
        std::cerr << "ERROR in LauChebychevPdf constructor: Too many coeffs - can only cope with order <=7 polynomial."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    for ( UInt_t i( 0 ); i < this->nParameters(); ++i ) {
        TString name  = "c";
        name         += i + 1;
        coeffs_[i]    = this->findParameter( name );
    }

    // Cache the normalisation factor
    this->calcNorm();
}

LauChebychevPdf::~LauChebychevPdf()
{
    // Destructor
}

//inline static double p0(double ,double a) {  return a; }
inline static double p1( double t, double a, double b )
{
    return a * t + b;
}

inline static double p2( double t, double a, double b, double c )
{
    return p1( t, p1( t, a, b ), c );
}

inline static double p3( double t, double a, double b, double c, double d )
{
    return p2( t, p1( t, a, b ), c, d );
}

//inline static double p4(double t,double a,double b,double c,double d,double e) {  return p3(t,p1(t,a,b),c,d,e); }

void LauChebychevPdf::calcLikelihoodInfo( const LauAbscissas& abscissas )
{
    // Check that the given abscissa is within the allowed range
    if ( ! this->checkRange( abscissas ) ) {
        gSystem->Exit( EXIT_FAILURE );
    }

    // Get our abscissa
    Double_t abscissa = abscissas[0];

    Double_t xmin = this->getMinAbscissa();
    Double_t xmax = this->getMaxAbscissa();
    Double_t x( -1.0 + 2.0 * ( abscissa - xmin ) / ( xmax - xmin ) );
    Double_t x2( x * x );
    Double_t sum( 0 );
    switch ( coeffs_.size() ) {
        case 7 :
            sum += coeffs_[6]->unblindValue() * x * p3( x2, 64, -112, 56, -7 );
            [[fallthrough]];
        case 6 :
            sum += coeffs_[5]->unblindValue() * p3( x2, 32, -48, 18, -1 );
            [[fallthrough]];
        case 5 :
            sum += coeffs_[4]->unblindValue() * x * p2( x2, 16, -20, 5 );
            [[fallthrough]];
        case 4 :
            sum += coeffs_[3]->unblindValue() * p2( x2, 8, -8, 1 );
            [[fallthrough]];
        case 3 :
            sum += coeffs_[2]->unblindValue() * x * p1( x2, 4, -3 );
            [[fallthrough]];
        case 2 :
            sum += coeffs_[1]->unblindValue() * p1( x2, 2, -1 );
            [[fallthrough]];
        case 1 :
            sum += coeffs_[0]->unblindValue() * x;
            [[fallthrough]];
        case 0 :
            sum += 1;
    }

    this->setUnNormPDFVal( sum );

    // if the parameters are floating then we
    // need to recalculate the normalisation
    if ( ! this->cachePDF() && ! this->withinNormCalc() && ! this->withinGeneration() ) {
        this->calcNorm();
    }
}

void LauChebychevPdf::calcNorm()
{
    Double_t xmin = this->getMinAbscissa();
    Double_t xmax = this->getMaxAbscissa();

    Double_t norm( 0 );
    switch ( coeffs_.size() ) {
        case 7 :
        case 6 :
            norm += coeffs_[5]->unblindValue() * ( -1 + 18. / 3. - 48. / 5. + 32. / 7. );
            [[fallthrough]];
        case 5 :
        case 4 :
            norm += coeffs_[3]->unblindValue() * ( 1 - 8. / 3. + 8. / 5. );
            [[fallthrough]];
        case 3 :
        case 2 :
            norm += coeffs_[1]->unblindValue() * ( -1 + 2. / 3. );
            [[fallthrough]];
        case 1 :
        case 0 :
            norm += 1;
    }
    norm *= xmax - xmin;

    this->setNorm( norm );
}

void LauChebychevPdf::calcPDFHeight( const LauKinematics* /*kinematics*/ )
{
    // TODO - this method can hopefully be improved
    // At present it scans through the range and then increases by a 20% safety factor
    // Maybe there's a better way?

    if ( this->heightUpToDate() ) {
        return;
    }

    // Calculate the PDF height
    LauAbscissas maxPoint( 1 );

    Double_t minAbs = this->getMinAbscissa();
    Double_t maxAbs = this->getMaxAbscissa();
    Double_t range  = maxAbs - minAbs;
    Double_t maxHeight( 0.0 );

    // Just scan through the range
    for ( Double_t point = minAbs; point <= maxAbs; point += range / 1000.0 ) {
        maxPoint[0] = point;
        this->calcLikelihoodInfo( maxPoint );
        Double_t heightAtPoint = this->getUnNormLikelihood();
        if ( heightAtPoint > maxHeight ) {
            maxHeight = heightAtPoint;
        }
    }

    // Mutliply by 120% to be on the safe side
    maxHeight *= 1.2;

    this->setMaxHeight( maxHeight );
}
