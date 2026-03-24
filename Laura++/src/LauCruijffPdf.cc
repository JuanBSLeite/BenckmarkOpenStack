
/*
Copyright 2008 University of Warwick

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

/*! \file LauCruijffPdf.cc
    \brief File containing implementation of LauCruijffPdf class.
*/

/*****************************************************************************
 * Class based on RooFit/RooCruijff.                                         *
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

#include <iostream>
#include <vector>
using std::cerr;
using std::cout;
using std::endl;

#include "LauConstants.hh"
#include "LauCruijffPdf.hh"

#include "TMath.h"
#include "TSystem.h"

LauCruijffPdf::LauCruijffPdf( const TString& theVarName,
                              const std::vector<LauAbsRValue*>& params,
                              Double_t minAbscissa,
                              Double_t maxAbscissa ) :
    LauAbsPdf( theVarName, params, minAbscissa, maxAbscissa ),
    mean_( 0 ),
    sigmaL_( 0 ),
    sigmaR_( 0 ),
    alphaL_( 0 ),
    alphaR_( 0 )
{
    // Constructor for the Cruijff PDF.
    //
    // The parameters in params are the mean, sigmaR, sigmaL, alphaR
    // and alphaL.
    // The last two arguments specify the range in which the PDF is defined, and the PDF
    // will be normalised w.r.t. these limits.

    mean_   = this->findParameter( "mean" );
    sigmaR_ = this->findParameter( "sigmaR" );
    sigmaL_ = this->findParameter( "sigmaL" );
    alphaR_ = this->findParameter( "alphaR" );
    alphaL_ = this->findParameter( "alphaL" );

    if ( ( this->nParameters() != 5 ) || ( mean_ == 0 ) || ( sigmaR_ == 0 ) || ( sigmaL_ == 0 ) ||
         ( alphaL_ == 0 ) || ( alphaR_ == 0 ) ) {
        cerr << "ERROR in LauCruijffPdf constructor: LauCruijffPdf requires 5 parameters: \"mean\", \"sigmaL\", \"sigmaR\", \"alphaR\" and \"alphaL\"."
             << endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Cache the normalisation factor
    this->calcNorm();
}

LauCruijffPdf::~LauCruijffPdf()
{
    // Destructor
}

void LauCruijffPdf::calcLikelihoodInfo( const LauAbscissas& abscissas )
{

    // Check that the given abscissa is within the allowed range
    if ( ! this->checkRange( abscissas ) ) {
        gSystem->Exit( EXIT_FAILURE );
    }

    // Get our abscissa
    Double_t abscissa = abscissas[0];

    // Get the up to date parameter values
    Double_t mean   = mean_->unblindValue();
    Double_t sigmaL = sigmaL_->unblindValue();
    Double_t sigmaR = sigmaR_->unblindValue();
    Double_t alphaL = alphaL_->unblindValue();
    Double_t alphaR = alphaR_->unblindValue();

    // Evaluate the LauCruijff PDF value

    Double_t arg = abscissa - mean;
    Double_t coef( 0.0 );
    Double_t value( 0.0 );

    if ( arg < 0.0 ) {
        if ( TMath::Abs( sigmaL ) > 1e-30 ) {
            coef = -1.0 / ( 2.0 * sigmaL * sigmaL + alphaL * arg * arg );
        }
    } else {
        if ( TMath::Abs( sigmaR ) > 1e-30 ) {
            coef = -1.0 / ( 2.0 * sigmaR * sigmaR + alphaR * arg * arg );
        }
    }
    value = TMath::Exp( coef * arg * arg );

    // if the parameters are floating then we
    // need to recalculate the normalisation
    if ( ! this->cachePDF() && ! this->withinNormCalc() && ! this->withinGeneration() ) {
        this->calcNorm();
    }

    this->setUnNormPDFVal( value );
}

void LauCruijffPdf::calcPDFHeight( const LauKinematics* /*kinematics*/ )
{
    if ( this->heightUpToDate() ) {
        return;
    }

    // Get the up to date parameter values
    Double_t mean = mean_->unblindValue();

    LauAbscissas maxPoint( 1 );
    maxPoint[0] = mean;

    // Calculate the PDF height

    if ( mean < this->getMinAbscissa() ) {
        maxPoint[0] = this->getMinAbscissa();
    } else if ( mean > this->getMaxAbscissa() ) {
        maxPoint[0] = this->getMaxAbscissa();
    }

    this->calcLikelihoodInfo( maxPoint );
    Double_t height = this->getUnNormLikelihood();

    // Multiply by a small factor to avoid problems from rounding errors
    height *= ( 1.0 + 1e-1 );

    this->setMaxHeight( height );
}
