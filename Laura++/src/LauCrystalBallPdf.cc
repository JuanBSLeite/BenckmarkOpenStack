
/*
Copyright 2006 University of Warwick

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

/*! \file LauCrystalBallPdf.cc
    \brief File containing implementation of LauCrystalBallPdf class.
*/

/*****************************************************************************
 * Class based on RooFit/RooCBShape.                                         *
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
#include "LauCrystalBallPdf.hh"

#include "TMath.h"
#include "TSystem.h"

LauCrystalBallPdf::LauCrystalBallPdf( const TString& theVarName,
                                      const std::vector<LauAbsRValue*>& params,
                                      Double_t minAbscissa,
                                      Double_t maxAbscissa ) :
    LauAbsPdf( theVarName, params, minAbscissa, maxAbscissa ),
    mean_( 0 ),
    sigma_( 0 ),
    alpha_( 0 ),
    n_( 0 )
{
    // Constructor for the Crystal Ball PDF, which is a gaussian and a decaying tail
    // smoothly matched up. The tail goes as a 1/x^n
    //
    // The parameters in params are the mean and the sigma (half the width) of the gaussian,
    // the distance from the mean in which the gaussian and the tail are matched up (which
    // can be negative or positive), and the power "n" for the tail.
    // The last two arguments specify the range in which the PDF is defined, and the PDF
    // will be normalised w.r.t. these limits.

    mean_  = this->findParameter( "mean" );
    sigma_ = this->findParameter( "sigma" );
    alpha_ = this->findParameter( "alpha" );
    n_     = this->findParameter( "order" );

    if ( ( this->nParameters() != 4 ) || ( mean_ == 0 ) || ( sigma_ == 0 ) || ( alpha_ == 0 ) ||
         ( n_ == 0 ) ) {
        cerr << "ERROR in LauCrystalBallPdf constructor: LauCrystalBallPdf requires 4 parameters: \"mean\", \"sigma\", \"alpha\" and \"order\"."
             << endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Cache the normalisation factor
    this->calcNorm();
}

LauCrystalBallPdf::~LauCrystalBallPdf()
{
    // Destructor
}

void LauCrystalBallPdf::calcLikelihoodInfo( const LauAbscissas& abscissas )
{
    // Check that the given abscissa is within the allowed range
    if ( ! this->checkRange( abscissas ) ) {
        gSystem->Exit( EXIT_FAILURE );
    }

    // Get our abscissa
    Double_t abscissa = abscissas[0];

    // Get the up to date parameter values
    Double_t mean  = mean_->unblindValue();
    Double_t sigma = sigma_->unblindValue();
    Double_t alpha = alpha_->unblindValue();
    Double_t n     = n_->unblindValue();

    Double_t result( 0.0 );
    Double_t t = ( abscissa - mean ) / sigma;
    if ( alpha < 0.0 ) {
        t = -t;
    }

    Double_t absAlpha = TMath::Abs( alpha );

    if ( t >= -absAlpha ) {
        result = TMath::Exp( -0.5 * t * t );
    } else {
        Double_t a = TMath::Power( n / absAlpha, n ) * TMath::Exp( -0.5 * absAlpha * absAlpha );
        Double_t b = n / absAlpha - absAlpha;

        result = a / TMath::Power( b - t, n );
    }

    this->setUnNormPDFVal( result );

    // if the parameters are floating then we
    // need to recalculate the normalisation
    if ( ! this->cachePDF() && ! this->withinNormCalc() && ! this->withinGeneration() ) {
        this->calcNorm();
    }
}

void LauCrystalBallPdf::calcNorm()
{
    // Get the up to date parameter values
    Double_t mean  = mean_->unblindValue();
    Double_t sigma = sigma_->unblindValue();
    Double_t alpha = alpha_->unblindValue();
    Double_t n     = n_->unblindValue();

    Double_t result = 0.0;
    Bool_t useLog   = kFALSE;

    if ( TMath::Abs( n - 1.0 ) < 1.0e-05 ) {
        useLog = kTRUE;
    }

    Double_t sig = TMath::Abs( sigma );

    Double_t tmin = ( this->getMinAbscissa() - mean ) / sig;
    Double_t tmax = ( this->getMaxAbscissa() - mean ) / sig;

    if ( alpha < 0 ) {
        Double_t tmp = tmin;
        tmin         = -tmax;
        tmax         = -tmp;
    }

    Double_t absAlpha = TMath::Abs( alpha );

    if ( tmin >= -absAlpha ) {
        result += sig * LauConstants::rootPiBy2 *
                  ( this->approxErf( tmax / LauConstants::root2 ) -
                    approxErf( tmin / LauConstants::root2 ) );
    } else if ( tmax <= -absAlpha ) {
        Double_t a = TMath::Power( n / absAlpha, n ) * TMath::Exp( -0.5 * absAlpha * absAlpha );
        Double_t b = n / absAlpha - absAlpha;

        if ( useLog == kTRUE ) {
            result += a * sig * ( TMath::Log( b - tmin ) - TMath::Log( b - tmax ) );
        } else {
            result += a * sig / ( 1.0 - n ) *
                      ( 1.0 / ( TMath::Power( b - tmin, n - 1.0 ) ) -
                        1.0 / ( TMath::Power( b - tmax, n - 1.0 ) ) );
        }
    } else {
        Double_t a = TMath::Power( n / absAlpha, n ) * TMath::Exp( -0.5 * absAlpha * absAlpha );
        Double_t b = n / absAlpha - absAlpha;

        Double_t term1 = 0.0;
        if ( useLog == kTRUE )
            term1 = a * sig * ( TMath::Log( b - tmin ) - TMath::Log( n / absAlpha ) );
        else
            term1 = a * sig / ( 1.0 - n ) *
                    ( 1.0 / ( TMath::Power( b - tmin, n - 1.0 ) ) -
                      1.0 / ( TMath::Power( n / absAlpha, n - 1.0 ) ) );

        Double_t term2 = sig * LauConstants::rootPiBy2 *
                         ( this->approxErf( tmax / LauConstants::root2 ) -
                           this->approxErf( -absAlpha / LauConstants::root2 ) );

        result += term1 + term2;
    }
    this->setNorm( result );
}

void LauCrystalBallPdf::calcPDFHeight( const LauKinematics* /*kinematics*/ )
{
    if ( this->heightUpToDate() ) {
        return;
    }

    // Get the up to date parameter values
    Double_t mean = mean_->unblindValue();

    // The Crystall Ball function is a Gaussian with an exponentially decaying tail
    // Therefore, calculate the PDF height for the Gaussian function.
    LauAbscissas abscissa( 1 );
    abscissa[0] = mean;
    this->calcLikelihoodInfo( abscissa );

    Double_t height = this->getUnNormLikelihood();
    this->setMaxHeight( height );
}

Double_t LauCrystalBallPdf::approxErf( Double_t arg ) const
{
    static const Double_t erflim = 5.0;
    if ( arg > erflim ) {
        return 1.0;
    }
    if ( arg < -erflim ) {
        return -1.0;
    }

    return TMath::Erf( arg );
}
