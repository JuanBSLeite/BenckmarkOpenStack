
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

/*! \file LauGaussPdf.cc
    \brief File containing implementation of LauGaussPdf class.
*/

#include <iostream>
#include <vector>
using std::cerr;
using std::cout;
using std::endl;

#include "LauConstants.hh"
#include "LauGaussPdf.hh"

#include "TMath.h"
#include "TSystem.h"

LauGaussPdf::LauGaussPdf( const TString& theVarName,
                          const std::vector<LauAbsRValue*>& params,
                          Double_t minAbscissa,
                          Double_t maxAbscissa ) :
    LauAbsPdf( theVarName, params, minAbscissa, maxAbscissa ),
    mean_( 0 ),
    sigma_( 0 )
{
    // Constructor for the Gaussian PDF.
    //
    // The parameters in params are the mean and the sigma (half the width) of the gaussian.
    // The last two arguments specify the range in which the PDF is defined, and the PDF
    // will be normalised w.r.t. these limits.

    mean_  = this->findParameter( "mean" );
    sigma_ = this->findParameter( "sigma" );

    if ( ( this->nParameters() != 2 ) || ( mean_ == 0 ) || ( sigma_ == 0 ) ) {
        cerr << "ERROR in LauGaussPdf constructor: LauGaussPdf requires 2 parameters: \"mean\" and \"sigma\"."
             << endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Cache the normalisation factor
    this->calcNorm();
}

LauGaussPdf::~LauGaussPdf()
{
    // Destructor
}

void LauGaussPdf::calcLikelihoodInfo( const LauAbscissas& abscissas )
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

    // Calculate the value of the Gaussian for the given value of the abscissa.
    Double_t arg = abscissa - mean;

    Double_t exponent( 0.0 );
    if ( TMath::Abs( sigma ) > 1e-10 ) {
        exponent = -0.5 * arg * arg / ( sigma * sigma );
    }

    Double_t value = TMath::Exp( exponent );
    this->setUnNormPDFVal( value );

    // if the parameters are floating then we
    // need to recalculate the normalisation
    if ( ! this->cachePDF() && ! this->withinNormCalc() && ! this->withinGeneration() ) {
        this->calcNorm();
    }
}

void LauGaussPdf::calcNorm()
{
    // Get the up to date parameter values
    Double_t mean  = mean_->unblindValue();
    Double_t sigma = sigma_->unblindValue();

    // Calculate the normalisation of the gaussian and cache it.
    Double_t scale = LauConstants::root2 * sigma;
    Double_t norm( 0.0 );

    if ( TMath::Abs( sigma ) > 1e-10 ) {
        norm = LauConstants::rootPiBy2 * sigma *
               ( TMath::Erf( ( this->getMaxAbscissa() - mean ) / scale ) -
                 TMath::Erf( ( this->getMinAbscissa() - mean ) / scale ) );
    }

    this->setNorm( norm );
}

void LauGaussPdf::calcPDFHeight( const LauKinematics* /*kinematics*/ )
{
    if ( this->heightUpToDate() ) {
        return;
    }

    // Get the up to date parameter values
    Double_t mean = mean_->unblindValue();

    LauAbscissas maxPoint( 1 );
    maxPoint[0] = mean;

    // Calculate the PDF height for the Gaussian function.
    if ( mean > this->getMaxAbscissa() ) {
        maxPoint[0] = this->getMaxAbscissa();
    } else if ( mean < this->getMinAbscissa() ) {
        maxPoint[0] = this->getMinAbscissa();
    }
    this->calcLikelihoodInfo( maxPoint );

    Double_t height = this->getUnNormLikelihood();
    this->setMaxHeight( height );
}
