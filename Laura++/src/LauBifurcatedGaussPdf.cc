
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

/*! \file LauBifurcatedGaussPdf.cc
    \brief File containing implementation of LauBifurcatedGaussPdf class.
*/

#include "LauBifurcatedGaussPdf.hh"

#include "LauConstants.hh"

#include "TMath.h"
#include "TSystem.h"

#include <iostream>
#include <vector>

LauBifurcatedGaussPdf::LauBifurcatedGaussPdf( const TString& theVarName,
                                              const std::vector<LauAbsRValue*>& params,
                                              Double_t minAbscissa,
                                              Double_t maxAbscissa ) :
    LauAbsPdf( theVarName, params, minAbscissa, maxAbscissa ),
    mean_( 0 ),
    sigmaL_( 0 ),
    sigmaR_( 0 )
{
    // Constructor for the bifurcated Gaussian PDF.
    // The bifurcated Gaussian combines the left half of a
    // Gaussian with resolution sigmaL with the right half
    // of a Gaussian with resolution sigmaR, both having
    // a common mean (or more correctly, peak position).

    // NB the parameters in params are the mean, sigmaL and sigmaR
    // The last two arguments specify the range in which the PDF is defined, and the PDF
    // will be normalised w.r.t. these limits.

    mean_   = this->findParameter( "mean" );
    sigmaL_ = this->findParameter( "sigmaL" );
    sigmaR_ = this->findParameter( "sigmaR" );

    if ( ( this->nParameters() != 3 ) || ( mean_ == 0 ) || ( sigmaL_ == 0 ) || ( sigmaR_ == 0 ) ) {
        std::cerr << "ERROR in LauBifurcatedGaussPdf constructor: LauBifurcatedGaussPdf requires 3 parameters: \"mean\", \"sigmaL\" and \"sigmaR\"."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Cache the normalisation factor
    this->calcNorm();
}

LauBifurcatedGaussPdf::~LauBifurcatedGaussPdf()
{
    // Destructor
}

void LauBifurcatedGaussPdf::calcLikelihoodInfo( const LauAbscissas& abscissas )
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

    // Evaluate the Birfucated Gaussian PDF value
    Double_t arg = abscissa - mean;
    Double_t coef( 0.0 );
    Double_t value( 0.0 );

    if ( arg < 0.0 ) {
        if ( TMath::Abs( sigmaL ) > 1e-30 ) {
            coef = -0.5 / ( sigmaL * sigmaL );
        }
    } else {
        if ( TMath::Abs( sigmaR ) > 1e-30 ) {
            coef = -0.5 / ( sigmaR * sigmaR );
        }
    }
    value = TMath::Exp( coef * arg * arg );

    // Calculate the norm
    Double_t xscaleL = LauConstants::root2 * sigmaL;
    Double_t xscaleR = LauConstants::root2 * sigmaR;

    Double_t integral( 0.0 );
    Double_t norm( 0.0 );
    Double_t result( 0.0 );

    if ( this->getMaxAbscissa() < mean ) {
        integral = sigmaL * ( TMath::Erf( ( this->getMaxAbscissa() - mean ) / xscaleL ) -
                              TMath::Erf( ( this->getMinAbscissa() - mean ) / xscaleL ) );
    } else if ( this->getMinAbscissa() > mean ) {
        integral = sigmaR * ( TMath::Erf( ( this->getMaxAbscissa() - mean ) / xscaleR ) -
                              TMath::Erf( ( this->getMinAbscissa() - mean ) / xscaleR ) );
    } else {
        integral = sigmaR * TMath::Erf( ( this->getMaxAbscissa() - mean ) / xscaleR ) -
                   sigmaL * TMath::Erf( ( this->getMinAbscissa() - mean ) / xscaleL );
    }

    norm = LauConstants::rootPiBy2 * integral;

    // the result
    result = value / norm;
    this->setUnNormPDFVal( result );
}

void LauBifurcatedGaussPdf::calcNorm()
{
    // Nothing to do here, since it already normalized
    this->setNorm( 1.0 );
}

void LauBifurcatedGaussPdf::calcPDFHeight( const LauKinematics* /*kinematics*/ )
{
    if ( this->heightUpToDate() ) {
        return;
    }

    // Get the up to date parameter values
    Double_t mean = mean_->unblindValue();

    LauAbscissas maxPoint( 1 );
    maxPoint[0] = mean;

    // Calculate the PDF height for the Bifurcated Gaussian function.
    if ( mean > this->getMaxAbscissa() ) {
        maxPoint[0] = this->getMaxAbscissa();
    } else if ( mean < this->getMinAbscissa() ) {
        maxPoint[0] = this->getMinAbscissa();
    }
    this->calcLikelihoodInfo( maxPoint );

    Double_t height = this->getUnNormLikelihood();
    this->setMaxHeight( height );
}
