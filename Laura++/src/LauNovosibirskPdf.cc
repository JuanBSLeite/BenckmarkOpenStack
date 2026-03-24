
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

/*! \file LauNovosibirskPdf.cc
    \brief File containing implementation of LauNovosibirskPdf class.
*/

#include <iostream>

using std::cerr;
using std::cout;
using std::endl;

#include "LauConstants.hh"
#include "LauNovosibirskPdf.hh"

#include "TMath.h"
#include "TSystem.h"

LauNovosibirskPdf::LauNovosibirskPdf( const TString& theVarName,
                                      const std::vector<LauAbsRValue*>& params,
                                      Double_t minAbscissa,
                                      Double_t maxAbscissa ) :
    LauAbsPdf( theVarName, params, minAbscissa, maxAbscissa ),
    mean_( 0 ),
    sigma_( 0 ),
    tail_( 0 )
{
    // Constructor for the Novosibirsk PDF.
    //
    // The parameters in params are the mean, sigma and tail.
    // The last two arguments specify the range in which the PDF is defined, and the PDF
    // will be normalised w.r.t. these limits.

    mean_  = this->findParameter( "mean" );
    sigma_ = this->findParameter( "sigma" );
    tail_  = this->findParameter( "tail" );

    if ( ( this->nParameters() != 3 ) || ( mean_ == 0 ) || ( sigma_ == 0 ) || ( tail_ == 0 ) ) {
        cerr << "ERROR in LauNovosibirskPdf constructor: LauNovosibirskPdf requires 3 parameters: \"mean\", \"sigma\"and \"tail\" "
             << endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Cache the normalisation factor
    this->calcNorm();
}

LauNovosibirskPdf::~LauNovosibirskPdf()
{
    // Destructor
}

void LauNovosibirskPdf::calcLikelihoodInfo( const LauAbscissas& abscissas )
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
    Double_t tail  = tail_->unblindValue();

    // Evaluate the Novosibirsk PDF value

    Double_t qa( 0.0 ), qb( 0.0 ), qx( 0.0 ), qy( 0.0 );
    Double_t arg( 0.0 );
    Double_t value( 0.0 );

    if ( TMath::Abs( tail ) < 1.e-7 )
        arg = 0.5 * ( ( abscissa - mean ) / sigma ) * ( ( abscissa - mean ) / sigma );
    else {
        qa = tail * TMath::Sqrt( LauConstants::log4 );
        qb = TMath::SinH( qa ) / qa;
        qx = ( abscissa - mean ) / sigma * qb;
        qy = 1.0 + tail * qx;

        //---- Cutting curve from right side

        if ( qy > 1.E-7 ) {
            arg = 0.5 * ( ( log( qy ) / tail ) * ( log( qy ) / tail ) + tail * tail );
        } else {
            arg = 15.0;
        }
    }

    value = TMath::Exp( -arg );

    // if the parameters are floating then we
    // need to recalculate the normalisation
    if ( ! this->cachePDF() && ! this->withinNormCalc() && ! this->withinGeneration() ) {
        this->calcNorm();
    }
    this->setUnNormPDFVal( value );
}

void LauNovosibirskPdf::calcPDFHeight( const LauKinematics* /*kinematics*/ )
{
    if ( this->heightUpToDate() ) {
        return;
    }

    // Get the up to date parameter values
    Double_t mean = mean_->unblindValue();

    LauAbscissas maxPoint( 1 );
    maxPoint[0] = mean;

    // Calculate the PDF height for the Bifurcated Gaussian function.

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
