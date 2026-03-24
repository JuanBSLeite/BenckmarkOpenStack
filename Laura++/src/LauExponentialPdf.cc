
/*
Copyright 2011 University of Warwick

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

/*! \file LauExponentialPdf.cc
    \brief File containing implementation of LauExponentialPdf class.
*/

#include <iostream>
#include <vector>
using std::cerr;
using std::cout;
using std::endl;

#include "LauConstants.hh"
#include "LauExponentialPdf.hh"

#include "TMath.h"
#include "TSystem.h"

LauExponentialPdf::LauExponentialPdf( const TString& theVarName,
                                      const std::vector<LauAbsRValue*>& params,
                                      Double_t minAbscissa,
                                      Double_t maxAbscissa ) :
    LauAbsPdf( theVarName, params, minAbscissa, maxAbscissa ),
    slope_( 0 )
{
    // Constructor for the Exponential PDF.
    //
    // The parameter in param is the slope of the exponential ie exp(slope*x).
    // The last two arguments specify the range in which the PDF is defined, and the PDF
    // will be normalised w.r.t. these limits.

    slope_ = this->findParameter( "slope" );

    if ( ( this->nParameters() != 1 ) || ( slope_ == 0 ) ) {
        cerr << "ERROR in LauExponentialPdf constructor: LauExponentialPdf requires 1 parameter: \"slope\"."
             << endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Cache the normalisation factor
    this->calcNorm();
}

LauExponentialPdf::~LauExponentialPdf()
{
    // Destructor
}

void LauExponentialPdf::calcLikelihoodInfo( const LauAbscissas& abscissas )
{
    // Check that the given abscissa is within the allowed range
    if ( ! this->checkRange( abscissas ) ) {
        gSystem->Exit( EXIT_FAILURE );
    }

    // Get our abscissa
    Double_t abscissa = abscissas[0];

    // Get the up to date parameter values
    Double_t slope = slope_->unblindValue();

    // Calculate the value of the Exponential for the given value of the abscissa.

    Double_t exponent( 0.0 );
    if ( slope != 0 ) {
        exponent = slope * abscissa;
    }

    Double_t value = TMath::Exp( exponent );
    this->setUnNormPDFVal( value );

    // if the parameters are floating then we
    // need to recalculate the normalisation
    if ( ! this->cachePDF() && ! this->withinNormCalc() && ! this->withinGeneration() ) {
        this->calcNorm();
    }
}

void LauExponentialPdf::calcNorm()
{
    // Get the up to date parameter values
    Double_t slope = slope_->unblindValue();

    // Calculate the normalisation of the exponential and cache it.
    Double_t norm( 0.0 );

    if ( slope != 0 ) {
        norm = ( 1 / slope ) * ( TMath::Exp( slope * this->getMaxAbscissa() ) -
                                 TMath::Exp( slope * this->getMinAbscissa() ) );
    }

    this->setNorm( norm );
}

void LauExponentialPdf::calcPDFHeight( const LauKinematics* /*kinematics*/ )
{
    if ( this->heightUpToDate() ) {
        return;
    }

    // Get the up to date parameter values
    Double_t slope = slope_->unblindValue();

    LauAbscissas maxPoint( 1 );
    maxPoint[0] = 0;

    // Calculate the PDF height for the Exponential function.
    if ( slope > 0 ) {
        maxPoint[0] = this->getMaxAbscissa();
    } else if ( slope < 0 ) {
        maxPoint[0] = this->getMinAbscissa();
    }
    this->calcLikelihoodInfo( maxPoint );

    Double_t height = this->getUnNormLikelihood();
    this->setMaxHeight( height );
}
