
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

/*! \file LauLinearPdf.cc
    \brief File containing implementation of LauLinearPdf class.
*/

#include <iostream>
#include <vector>
using std::cerr;
using std::cout;
using std::endl;

#include "LauConstants.hh"
#include "LauLinearPdf.hh"

#include "TMath.h"
#include "TSystem.h"

LauLinearPdf::LauLinearPdf( const TString& theVarName,
                            const std::vector<LauAbsRValue*>& params,
                            Double_t minAbscissa,
                            Double_t maxAbscissa ) :
    LauAbsPdf( theVarName, params, minAbscissa, maxAbscissa ),
    slope_( 0 ),
    posflag_( kTRUE )
{
    // Constructor for the linear PDF.
    //
    // The parameters in params are the slope and y-intercept of the line.
    // The last two arguments specify the range in which the PDF is defined, and the PDF
    // will be normalised w.r.t. these limits.

    slope_ = this->findParameter( "slope" );

    if ( ( this->nParameters() != 1 ) || ( slope_ == 0 ) ) {
        cerr << "Warning. LauLinearPdf requires 1 parameter: \"slope\"." << endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Cache the normalisation factor.
    this->calcNorm();
}

LauLinearPdf::~LauLinearPdf()
{
    // Destructor
}

void LauLinearPdf::calcLikelihoodInfo( const LauAbscissas& abscissas )
{
    // Check that the given abscissa is within the allowed range
    if ( ! this->checkRange( abscissas ) ) {
        gSystem->Exit( EXIT_FAILURE );
    }

    // Get our abscissa
    Double_t abscissa = abscissas[0];

    // Get the up to date parameter values
    Double_t slope = slope_->unblindValue();

    // Calculate the value of the straight line for the given value of the abscissa.
    Double_t constVal  = 1.0 / ( this->getMaxAbscissa() - this->getMinAbscissa() );
    constVal          -= slope * ( this->getMaxAbscissa() + this->getMinAbscissa() ) * 0.5;

    Double_t value = slope * abscissa + constVal;

    // Make sure the PDF doesn't go negative
    if ( value < 0.0 ) {
        if ( posflag_ ) {
            std::cerr << "WARNING in LauLinearPdf::calcLikelihoodInfo : The PDF is negative, setting to zero"
                      << std::endl;
            posflag_ = kFALSE;
        }
        value = 0.0;
    }

    this->setUnNormPDFVal( value );
}

void LauLinearPdf::calcNorm()
{
    // Nothing to calculate here since the PDF is already normalised to 1
    this->setNorm( 1.0 );
}

void LauLinearPdf::calcPDFHeight( const LauKinematics* /*kinematics*/ )
{
    if ( this->heightUpToDate() ) {
        return;
    }

    // Get the up to date parameter values
    Double_t slope = slope_->unblindValue();

    // Calculate the PDF height for the straight line
    LauAbscissas maxPoint( 1 );
    if ( slope > 0.0 ) {
        maxPoint[0] = this->getMaxAbscissa();
    } else {
        maxPoint[0] = this->getMinAbscissa();
    }
    this->calcLikelihoodInfo( maxPoint );

    Double_t height = this->getUnNormLikelihood();

    this->setMaxHeight( height );
}
