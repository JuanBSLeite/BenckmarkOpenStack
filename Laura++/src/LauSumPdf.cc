
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

/*! \file LauSumPdf.cc
    \brief File containing implementation of LauSumPdf class.
*/

#include "LauSumPdf.hh"

#include "LauConstants.hh"
#include "LauParameter.hh"

#include "TMath.h"
#include "TSystem.h"

#include <iostream>
#include <vector>

LauSumPdf::LauSumPdf( LauAbsPdf* pdf1, LauAbsPdf* pdf2, LauParameter* frac ) :
    LauAbsPdf( pdf1 ? pdf1->varNames() : std::vector<TString>(),
               std::vector<LauAbsRValue*>(),
               pdf1 ? pdf1->getMinAbscissas() : LauFitData(),
               pdf1 ? pdf1->getMaxAbscissas() : LauFitData() ),
    pdf1_( pdf1 ),
    pdf2_( pdf2 ),
    frac_( frac )
{
    // Constructor for the sum PDF.
    // We are defining the sum as:
    // f x (PDF1/S(PDF1)) + (1-f) x (PDF2/S(PDF2))
    // where f is the fraction, x is multiplication, PDFi is the i'th PDF,
    // and S(PDFi) is the integral of the i'th PDF.
    //
    // The last two arguments specify the range in which the PDF is defined, and the PDF
    // will be normalised w.r.t. these limits.

    // So the first thing we have to do is check the pointers are all valid.
    if ( ! frac || ! pdf1 || ! pdf2 ) {
        std::cerr << "ERROR in LauSumPdf constructor : one of the 3 pointers is null." << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Next check that the abscissa ranges are the same for each PDF
    if ( pdf1->getMinAbscissa() != pdf2->getMinAbscissa() ) {
        std::cerr << "ERROR in LauSumPdf constructor : minimum abscissa values not the same for the two PDFs."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }
    if ( pdf1->getMaxAbscissa() != pdf2->getMaxAbscissa() ) {
        std::cerr << "ERROR in LauSumPdf constructor : maximum abscissa values not the same for the two PDFs."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Also check that both PDFs are expecting the same number of input variables
    if ( pdf1->nInputVars() != pdf2->nInputVars() ) {
        std::cerr << "ERROR in LauSumPdf constructor : number of input variables not the same for the two PDFs."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Also check that both PDFs are expecting the same variable name
    if ( pdf1->varNames() != pdf2->varNames() ) {
        std::cerr << "ERROR in LauSumPdf constructor : variable name not the same for the two PDFs."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Then we need to grab all the parameters and pass them to the base class.
    // This is so that when we are asked for them they can be put into the fit.
    // The number of parameters is the number in PDF1 + the number in PDF2 + 1 for the fraction.
    UInt_t nPar( pdf1->nParameters() + pdf2->nParameters() + 1 );
    std::vector<LauAbsRValue*> params;
    params.reserve( nPar );
    std::vector<LauAbsRValue*>& pdf1pars = pdf1->getParameters();
    std::vector<LauAbsRValue*>& pdf2pars = pdf2->getParameters();
    params.push_back( frac );
    for ( std::vector<LauAbsRValue*>::iterator iter = pdf1pars.begin(); iter != pdf1pars.end();
          ++iter ) {
        params.push_back( *iter );
    }
    for ( std::vector<LauAbsRValue*>::iterator iter = pdf2pars.begin(); iter != pdf2pars.end();
          ++iter ) {
        params.push_back( *iter );
    }
    this->addParameters( params );

    // Now check that we can find the fraction parameter ok
    frac_ = this->findParameter( "frac" );
    if ( frac_ == 0 ) {
        std::cerr << "ERROR in LauSumPdf constructor : parameter \"frac\" not found." << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Cache the normalisation factor
    this->calcNorm();
}

LauSumPdf::~LauSumPdf()
{
    // Destructor
}

void LauSumPdf::calcLikelihoodInfo( const LauAbscissas& abscissas )
{
    // Check that the given abscissa is within the allowed range
    if ( ! this->checkRange( abscissas ) ) {
        gSystem->Exit( EXIT_FAILURE );
    }

    // Get the up to date parameter values
    Double_t frac = frac_->unblindValue();

    // Evaluate the normalised PDF values
    pdf1_->calcLikelihoodInfo( abscissas );
    pdf2_->calcLikelihoodInfo( abscissas );
    Double_t result1 = pdf1_->getLikelihood();
    Double_t result2 = pdf2_->getLikelihood();

    // Add them together
    Double_t result = frac * result1 + ( 1.0 - frac ) * result2;
    this->setUnNormPDFVal( result );
}

void LauSumPdf::calcNorm()
{
    // Nothing to do here, since it is already normalized
    this->setNorm( 1.0 );
}

void LauSumPdf::calcPDFHeight( const LauKinematics* kinematics )
{
    if ( this->heightUpToDate() && ! this->isDPDependent() ) {
        return;
    }

    // This method gives you the maximum possible height of the PDF.
    // It combines the maximum heights of the two individual PDFs.
    // So it would give the true maximum if the two individual maxima coincided.
    // It is guaranteed to always be >= the true maximum.

    // Update the heights of the individual PDFs
    pdf1_->calcPDFHeight( kinematics );
    pdf2_->calcPDFHeight( kinematics );

    // Get the up to date parameter values
    Double_t frac = frac_->unblindValue();

    // Find the (un-normalised) individual PDF maxima
    Double_t height1 = pdf1_->getMaxHeight();
    Double_t height2 = pdf2_->getMaxHeight();

    // Get the individual PDF normalisation factors
    Double_t norm1 = pdf1_->getNorm();
    Double_t norm2 = pdf2_->getNorm();

    // Calculate the normalised individual PDF maxima
    height1 /= norm1;
    height2 /= norm2;

    // Combine these heights together
    Double_t height = frac * height1 + ( 1 - frac ) * height2;
    this->setMaxHeight( height );
}

// Override the base class methods for cacheInfo and calcLikelihoodInfo(UInt_t iEvt).
// For both of these we delegate to the two constituent PDFs.

void LauSumPdf::cacheInfo( const LauFitDataTree& inputData )
{
    pdf1_->cacheInfo( inputData );
    pdf2_->cacheInfo( inputData );

    // determine whether we are caching our PDF value
    Bool_t doCaching( this->nFixedParameters() == this->nParameters() );
    this->cachePDF( doCaching );

    // Set whether the height is up-to-date only if it is true for both of the sub-PDFs
    Bool_t hutd = pdf1_->heightUpToDate() && pdf2_->heightUpToDate();
    this->heightUpToDate( hutd );
}

void LauSumPdf::calcLikelihoodInfo( UInt_t iEvt )
{
    // Get the up to date parameter values
    Double_t frac = frac_->unblindValue();

    // Evaluate the normalised PDF values
    pdf1_->calcLikelihoodInfo( iEvt );
    pdf2_->calcLikelihoodInfo( iEvt );
    Double_t result1 = pdf1_->getLikelihood();
    Double_t result2 = pdf2_->getLikelihood();

    // Add them together
    Double_t result = frac * result1 + ( 1.0 - frac ) * result2;
    this->setUnNormPDFVal( result );
}
