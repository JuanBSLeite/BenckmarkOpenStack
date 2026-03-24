
/*
Copyright 2013 University of Warwick

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

/*! \file LauDPDepMapPdf.cc
    \brief File containing implementation of LauDPDepMapPdf class.
*/

#include <iostream>
#include <vector>
using std::cerr;
using std::cout;
using std::endl;

#include "Lau2DHistDP.hh"
#include "LauConstants.hh"
#include "LauDPDepMapPdf.hh"
#include "LauDaughters.hh"
#include "LauParameter.hh"

#include "TH1.h"
#include "TMath.h"
#include "TSystem.h"

LauDPDepMapPdf::LauDPDepMapPdf( const std::vector<LauAbsPdf*>& pdfs,
                                const LauDaughters* daughters,
                                const TH2* dpHisto,
                                Bool_t upperHalf ) :
    LauAbsPdf( ( ! pdfs.empty() && pdfs[0] ) ? pdfs[0]->varNames() : std::vector<TString>(),
               std::vector<LauAbsRValue*>(),
               ( ! pdfs.empty() && pdfs[0] ) ? pdfs[0]->getMinAbscissas() : LauFitData(),
               ( ! pdfs.empty() && pdfs[0] ) ? pdfs[0]->getMaxAbscissas() : LauFitData() ),
    daughters_( new LauDaughters( *daughters ) ),
    pdfs_( pdfs ),
    dpDependence_(
        new Lau2DHistDP( dpHisto, daughters, kFALSE, kFALSE, 0, 0, upperHalf, daughters->squareDP() ) ),
    dpAxisDependence_( 0 ),
    dpAxis_( CentreDist ),
    indices_( 0 )
{
    // Constructor for the PDF map.
    // The index into the PDF collection is read from the DP histogram.

    // So the first thing we have to do is check the pointers are all valid.
    for ( std::vector<LauAbsPdf*>::const_iterator iter = pdfs_.begin(); iter != pdfs_.end(); ++iter ) {

        if ( ! ( *iter ) ) {
            cerr << "ERROR in LauDPDepMapPdf constructor: one of the PDF pointers is null." << endl;
            gSystem->Exit( EXIT_FAILURE );
        }

        // Next check that the abscissa ranges are the same for each PDF
        if ( pdfs_[0]->getMinAbscissa() != ( *iter )->getMinAbscissa() ) {
            cerr << "ERROR in LauDPDepMapPdf constructor: minimum abscissa values not the same for two PDFs."
                 << endl;
            gSystem->Exit( EXIT_FAILURE );
        }
        if ( pdfs_[0]->getMaxAbscissa() != ( *iter )->getMaxAbscissa() ) {
            cerr << "ERROR in LauDPDepMapPdf constructor: maximum abscissa values not the same for two PDFs."
                 << endl;
            gSystem->Exit( EXIT_FAILURE );
        }

        // Also check that both PDFs are expecting the same number of input variables
        if ( pdfs_[0]->nInputVars() != ( *iter )->nInputVars() ) {
            cerr << "ERROR in LauDPDepMapPdf constructor: number of input variables not the same for two PDFs."
                 << endl;
            gSystem->Exit( EXIT_FAILURE );
        }

        // Also check that both PDFs are expecting the same variable name(s)
        if ( pdfs_[0]->varNames() != ( *iter )->varNames() ) {
            cerr << "ERROR in LauDPDepMapPdf constructor: variable name(s) not the same for two PDFs."
                 << endl;
            gSystem->Exit( EXIT_FAILURE );
        }
    }

    // Then we need to grab all the parameters and pass them to the base class.
    // This is so that when we are asked for them they can be put into the fit.

    for ( std::vector<LauAbsPdf*>::const_iterator iter = pdfs_.begin(); iter != pdfs_.end(); ++iter ) {
        std::vector<LauAbsRValue*>& pdfpars = ( *iter )->getParameters();
        this->addParameters( pdfpars );
    }

    // Cache the normalisation factor
    this->calcNorm();
}

LauDPDepMapPdf::LauDPDepMapPdf( const std::vector<LauAbsPdf*>& pdfs,
                                const LauDaughters* daughters,
                                const TH1* dpAxisHisto,
                                DPAxis dpAxis ) :
    LauAbsPdf( ( ! pdfs.empty() && pdfs[0] ) ? pdfs[0]->varNames() : std::vector<TString>(),
               std::vector<LauAbsRValue*>(),
               ( ! pdfs.empty() && pdfs[0] ) ? pdfs[0]->getMinAbscissas() : LauFitData(),
               ( ! pdfs.empty() && pdfs[0] ) ? pdfs[0]->getMaxAbscissas() : LauFitData() ),
    daughters_( new LauDaughters( *daughters ) ),
    pdfs_( pdfs ),
    dpDependence_( 0 ),
    dpAxisDependence_( dpAxisHisto ? dynamic_cast<TH1*>( dpAxisHisto->Clone() ) : 0 ),
    dpAxis_( dpAxis ),
    indices_( 0 )
{
    // Constructor for the PDFs map.
    // The index into the PDF collection is read from one of
    // the DP axes.

    if ( dpAxisDependence_ ) {
        dpAxisDependence_->SetDirectory( 0 );
    }

    // So the first thing we have to do is check the pointers are all valid.
    for ( std::vector<LauAbsPdf*>::const_iterator iter = pdfs_.begin(); iter != pdfs_.end(); ++iter ) {

        if ( ! ( *iter ) ) {
            cerr << "ERROR in LauDPDepMapPdf constructor: one of the PDF pointers is null." << endl;
            gSystem->Exit( EXIT_FAILURE );
        }

        // Next check that the abscissa ranges are the same for each PDF
        if ( pdfs_[0]->getMinAbscissa() != ( *iter )->getMinAbscissa() ) {
            cerr << "ERROR in LauDPDepMapPdf constructor: minimum abscissa values not the same for two PDFs."
                 << endl;
            gSystem->Exit( EXIT_FAILURE );
        }
        if ( pdfs_[0]->getMaxAbscissa() != ( *iter )->getMaxAbscissa() ) {
            cerr << "ERROR in LauDPDepMapPdf constructor: maximum abscissa values not the same for two PDFs."
                 << endl;
            gSystem->Exit( EXIT_FAILURE );
        }

        // Also check that both PDFs are expecting the same number of input variables
        if ( pdfs_[0]->nInputVars() != ( *iter )->nInputVars() ) {
            cerr << "ERROR in LauDPDepMapPdf constructor: number of input variables not the same for two PDFs."
                 << endl;
            gSystem->Exit( EXIT_FAILURE );
        }

        // Also check that both PDFs are expecting the same variable name(s)
        if ( pdfs_[0]->varNames() != ( *iter )->varNames() ) {
            cerr << "ERROR in LauDPDepMapPdf constructor: variable name(s) not the same for two PDFs."
                 << endl;
            gSystem->Exit( EXIT_FAILURE );
        }
    }

    // Then we need to grab all the parameters and pass them to the base class.
    // This is so that when we are asked for them they can be put into the fit.

    for ( std::vector<LauAbsPdf*>::const_iterator iter = pdfs_.begin(); iter != pdfs_.end(); ++iter ) {
        std::vector<LauAbsRValue*>& pdfpars = ( *iter )->getParameters();
        this->addParameters( pdfpars );
    }

    // Cache the normalisation factor
    this->calcNorm();
}

LauDPDepMapPdf::~LauDPDepMapPdf()
{
    // Destructor
    delete daughters_;
    daughters_ = 0;
    delete dpDependence_;
    dpDependence_ = 0;
    delete dpAxisDependence_;
    dpAxisDependence_ = 0;
}

UInt_t LauDPDepMapPdf::determineDPRegion( Double_t m13Sq, Double_t m23Sq ) const
{
    UInt_t regionIndex( 0 );

    LauKinematics* kinematics = daughters_->getKinematics();
    kinematics->updateKinematics( m13Sq, m23Sq );

    if ( dpDependence_ ) {
        if ( daughters_->squareDP() ) {
            Double_t mprime  = kinematics->getmPrime();
            Double_t thprime = kinematics->getThetaPrime();
            regionIndex = static_cast<UInt_t>( dpDependence_->interpolateXY( mprime, thprime ) );
        } else {
            regionIndex = static_cast<UInt_t>( dpDependence_->interpolateXY( m13Sq, m23Sq ) );
        }
    } else if ( dpAxisDependence_ ) {
        Double_t dpPos( 0.0 );
        if ( dpAxis_ == M12 ) {
            dpPos = kinematics->calcThirdMassSq( m13Sq, m23Sq );
        } else if ( dpAxis_ == M13 ) {
            dpPos = m13Sq;
        } else if ( dpAxis_ == M23 ) {
            dpPos = m23Sq;
        } else if ( dpAxis_ == MMIN ) {
            dpPos = TMath::Min( m13Sq, m23Sq );
        } else if ( dpAxis_ == MMAX ) {
            dpPos = TMath::Max( m13Sq, m23Sq );
        } else {
            dpPos = kinematics->distanceFromDPCentre( m13Sq, m23Sq );
        }
        Int_t bin   = dpAxisDependence_->FindFixBin( dpPos );
        regionIndex = static_cast<UInt_t>( dpAxisDependence_->GetBinContent( bin ) );
    } else {
        // This should never happen
        cerr << "ERROR in LauDPDepMapPdf::determineDPRegion : No means of determining the region!"
             << endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    return regionIndex;
}

void LauDPDepMapPdf::calcLikelihoodInfo( const LauAbscissas& abscissas )
{
    // Check that the given abscissa is within the allowed range
    if ( ! this->checkRange( abscissas ) ) {
        gSystem->Exit( EXIT_FAILURE );
    }

    // Create a set of absissas that doesn't contain the DP variables
    UInt_t nVars = this->nInputVars();
    LauAbscissas noDPVars( nVars );
    for ( UInt_t i( 0 ); i < nVars; ++i ) {
        noDPVars[i] = abscissas[i];
    }

    // Find which region of the DP we're in
    // The DP variables will be abscissas[nInputVars] and
    // abscissas[nInputVars+1] (if present).
    UInt_t regionIndex( 0 );
    if ( abscissas.size() == nVars + 2 ) {
        Double_t m13Sq = abscissas[nVars];
        Double_t m23Sq = abscissas[nVars + 1];
        regionIndex    = this->determineDPRegion( m13Sq, m23Sq );
    } else {
        // This should never happen
        cerr << "ERROR in LauDPDepMapPdf::calcLikelihoodInfo : DP vars not supplied, no means of determining the region!"
             << endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Check that the region index is valid
    if ( regionIndex >= pdfs_.size() ) {
        cerr << "ERROR in LauDPDepMapPdf::calcLikelihoodInfo : No PDF supplied for region "
             << regionIndex << endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Evaluate the normalised PDF values
    LauAbsPdf* pdf = pdfs_[regionIndex];

    if ( pdf->isDPDependent() ) {
        pdf->calcLikelihoodInfo( abscissas );
    } else {
        pdf->calcLikelihoodInfo( noDPVars );
    }

    Double_t result = pdf->getUnNormLikelihood();
    this->setUnNormPDFVal( result );

    Double_t norm = pdf->getNorm();
    this->setNorm( norm );
}

void LauDPDepMapPdf::calcNorm()
{
    // Nothing to do here, since it is already normalized
    this->setNorm( 1.0 );
}

void LauDPDepMapPdf::calcPDFHeight( const LauKinematics* kinematics )
{
    // This method gives you the maximum possible height of the PDF.
    // It combines the maximum heights of the two individual PDFs.
    // So it would give the true maximum if the two individual maxima coincided.
    // It is guaranteed to always be >= the true maximum.

    Double_t m13Sq = kinematics->getm13Sq();
    Double_t m23Sq = kinematics->getm23Sq();

    // Find which region of the DP we're in
    UInt_t regionIndex = this->determineDPRegion( m13Sq, m23Sq );

    // Check that the region index is valid
    if ( regionIndex >= pdfs_.size() ) {
        cerr << "ERROR in LauDPDepMapPdf::calcPDFHeight : No PDF supplied for region "
             << regionIndex << endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Evaluate the normalised PDF values
    LauAbsPdf* pdf = pdfs_[regionIndex];

    // Update the heights of the individual PDFs
    pdf->calcPDFHeight( kinematics );

    // Find the PDF maximum
    Double_t height = pdf->getMaxHeight();
    this->setMaxHeight( height );
}

// Override the base class methods for cacheInfo and calcLikelihoodInfo(UInt_t iEvt).
// For both of these we delegate some functionality to the constituent PDFs.

void LauDPDepMapPdf::cacheInfo( const LauFitDataTree& inputData )
{
    // delegate to the sub-PDFs to cache their information
    for ( std::vector<LauAbsPdf*>::const_iterator iter = pdfs_.begin(); iter != pdfs_.end(); ++iter ) {
        ( *iter )->cacheInfo( inputData );
    }

    // now check that the DP variables are included in the data
    Bool_t hasBranch  = inputData.haveBranch( "m13Sq" );
    hasBranch        &= inputData.haveBranch( "m23Sq" );
    if ( ! hasBranch ) {
        cerr << "ERROR in LauDPDepMapPdf::cacheInfo : Input data does not contain Dalitz plot variables m13Sq and/or m23Sq."
             << endl;
        return;
    }

    // determine whether we are caching our PDF value
    Bool_t doCaching( this->nFixedParameters() == this->nParameters() );
    this->cachePDF( doCaching );

    if ( ! doCaching ) {
        // in this case we seem to be doing a fit where the parameters are floating
        // so need to mark that the PDF height is no longer up to date
        this->heightUpToDate( kFALSE );
    }

    // clear the vector and reserve enough space
    UInt_t nEvents = inputData.nEvents();
    indices_.clear();
    indices_.reserve( nEvents );

    // loop through the events, determine the fraction and store
    for ( UInt_t iEvt = 0; iEvt < nEvents; ++iEvt ) {

        const LauFitData& dataValues = inputData.getData( iEvt );

        Double_t m13Sq = dataValues.find( "m13Sq" )->second;
        Double_t m23Sq = dataValues.find( "m23Sq" )->second;

        UInt_t regionIndex = this->determineDPRegion( m13Sq, m23Sq );
        indices_.push_back( regionIndex );
    }
}

void LauDPDepMapPdf::calcLikelihoodInfo( UInt_t iEvt )
{
    // Get the fraction value for this event
    UInt_t regionIndex = indices_[iEvt];

    // Evaluate the normalised PDF value
    LauAbsPdf* pdf = pdfs_[regionIndex];
    pdf->calcLikelihoodInfo( iEvt );

    Double_t result = pdf->getUnNormLikelihood();
    this->setUnNormPDFVal( result );

    Double_t norm = pdf->getNorm();
    this->setNorm( norm );
}
