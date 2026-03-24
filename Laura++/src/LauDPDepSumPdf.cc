
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

/*! \file LauDPDepSumPdf.cc
    \brief File containing implementation of LauDPDepSumPdf class.
*/

#include <iostream>
#include <vector>
using std::cerr;
using std::cout;
using std::endl;

#include "LauConstants.hh"
#include "LauDPDepSumPdf.hh"
#include "LauDaughters.hh"
#include "LauEffModel.hh"
#include "LauParameter.hh"

#include "TMath.h"
#include "TSystem.h"

LauDPDepSumPdf::LauDPDepSumPdf( LauAbsPdf* pdf1,
                                LauAbsPdf* pdf2,
                                const LauDaughters* daughters,
                                const TH2* dpHisto,
                                Bool_t upperHalf,
                                Bool_t useSpline ) :
    LauAbsPdf( pdf1 ? pdf1->varNames() : std::vector<TString>(),
               std::vector<LauAbsRValue*>(),
               pdf1 ? pdf1->getMinAbscissas() : LauFitData(),
               pdf1 ? pdf1->getMaxAbscissas() : LauFitData() ),
    daughters_( new LauDaughters( *daughters ) ),
    pdf1_( pdf1 ),
    pdf2_( pdf2 ),
    frac_( 0 ),
    fracVal_( 0.5 ),
    dpDependence_( new LauEffModel( daughters, 0 ) ),
    dpAxis_( CentreDist )
{
    // Constructor for the sum PDF.
    // We are defining the sum as:
    // f x (PDF1/S(PDF1)) + (1-f) x (PDF2/S(PDF2))
    // where f is the fraction, x is multiplication, PDFi is the i'th PDF,
    // and S(PDFi) is the integral of the i'th PDF.
    // The value of the fraction is read from the DP histogram.

    if ( useSpline ) {
        dpDependence_->setEffSpline( dpHisto, kFALSE, 0.0, 0.0, upperHalf, daughters->squareDP() );
    } else {
        dpDependence_->setEffHisto( dpHisto, kTRUE, kFALSE, 0.0, 0.0, upperHalf, daughters->squareDP() );
    }

    // So the first thing we have to do is check the pointers are all valid.
    if ( ! pdf1 || ! pdf2 ) {
        cerr << "ERROR in LauDPDepSumPdf constructor: one of the 2 PDF pointers is null." << endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Next check that the abscissa ranges are the same for each PDF
    if ( pdf1->getMinAbscissa() != pdf2->getMinAbscissa() ) {
        cerr << "ERROR in LauDPDepSumPdf constructor: minimum abscissa values not the same for the two PDFs."
             << endl;
        gSystem->Exit( EXIT_FAILURE );
    }
    if ( pdf1->getMaxAbscissa() != pdf2->getMaxAbscissa() ) {
        cerr << "ERROR in LauDPDepSumPdf constructor: maximum abscissa values not the same for the two PDFs."
             << endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Also check that both PDFs are expecting the same number of input variables
    if ( pdf1->nInputVars() != pdf2->nInputVars() ) {
        cerr << "ERROR in LauDPDepSumPdf constructor: number of input variables not the same for the two PDFs."
             << endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Also check that both PDFs are expecting the same variable name(s)
    if ( pdf1->varNames() != pdf2->varNames() ) {
        cerr << "ERROR in LauDPDepSumPdf constructor: variable name(s) not the same for the two PDFs."
             << endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Then we need to grab all the parameters and pass them to the base class.
    // This is so that when we are asked for them they can be put into the fit.
    // The number of parameters is the number in PDF1 + the number in PDF2.
    UInt_t nPar( pdf1->nParameters() + pdf2->nParameters() );
    std::vector<LauAbsRValue*> params;
    params.reserve( nPar );
    std::vector<LauAbsRValue*>& pdf1pars = pdf1->getParameters();
    std::vector<LauAbsRValue*>& pdf2pars = pdf2->getParameters();
    for ( std::vector<LauAbsRValue*>::iterator iter = pdf1pars.begin(); iter != pdf1pars.end();
          ++iter ) {
        params.push_back( *iter );
    }
    for ( std::vector<LauAbsRValue*>::iterator iter = pdf2pars.begin(); iter != pdf2pars.end();
          ++iter ) {
        params.push_back( *iter );
    }
    this->addParameters( params );

    // Cache the normalisation factor
    this->calcNorm();
}

LauDPDepSumPdf::LauDPDepSumPdf( LauAbsPdf* pdf1,
                                LauAbsPdf* pdf2,
                                LauParameter* frac,
                                const LauDaughters* daughters,
                                const std::vector<Double_t>& fracCoeffs,
                                DPAxis dpAxis ) :
    LauAbsPdf( pdf1 ? pdf1->varNames() : std::vector<TString>(),
               std::vector<LauAbsRValue*>(),
               pdf1 ? pdf1->getMinAbscissas() : LauFitData(),
               pdf1 ? pdf1->getMaxAbscissas() : LauFitData() ),
    daughters_( new LauDaughters( *daughters ) ),
    pdf1_( pdf1 ),
    pdf2_( pdf2 ),
    frac_( frac ),
    fracVal_( frac ? frac->unblindValue() : 0.0 ),
    dpDependence_( 0 ),
    fracCoeffs_( fracCoeffs ),
    dpAxis_( dpAxis )
{
    // Constructor for the sum PDF.
    // We are defining the sum as:
    // f x (PDF1/S(PDF1)) + (1-f) x (PDF2/S(PDF2))
    // where f is the fraction, x is multiplication, PDFi is the i'th PDF,
    // and S(PDFi) is the integral of the i'th PDF.
    // The value of the fraction has a polynomial dependence on one of
    // the DP axes.

    // So the first thing we have to do is check the pointers are all valid.
    if ( ! pdf1 || ! pdf2 ) {
        cerr << "ERROR in LauDPDepSumPdf constructor: one of the 2 PDF pointers is null." << endl;
        gSystem->Exit( EXIT_FAILURE );
    }
    if ( ! frac ) {
        cerr << "ERROR in LauDPDepSumPdf constructor: the fraction parameter pointer is null."
             << endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Next check that the abscissa ranges are the same for each PDF
    if ( pdf1->getMinAbscissa() != pdf2->getMinAbscissa() ) {
        cerr << "ERROR in LauDPDepSumPdf constructor: minimum abscissa values not the same for the two PDFs."
             << endl;
        gSystem->Exit( EXIT_FAILURE );
    }
    if ( pdf1->getMaxAbscissa() != pdf2->getMaxAbscissa() ) {
        cerr << "ERROR in LauDPDepSumPdf constructor: maximum abscissa values not the same for the two PDFs."
             << endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Also check that both PDFs are expecting the same number of input variables
    if ( pdf1->nInputVars() != pdf2->nInputVars() ) {
        cerr << "ERROR in LauDPDepSumPdf constructor: number of input variables not the same for the two PDFs."
             << endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Also check that both PDFs are expecting the same variable name(s)
    if ( pdf1->varNames() != pdf2->varNames() ) {
        cerr << "ERROR in LauDPDepSumPdf constructor: variable name(s) not the same for the two PDFs."
             << endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Then we need to grab all the parameters and pass them to the base class.
    // This is so that when we are asked for them they can be put into the fit.
    // The number of parameters is the number in PDF1 + the number in PDF2.
    UInt_t nPar( pdf1->nParameters() + pdf2->nParameters() + 1 );
    std::vector<LauAbsRValue*> params;
    params.reserve( nPar );
    params.push_back( frac );
    std::vector<LauAbsRValue*>& pdf1pars = pdf1->getParameters();
    std::vector<LauAbsRValue*>& pdf2pars = pdf2->getParameters();
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
        cerr << "ERROR in LauDPDepSumPdf constructor: parameter \"frac\" not found." << endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Cache the normalisation factor
    this->calcNorm();
}

LauDPDepSumPdf::~LauDPDepSumPdf()
{
    // Destructor
    delete daughters_;
    daughters_ = 0;
    delete dpDependence_;
    dpDependence_ = 0;
}

void LauDPDepSumPdf::calcLikelihoodInfo( const LauAbscissas& abscissas )
{
    // Check that the given abscissa is within the allowed range
    if ( ! this->checkRange( abscissas ) ) {
        gSystem->Exit( EXIT_FAILURE );
    }

    LauAbscissas noDPVars( 1 );
    noDPVars[0] = abscissas[0];

    // Evaluate the normalised PDF values
    if ( pdf1_->isDPDependent() ) {
        pdf1_->calcLikelihoodInfo( abscissas );
    } else {
        pdf1_->calcLikelihoodInfo( noDPVars );
    }
    if ( pdf2_->isDPDependent() ) {
        pdf2_->calcLikelihoodInfo( abscissas );
    } else {
        pdf2_->calcLikelihoodInfo( noDPVars );
    }
    Double_t result1 = pdf1_->getLikelihood();
    Double_t result2 = pdf2_->getLikelihood();

    // Determine the fraction
    // The DP variables will be abscissas[nInputVars] and
    // abscissas[nInputVars+1] (if present).
    UInt_t nVars = this->nInputVars();
    if ( abscissas.size() == nVars + 2 ) {
        Double_t m13Sq            = abscissas[nVars];
        Double_t m23Sq            = abscissas[nVars + 1];
        LauKinematics* kinematics = daughters_->getKinematics();
        if ( dpDependence_ ) {
            kinematics->updateKinematics( m13Sq, m23Sq );
            fracVal_ = dpDependence_->calcEfficiency( kinematics );
        } else {
            fracVal_ = frac_->unblindValue();
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
            this->scaleFrac( dpPos );
        }
    }

    // Add them together
    Double_t result = fracVal_ * result1 + ( 1.0 - fracVal_ ) * result2;
    this->setUnNormPDFVal( result );
}

void LauDPDepSumPdf::scaleFrac( Double_t dpPos )
{
    Int_t power = 1;
    for ( std::vector<Double_t>::const_iterator iter = fracCoeffs_.begin(); iter != fracCoeffs_.end();
          ++iter ) {
        Double_t coeff  = ( *iter );
        fracVal_       += coeff * TMath::Power( dpPos, power );
        ++power;
    }
}

void LauDPDepSumPdf::calcNorm()
{
    // Nothing to do here, since it is already normalized
    this->setNorm( 1.0 );
}

void LauDPDepSumPdf::calcPDFHeight( const LauKinematics* kinematics )
{
    // This method gives you the maximum possible height of the PDF.
    // It combines the maximum heights of the two individual PDFs.
    // So it would give the true maximum if the two individual maxima coincided.
    // It is guaranteed to always be >= the true maximum.

    // Update the heights of the individual PDFs
    pdf1_->calcPDFHeight( kinematics );
    pdf2_->calcPDFHeight( kinematics );

    // Get the up to date parameter values
    if ( dpDependence_ ) {
        fracVal_ = dpDependence_->calcEfficiency( kinematics );
    } else {
        fracVal_ = frac_->unblindValue();
        Double_t dpPos( 0.0 );
        if ( dpAxis_ == M12 ) {
            dpPos = kinematics->getm12Sq();
        } else if ( dpAxis_ == M13 ) {
            dpPos = kinematics->getm13Sq();
        } else if ( dpAxis_ == M23 ) {
            dpPos = kinematics->getm23Sq();
        } else if ( dpAxis_ == MMIN ) {
            Double_t m13Sq = kinematics->getm13Sq();
            Double_t m23Sq = kinematics->getm23Sq();
            dpPos          = TMath::Min( m13Sq, m23Sq );
        } else if ( dpAxis_ == MMAX ) {
            Double_t m13Sq = kinematics->getm13Sq();
            Double_t m23Sq = kinematics->getm23Sq();
            dpPos          = TMath::Max( m13Sq, m23Sq );
        } else {
            dpPos = kinematics->distanceFromDPCentre();
        }
        this->scaleFrac( dpPos );
    }

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
    Double_t height = fracVal_ * height1 + ( 1 - fracVal_ ) * height2;
    this->setMaxHeight( height );
}

// Override the base class methods for cacheInfo and calcLikelihoodInfo(UInt_t iEvt).
// For both of these we delegate to the two constituent PDFs.

void LauDPDepSumPdf::cacheInfo( const LauFitDataTree& inputData )
{
    // delegate to the two sub-PDFs to cache their information
    pdf1_->cacheInfo( inputData );
    pdf2_->cacheInfo( inputData );

    // now check that the DP variables are included in the data
    Bool_t hasBranch  = inputData.haveBranch( "m13Sq" );
    hasBranch        &= inputData.haveBranch( "m23Sq" );
    if ( ! hasBranch ) {
        cerr << "ERROR in LauDPDepSumPdf::cacheInfo : Input data does not contain Dalitz plot variables m13Sq and/or m23Sq."
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
    fractions_.clear();
    fractions_.reserve( nEvents );

    // loop through the events, determine the fraction and store
    for ( UInt_t iEvt = 0; iEvt < nEvents; iEvt++ ) {

        const LauFitData& dataValues = inputData.getData( iEvt );

        Double_t m13Sq = dataValues.find( "m13Sq" )->second;
        Double_t m23Sq = dataValues.find( "m23Sq" )->second;

        LauKinematics* kinematics = daughters_->getKinematics();
        if ( dpDependence_ ) {
            // if we're using the histogram then just
            // determine the fraction and store
            kinematics->updateKinematics( m13Sq, m23Sq );
            fracVal_ = dpDependence_->calcEfficiency( kinematics );
        } else {
            // if we're scaling the fraction parameter then we
            // just store the scaling info since the parameter
            // might be floating
            fracVal_ = frac_->unblindValue();
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
            this->scaleFrac( dpPos );
            fracVal_ -= frac_->unblindValue();
        }

        fractions_.push_back( fracVal_ );
    }
}

void LauDPDepSumPdf::calcLikelihoodInfo( UInt_t iEvt )
{
    // Get the fraction value for this event
    fracVal_ = fractions_[iEvt];
    if ( frac_ ) {
        // if we're scaling the parameter then need to add the
        // current value of the parameter
        fracVal_ += frac_->unblindValue();
    }

    // Evaluate the normalised PDF values
    pdf1_->calcLikelihoodInfo( iEvt );
    pdf2_->calcLikelihoodInfo( iEvt );
    Double_t result1 = pdf1_->getLikelihood();
    Double_t result2 = pdf2_->getLikelihood();

    // Add them together
    Double_t result = fracVal_ * result1 + ( 1 - fracVal_ ) * result2;
    this->setUnNormPDFVal( result );
}
