
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

/*! \file LauDPDepGaussPdf.cc
    \brief File containing implementation of LauDPDepGaussPdf class.
*/

#include <iostream>
#include <vector>
using std::cerr;
using std::cout;
using std::endl;

#include "LauConstants.hh"
#include "LauDPDepGaussPdf.hh"
#include "LauDaughters.hh"
#include "LauKinematics.hh"

#include "TMath.h"
#include "TSystem.h"

LauDPDepGaussPdf::LauDPDepGaussPdf( const TString& theVarName,
                                    const std::vector<LauAbsRValue*>& params,
                                    Double_t minAbscissa,
                                    Double_t maxAbscissa,
                                    const LauDaughters* daughters,
                                    const std::vector<Double_t>& meanCoeffs,
                                    const std::vector<Double_t>& sigmaCoeffs,
                                    DPAxis dpAxis ) :
    LauAbsPdf( theVarName, params, minAbscissa, maxAbscissa ),
    kinematics_( daughters ? daughters->getKinematics() : 0 ),
    mean_( 0 ),
    sigma_( 0 ),
    meanVal_( 0.0 ),
    sigmaVal_( 0.0 ),
    meanCoeffs_( meanCoeffs ),
    sigmaCoeffs_( sigmaCoeffs ),
    dpAxis_( dpAxis )
{
    // Constructor for the Gaussian PDF.

    // Check we have a valid kinematics object
    if ( ! kinematics_ ) {
        cerr << "ERROR in LauDPDepGaussPdf::LauDPDepGaussPdf : Have not been provided with a valid DP kinematics object."
             << endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // The parameters in params are the mean and the sigma (half the width) of the gaussian.
    // The last two arguments specify the range in which the PDF is defined, and the PDF
    // will be normalised w.r.t. these limits.

    mean_  = this->findParameter( "mean" );
    sigma_ = this->findParameter( "sigma" );

    if ( ( this->nParameters() != 2 ) || ( mean_ == 0 ) || ( sigma_ == 0 ) ) {
        cerr << "ERROR in LauDPDepGaussPdf constructor: LauDPDepGaussPdf requires 2 parameters: \"mean\" and \"sigma\"."
             << endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Cache the normalisation factor
    this->calcNorm();
}

LauDPDepGaussPdf::~LauDPDepGaussPdf()
{
    // Destructor
}

void LauDPDepGaussPdf::calcNorm()
{
    this->setNorm( 1.0 );
}

void LauDPDepGaussPdf::calcLikelihoodInfo( const LauAbscissas& abscissas )
{
    // Check that the given abscissa is within the allowed range
    if ( ! this->checkRange( abscissas ) ) {
        gSystem->Exit( EXIT_FAILURE );
    }

    // Get our abscissa
    Double_t abscissa = abscissas[0];

    // Get the up to date parameter values
    meanVal_  = mean_->unblindValue();
    sigmaVal_ = sigma_->unblindValue();

    // Find out the DP position
    Double_t dpPos( 0.0 );
    UInt_t nVars = this->nInputVars();
    if ( abscissas.size() == nVars + 2 ) {
        Double_t m13Sq = abscissas[nVars];
        Double_t m23Sq = abscissas[nVars + 1];

        if ( dpAxis_ == M12 ) {
            dpPos = kinematics_->calcThirdMassSq( m13Sq, m23Sq );
        } else if ( dpAxis_ == M13 ) {
            dpPos = m13Sq;
        } else if ( dpAxis_ == M23 ) {
            dpPos = m23Sq;
        } else if ( dpAxis_ == MMIN ) {
            dpPos = TMath::Min( m13Sq, m23Sq );
        } else if ( dpAxis_ == MMAX ) {
            dpPos = TMath::Max( m13Sq, m23Sq );
        } else {
            dpPos = kinematics_->distanceFromDPCentre( m13Sq, m23Sq );
        }
    }

    // Scale the parameters according to the DP position
    this->scalePars( dpPos );

    // Calculate the value of the Gaussian for the given value of the abscissa.
    Double_t arg = abscissa - meanVal_;

    Double_t exponent( 0.0 );
    if ( TMath::Abs( sigmaVal_ ) > 1e-10 ) {
        exponent = -0.5 * arg * arg / ( sigmaVal_ * sigmaVal_ );
    }

    Double_t value = TMath::Exp( exponent );

    Double_t norm( 0.0 );
    Double_t scale = LauConstants::root2 * sigmaVal_;
    if ( TMath::Abs( sigmaVal_ ) > 1e-10 ) {
        norm = LauConstants::rootPiBy2 * sigmaVal_ *
               ( TMath::Erf( ( this->getMaxAbscissa() - meanVal_ ) / scale ) -
                 TMath::Erf( ( this->getMinAbscissa() - meanVal_ ) / scale ) );
    }

    value /= norm;

    this->setUnNormPDFVal( value );
}

void LauDPDepGaussPdf::scalePars( Double_t dpPos )
{
    Int_t power = 1;
    for ( std::vector<Double_t>::const_iterator iter = meanCoeffs_.begin(); iter != meanCoeffs_.end();
          ++iter ) {
        Double_t coeff  = ( *iter );
        meanVal_       += coeff * TMath::Power( dpPos, power );
        ++power;
    }

    power = 1;
    for ( std::vector<Double_t>::const_iterator iter = sigmaCoeffs_.begin();
          iter != sigmaCoeffs_.end();
          ++iter ) {
        Double_t coeff  = ( *iter );
        sigmaVal_      += coeff * TMath::Power( dpPos, power );
        ++power;
    }
}

void LauDPDepGaussPdf::calcPDFHeight( const LauKinematics* kinematics )
{
    // Get the up to date parameter values
    meanVal_ = mean_->unblindValue();

    // Find out the DP position
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

    // Scale the parameters according to the DP position
    this->scalePars( dpPos );

    LauAbscissas maxPoint( 3 );
    maxPoint[0] = meanVal_;
    maxPoint[1] = kinematics->getm13Sq();
    maxPoint[2] = kinematics->getm23Sq();

    // Calculate the PDF height for the Gaussian function.
    if ( meanVal_ > this->getMaxAbscissa() ) {
        maxPoint[0] = this->getMaxAbscissa();
    } else if ( meanVal_ < this->getMinAbscissa() ) {
        maxPoint[0] = this->getMinAbscissa();
    }
    this->calcLikelihoodInfo( maxPoint );

    Double_t height = this->getUnNormLikelihood();
    this->setMaxHeight( height );
}
