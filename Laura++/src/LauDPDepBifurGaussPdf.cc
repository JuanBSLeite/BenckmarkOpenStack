
/*
Copyright 2007 University of Warwick

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

/*! \file LauDPDepBifurGaussPdf.cc
    \brief File containing implementation of LauDPDepBifurGaussPdf class.
*/

#include <iostream>
#include <vector>
using std::cerr;
using std::cout;
using std::endl;

#include "Lau1DHistPdf.hh"
#include "LauComplex.hh"
#include "LauConstants.hh"
#include "LauDPDepBifurGaussPdf.hh"
#include "LauDaughters.hh"
#include "LauFitDataTree.hh"
#include "LauKinematics.hh"

#include "TMath.h"
#include "TRandom.h"
#include "TSystem.h"

LauDPDepBifurGaussPdf::LauDPDepBifurGaussPdf( const TString& theVarName,
                                              const std::vector<LauAbsRValue*>& params,
                                              Double_t minAbscissa,
                                              Double_t maxAbscissa,
                                              const LauDaughters* daughters,
                                              const std::vector<Double_t>& meanCoeffs,
                                              const std::vector<Double_t>& sigmaLCoeffs,
                                              const std::vector<Double_t>& sigmaRCoeffs,
                                              DPAxis dpAxis ) :
    LauAbsPdf( theVarName, params, minAbscissa, maxAbscissa ),
    kinematics_( daughters ? daughters->getKinematics() : 0 ),
    mean_( 0 ),
    sigmaL_( 0 ),
    sigmaR_( 0 ),
    meanVal_( 0.0 ),
    sigmaLVal_( 0.0 ),
    sigmaRVal_( 0.0 ),
    meanCoeffs_( meanCoeffs ),
    sigmaLCoeffs_( sigmaLCoeffs ),
    sigmaRCoeffs_( sigmaRCoeffs ),
    dpAxis_( dpAxis ),
    scaleMethod_( poly )
{
    // Check we have a valid kinematics object
    if ( ! kinematics_ ) {
        cerr << "ERROR in LauDPDepBifurGaussPdf::LauDPDepBifurGaussPdf : Have not been provided with a valid DP kinematics object."
             << endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // The parameters are:
    // - the mean, the sigmaL and the sigmaR of the Bifurcated Gaussian
    //
    // The next two arguments specify the range in which the PDF is defined,
    // and the PDF will be normalised w.r.t. these limits.
    //
    // The final three argument define whether the BF Gaussian parameters should be scaled or not

    mean_   = this->findParameter( "mean" );
    sigmaL_ = this->findParameter( "sigmaL" );
    sigmaR_ = this->findParameter( "sigmaR" );

    if ( ( this->nParameters() != 3 ) || ( mean_ == 0 ) || ( sigmaL_ == 0 ) || ( sigmaR_ == 0 ) ) {
        cerr << "ERROR in LauDPDepBifurGaussPdf constructor: LauDPDepBifurGaussPdf requires 3 parameters:"
             << " \"mean\", \"sigmaL\" and \"sigmaR\"." << endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Cache the normalisation factor
    this->calcNorm();
}

LauDPDepBifurGaussPdf::~LauDPDepBifurGaussPdf()
{
    // Destructor
}

void LauDPDepBifurGaussPdf::calcLikelihoodInfo( const LauAbscissas& abscissas )
{
    // Check that the given abscissas are within the allowed range
    if ( ! this->checkRange( abscissas ) ) {
        gSystem->Exit( EXIT_FAILURE );
    }

    // Get our abscissa
    Double_t abscissa = abscissas[0];

    // Get the up to date parameter values
    meanVal_   = mean_->unblindValue();
    sigmaLVal_ = sigmaL_->unblindValue();
    sigmaRVal_ = sigmaR_->unblindValue();

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

    // Scale the gaussian parameters by the dpPos (if appropriate)
    ScaleMethod scale = this->scaleMethod();
    if ( scale == poly ) {
        this->scalePars_poly( dpPos );
    } else if ( scale == polyNegPower ) {
        this->scalePars_polyNegPower( dpPos );
    } else {
        cerr << "Scaling method unknown!  Methods available: <poly> and <polyNegPower>." << endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Evaluate the Birfucated Gaussian PDF value
    Double_t arg = abscissa - meanVal_;
    Double_t coef( 0.0 );
    Double_t value( 0.0 );

    if ( arg < 0.0 ) {
        if ( TMath::Abs( sigmaLVal_ ) > 1e-30 ) {
            coef = -0.5 / ( sigmaLVal_ * sigmaLVal_ );
        }
    } else {
        if ( TMath::Abs( sigmaRVal_ ) > 1e-30 ) {
            coef = -0.5 / ( sigmaRVal_ * sigmaRVal_ );
        }
    }
    value = TMath::Exp( coef * arg * arg );

    // Calculate the norm
    Double_t xscaleL = LauConstants::root2 * sigmaLVal_;
    Double_t xscaleR = LauConstants::root2 * sigmaRVal_;

    Double_t integral( 0.0 );
    Double_t norm( 0.0 );
    Double_t result( 0.0 );

    if ( this->getMaxAbscissa() < meanVal_ ) {
        integral = sigmaLVal_ * ( TMath::Erf( ( this->getMaxAbscissa() - meanVal_ ) / xscaleL ) -
                                  TMath::Erf( ( this->getMinAbscissa() - meanVal_ ) / xscaleL ) );
    } else if ( this->getMinAbscissa() > meanVal_ ) {
        integral = sigmaRVal_ * ( TMath::Erf( ( this->getMaxAbscissa() - meanVal_ ) / xscaleR ) -
                                  TMath::Erf( ( this->getMinAbscissa() - meanVal_ ) / xscaleR ) );
    } else {
        integral = sigmaRVal_ * TMath::Erf( ( this->getMaxAbscissa() - meanVal_ ) / xscaleR ) -
                   sigmaLVal_ * TMath::Erf( ( this->getMinAbscissa() - meanVal_ ) / xscaleL );
    }

    norm = LauConstants::rootPiBy2 * integral;

    // the result
    result = value / norm;
    this->setUnNormPDFVal( result );
}

void LauDPDepBifurGaussPdf::scalePars_poly( Double_t dpPos )
{
    Int_t power = 1;
    for ( std::vector<Double_t>::const_iterator iter = meanCoeffs_.begin(); iter != meanCoeffs_.end();
          ++iter ) {
        Double_t coeff  = ( *iter );
        meanVal_       += coeff * TMath::Power( dpPos, power );
        ++power;
    }

    power = 1;
    for ( std::vector<Double_t>::const_iterator iter = sigmaLCoeffs_.begin();
          iter != sigmaLCoeffs_.end();
          ++iter ) {
        Double_t coeff  = ( *iter );
        sigmaLVal_     += coeff * TMath::Power( dpPos, power );
        ++power;
    }

    power = 1;
    for ( std::vector<Double_t>::const_iterator iter = sigmaRCoeffs_.begin();
          iter != sigmaRCoeffs_.end();
          ++iter ) {
        Double_t coeff  = ( *iter );
        sigmaRVal_     += coeff * TMath::Power( dpPos, power );
        ++power;
    }
}

void LauDPDepBifurGaussPdf::scalePars_polyNegPower( Double_t dpPos )
{
    Int_t power = -1;
    for ( std::vector<Double_t>::const_iterator iter = meanCoeffs_.begin(); iter != meanCoeffs_.end();
          ++iter ) {
        Double_t coeff  = ( *iter );
        meanVal_       += coeff * TMath::Power( dpPos, power );
        --power;
    }

    power = -1;
    for ( std::vector<Double_t>::const_iterator iter = sigmaLCoeffs_.begin();
          iter != sigmaLCoeffs_.end();
          ++iter ) {
        Double_t coeff  = ( *iter );
        sigmaLVal_     += coeff * TMath::Power( dpPos, power );
        --power;
    }

    power = -1;
    for ( std::vector<Double_t>::const_iterator iter = sigmaRCoeffs_.begin();
          iter != sigmaRCoeffs_.end();
          ++iter ) {
        Double_t coeff  = ( *iter );
        sigmaRVal_     += coeff * TMath::Power( dpPos, power );
        --power;
    }
}

void LauDPDepBifurGaussPdf::calcNorm()
{
    this->setNorm( 1.0 );
}

void LauDPDepBifurGaussPdf::calcPDFHeight( const LauKinematics* kinematics )
{
    // Get the up to date parameter values
    meanVal_   = mean_->unblindValue();
    sigmaLVal_ = sigmaL_->unblindValue();
    sigmaRVal_ = sigmaR_->unblindValue();

    // Scale the gaussian parameters by the dpCentreDist (if appropriate)
    Double_t dpCentreDist = kinematics->distanceFromDPCentre();
    ScaleMethod scale     = this->scaleMethod();
    if ( scale == poly ) {
        this->scalePars_poly( dpCentreDist );
    } else if ( scale == polyNegPower ) {
        this->scalePars_polyNegPower( dpCentreDist );
    } else {
        cerr << "Scaling method unknown!  Methods available: <poly> and <polyNegPower>." << endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Calculate the PDF height for the Bifurcated Gaussian function.
    LauAbscissas maxPoint( 3 );
    maxPoint[0] = meanVal_;
    maxPoint[1] = kinematics->getm13Sq();
    maxPoint[2] = kinematics->getm23Sq();
    if ( meanVal_ > this->getMaxAbscissa() ) {
        maxPoint[0] = this->getMaxAbscissa();
    } else if ( meanVal_ < this->getMinAbscissa() ) {
        maxPoint[0] = this->getMinAbscissa();
    }
    this->calcLikelihoodInfo( maxPoint );
    Double_t height = this->getUnNormLikelihood();

    // Multiply by a small factor to avoid problems from rounding errors
    height *= 1.001;

    this->setMaxHeight( height );
}
