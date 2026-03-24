
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

/*! \file LauMinuit.cc
    \brief File containing implementation of LauMinuit methods.
*/

#include "LauMinuit.hh"

#include "LauFitObject.hh"
#include "LauFitter.hh"
#include "LauParamFixed.hh"
#include "LauParameter.hh"

#include "TMath.h"
#include "TMatrixD.h"
#include "TVirtualFitter.h"

#include <algorithm>
#include <array>
#include <iostream>

// It's necessary to define an external function that specifies the address of the function
// that Minuit needs to minimise. Minuit doesn't know about any classes - therefore
// use gMinuit->SetFCN(external_function), gMinuit->SetObjectFit(this).
// Here, we use TVirtualFitter* fitter instead of gMinuit, defined below.
// Then, within the external function, invoke an object from this class (LauAllModel),
// and use the member functions to access the parameters/variables.
extern void logLikeFun( Int_t& npar, Double_t* gin, Double_t& f, Double_t* par, Int_t iflag );

LauMinuit::LauMinuit( const UInt_t maxPar, const LauOutputLevel verbosity ) :
    LauAbsFitter(),
    maxPar_ { maxPar },
    outputLevel_ { verbosity }
{
    TVirtualFitter::SetDefaultFitter( "Minuit" );
    minuit_ = TVirtualFitter::Fitter( nullptr, maxPar_ );

    // Set the printout level
    std::array<Double_t, 1> argL { static_cast<Double_t>( outputLevel_ ) };
    minuit_->ExecuteCommand( "SET PRINT", argL.data(), argL.size() );
    if ( outputLevel_ == LauOutputLevel::None ) {
        minuit_->ExecuteCommand( "SET NOWARNINGS", argL.data(), 0 );
    }
}

void LauMinuit::initialise( LauFitObject* fitObj, const std::vector<LauParameter*>& parameters )
{
    // Check whether we're going to use asymmetric errors
    if ( outputLevel_ > LauOutputLevel::Quiet ) {
        if ( useAsymmFitErrors_ == kTRUE ) {
            std::cout << "INFO in LauMinuit::fit : We are going to calculate the asymmetric fit errors."
                      << std::endl;
            std::cout << "                       : This will, in general, significantly increase the CPU time required for fitting."
                      << std::endl;
        } else {
            std::cout << "INFO in LauMinuit::fit : Only parabolic errors will be calculated."
                      << std::endl;
        }
    }

    // Store the parameters
    params_ = parameters;

    // Hook the external likelihood function to this LauFitter::fitter() and this class.
    minuit_->SetFCN( logLikeFun );
    minuit_->SetObjectFit( fitObj );

    // Clear any stored parameters etc... before using
    minuit_->Clear();

    nParams_ = params_.size();
    if ( outputLevel_ > LauOutputLevel::Quiet ) {
        std::cout << "INFO in LauMinuit::initialise : Setting fit parameters" << std::endl;
        std::cout << "                              : Total number of parameters = " << nParams_
                  << std::endl;
    }

    // Define the default relative error
    const Double_t defaultError( 0.01 );

    // Set-up the parameters
    for ( UInt_t i = 0; i < nParams_; ++i ) {
        TString name     = params_[i]->name();
        Double_t initVal = params_[i]->initValue();
        Double_t initErr = params_[i]->error();
        // If we do not have a supplied estimate of the error, we should make a reasonable guess
        if ( initErr == 0.0 ) {
            if ( initVal == 0.0 ) {
                initErr = defaultError;
            } else if ( TMath::Abs( initErr / initVal ) < 1e-6 ) {
                initErr = TMath::Abs( defaultError * initVal );
            }
        }
        Double_t minVal    = params_[i]->minValue();
        Double_t maxVal    = params_[i]->maxValue();
        Bool_t secondStage = params_[i]->secondStage();
        if ( this->twoStageFit() && secondStage == kTRUE ) {
            params_[i]->fixed( kTRUE );
        }
        Bool_t fixVar = params_[i]->fixed();

        if ( outputLevel_ > LauOutputLevel::Quiet ) {
            std::cout << "                              : Setting parameter " << i << " called "
                      << name << " to have initial value " << initVal << ", error " << initErr
                      << " and range " << minVal << " to " << maxVal << std::endl;
        }
        minuit_->SetParameter( i, name, initVal, initErr, minVal, maxVal );

        // Fix parameter if required
        if ( fixVar == kTRUE ) {
            if ( outputLevel_ > LauOutputLevel::Quiet ) {
                std::cout << "                              : Fixing parameter " << i << std::endl;
            }
            minuit_->FixParameter( i );
        }
    }

    LauParamFixed pred;
    nFreeParams_ = nParams_ - std::count_if( params_.begin(), params_.end(), pred );

    // Need to set the "SET ERR" command to +0.5 for +/-1 sigma errors
    // for maximum likelihood fit. Very important command, otherwise all
    // extracted errors will be too big, and pull distributions will be too narrow!
    // TODO - The alternative to this is to make FCN = -2log(L) rather than -log(L)

    std::array<Double_t, 1> argL { 0.5 };
    fitStatus_.status = minuit_->ExecuteCommand( "SET ERR", argL.data(), argL.size() );

    //argL[0] = 0;
    //fitStatus_.status = minuit_->ExecuteCommand("SET STRATEGY", argL.data(), argL.size());
}

LauFitObject* LauMinuit::getFitObject()
{
    return ( minuit_ != nullptr ) ? dynamic_cast<LauFitObject*>( minuit_->GetObjectFit() ) : nullptr;
}

const LauAbsFitter::FitStatus& LauMinuit::minimise()
{
    std::array<Double_t, 2> arglist {
        1000.0 * nParams_,    // maximum iterations
        0.05                  // tolerance -> min EDM = 0.001*tolerance (0.05)
    };
    fitStatus_.status = minuit_->ExecuteCommand( "MIGRAD", arglist.data(), arglist.size() );

    // Dummy variables - need to feed them to the function
    // used for getting NLL, EDM and error matrix status
    Double_t errdef;
    Int_t nvpar, nparx;

    if ( fitStatus_.status != 0 ) {

        if ( outputLevel_ > LauOutputLevel::None ) {
            std::cerr << "ERROR in LauMinuit::minimise : Error in minimising loglike." << std::endl;
        }

    } else {

        // Check that the error matrix is ok
        fitStatus_.status = minuit_->GetStats( fitStatus_.NLL, fitStatus_.EDM, errdef, nvpar, nparx );
        if ( outputLevel_ > LauOutputLevel::Quiet ) {
            std::cout << "INFO in LauMinuit::minimise : Error matrix status after MIGRAD is: "
                      << fitStatus_.status << std::endl;
        }
        // 0= not calculated at all
        // 1= approximation only, not accurate
        // 2= full matrix, but forced positive-definite
        // 3= full accurate covariance matrix

        // Fit result was OK. Now get the more precise errors.
        fitStatus_.status = minuit_->ExecuteCommand( "HESSE", arglist.data(), 1 );

        if ( fitStatus_.status != 0 ) {

            if ( outputLevel_ > LauOutputLevel::None ) {
                std::cerr << "ERROR in LauMinuit::minimise : Error in HESSE routine." << std::endl;
            }

        } else {

            // Check that the error matrix is ok
            fitStatus_.status =
                minuit_->GetStats( fitStatus_.NLL, fitStatus_.EDM, errdef, nvpar, nparx );
            if ( outputLevel_ > LauOutputLevel::Quiet ) {
                std::cout << "INFO in LauMinuit::minimise : Error matrix status after HESSE is: "
                          << fitStatus_.status << std::endl;
            }
            // 0= not calculated at all
            // 1= approximation only, not accurate
            // 2= full matrix, but forced positive-definite
            // 3= full accurate covariance matrix

            // Symmetric errors and eror matrix were OK.
            // Get asymmetric errors if asked for.
            if ( useAsymmFitErrors_ == kTRUE ) {
                LauFitObject* fitObj = this->getFitObject();
                fitObj->withinAsymErrorCalc( kTRUE );
                fitStatus_.status = minuit_->ExecuteCommand( "MINOS", arglist.data(), 1 );
                fitObj->withinAsymErrorCalc( kFALSE );
                if ( fitStatus_.status != 0 ) {
                    std::cerr << "ERROR in LauMinuit::minimise : Error in MINOS routine."
                              << std::endl;
                }
            }
        }
    }

    // Print results
    fitStatus_.status = minuit_->GetStats( fitStatus_.NLL, fitStatus_.EDM, errdef, nvpar, nparx );
    if ( outputLevel_ > LauOutputLevel::None ) {
        std::cout << "INFO in LauMinuit::minimise : Final error matrix status is: "
                  << fitStatus_.status << std::endl;
    }
    // 0= not calculated at all
    // 1= approximation only, not accurate
    // 2= full matrix, but forced positive-definite
    // 3= full accurate covariance matrix

    if ( outputLevel_ > LauOutputLevel::None ) {
        minuit_->PrintResults( 3, fitStatus_.NLL );
    }

    // Retrieve the covariance matrix from the fitter
    // For some reason the array returned is as if the matrix is of dimension nParams_ x nParams_
    // but only the elements within the sub-matrix nFreeParams_ x nFreeParams_ have values,
    // the "trailing" elements are zero, so we trim them off.
    Double_t* covMatrix = minuit_->GetCovarianceMatrix();
    covMatrix_.Clear();
    covMatrix_.ResizeTo( nParams_, nParams_ );
    covMatrix_.SetMatrixArray( covMatrix );
    covMatrix_.ResizeTo( nFreeParams_, nFreeParams_ );

    return fitStatus_;
}

void LauMinuit::fixSecondStageParameters()
{
    for ( UInt_t i { 0 }; i < nParams_; ++i ) {
        if ( params_[i]->secondStage() ) {
            params_[i]->fixed( kTRUE );
            minuit_->FixParameter( i );
        }
    }
    LauParamFixed pred;
    nFreeParams_ = nParams_ - std::count_if( params_.begin(), params_.end(), pred );
}

void LauMinuit::releaseSecondStageParameters()
{
    for ( UInt_t i { 0 }; i < nParams_; ++i ) {
        if ( params_[i]->secondStage() ) {
            params_[i]->fixed( kFALSE );
            minuit_->ReleaseParameter( i );
        }
    }
    LauParamFixed pred;
    nFreeParams_ = nParams_ - std::count_if( params_.begin(), params_.end(), pred );
}

void LauMinuit::updateParameters()
{
    for ( UInt_t i { 0 }; i < nParams_; ++i ) {
        // Get the value and errors from MINUIT
        Double_t value { minuit_->GetParameter( i ) };
        Double_t error { 0.0 };
        Double_t negError { 0.0 };
        Double_t posError { 0.0 };
        Double_t globalcc { 0.0 };
        minuit_->GetErrors( i, posError, negError, error, globalcc );
        params_[i]->valueAndErrors( value, error, negError, posError );
        params_[i]->globalCorrelationCoeff( globalcc );
    }
}

// Definition of the fitting function for Minuit
void logLikeFun( Int_t& npar, Double_t* /*first_derivatives*/, Double_t& f, Double_t* par, Int_t /*iflag*/ )
{
    // Routine that specifies the negative log-likelihood function for the fit.
    // Used by the MINUIT minimising code.

    LauFitObject* theModel = LauFitter::fitter().getFitObject();

    // Set the internal parameters for this model using parameters from Minuit (pars):
    theModel->setParsFromMinuit( par, npar );

    // Set the value of f to be the total negative log-likelihood for the data sample.
    f = theModel->getTotNegLogLikelihood();
}
