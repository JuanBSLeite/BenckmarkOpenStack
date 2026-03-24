
/*
Copyright 2004 University of Warwick

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

/*! \file LauFitNtuple.cc
    \brief File containing implementation of LauFitNtuple class.
*/

#include "LauFitNtuple.hh"

#include "LauAbsFitter.hh"
#include "LauFitter.hh"
#include "LauParamFixed.hh"
#include "LauParameter.hh"

#include "TFile.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TSystem.h"
#include "TTree.h"

#include <iostream>

LauFitNtuple::LauFitNtuple( const TString& fileName, Bool_t storeAsymErrors ) :
    rootFileName_( fileName ),
    rootFile_( 0 ),
    fitResults_( 0 ),
    definedFitTree_( kFALSE ),
    storeAsymErrors_( storeAsymErrors ),
    fitStatus_( { -1, 0.0, 0.0 } ),
    nFitPars_( 0 ),
    nFreePars_( 0 ),
    nExtraPars_( 0 ),
    iExpt_( 0 )
{
    rootFile_ = TFile::Open( rootFileName_, "recreate" );
    rootFile_->cd();
    fitResults_ = new TTree( "fitResults", "fitResults" );
    fitResults_->SetDirectory( rootFile_ );

    fitVars_.clear();
    constrainedVars_.clear();
    extraVars_.clear();
}

LauFitNtuple::~LauFitNtuple()
{
    if ( rootFile_ && rootFile_->IsOpen() ) {
        delete fitResults_;
        fitResults_ = 0;
    }
    delete rootFile_;
    rootFile_ = 0;
}

void LauFitNtuple::storeCorrMatrix( const UInt_t iExpt,
                                    const LauAbsFitter::FitStatus& fitStatus,
                                    const TMatrixD& covMatrix )
{
    // store the minimised NLL value, correlation matrix status and experiment number
    iExpt_     = iExpt;
    fitStatus_ = fitStatus;

    // make the correlation matrix the correct dimensions
    if ( definedFitTree_ == kFALSE ) {
        corrMatrix_.clear();
        corrMatrix_.resize( nFitPars_ );
        for ( UInt_t i = 0; i < nFitPars_; ++i ) {
            corrMatrix_[i].resize( nFitPars_ );
        }
    }

    // under certain circumstances, e.g. if the fit has failed in the first
    // stage of a two-stage fit, the covariance matrix might not have the
    // expected dimensions, or it might even be empty
    Bool_t needsPadding    = kFALSE;
    const UInt_t nElements = covMatrix.GetNoElements();
    if ( nElements == 0 ) {

        // if it's empty we can just make a diagonal matrix
        std::cerr << "WARNING in LauFitNtuple::storeCorrMatrix : received empty covariance matrix - will store diagonal correlation matrix"
                  << std::endl;
        for ( UInt_t i = 0; i < nFitPars_; ++i ) {
            for ( UInt_t j = 0; j < nFitPars_; ++j ) {
                if ( i == j ) {
                    corrMatrix_[i][j] = 1.0;
                } else {
                    corrMatrix_[i][j] = 0.0;
                }
            }
        }
        return;

    } else if ( nElements != nFreePars_ * nFreePars_ ) {

        UInt_t dimension    = covMatrix.GetNrows();
        UInt_t nSecondStage = 0;
        for ( UInt_t i = 0; i < nFitPars_; ++i ) {
            if ( fitVars_[i]->secondStage() ) {
                ++nSecondStage;
            }
        }

        if ( ( dimension + nSecondStage ) == nFreePars_ ) {
            needsPadding = kTRUE;
            std::cerr << "WARNING in LauFitNtuple::storeCorrMatrix : received smaller covariance matrix than expected, likely due to failure of first stage fit - will pad the correlation matrix"
                      << std::endl;
            std::cerr << "                                         : nFitPars_      = " << nFitPars_
                      << std::endl;
            std::cerr << "                                         : nFreePars_     = " << nFreePars_
                      << std::endl;
            std::cerr << "                                         : nSecondStage   = "
                      << nSecondStage << std::endl;
            std::cerr << "                                         : covMatrix size = " << dimension
                      << std::endl;
        } else {
            std::cerr << "WARNING in LauFitNtuple::storeCorrMatrix : received smaller covariance matrix than expected, for unknown reasons - will store diagonal correlation matrix"
                      << std::endl;
            std::cerr << "                                         : nFitPars_      = " << nFitPars_
                      << std::endl;
            std::cerr << "                                         : nFreePars_     = " << nFreePars_
                      << std::endl;
            std::cerr << "                                         : nSecondStage   = "
                      << nSecondStage << std::endl;
            std::cerr << "                                         : covMatrix size = " << dimension
                      << std::endl;
            for ( UInt_t i = 0; i < nFitPars_; ++i ) {
                for ( UInt_t j = 0; j < nFitPars_; ++j ) {
                    if ( i == j ) {
                        corrMatrix_[i][j] = 1.0;
                    } else {
                        corrMatrix_[i][j] = 0.0;
                    }
                }
            }
            return;
        }
    }

    // calculate the  correlation matrix information from the fitter covariance matrix
    Bool_t iFixed( kFALSE );
    Bool_t jFixed( kFALSE );
    Bool_t iSecondStage( kFALSE );
    Bool_t jSecondStage( kFALSE );
    UInt_t iFree( 0 );
    UInt_t jFree( 0 );

    for ( UInt_t i = 0; i < nFitPars_; ++i ) {

        iFixed       = fitVars_[i]->fixed();
        iSecondStage = fitVars_[i]->secondStage();

        // reset the "j" free parameter counter
        jFree = 0;

        // NB the supplied covariance matrix is of dimension nFreePars_ x nFreePars_
        for ( UInt_t j = 0; j < nFitPars_; ++j ) {

            jFixed       = fitVars_[j]->fixed();
            jSecondStage = fitVars_[j]->secondStage();

            if ( i == j ) {
                corrMatrix_[i][j] = 1.0;
            } else if ( iFixed == kTRUE || jFixed == kTRUE ) {
                corrMatrix_[i][j] = 0.0;
            } else if ( needsPadding && ( iSecondStage || jSecondStage ) ) {
                corrMatrix_[i][j] = 0.0;
            } else {
                Double_t r_ij  = covMatrix( iFree, jFree );
                Double_t r_ii  = covMatrix( iFree, iFree );
                Double_t r_jj  = covMatrix( jFree, jFree );
                Double_t denom = r_ii * r_jj;
                if ( denom < 0.0 ) {
                    r_ij  = 0.0;
                    denom = 1.0;
                }
                denom = TMath::Sqrt( denom );
                if ( denom < 1e-30 ) {
                    r_ij  = 0.0;
                    denom = 1.0;
                }
                corrMatrix_[i][j] = r_ij / denom;
            }

            if ( ! jFixed && ! ( needsPadding && jSecondStage ) ) {
                ++jFree;
            }
        }

        if ( ! iFixed && ! ( needsPadding && iSecondStage ) ) {
            ++iFree;
        }
    }
}

void LauFitNtuple::storeParsAndErrors( const std::vector<LauParameter*>& fitVars,
                                       const std::set<TString>& constrainedVars,
                                       const std::vector<LauParameter>& extraVars )
{
    fitVars_        = fitVars;
    UInt_t nFitPars = fitVars_.size();

    // the number of parameters being given to us should be the same as the number from the last fit
    // OR it's the first time so the "last" number is zero
    if ( nFitPars_ != 0 && nFitPars_ != nFitPars ) {
        std::cerr << "ERROR in LauFitNtuple::storeParsAndErrors : expected total number of parameters ("
                  << nFitPars_ << ") not the same as the number provided (" << nFitPars << ")."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    LauParamFixed pred;
    UInt_t nFreePars = nFitPars - std::count_if( fitVars_.begin(), fitVars_.end(), pred );

    // the number of free parameters being given to us should be the same as the number from the last fit
    // OR it's the first time so the "last" number is zero
    // (NB we check whether nFitPars_ is zero for this since it is possible to have zero free parameters, albeit rather daft)
    if ( nFitPars_ != 0 && nFreePars_ != nFreePars ) {
        std::cerr << "ERROR in LauFitNtuple::storeParsAndErrors : expected number of free parameters ("
                  << nFreePars_ << ") not the same as the number provided (" << nFreePars << ")."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    nFitPars_  = nFitPars;
    nFreePars_ = nFreePars;

    constrainedVars_ = constrainedVars;

    extraVars_  = extraVars;
    nExtraPars_ = extraVars_.size();
}

void LauFitNtuple::updateFitNtuple()
{
    // Now create and fill the stored fit results into an ntuple (TTree)
    if ( definedFitTree_ == kFALSE ) {

        std::cout << "INFO in LauFitNtuple::updateFitNtuple : totNoPars = " << nFitPars_ << std::endl;

        // Add experiment number as a branch
        fitResults_->Branch( "iExpt", &iExpt_, "iExpt/I" );
        fitResults_->Branch( "fitStatus", &fitStatus_.status, "fitStatus/I" );

        // Add NLL (negative log-likelihood) and EDM values from fit
        fitResults_->Branch( "NLL", &fitStatus_.NLL, "NLL/D" );
        fitResults_->Branch( "EDM", &fitStatus_.EDM, "EDM/D" );

        for ( UInt_t i = 0; i < nFitPars_; i++ ) {

            TString parName = fitVars_[i]->name();
            TString parNameD( parName );
            parNameD += "/D";
            fitResults_->Branch( parName.Data(), &fitVars_[i]->value_, parNameD.Data() );

            TString parInitName( parName );
            parInitName += "_True";
            TString parInitNameD( parInitName );
            parInitNameD += "/D";
            fitResults_->Branch( parInitName.Data(), &fitVars_[i]->genValue_, parInitNameD.Data() );

            if ( ! fitVars_[i]->fixed() ) {
                TString parErrName( parName );
                parErrName += "_Error";
                TString parErrNameD( parErrName );
                parErrNameD += "/D";
                fitResults_->Branch( parErrName.Data(), &fitVars_[i]->error_, parErrNameD.Data() );

                if ( storeAsymErrors_ ) {
                    TString parNegErrName( parName );
                    parNegErrName += "_NegError";
                    TString parNegErrNameD( parNegErrName );
                    parNegErrNameD += "/D";
                    fitResults_->Branch( parNegErrName.Data(),
                                         &fitVars_[i]->negError_,
                                         parNegErrNameD.Data() );

                    TString parPosErrName( parName );
                    parPosErrName += "_PosError";
                    TString parPosErrNameD( parPosErrName );
                    parPosErrNameD += "/D";
                    fitResults_->Branch( parPosErrName.Data(),
                                         &fitVars_[i]->posError_,
                                         parPosErrNameD.Data() );
                }

                TString parPullName( parName );
                parPullName += "_Pull";
                TString parPullNameD( parPullName );
                parPullNameD += "/D";
                fitResults_->Branch( parPullName.Data(), &fitVars_[i]->pull_, parPullNameD.Data() );

                if ( fitVars_[i]->gaussConstraint() ||
                     constrainedVars_.find( parName ) != constrainedVars_.end() ) {
                    TString parConstMeanName( parName );
                    parConstMeanName += "_ConstraintMean";
                    TString parConstMeanNameD( parConstMeanName );
                    parConstMeanNameD += "/D";
                    fitResults_->Branch( parConstMeanName.Data(),
                                         &fitVars_[i]->constraintMean_,
                                         parConstMeanNameD.Data() );
                }
            }

            // Now add in the correlation matrix values (only for floating parameters)
            if ( ! fitVars_[i]->fixed() ) {
                // First the global correlation coeffs
                TString parGCCName( parName );
                parGCCName += "_GCC";
                TString parGCCNameD( parGCCName );
                parGCCNameD += "/D";
                fitResults_->Branch( parGCCName.Data(), &fitVars_[i]->gcc_, parGCCNameD.Data() );

                if ( ! corrMatrix_.empty() ) {
                    // Then the rest
                    for ( UInt_t j = 0; j < nFitPars_; j++ ) {
                        if ( ! fitVars_[j]->fixed() && i != j ) {
                            TString parName2 = fitVars_[j]->name();
                            TString corrName( "corr__" );
                            corrName += parName;
                            corrName += "__";
                            corrName += parName2;

                            TString corrNameD( corrName );
                            corrNameD += "/D";
                            fitResults_->Branch( corrName.Data(),
                                                 &corrMatrix_[i][j],
                                                 corrNameD.Data() );
                        }
                    }
                }
            }
        }

        // Update extra parameter values...
        for ( UInt_t i = 0; i < nExtraPars_; i++ ) {

            TString parName = extraVars_[i].name();
            TString parNameD( parName );
            parNameD += "/D";
            fitResults_->Branch( parName.Data(), &extraVars_[i].value_, parNameD.Data() );

            TString parInitName( parName );
            parInitName += "_True";
            TString parInitNameD( parInitName );
            parInitNameD += "/D";
            fitResults_->Branch( parInitName.Data(), &extraVars_[i].genValue_, parInitNameD.Data() );

            //TString parErrName(parName); parErrName += "_Error";
            //TString parErrNameD(parErrName); parErrNameD += "/D";
            //fitResults_->Branch(parErrName.Data(), &extraVars_[i].error_, parErrNameD.Data());

            // Also find the fit fraction pull and store it
            //TString pullName(parName); pullName += "_Pull";
            //TString pullNameD(pullName); pullNameD += "/D";
            //fitResults_->Branch(pullName.Data(), &extraVars_[i].pull_, pullNameD.Data());
        }

        definedFitTree_ = kTRUE;
    }

    std::cout << "INFO in LauFitNtuple::updateFitNtuple : Stored fit values for experiment "
              << iExpt_ << "\n";
    std::cout << "                                      : fitStatus = " << fitStatus_.status
              << std::endl;
    std::cout << "                                      : NLL       = " << fitStatus_.NLL
              << std::endl;
    std::cout << "                                      : EDM       = " << fitStatus_.EDM
              << std::endl;

    fitResults_->Fill();
}

void LauFitNtuple::writeOutFitResults()
{
    // Write out the fit ntuple to the appropriate root file
    rootFile_->cd();
    fitResults_->Write( "", TObject::kOverwrite );
    rootFile_->Close();
    delete rootFile_;
    rootFile_ = 0;
}
