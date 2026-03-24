
/*
Copyright 2017 University of Warwick

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

/*! \file LauFitObject.cc
    \brief File containing implementation of LauFitObject class.
*/

#include "LauFitObject.hh"

#include "LauRandom.hh"

#include "TDecompChol.h"
#include "TMatrixD.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TVectorD.h"

#include <iostream>

LauFitObject::LauFitObject() :
    TObject(),
    twoStageFit_( kFALSE ),
    useAsymmFitErrors_( kFALSE ),
    nParams_( 0 ),
    nFreeParams_( 0 ),
    withinAsymErrorCalc_( kFALSE ),
    toyExpts_( kFALSE ),
    firstExpt_( 0 ),
    nExpt_( 1 ),
    iExpt_( 0 ),
    evtsPerExpt_( 0 ),
    fitStatus_( { -1, 0.0, 0.0 } ),
    worstLogLike_( std::numeric_limits<Double_t>::max() ),
    covMatrix_(),
    numberOKFits_( 0 ),
    numberBadFits_( 0 )
{
}

void LauFitObject::setNExpts( UInt_t nExperiments, UInt_t firstExperiment, Bool_t toyExpts )
{
    nExpt_     = nExperiments;
    firstExpt_ = firstExperiment;
    toyExpts_  = toyExpts;

    if ( ! toyExpts && ( nExperiments > 1 || firstExperiment > 0 ) ) {
        std::cerr << "WARNING in LauFitObject::setNExpts : toyExpts is set to kFALSE but the values of nExperiments and firstExperiment indicate otherwise, please check"
                  << std::endl;
    } else if ( toyExpts && nExperiments == 1 && firstExperiment == 0 ) {
        std::cerr << "WARNING in LauFitObject::setNExpts : toyExpts is set to kTRUE but the values of nExperiments and firstExperiment perhaps indicate otherwise, please check"
                  << std::endl;
    }
}

void LauFitObject::resetFitCounters()
{
    numberOKFits_  = 0;
    numberBadFits_ = 0;
    fitStatus_     = { -1, 0.0, 0.0 };
}

void LauFitObject::startNewFit( const UInt_t nPars, const UInt_t nFreePars )
{
    // Reset the worst likelihood found to its catch-all value
    worstLogLike_ = std::numeric_limits<Double_t>::max();

    // Store the number of fit parameters (total and floating)
    nParams_     = nPars;
    nFreeParams_ = nFreePars;
}

void LauFitObject::storeFitStatus( const LauAbsFitter::FitStatus& status, const TMatrixD& covMatrix )
{
    fitStatus_ = status;

    covMatrix_.Clear();
    covMatrix_.ResizeTo( covMatrix.GetNrows(), covMatrix.GetNcols() );
    covMatrix_.SetMatrixArray( covMatrix.GetMatrixArray() );

    // Keep track of how many fits worked or failed
    // NB values of fitStatus_ indicate the status of the error matrix:
    // 0= not calculated at all
    // 1= approximation only, not accurate
    // 2= full matrix, but forced positive-definite
    // 3= full accurate covariance matrix
    if ( fitStatus_.status == 3 ) {
        ++numberOKFits_;
    } else {
        ++numberBadFits_;
    }
}

void LauFitObject::addConstraint( const TString& formula,
                                  const std::vector<TString>& pars,
                                  const Double_t mean,
                                  const Double_t width )
{
    std::cerr << "WARNING in LauFitObject::addConstraint : This function is deprecated, please switch to addFormulaConstraint!"
              << std::endl;
    this->addFormulaConstraint( formula, pars, mean, width );
}

void LauFitObject::addFormulaConstraint( const TString& formula,
                                         const std::vector<TString>& pars,
                                         const Double_t mean,
                                         const Double_t width )
{
    if ( ! this->checkRepetition( pars, ConstraintType::Formula ) ) {
        std::cerr << "WARNING in LauFitObject::addFormulaConstraint : Parameter(s) added to multiple constraints!"
                  << std::endl;
    }

    formulaConstraints_.emplace_back( FormulaConstraint { formula, pars, mean, width, nullptr } );

    std::cout << "INFO in LauFitObject::addFormulaConstraint : Added formula constraint" << std::endl;
}

void LauFitObject::addMultiDimConstraint( const std::vector<TString>& pars,
                                          const TVectorD& means,
                                          const TMatrixD& covMat )
{
    if ( ! this->checkRepetition( pars, ConstraintType::MultDim ) ) {
        std::cerr << "WARNING in LauFitObject::addMultiDimConstraint : Parameter(s) added to multiple constraints!"
                  << std::endl;
    }

    multiDimConstraints_.emplace_back( pars, means, covMat );

    std::cout << "INFO in LauFitObject::addMultiDimConstraint : Added multi-dimensional constraint"
              << std::endl;
}

void LauFitObject::generateConstraintMeans( std::vector<LauAbsRValue*>& conVars )
{
    if ( ! this->toyExpts() ) {
        return;
    }

    // For reproducibility, set a seed based on the experiment number
    // First, store the current seed, so that it can be restored afterwards
    const UInt_t oldSeed { LauRandom::randomFun()->GetSeed() };
    LauRandom::randomFun()->SetSeed( 827375 + this->iExpt() );

    for ( LauAbsRValue* par : conVars ) {
        par->generateConstraintMean();
    }

    for ( auto& constraint : multiDimConstraints_ ) {
        constraint.generateConstraintMeans();
    }

    // Restore the old random seed
    LauRandom::randomFun()->SetSeed( oldSeed );
}

Bool_t LauFitObject::checkRepetition( const std::vector<TString>& names, const ConstraintType conType )
{
    Bool_t allOK { kTRUE };

    for ( auto& newname : names ) {
        // Check in formula constraints
        if ( formulaConstrainedPars_.find( newname ) != formulaConstrainedPars_.end() ) {
            std::cerr << "WARNING in LauFitObject::checkRepetition: named parameter " << newname
                      << " already used in a formula constraint" << std::endl;
            allOK = kFALSE;
        }

        // Check in ND constraints
        if ( multiDimConstrainedPars_.find( newname ) != multiDimConstrainedPars_.end() ) {
            std::cerr << "WARNING in LauFitObject::checkRepetition: named parameter " << newname
                      << " already used in a multi-dimensional constraint" << std::endl;
            allOK = kFALSE;
        }

        // Add the new names to the appropriate set
        switch ( conType ) {
            case ConstraintType::Formula :
                formulaConstrainedPars_.insert( newname );
                break;
            case ConstraintType::MultDim :
                multiDimConstrainedPars_.insert( newname );
                break;
        }
    }

    return allOK;
}

LauFitObject::MultiDimConstraint::MultiDimConstraint( const std::vector<TString>& parNames,
                                                      const TVectorD& means,
                                                      const TMatrixD& covMat ) :
    conPars_ { parNames },
    trueMeans_ { means },
    means_ { means },
    invCovMat_ { covMat.GetNrows(), covMat.GetNcols() },
    sqrtCovMat_ { covMat.GetNrows(), covMat.GetNcols() }
{
    if ( covMat.GetNcols() != covMat.GetNrows() ) {
        std::cerr << "ERROR in LauFitObject::MultiDimConstraint : Covariance matrix is not square!"
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    if ( ( parNames.size() != static_cast<std::size_t>( means.GetNrows() ) ) ||
         ( parNames.size() != static_cast<std::size_t>( covMat.GetNrows() ) ) ) {
        std::cerr << "ERROR in LauFitObject::MultiDimConstraint : Different number of elements in vectors/covariance matrix!"
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Check invertion of the covariance matrix was successful
    TMatrixD invCovMat { TMatrixD::kInverted, covMat };
    if ( invCovMat == covMat ) {
        std::cerr << "ERROR in LauFitObject::MultiDimConstraint : covariance matrix inversion failed, check your input!"
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }
    invCovMat_ = invCovMat;

    // Check invertion of the covariance matrix was successful
    TDecompChol cholDecomp { covMat };
    if ( ! cholDecomp.Decompose() ) {
        std::cerr << "ERROR in LauFitObject::MultiDimConstraint : covariance matrix decomposition failed, check your input!"
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }
    sqrtCovMat_ = TMatrixD { TMatrixD::kTransposed, cholDecomp.GetU() };
}

Double_t LauFitObject::MultiDimConstraint::constraintPenalty() const
{
    TVectorD diff { means_.GetNrows() };
    for ( ULong_t j { 0 }; j < conLauPars_.size(); ++j ) {
        LauParameter* param = conLauPars_[j];
        diff[j]             = param->unblindValue();
    }
    diff -= means_;
    return 0.5 * invCovMat_.Similarity( diff );
}

void LauFitObject::MultiDimConstraint::generateConstraintMeans()
{
    TRandom* random = LauRandom::randomFun();

    for ( Int_t j { 0 }; j < trueMeans_.GetNrows(); ++j ) {
        means_[j] = random->Gaus( 0.0, 1.0 );
    }

    means_ *= sqrtCovMat_;
    means_ += trueMeans_;

    // Store the new mean in the parameter itself to enable writing it out to the results ntuple
    // (this is safe because a parameter cannot have both a 1D constraint and an ND constraint on it at the same time)
    for ( Int_t j { 0 }; j < trueMeans_.GetNrows(); ++j ) {
        conLauPars_[j]->constraintMean_ = means_[j];
        std::cout << "INFO in LauFitObject::MultiDimConstraint::generateConstraintMeans : set constraint mean for parameter \""
                  << conLauPars_[j]->name() << "\" to " << conLauPars_[j]->constraintMean_
                  << std::endl;
    }
}
