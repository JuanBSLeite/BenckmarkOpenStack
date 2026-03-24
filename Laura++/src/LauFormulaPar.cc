
/*
Copyright 2014 University of Warwick

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

/*! \file LauFormulaPar.cc
    \brief File containing implementation of LauFormulaPar class.
*/

#include "LauFormulaPar.hh"

#include "LauParameter.hh"
#include "LauRandom.hh"

#include "TMessage.h"
#include "TRandom.h"
#include "TSystem.h"

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

LauFormulaPar::LauFormulaPar( const TString& forName,
                              const TString& formula,
                              const std::vector<LauParameter*>& params ) :
    name_( forName ),
    formula_( forName, formula ),
    paramVec_( params ),
    paramArray_( nullptr ),
    gaussConstraint_( kFALSE ),
    constraintTrueMean_( 0.0 ),
    constraintMean_( 0.0 ),
    constraintWidth_( 0.0 )
{
    // Check length of vector matches number of parameter in the formula
    Int_t nPars = paramVec_.size();
    if ( formula_.GetNpar() != nPars ) {
        std::cerr << "ERROR in LauFormulaPar::evaluate : Number of parameters in the formula is : "
                  << formula_.GetNpar() << " and the number of LauParameters is : " << nPars
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    if ( formula_.GetNdim() != 0 ) {
        std::cerr << "ERROR in LauFormulaPar::evaluate : Given formula of dimension: "
                  << formula_.GetNdim() << " and not 0" << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Array of input parameters
    paramArray_ = new Double_t[nPars];
}

LauFormulaPar::~LauFormulaPar()
{
    delete[] paramArray_;
}

LauFormulaPar::LauFormulaPar( const LauFormulaPar& rhs ) :
    LauAbsRValue( rhs ),
    name_( rhs.name_ ),
    formula_( rhs.formula_ ),
    paramVec_( rhs.paramVec_ ),
    paramArray_( nullptr ),
    gaussConstraint_( rhs.gaussConstraint_ ),
    constraintTrueMean_( rhs.constraintTrueMean_ ),
    constraintMean_( rhs.constraintMean_ ),
    constraintWidth_( rhs.constraintWidth_ )
{
    // Check length of vector matches number of parameter in the formula
    Int_t nPars = paramVec_.size();
    if ( formula_.GetNpar() != nPars ) {
        std::cerr << "ERROR in LauFormulaPar::evaluate : Number of parameters in the formula is : "
                  << formula_.GetNpar() << " and the number of LauParameters is : " << nPars
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    if ( formula_.GetNdim() != 0 ) {
        std::cerr << "ERROR in LauFormulaPar::evaluate : Given formula of dimension: "
                  << formula_.GetNdim() << " and not 0" << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Array of input parameters
    paramArray_ = new Double_t[nPars];
}

LauFormulaPar& LauFormulaPar::operator=( const LauFormulaPar& rhs )
{
    if ( &rhs != this ) {
        name_    = rhs.name_;
        formula_ = rhs.formula_;

        Int_t nOldPars = paramVec_.size();
        Int_t nNewPars = rhs.paramVec_.size();

        paramVec_ = rhs.paramVec_;
        if ( nOldPars != nNewPars ) {
            delete[] paramArray_;
            paramArray_ = new Double_t[nNewPars];
        }

        gaussConstraint_    = rhs.gaussConstraint_;
        constraintTrueMean_ = rhs.constraintTrueMean_;
        constraintMean_     = rhs.constraintMean_;
        constraintWidth_    = rhs.constraintWidth_;
    }
    return *this;
}

Double_t LauFormulaPar::value() const
{
    //Assign vector values to array
    Int_t nPars = paramVec_.size();

    for ( Int_t i = 0; i < nPars; ++i ) {
        paramArray_[i] = paramVec_[i]->value();
    }

    return formula_.EvalPar( nullptr, paramArray_ );
}

Double_t LauFormulaPar::unblindValue() const
{
    //Assign vector values to array
    Int_t nPars = paramVec_.size();

    for ( Int_t i = 0; i < nPars; ++i ) {
        paramArray_[i] = paramVec_[i]->unblindValue();
    }

    return formula_.EvalPar( nullptr, paramArray_ );
}

Double_t LauFormulaPar::genValue() const
{
    //Assign vector values to array
    Int_t nPars = paramVec_.size();

    for ( Int_t i = 0; i < nPars; ++i ) {
        paramArray_[i] = paramVec_[i]->genValue();
    }

    return formula_.EvalPar( nullptr, paramArray_ );
}

Double_t LauFormulaPar::initValue() const
{
    //Assign vector values to array
    Int_t nPars = paramVec_.size();

    for ( Int_t i = 0; i < nPars; ++i ) {
        paramArray_[i] = paramVec_[i]->initValue();
    }

    return formula_.EvalPar( nullptr, paramArray_ );
}

Bool_t LauFormulaPar::fixed() const
{
    for ( std::vector<LauParameter*>::const_iterator iter = paramVec_.begin();
          iter != paramVec_.end();
          ++iter ) {
        if ( ! ( *iter )->fixed() ) {
            return kFALSE;
        }
    }
    return kTRUE;
}

Bool_t LauFormulaPar::blind() const
{
    for ( std::vector<LauParameter*>::const_iterator iter = paramVec_.begin();
          iter != paramVec_.end();
          ++iter ) {
        if ( ( *iter )->blind() ) {
            return kTRUE;
        }
    }
    return kFALSE;
}

void LauFormulaPar::addGaussianConstraint( Double_t newGaussMean, Double_t newGaussWidth )
{
    gaussConstraint_    = kTRUE;
    constraintTrueMean_ = newGaussMean;
    constraintMean_     = newGaussMean;
    constraintWidth_    = newGaussWidth;
}

void LauFormulaPar::removeGaussianConstraint()
{
    gaussConstraint_ = kFALSE;
}

void LauFormulaPar::generateConstraintMean()
{
    constraintMean_ = LauRandom::randomFun()->Gaus( constraintTrueMean_, constraintWidth_ );
}

Double_t LauFormulaPar::constraintPenalty() const
{
    const Double_t val { this->unblindValue() };
    const Double_t diff { val - constraintMean_ };
    const Double_t term { diff * diff };

    return term / ( 2.0 * constraintWidth_ * constraintWidth_ );
}
