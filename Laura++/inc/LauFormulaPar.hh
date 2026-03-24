
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

/*! \file LauFormulaPar.hh
    \brief File containing declaration of LauFormulaPar class.
*/

/*! \class LauFormulaPar
    \brief Class for defining combinations of fit parameter objects.

    Allows for combinations of LauParameters to be passed to the fit models.
    Equations of the form [0]*2 + 3*[1] are accepted, with a vector of LauParameters to define
    inputs [0] and [1]. The parameter for [0] must be stored first in the vector, and so on.
*/

#ifndef LAU_FORMULAPAR
#define LAU_FORMULAPAR

#include "LauAbsRValue.hh"
#include "LauParameter.hh"

#include "TFormula.h"
#include "TString.h"

#include <map>
#include <vector>

class LauFormulaPar : public LauAbsRValue {

  public:
    //! Constructor
    /*!
        \param [in] forName the name of the formula
        \param [in] formula the desired expression, using TFormula syntax
        \param [in] params a vector of LauParameters used in the formula
    */
    LauFormulaPar( const TString& forName,
                   const TString& formula,
                   const std::vector<LauParameter*>& params );

    // Destructor
    virtual ~LauFormulaPar();

    //! Copy constructor
    LauFormulaPar( const LauFormulaPar& rhs );

    //! Copy assignment operator
    LauFormulaPar& operator=( const LauFormulaPar& rhs );

    //! Return the value of the LauFormalaPar
    /*!
        \return the value of the formula
    */
    Double_t value() const;

    //! The unblinded value of the parameter
    /*!
        \return the unblinded value of the parameter
    */
    Double_t unblindValue() const;

    //! The value generated for the parameter
    /*!
        \return the value generated for the parameter
    */
    Double_t genValue() const;

    //! The initial value of the parameter
    /*!
        \return the initial value of the parameter given to the fitter
    */
    Double_t initValue() const;

    //! The parameter name
    /*!
        \return the name of the parameter
    */
    inline const TString& name() const { return name_; }

    //! Set the parameter name
    /*!
        \param [in] newName the name of the parameter
    */
    inline void name( const TString& newName ) { name_ = newName; };

    //! Get the LauParameters used in LauFormulaPar
    /*!
        \return the list of LauParameters
    */
    std::vector<LauParameter*> getPars() { return paramVec_; }

    //! Boolean to say it is not an L value
    /*!
        \return kFALSE, LauFormulaPars are not L values
    */
    inline Bool_t isLValue() const { return kFALSE; }

    //! Boolean to say if the LauFormulaPar is fixed
    /*!
        \return kFALSE unless all LauParameters in the formula are fixed
    */
    Bool_t fixed() const;

    //! The blinding state
    /*!
        \return kTRUE if any of the LauParameters in the formula are blinded, kFALSE otherwise
    */
    Bool_t blind() const;

    //! Check whether a Gaussian constraints is applied
    /*!
        \return the boolean flag true/false whether a Gaussian constraint is applied
    */
    inline Bool_t gaussConstraint() const { return gaussConstraint_; }

    //! The penalty term from the Gaussian constraint
    /*!
        \return the penalty term from the Gaussian constraint
    */
    Double_t constraintPenalty() const;

    //! Add a Gaussian constraint (or modify an existing one)
    /*!
        \param [in] newGaussMean the new value of the Gaussian constraint mean
        \param [in] newGaussWidth the new value of the Gaussian constraint width
    */
    void addGaussianConstraint( Double_t newGaussMean, Double_t newGaussWidth );

    //! Remove the Gaussian constraint
    void removeGaussianConstraint();

    //! Generate per-experiment constraint mean
    void generateConstraintMean();

  protected:

  private:
    //! The parameter name
    TString name_;

    //! The formula
    mutable TFormula formula_;

    //! Vector of LauParameters in the formula
    std::vector<LauParameter*> paramVec_;

    //! Array to hold parameter values to pass to formula
    Double_t* paramArray_;

    //! Choice to use Gaussian constraint
    Bool_t gaussConstraint_;
    //! True mean of the Gaussian constraint
    Double_t constraintTrueMean_;
    //! Mean value of the Gaussian constraint
    Double_t constraintMean_;
    //! Width of the Gaussian constraint
    Double_t constraintWidth_;

    ClassDef( LauFormulaPar, 0 )
};

#endif
