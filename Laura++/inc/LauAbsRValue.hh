
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

/*! \file LauAbsRValue.hh
    \brief File containing declaration of LauAbsRValue class.
*/

/*! \class LauAbsRValue
    \brief Pure abstract base class for defining a parameter containing an R value

    Pure abstract base class for defining a parameter containing an R value,
    either a LauParameter or a LauFormulaPar
*/

#ifndef LAU_ABSRVALUE
#define LAU_ABSRVALUE

#include "TString.h"

#include <vector>

class LauParameter;

class LauAbsRValue {

  public:
    //! Constructor
    LauAbsRValue() {}

    //! Destructor
    virtual ~LauAbsRValue() {}

    //! Copy constructor
    LauAbsRValue( const LauAbsRValue& /*rhs*/ ) {}

    //! Copy assignment operator
    LauAbsRValue& operator=( const LauAbsRValue& /*rhs*/ ) { return *this; }

    //! Return the name of the parameter
    /*!
        \return the name of the parameter
    */
    virtual const TString& name() const = 0;

    //! Set the parameter name
    /*!
        \param [in] newName the name of the parameter
    */
    virtual void name( const TString& newName ) = 0;

    //! Return the value of the parameter
    /*!
        \return the value of the parameter
    */
    virtual Double_t value() const = 0;

    //! The unblinded value of the parameter
    /*!
        \return the unblinded value of the parameter
    */
    virtual Double_t unblindValue() const = 0;

    //! The value generated for the parameter
    /*!
        \return the value generated for the parameter
    */
    virtual Double_t genValue() const = 0;

    //! The initial value of the parameter
    /*!
        \return the initial value of the parameter given to the fitter
    */
    virtual Double_t initValue() const = 0;

    //! Check whether a Gaussian constraints is applied
    /*!
        \return the boolean flag true/false whether a Gaussian constraint is applied
    */
    virtual Bool_t gaussConstraint() const = 0;

    //! The penalty term from the Gaussian constraint
    /*!
        \return the penalty term from the Gaussian constraint
    */
    virtual Double_t constraintPenalty() const = 0;

    //! Generate per-experiment constraint mean
    virtual void generateConstraintMean() = 0;

    //! Return the list of LauParameters on which the LauAbsRValue depends
    /*!
        \return the list of LauParameters
    */
    virtual std::vector<LauParameter*> getPars() = 0;

    //! Is the parameter also an L value or not
    /*!
        \return whether the parameter is also an L value
    */
    virtual Bool_t isLValue() const = 0;

    //! Check is the parameter is fixed or floated
    /*!
        \return the boolean flag whether the parameter is fixed
    */
    virtual Bool_t fixed() const = 0;

    //! The blinding state
    /*!
        \return the blinding state: kTRUE means that it is blinded, kFALSE that it is not blinded
    */
    virtual Bool_t blind() const = 0;

  private:
    ClassDef( LauAbsRValue, 1 )    // Abstract base class for R parameters
};

#endif
