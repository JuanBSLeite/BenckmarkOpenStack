
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

/*! \file LauAsymmCalc.hh
    \brief File containing declaration of LauAsymmCalc class.
*/

/*! \class LauAsymmCalc
    \brief Class for calculating the asymmetry between two variables.
*/

#ifndef LAU_ASYMM_CALC
#define LAU_ASYMM_CALC

#include "Rtypes.h"

class LauAsymmCalc {

  public:
    //! Constructor
    /*!
        \param [in] negValue the negative value
        \param [in] posValue the positive value
    */
    LauAsymmCalc( Double_t negValue, Double_t posValue );

    //! Destructor
    virtual ~LauAsymmCalc();

    //! Copy constructor
    LauAsymmCalc( const LauAsymmCalc& rhs );

    //! Copy assignment operator
    LauAsymmCalc& operator=( const LauAsymmCalc& rhs );

    //! Retrieve the asymmetry
    /*!
        \return the asymmetry
    */
    Double_t getAsymmetry() { return asymm_; }

  protected:
    //! Calculate the asymmetry
    /*!
        \return the asymmetry
    */
    Double_t calcAsymmetry();

  private:
    //! The negative value
    Double_t negValue_;

    //! The positive value
    Double_t posValue_;

    //! The asymmetry
    Double_t asymm_;

    ClassDef( LauAsymmCalc, 0 )    // Calculate asymmetry between two variables
};

#endif
