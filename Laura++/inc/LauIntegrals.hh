
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

/*! \file LauIntegrals.hh
    \brief File containing declaration of LauIntegrals class.
*/

/*! \class LauIntegrals
    \brief Class for performing numerical integration routines.

    Class for calculating the Gauss-Legendre or Gauss-Hermite weights
    required for numerical integration of DP and other PDFs.
*/

#ifndef LAU_INTEGRALS
#define LAU_INTEGRALS

#include "Rtypes.h"

#include <vector>

class LauIntegrals {

  public:
    //! Constructor
    /*!
        \param [in] weightsPrecision the desired precision to which the weights should be calculated
    */
    LauIntegrals( Double_t weightsPrecision = 1.0e-6 );

    //! Destructor
    virtual ~LauIntegrals();

    //! Copy constructor
    LauIntegrals( const LauIntegrals& rhs );

    //! Copy assignment operator
    LauIntegrals& operator=( const LauIntegrals& rhs );

    //! Calculate the Gauss-Legendre weights
    /*!
        \param [in] numPoints the number of points in the integration
        \param [out] abscissas the values of the abscissa at each integration point
        \param [out] weights the weight at each integration point
    */
    void calcGaussLegendreWeights( const Int_t numPoints,
                                   std::vector<Double_t>& abscissas,
                                   std::vector<Double_t>& weights );

    //! Calculate the Gauss-Hermite weights
    /*!
        \param [in] numPoints the number of points in the integration
        \param [out] abscissas the values of the abscissa at each integration point
        \param [out] weights the weight at each integration point
    */
    void calcGaussHermiteWeights( const Int_t numPoints,
                                  std::vector<Double_t>& abscissas,
                                  std::vector<Double_t>& weights );

  private:
    //! The precision to which the weights should be calculated
    Double_t weightsPrecision_;

    ClassDef( LauIntegrals, 0 )
};

#endif
