
/*
Copyright 2006 University of Warwick

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

/*! \file LauBifurcatedGaussPdf.hh
    \brief File containing declaration of LauBifurcatedGaussPdf class.
*/

/*! \class LauBifurcatedGaussPdf
    \brief Class for defining a bifurcated Gaussian PDF.

    Class that allows the user to define a bifurcated Gaussian PDF,
    a Gaussian that has different widths either side of the peak position.
*/

#ifndef LAU_BIFURCATED_GAUSS_PDF
#define LAU_BIFURCATED_GAUSS_PDF

#include "LauAbsPdf.hh"

#include "TString.h"

#include <vector>

class LauParameter;

class LauBifurcatedGaussPdf : public LauAbsPdf {

  public:
    //! Constructor
    /*!
        \param [in] theVarName the name of the abscissa variable
        \param [in] params the PDF parameters - mean, sigmaL and sigmaR
        \param [in] minAbscissa the minimum value of the abscissa
        \param [in] maxAbscissa the maximum value of the abscissa
    */
    LauBifurcatedGaussPdf( const TString& theVarName,
                           const std::vector<LauAbsRValue*>& params,
                           Double_t minAbscissa,
                           Double_t maxAbscissa );

    //! Destructor
    virtual ~LauBifurcatedGaussPdf();

    //! Calculate the likelihood (and intermediate info) for a given abscissa
    /*!
        \param [in] abscissas the values of the abscissa(s)
    */
    virtual void calcLikelihoodInfo( const LauAbscissas& abscissas );

    using LauAbsPdf::calcLikelihoodInfo;

    //! Calculate the normalisation
    virtual void calcNorm();

    //! Calculate the PDF height
    /*!
        \param [in] kinematics the current DP kinematics
    */
    virtual void calcPDFHeight( const LauKinematics* kinematics );

  protected:

  private:
    //! Copy constructor (not implemented)
    LauBifurcatedGaussPdf( const LauBifurcatedGaussPdf& other );

    //! Copy assignment operator (not implemented)
    LauBifurcatedGaussPdf& operator=( const LauBifurcatedGaussPdf& other );

    //! Gaussian mean
    LauAbsRValue* mean_;

    //! Sigma of left Gaussian
    LauAbsRValue* sigmaL_;

    //! Sigma of right Gaussian
    LauAbsRValue* sigmaR_;

    ClassDef( LauBifurcatedGaussPdf, 0 )    // Define the Bifurcated Gaussian PDF
};

#endif
