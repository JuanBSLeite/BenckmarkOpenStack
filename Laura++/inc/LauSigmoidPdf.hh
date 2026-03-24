
/*
Copyright 2012 University of Warwick

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

/*! \file LauSigmoidPdf.hh
    \brief File containing declaration of LauSigmoidPdf class.
*/

/*! \class LauSigmoidPdf
    \brief Class for defining a generalised sigmoid PDF.

    Class for defining a generalised sigmoid PDF:
    f(x) = 1 / ( 1 + e^(b - ax) )
    The standard sigmoid is obtained when a = 1 and b = 0.
*/

#ifndef LAU_SIGMOID_PDF
#define LAU_SIGMOID_PDF

#include "LauAbsPdf.hh"

#include "TString.h"

#include <vector>

class LauParameter;

class LauSigmoidPdf : public LauAbsPdf {

  public:
    //! Constructor
    /*!
        \param [in] theVarName the name of the abscissa variable
        \param [in] params the PDF parameters - a and b
        \param [in] minAbscissa the minimum value of the abscissa
        \param [in] maxAbscissa the maximum value of the abscissa
    */
    LauSigmoidPdf( const TString& theVarName,
                   const std::vector<LauAbsRValue*>& params,
                   Double_t minAbscissa,
                   Double_t maxAbscissa );

    //! Destructor
    virtual ~LauSigmoidPdf();

    //! Copy constructor
    LauSigmoidPdf( const LauSigmoidPdf& other );

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
    //! Parameter a
    LauAbsRValue* a_;
    //! Parameter b
    LauAbsRValue* b_;

    ClassDef( LauSigmoidPdf, 0 )    // Define the Sigmoid PDF
};

#endif
