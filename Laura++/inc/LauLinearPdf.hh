
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

/*! \file LauLinearPdf.hh
    \brief File containing declaration of LauLinearPdf class.
*/

/*! \class LauLinearPdf
    \brief Class for defining a straight line PDF.

    Class that allows the user to define a straight line PDF.
*/

#ifndef LAU_LINEAR_PDF
#define LAU_LINEAR_PDF

#include "LauAbsPdf.hh"

#include "TString.h"

#include <vector>

class LauLinearPdf : public LauAbsPdf {

  public:
    //! Constructor
    /*!
        \param [in] theVarName the name of the abscissa variable
        \param [in] params the PDF parameter - slope
        \param [in] minAbscissa the minimum value of the abscissa
        \param [in] maxAbscissa the maximum value of the abscissa
    */
    LauLinearPdf( const TString& theVarName,
                  const std::vector<LauAbsRValue*>& params,
                  Double_t minAbscissa,
                  Double_t maxAbscissa );

    //! Destructor
    virtual ~LauLinearPdf();

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
    LauLinearPdf( const LauLinearPdf& rhs );

    //! Copy assignment operator (not implemented)
    LauLinearPdf& operator=( const LauLinearPdf& rhs );

    //! Line slope
    LauAbsRValue* slope_;

    //! Flag to control printing of warnings about the PDF going negative
    Bool_t posflag_;

    ClassDef( LauLinearPdf, 0 )    // Define the Linear PDF
};

#endif
