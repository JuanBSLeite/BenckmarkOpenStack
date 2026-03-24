
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

/*! \file LauArgusPdf.hh
    \brief File containing declaration of LauArgusPdf class.
*/

/*! \class LauArgusPdf
    \brief Class for defining an ARGUS PDF.

    Class that allows the user to define an ARGUS PDF:
    Phys. Lett. B 241, 278 (1990)
*/

#ifndef LAU_ARGUS_PDF
#define LAU_ARGUS_PDF

#include "LauAbsPdf.hh"

#include "TString.h"

#include <vector>

class LauParameter;

class LauArgusPdf : public LauAbsPdf {

  public:
    //! Constructor
    /*!
        \param [in] theVarName the name of the abscissa variable
        \param [in] params the curve parameters - the shape, xi, and the end-point, m0, of the curve
        \param [in] minAbscissa the minimum value of the abscissa
        \param [in] maxAbscissa the maximum value of the abscissa
    */
    LauArgusPdf( const TString& theVarName,
                 const std::vector<LauAbsRValue*>& params,
                 Double_t minAbscissa,
                 Double_t maxAbscissa );

    //! Destructor
    virtual ~LauArgusPdf();

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
    LauArgusPdf( const LauArgusPdf& other );

    //! Copy assignment operator (not implemented)
    LauArgusPdf& operator=( const LauArgusPdf& other );

    //! Shape of curve
    LauAbsRValue* xi_;

    //! Endpoint of curve
    LauAbsRValue* m0_;

    ClassDef( LauArgusPdf, 0 )    // Define the ARGUS PDF
};

#endif
