
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

/*! \file LauSumPdf.hh
    \brief File containing declaration of LauSumPdf class.
*/

/*! \class LauSumPdf
    \brief Class for defining a PDF that is the sum of two other PDFs.

    Class that allows the user to define PDF that is the sum of two other PDFs.
    The relative fraction of the two PDFs is user-defined and can be a free parameter.
*/

#ifndef LAU_SUM_PDF
#define LAU_SUM_PDF

#include "LauAbsPdf.hh"

#include "Rtypes.h"

class LauParameter;

class LauSumPdf : public LauAbsPdf {

  public:
    //! Constructor
    /*!
        \param [in] pdf1 the first PDF
        \param [in] pdf2 the second PDF
        \param [in] frac1 the fractional contribution of the first PDF to the final PDF
    */
    LauSumPdf( LauAbsPdf* pdf1, LauAbsPdf* pdf2, LauParameter* frac1 );

    //! Destructor
    virtual ~LauSumPdf();

    //! Returns the number of input variables
    /*!
        \return the number of input variables
    */
    UInt_t nInputVars() const { return pdf1_->nInputVars(); }

    //! Boolean for the DP dependence of PDFs in the sum
    /*!
        \return DP dependence of PDF
    */
    virtual Bool_t isDPDependent() const
    {
        return ( pdf1_->isDPDependent() || pdf2_->isDPDependent() );
    }

    //! Cache information from data
    /*!
        \param [in] inputData the input data
    */
    virtual void cacheInfo( const LauFitDataTree& inputData );

    //! Calculate the likelihood (and intermediate info) for a given abscissa
    /*!
        \param [in] abscissas the values of the abscissa(s)
    */
    virtual void calcLikelihoodInfo( const LauAbscissas& abscissas );

    //! Calculate the likelihood (and intermediate info) for a given event number
    /*!
        \param [in] iEvt event number
    */
    virtual void calcLikelihoodInfo( UInt_t iEvt );

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
    LauSumPdf( const LauSumPdf& rhs );

    //! Copy assignment operator (not implemented)
    LauSumPdf& operator=( const LauSumPdf& rhs );

    //! First PDF
    LauAbsPdf* pdf1_;

    //! Second PDF
    LauAbsPdf* pdf2_;

    //! Fractional contribution of first PDF to the final PDF
    LauAbsRValue* frac_;

    ClassDef( LauSumPdf, 0 )    // Define the sum PDF
};

#endif
