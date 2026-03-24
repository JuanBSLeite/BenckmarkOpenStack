
/*
Copyright 2009 University of Warwick

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

/*! \file LauDPDepSumPdf.hh
    \brief File containing declaration of LauDPDepSumPdf class.
*/

/*! \class LauDPDepSumPdf
    \brief Class for defining a PDF that is the DP-dependent sum of two other PDFs.

    Class for defining a PDF that is the sum of two other PDFs, where the relative fraction of the two PDFs depends on the Dalitz-plot position.
*/

#ifndef LAU_DPDEPSUM_PDF
#define LAU_DPDEPSUM_PDF

#include "LauAbsPdf.hh"

#include "Rtypes.h"

class LauDaughters;
class LauEffModel;
class LauParameter;
class TH2;

class LauDPDepSumPdf : public LauAbsPdf {

  public:
    //! Define possibilties for the DP axes
    enum DPAxis {
        M12,        /*!< dependence is on m^2_12 */
        M13,        /*!< dependence is on m^2_13 */
        M23,        /*!< dependence is on m^2_23 */
        CentreDist, /*!< dependence is on the distance from the DP centre */
        MMIN,       /*!< dependence is on the minimum of m^2_13 and m^2_23 */
        MMAX        /*!< dependence is on the maximum of m^2_13 and m^2_23 */
    };

    //! Constructor - fraction determined by 2D histogram
    /*!
        \param [in] pdf1 the first PDF
        \param [in] pdf2 the second PDF
        \param [in] daughters the daughter particles
        \param [in] dpHisto the 2D histogram
        \param [in] upperHalf boolean flag to specify that the supplied histogram contains only the upper half of a symmetric DP (or lower half if using square DP coordinates)
        \param [in] useSpline specifies whether a spline is used to interpolate the histogram
    */
    LauDPDepSumPdf( LauAbsPdf* pdf1,
                    LauAbsPdf* pdf2,
                    const LauDaughters* daughters,
                    const TH2* dpHisto,
                    Bool_t upperHalf = kFALSE,
                    Bool_t useSpline = kFALSE );

    //! Constructor - fraction determined by a polynomial of a DP "axis"
    /*!
        \param [in] pdf1 the first PDF
        \param [in] pdf2 the second PDF
        \param [in] frac the fractional contribution of the first PDF to the final PDF
        \param [in] daughters the daughter particles
        \param [in] fracCoeffs the coefficients of the DP dependence of the PDF fraction
        \param [in] dpAxis the DP axis that defines the parameter dependence
    */
    LauDPDepSumPdf( LauAbsPdf* pdf1,
                    LauAbsPdf* pdf2,
                    LauParameter* frac,
                    const LauDaughters* daughters,
                    const std::vector<Double_t>& fracCoeffs,
                    DPAxis dpAxis );

    //! Destructor
    virtual ~LauDPDepSumPdf();

    //! Specifies whether or not the PDF is DP dependent.
    /*!
        \return true if the PDF is DP-dependent (the default)
    */
    virtual Bool_t isDPDependent() const { return kTRUE; }

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
    //! Scale fraction according to DP position
    /*!
        \param [in] dpPos the DP position
    */
    void scaleFrac( Double_t dpPos );

  private:
    //! Copy constructor - not implemented
    LauDPDepSumPdf( const LauDPDepSumPdf& other );

    //! Copy assignment operator - not implemented
    LauDPDepSumPdf& operator=( const LauDPDepSumPdf& other );

    //! Daughter particles
    LauDaughters* daughters_;

    //! First PDF
    LauAbsPdf* pdf1_;
    //! Second PDF
    LauAbsPdf* pdf2_;

    //! Fractional contribution of first PDF to the final PDF
    LauAbsRValue* frac_;
    //! Fractional contribution of first PDF to the final PDF
    Double_t fracVal_;
    //! DP dependence
    LauEffModel* dpDependence_;

    //! Polynomial used to scale fraction
    const std::vector<Double_t> fracCoeffs_;

    //! The DP axis we depend on
    DPAxis dpAxis_;

    //! Cached values of fractions
    std::vector<Double_t> fractions_;

    ClassDef( LauDPDepSumPdf, 0 )    // Define the sum PDF
};

#endif
