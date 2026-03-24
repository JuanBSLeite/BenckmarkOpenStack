
/*
Copyright 2013 University of Warwick

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

/*! \file LauDPDepMapPdf.hh
    \brief File containing declaration of LauDPDepMapPdf class.
*/

/*! \class LauDPDepMapPdf
    \brief Class to allow having different PDFs in different regions of the DP.

    Class to allow having different PDFs in different regions of the DP.
    The DP regions are specified in a 2D histogram or in a 1D histogram of some projection of the DP.
    Each PDF is then associated with that region.
*/

#ifndef LAU_DPDEPMAP_PDF
#define LAU_DPDEPMAP_PDF

#include "LauAbsPdf.hh"

#include "Rtypes.h"

#include <vector>

class LauDaughters;
class Lau2DAbsDP;
class LauParameter;
class TH2;
class TH1;

class LauDPDepMapPdf : public LauAbsPdf {

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

    //! Constructor - map described by 2D histogram
    /*!
        \param [in] pdfs the PDFs
        \param [in] daughters the daughter particles
        \param [in] dpHisto the 2D histogram
        \param [in] upperHalf boolean flag to specify that the supplied histogram contains only the upper half of a symmetric DP (or lower half if using square DP coordinates)
    */
    LauDPDepMapPdf( const std::vector<LauAbsPdf*>& pdfs,
                    const LauDaughters* daughters,
                    const TH2* dpHisto,
                    Bool_t upperHalf = kFALSE );

    //! Constructor - map described by 1D histogram of a DP "axis"
    /*!
        \param [in] pdfs the PDFs
        \param [in] daughters the daughter particles
        \param [in] dpAxisHisto the 1D histogram
        \param [in] dpAxis the DP axis
    */
    LauDPDepMapPdf( const std::vector<LauAbsPdf*>& pdfs,
                    const LauDaughters* daughters,
                    const TH1* dpAxisHisto,
                    DPAxis dpAxis );

    //! Destructor
    virtual ~LauDPDepMapPdf();

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
    //! Determine the DP region
    /*!
        \param [in] m13Sq the invariant mass squared of daughters 1 and 3
        \param [in] m23Sq the invariant mass squared of daughters 2 and 3
    */
    UInt_t determineDPRegion( Double_t m13Sq, Double_t m23Sq ) const;

  private:
    //! Copy constructor - not implemented
    LauDPDepMapPdf( const LauDPDepMapPdf& other );

    //! Copy assignment operator - not implemented
    LauDPDepMapPdf& operator=( const LauDPDepMapPdf& other );

    //! Daughter particles
    LauDaughters* daughters_;

    //! The PDFs
    std::vector<LauAbsPdf*> pdfs_;

    //! 2D histogram - DP
    Lau2DAbsDP* dpDependence_;
    //! 1D histogram - DP axis
    TH1* dpAxisDependence_;

    //! The DP axis we depend on
    DPAxis dpAxis_;

    //! Cached indices values
    std::vector<UInt_t> indices_;

    ClassDef( LauDPDepMapPdf, 0 )    // Define the sum PDF
};

#endif
