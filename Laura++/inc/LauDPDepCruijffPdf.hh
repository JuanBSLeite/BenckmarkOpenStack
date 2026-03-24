
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

/*! \file LauDPDepCruijffPdf.hh
    \brief File containing declaration of LauDPDepCruijffPdf class.
*/

/*! \class LauDPDepCruijffPdf
    \brief Class for defining a Cruijff PDF (with DP dependence).

    Class that allows the user to define a Cruijff PDF where one or more of the parameters have a polynomial dependence on the DP position.
*/

#ifndef LAU_DPDEP_CRUIJFF_PDF
#define LAU_DPDEP_CRUIJFF_PDF

#include "LauAbsPdf.hh"

#include "TString.h"

#include <vector>

class LauDaughters;
class LauKinematics;
class LauParameter;

class LauDPDepCruijffPdf : public LauAbsPdf {

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

    //! Constructor
    /*!
        \param [in] theVarName the name of the abscissa variable
        \param [in] params the PDF parameters - mean, sigmaR, sigmaL, alphaR and alphaL
        \param [in] minAbscissa the minimum value of the abscissa
        \param [in] maxAbscissa the maximum value of the abscissa
        \param [in] daughters the daughter particles
        \param [in] meanCoeffs the coefficients of the DP dependence of the Gaussian mean
        \param [in] sigmaLCoeffs the coefficients of the DP dependence of the sigma for the left Gaussian
        \param [in] sigmaRCoeffs the coefficients of the DP dependence of the sigma for the right Gaussian
        \param [in] alphaLCoeffs the coefficients of the DP dependence of the alpha for the left Gaussian
        \param [in] alphaRCoeffs the coefficients of the DP dependence of the alpha for the right Gaussian
        \param [in] dpAxis the DP axis that defines the parameter dependence
    */
    LauDPDepCruijffPdf( const TString& theVarName,
                        const std::vector<LauAbsRValue*>& params,
                        Double_t minAbscissa,
                        Double_t maxAbscissa,
                        const LauDaughters* daughters,
                        const std::vector<Double_t>& meanCoeffs,
                        const std::vector<Double_t>& sigmaLCoeffs,
                        const std::vector<Double_t>& sigmaRCoeffs,
                        const std::vector<Double_t>& alphaLCoeffs,
                        const std::vector<Double_t>& alphaRCoeffs,
                        DPAxis dpAxis );

    //! Destructor
    virtual ~LauDPDepCruijffPdf();

    //! Specifies whether or not the PDF is DP dependent.
    /*!
        \return true if the PDF is DP-dependent (the default)
    */
    virtual Bool_t isDPDependent() const { return kTRUE; }

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
    //! Scale parameters by their dependence on the DP position
    /*!
        \param [in] dpPos the DP position
    */
    void scalePars( Double_t dpPos );

    //! Current PDF value
    /*!
        \param [in] abscissa the values of the abscissa(s)
    */
    Double_t currentPDFValue( Double_t abscissa ) const;

    //! Integrate the PDF using the Gauss-Legendre method
    /*!
        \return the integral of the PDF
    */
    virtual Double_t integrGaussLegendre();

    //! Integrate the PDF using the simple trapezoid method
    /*!
        \return the integral of the PDF
    */
    virtual Double_t integTrapezoid();

  private:
    //! Copy constructor (not implemented)
    LauDPDepCruijffPdf( const LauDPDepCruijffPdf& other );

    //! Copy assignment operator (not implemented)
    LauDPDepCruijffPdf& operator=( const LauDPDepCruijffPdf& other );

    //! The current DP kinematics
    const LauKinematics* kinematics_;

    //! Gaussian mean
    LauAbsRValue* mean_;
    //! Sigma of left Gaussian
    LauAbsRValue* sigmaL_;
    //! Sigma of right Gaussian
    LauAbsRValue* sigmaR_;
    //! Alpha of left Gaussian
    LauAbsRValue* alphaL_;
    //! Alpha of right Gaussian
    LauAbsRValue* alphaR_;

    //! Gaussian mean
    Double_t meanVal_;
    //! Sigma of left Gaussian
    Double_t sigmaLVal_;
    //! Sigma of right Gaussian
    Double_t sigmaRVal_;
    //! Alpha of left Gaussian
    Double_t alphaLVal_;
    //! Alpha of right Gaussian
    Double_t alphaRVal_;

    //! Coefficients of Gaussian mean
    const std::vector<Double_t> meanCoeffs_;
    //! Coefficients of sigma for the left Gaussian
    const std::vector<Double_t> sigmaLCoeffs_;
    //! Coefficients of sigma for the right Gaussian
    const std::vector<Double_t> sigmaRCoeffs_;
    //! Coefficients of alpha for the left Gaussian
    const std::vector<Double_t> alphaLCoeffs_;
    //! Coefficients of alpha for the right Gaussian
    const std::vector<Double_t> alphaRCoeffs_;

    //! The DP axis we depend on
    DPAxis dpAxis_;

    ClassDef( LauDPDepCruijffPdf, 0 )    // Define the Cruijff PDF
};

#endif
