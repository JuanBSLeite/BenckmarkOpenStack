
/*
Copyright 2007 University of Warwick

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

/*! \file LauDPDepBifurGaussPdf.hh
    \brief File containing declaration of LauDPDepBifurGaussPdf class.
*/

/*! \class LauDPDepBifurGaussPdf
    \brief Class for defining a Bifurcated Gaussian PDF (DP dependent).

    Class that allows the user to define a Bifurcated Gaussian PDF where one or more of the parameters have a dependence on the DP position.
*/

#ifndef LAU_DPDEP_BIFURGAUSS_PDF
#define LAU_DPDEP_BIFURGAUSS_PDF

#include "LauAbsPdf.hh"

#include "TString.h"

#include <vector>

class LauDaughters;

class LauDPDepBifurGaussPdf : public LauAbsPdf {

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

    //! Define possibilties for the scaling method
    enum ScaleMethod {
        poly,        /*!< dependence is polynomial */
        polyNegPower /*!< dependence is a Laurent polynomial with only negative powers */
    };

    //! Constructor
    /*!
        \param [in] theVarName the name of the abscissa variable
        \param [in] params the PDF parameters - mean and sigma
        \param [in] minAbscissa the minimum value of the abscissa
        \param [in] maxAbscissa the maximum value of the abscissa
        \param [in] daughters the daughter particles
        \param [in] meanCoeffs the coefficients of the DP dependence of the Gaussian mean
        \param [in] sigmaLCoeffs the coefficients of the DP dependence of the sigma for the left Gaussian
        \param [in] sigmaRCoeffs the coefficients of the DP dependence of the sigma for the right Gaussian
        \param [in] dpAxis the DP axis that defines the parameter dependence
    */
    LauDPDepBifurGaussPdf( const TString& theVarName,
                           const std::vector<LauAbsRValue*>& params,
                           Double_t minAbscissa,
                           Double_t maxAbscissa,
                           const LauDaughters* daughters,
                           const std::vector<Double_t>& meanCoeffs,
                           const std::vector<Double_t>& sigmaLCoeffs,
                           const std::vector<Double_t>& sigmaRCoeffs,
                           DPAxis dpAxis );

    //! Destructor
    virtual ~LauDPDepBifurGaussPdf();

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

    //! Retrieve scaling method information
    /*!
        \return scaling method
    */
    ScaleMethod scaleMethod() const { return scaleMethod_; }

    //! Set the scaling method
    /*!
        \param [in] method the scaling method
    */
    void scaleMethod( ScaleMethod method ) { scaleMethod_ = method; }

  protected:
    //! Scale the Gaussian parameters with polynomial method
    /*!
        \param [in] perEventDist the event distribution
    */
    void scalePars_poly( Double_t perEventDist );

    //! Scale the Gaussian parameters with negative power polynomial method
    /*!
        \param [in] perEventDist the event distribution
    */
    void scalePars_polyNegPower( Double_t perEventDist );

  private:
    //! Copy constructor (not implemented)
    LauDPDepBifurGaussPdf( const LauDPDepBifurGaussPdf& other );

    //! Copy assignment operator (not implemented)
    LauDPDepBifurGaussPdf& operator=( const LauDPDepBifurGaussPdf& other );

    //! The current DP kinematics
    const LauKinematics* kinematics_;

    //! Gaussian mean
    LauAbsRValue* mean_;
    //! Left Gaussian sigma
    LauAbsRValue* sigmaL_;
    //! Right Gaussian sigma
    LauAbsRValue* sigmaR_;

    //! Gaussian mean
    Double_t meanVal_;
    //! Left Gaussian sigma
    Double_t sigmaLVal_;
    //! Right Gaussian sigma
    Double_t sigmaRVal_;

    //! Coefficients of Gaussian mean
    const std::vector<Double_t> meanCoeffs_;
    //! Coefficients of left Gaussian sigma
    const std::vector<Double_t> sigmaLCoeffs_;
    //! Coefficients of right Gaussian sigma
    const std::vector<Double_t> sigmaRCoeffs_;

    //! The DP axis we depend on
    DPAxis dpAxis_;

    //! Scaling method information
    ScaleMethod scaleMethod_;

    ClassDef( LauDPDepBifurGaussPdf, 0 )    // Define the MVA PDF
};

#endif
