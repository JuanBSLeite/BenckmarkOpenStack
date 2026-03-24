
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

/*! \file LauDPDepGaussPdf.hh
    \brief File containing declaration of LauDPDepGaussPdf class.
*/

/*! \class LauDPDepGaussPdf
    \brief Class for defining a Gaussian PDF (DP dependent).

    Class that allows the user to define a Gaussian PDF where one or more of the parameters have a polynomial dependence on the DP position.
*/

#ifndef LAU_DPDEP_GAUSS_PDF
#define LAU_DPDEP_GAUSS_PDF

#include "LauAbsPdf.hh"

#include "TString.h"

#include <vector>

class LauDaughters;
class LauKinematics;
class LauParameter;

class LauDPDepGaussPdf : public LauAbsPdf {

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
        \param [in] params the PDF parameters - mean and sigma
        \param [in] minAbscissa the minimum value of the abscissa
        \param [in] maxAbscissa the maximum value of the abscissa
        \param [in] daughters the daughter particles
        \param [in] meanCoeffs the coefficients of the DP dependence of the Gaussian mean
        \param [in] sigmaCoeffs the coefficients of the DP dependence of the Gaussian sigma
        \param [in] dpAxis the DP axis that defines the parameter dependence
    */
    LauDPDepGaussPdf( const TString& theVarName,
                      const std::vector<LauAbsRValue*>& params,
                      Double_t minAbscissa,
                      Double_t maxAbscissa,
                      const LauDaughters* daughters,
                      const std::vector<Double_t>& meanCoeffs,
                      const std::vector<Double_t>& sigmaCoeffs,
                      DPAxis dpAxis );

    //! Destructor
    virtual ~LauDPDepGaussPdf();

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

  private:
    //! Copy constructor (not implemented)
    LauDPDepGaussPdf( const LauDPDepGaussPdf& other );

    //! Copy assignment operator (not implemented)
    LauDPDepGaussPdf& operator=( const LauDPDepGaussPdf& other );

    //! The current DP kinematics
    const LauKinematics* kinematics_;

    //! Gaussian mean
    LauAbsRValue* mean_;
    //! Gaussian sigma
    LauAbsRValue* sigma_;

    //! Gaussian mean
    Double_t meanVal_;
    //! Gaussian sigma
    Double_t sigmaVal_;

    //! Coefficients of Gaussian mean
    const std::vector<Double_t> meanCoeffs_;
    //! Coefficients of Gaussian sigma
    const std::vector<Double_t> sigmaCoeffs_;

    //! The DP axis we depend on
    DPAxis dpAxis_;

    ClassDef( LauDPDepGaussPdf, 0 )    // Define the Gaussian PDF
};

#endif
