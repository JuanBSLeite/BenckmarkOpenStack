
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

/*! \file LauMinuit.hh
    \brief File containing declaration of LauMinuit class.
*/

/*! \class LauMinuit
    \brief The interface to the Minuit fitter.

    The singleton interface to the Minuit fitter.
*/

#ifndef LAU_MINUIT
#define LAU_MINUIT

#include "LauAbsFitter.hh"
#include "LauPrint.hh"

#include "Rtypes.h"
#include "TMatrixD.h"

#include <utility>
#include <vector>

class LauParameter;
class TVirtualFitter;

class LauMinuit : public LauAbsFitter {

  public:
    //! Destructor
    virtual ~LauMinuit() = default;

    //! Initialise the fitter, setting the information on the parameters
    /*!
        \param [in] fitObj the object that controls the likelihood calculation
        \param [in] parameters the list of parameters of the fit
    */
    virtual void initialise( LauFitObject* fitObj, const std::vector<LauParameter*>& parameters );

    //! Get the object that controls the calculation of the likelihood
    virtual LauFitObject* getFitObject();

    //! Get the total number of fit parameters
    virtual UInt_t nParameters() const { return nParams_; }

    //! Get the number of floating fit parameters
    virtual UInt_t nFreeParameters() const { return nFreeParams_; }

    //! Determine whether the two-stage fit is enabled
    virtual Bool_t twoStageFit() const { return twoStageFit_; }

    //! Turn on or off the two stage fit
    /*!
        The two-stage fit allows certain parameters to be fixed
        in one stage and floated in another stage of the fit.
        Can be used, for example, in a CP fit where the
        CP-parameters are fixed to zero in the first stage
        (while the CP-average parameters are determined), then
        floated in the second.

        \param [in] doTwoStageFit boolean specifying whether or not the two-stage fit should be enabled
    */
    virtual void twoStageFit( Bool_t doTwoStageFit ) { twoStageFit_ = doTwoStageFit; }

    //! Determine whether calculation of asymmetric errors is enabled
    virtual Bool_t useAsymmFitErrors() const { return useAsymmFitErrors_; }

    //! Turn on or off the computation of asymmetric errors (i.e. the MINOS routine)
    /*!
        \param [in] useAsymmErrors boolean specifying whether or not the computation of asymmetric errors is enabled
    */
    virtual void useAsymmFitErrors( Bool_t useAsymmErrors ) { useAsymmFitErrors_ = useAsymmErrors; }

    //! Perform the minimisation of the fit function
    /*!
        \return the status code of the fit and the minimised value
    */
    virtual const FitStatus& minimise();

    //! Fix parameters marked as "second stage"
    virtual void fixSecondStageParameters();

    //! Release parameters marked as "second stage"
    virtual void releaseSecondStageParameters();

    //! Update the values and errors of the parameters based on the fit minimum
    virtual void updateParameters();

    //! Retrieve the fit covariance matrix
    virtual const TMatrixD& covarianceMatrix() const { return covMatrix_; }

  private:
    //! Allow the factory class to access private methods
    friend class LauFitter;

    //! Constructor
    LauMinuit( const UInt_t maxPar = 100, const LauOutputLevel verbosity = LauOutputLevel::Standard );

    //! Copy constructor - private and not implemented
    LauMinuit( const LauMinuit& ) = delete;

    //! Move constructor - private and not implemented
    LauMinuit( LauMinuit&& ) = delete;

    //! Copy assignment operator - private and not implemented
    LauMinuit& operator=( const LauMinuit& ) = delete;

    //! Move assignment operator - private and not implemented
    LauMinuit& operator=( LauMinuit&& ) = delete;

    //! The interface to Minuit
    TVirtualFitter* minuit_ { nullptr };

    //! The maximum number of parameters
    const UInt_t maxPar_ { 100 };

    //! The verbosity level of the fitter
    const LauOutputLevel outputLevel_ { LauOutputLevel::Standard };

    //! The fit parameters
    std::vector<LauParameter*> params_;

    //! The total number of parameters
    UInt_t nParams_ { 0 };

    //! The number of free parameters
    UInt_t nFreeParams_ { 0 };

    //! Option to perform a two stage fit
    Bool_t twoStageFit_ { kFALSE };

    //! Option to use asymmetric errors
    Bool_t useAsymmFitErrors_ { kFALSE };

    //! The status of the fit
    FitStatus fitStatus_ { -1, 0.0, 0.0 };

    //! The covariance matrix
    TMatrixD covMatrix_;

    ClassDef( LauMinuit, 0 );
};

#endif
