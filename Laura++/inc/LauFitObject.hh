
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

/*! \file LauFitObject.hh
    \brief File containing declaration of LauFitObject class.
*/

#ifndef LAU_FIT_OBJECT
#define LAU_FIT_OBJECT

#include "LauAbsFitter.hh"
#include "LauFormulaPar.hh"

#include "TMatrixD.h"
#include "TObject.h"
#include "TString.h"
#include "TVectorD.h"

#include <set>
#include <vector>

/*! \class LauFitObject
    \brief The abstract interface for the objects that control the calculation of the likelihood.
*/

class LauFitObject : public TObject {

  public:
    //! Destructor
    virtual ~LauFitObject() = default;

    //! Turn on or off the computation of asymmetric errors (e.g. MINOS routine in Minuit)
    /*!
        \param [in] useAsymmErrors boolean specifying whether or not the computation of asymmetric errors is enabled
    */
    void useAsymmFitErrors( Bool_t useAsymmErrors ) { useAsymmFitErrors_ = useAsymmErrors; }

    //! Report whether or not calculation of asymmetric errors is enabled
    Bool_t useAsymmFitErrors() const { return useAsymmFitErrors_; }

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
    void twoStageFit( Bool_t doTwoStageFit ) { twoStageFit_ = doTwoStageFit; }

    //! Report whether the two-stage fit is enabled
    Bool_t twoStageFit() const { return twoStageFit_; }

    //! Mark that the fit is calculating asymmetric errors
    /*!
        This is called by the fitter interface to mark when
        entering and exiting the asymmetric error calculation.

        \param [in] inAsymErrCalc boolean marking that the fit is calculating the asymmetric errors
    */
    virtual void withinAsymErrorCalc( const Bool_t inAsymErrCalc )
    {
        withinAsymErrorCalc_ = inAsymErrCalc;
    }

    //! Query whether the fit is calculating the asymmetric errors
    /*!
        \return kTRUE if the fit is calculating the asymmetric errors, kFALSE otherwise
    */
    Bool_t withinAsymErrorCalc() const { return withinAsymErrorCalc_; }

    //! Set the number of experiments, the first experiment, and whether this is toy
    /*!
        The default settings are nExperiments = 1, firstExperiment = 0, toyExpts = kFALSE,
        i.e. the settings for fitting a single data sample.
        As such, this function only needs to be called if generating/fitting toy samples.

        \param [in] nExperiments the number of experiments
        \param [in] firstExperiment the number of the first experiment
        \param [in] toyExpts whether this is toy - determines whether to generate per-experiment means for each Gaussian constraint, as per arXiv:1210.7141
    */
    void setNExpts( UInt_t nExperiments, UInt_t firstExperiment, Bool_t toyExpts );

    //! Obtain the total number of events in the current experiment
    UInt_t eventsPerExpt() const { return evtsPerExpt_; }

    //! Obtain the number of experiments
    UInt_t nExpt() const { return nExpt_; }

    //! Obtain the number of the first experiment
    UInt_t firstExpt() const { return firstExpt_; }

    //! Obtain the number of the current experiment
    UInt_t iExpt() const { return iExpt_; }

    //! Obtain whether this is toy
    Bool_t toyExpts() const { return toyExpts_; }

    //! This function sets the parameter values from Minuit
    /*!
        This function has to be public since it is called from the global FCN.
        It should not be called otherwise!

        \param [in] par an array storing the various parameter values
        \param [in] npar the number of free parameters
    */
    virtual void setParsFromMinuit( Double_t* par, Int_t npar ) = 0;

    //! Calculate the new value of the negative log likelihood
    /*!
        This function has to be public since it is called from the global FCN.
        It should not be called otherwise!
    */
    virtual Double_t getTotNegLogLikelihood() = 0;

    //! Store constraint information for fit parameters
    /*!
        \deprecated Renamed to addFormulaConstraint, please switch to use this.  Will be dropped in next major release.

        \param [in] formula the formula to be used in the LauFormulaPar
        \param [in] pars a vector of LauParameter names to be used in the Formula, in the order specified by the formula
        \param [in] mean the value of the mean of the Gaussian constraint
        \param [in] width the value of the width of the Gaussian constraint
    */
    void addConstraint( const TString& formula,
                        const std::vector<TString>& pars,
                        const Double_t mean,
                        const Double_t width );

    //! Store constraint information for fit parameters
    /*!
        \param [in] formula the formula to be used in the LauFormulaPar
        \param [in] pars a vector of LauParameter names to be used in the Formula, in the order specified by the formula
        \param [in] mean the value of the mean of the Gaussian constraint
        \param [in] width the value of the width of the Gaussian constraint
    */
    void addFormulaConstraint( const TString& formula,
                               const std::vector<TString>& pars,
                               const Double_t mean,
                               const Double_t width );

    //! Store n-dimensional constraint information for fit parameters
    /*!
        \param [in] pars a vector of LauParameter names to be used in the constraint
        \param [in] means the values of the means of the Gaussian constraint
        \param [in] covMat the covariance matrix of the parameters of the Gaussian constraint
    */
    void addMultiDimConstraint( const std::vector<TString>& pars,
                                const TVectorD& means,
                                const TMatrixD& covMat );

  protected:
    //! Constructor
    LauFitObject();

    /*!
        \enum ConstraintType
        \brief Enumeration of the different types of constraint
    */
    enum class ConstraintType {
        Formula,    //!< Formula-based constraint on a combination of parameters
        MultDim     //!< Multi-dimensional constraint on several parameters
    };

    // Setup a struct to store information on constrained fit parameters
    /*!
        \struct FormulaConstraint
        \brief Struct to store formula-based constraint information
    */
    struct FormulaConstraint {
        //! The formula to be used in the LauFormulaPar
        TString formula_;
        //! The list of LauParameter names to be used in the LauFormulaPar
        std::vector<TString> conPars_;
        //! The mean value of the Gaussian constraint to be applied
        Double_t mean_;
        //! The width of the Gaussian constraint to be applied
        Double_t width_;
        //! The LauFormulaPar pointer
        std::unique_ptr<LauFormulaPar> formulaPar_;
    };

    // Setup a class to store information on n-dimensional constrained fit parameters
    /*!
        \class MultiDimConstraint
        \brief Class to store n-dimensional constraint information
    */
    class MultiDimConstraint {
      public:
        //! Default constructor
        MultiDimConstraint() = default;

        //! Constructor
        MultiDimConstraint( const std::vector<TString>& parNames,
                            const TVectorD& means,
                            const TMatrixD& covMat );

        //! Get the penalty term
        Double_t constraintPenalty() const;

        //! Generate per-experiment constraint means
        void generateConstraintMeans();

        //! The list of LauParameter names to be used in the constraint
        std::vector<TString> conPars_;
        //! The true mean values of the constraint
        TVectorD trueMeans_;
        //! The per-experiment mean values of the constraint
        TVectorD means_;
        //! The inverse covariance matrix of the parameters
        TMatrixD invCovMat_;
        //! The Cholesky Decomposition of the covariance matrix of the parameters
        TMatrixD sqrtCovMat_;
        //! The LauParameters used in the constraints
        std::vector<LauParameter*> conLauPars_;
    };

    //! Check if parameters names for constraints have already been used elsewhere
    /*!
        Also adds the names to the appropriate set

        \param [in] names a vector of parameter names
        \param [in] conType the type of constraint being added - determines which set to which to add the parmaeters
        \return kTRUE if no repetitions found, kFALSE if one or more repetitions found
    */
    Bool_t checkRepetition( const std::vector<TString>& names, const ConstraintType conType );

    //! Generate per-experiment mean for each Gaussian constraint
    /*!
        Generates a new mean for all Gaussian constraints.
        The constraints on a single fit parameter and the
        formula-based constraints are provided as the argument,
        while the multi-dimensional constraints are already a
        member variable.

        \param [in,out] conVars the fit parameter and formula-based constraints
    */
    void generateConstraintMeans( std::vector<LauAbsRValue*>& conVars );

    //! Const access to the formula constraints store
    const std::vector<FormulaConstraint>& formulaConstraints() const { return formulaConstraints_; }

    //! Access to the formula constraints store
    std::vector<FormulaConstraint>& formulaConstraints() { return formulaConstraints_; }

    //! Const access to the ND constraints store
    const std::vector<MultiDimConstraint>& multiDimConstraints() const
    {
        return multiDimConstraints_;
    }

    //! Access to the ND constraints store
    std::vector<MultiDimConstraint>& multiDimConstraints() { return multiDimConstraints_; }

    //! Const access to the parameter names used in formula constraints
    const std::set<TString>& formulaConstrainedPars() const { return formulaConstrainedPars_; }

    //! Access to the parameter names used in formula constraints
    std::set<TString>& formulaConstrainedPars() { return formulaConstrainedPars_; }

    //! Const access to the parameter names used in ND constraints
    const std::set<TString>& multiDimConstrainedPars() const { return multiDimConstrainedPars_; }

    //! Access to the parameter names used in ND constraints
    std::set<TString>& multiDimConstrainedPars() { return multiDimConstrainedPars_; }

    //! Reset the good/bad fit counters
    void resetFitCounters();

    //! Set the ID of the current experiment
    /*!
        \param [in] curExpt the experiment number
    */
    void setCurrentExperiment( const UInt_t curExpt ) { iExpt_ = curExpt; }

    //! Indicate the start of a new fit
    /*!
        \param [in] nPars the total number of fit parameters
        \param [in] nFreePars the number of free fit parameters
    */
    void startNewFit( const UInt_t nPars, const UInt_t nFreePars );

    //! Set the number of events in the current experiment
    void eventsPerExpt( UInt_t nEvents ) { evtsPerExpt_ = nEvents; }

    //! Access the worst log likelihood found so far
    Double_t worstLogLike() const { return worstLogLike_; }

    //! Set a new value for the worst log likelihood
    /*!
        \param [in] newWorstLogLike the new value of the worst log likelihood
    */
    void worstLogLike( const Double_t newWorstLogLike ) { worstLogLike_ = newWorstLogLike; }

    //! Store fit status information
    /*!
        \param [in] status the status information of the fit
        \param [in] covMatrix the fit covariance matrix
    */
    void storeFitStatus( const LauAbsFitter::FitStatus& status, const TMatrixD& covMatrix );

    //! Access the total number of fit parameters
    UInt_t nTotParams() const { return nParams_; }

    //! Access the total number of fit parameters
    UInt_t nFreeParams() const { return nFreeParams_; }

    //! Access the fit status information
    const LauAbsFitter::FitStatus& fitStatus() const { return fitStatus_; }

    //! Access the current NLL value
    Double_t nll() const { return fitStatus_.NLL; }

    //! Access the current EDM value
    Double_t edm() const { return fitStatus_.EDM; }

    //! Access the fit status code
    Int_t statusCode() const { return fitStatus_.status; }

    //! Access the fit covariance matrix
    const TMatrixD& covarianceMatrix() const { return covMatrix_; }

    //! Access the number of successful fits
    UInt_t numberOKFits() const { return numberOKFits_; }

    //! Access the number of failed fits
    UInt_t numberBadFits() const { return numberBadFits_; }

  private:
    //! Copy constructor (not implemented)
    LauFitObject( const LauFitObject& rhs );

    //! Copy assignment operator (not implemented)
    LauFitObject& operator=( const LauFitObject& rhs );

    //! Store the constraints for fit parameters until initialisation is complete
    std::vector<FormulaConstraint> formulaConstraints_;

    //! Store the ND constraints for fit parameters until initialisation is complete
    std::vector<MultiDimConstraint> multiDimConstraints_;

    //! Store the names of all parameters used in all formula constraints
    std::set<TString> formulaConstrainedPars_;

    //! Store the names of all parameters used in all multi-dimensional constraints
    std::set<TString> multiDimConstrainedPars_;

    //! Option to perform a two stage fit
    Bool_t twoStageFit_;

    //! Option to use asymmetric errors
    Bool_t useAsymmFitErrors_;

    //! The number of fit parameters
    UInt_t nParams_;

    //! The number of free fit parameters
    UInt_t nFreeParams_;

    //! Flag to indicate if the asymmetric error calculation (e.g. MINOS) is currently running
    Bool_t withinAsymErrorCalc_;

    //! Flag to indicate whether this is toy
    /*! Determines whether to generate per-experiment means for each Gaussian constraint, as per arXiv:1210.7141 */
    Bool_t toyExpts_;

    //! The number of the first experiment to consider
    UInt_t firstExpt_;

    //! The number of experiments to consider
    UInt_t nExpt_;

    //! The current experiment number
    UInt_t iExpt_;

    //! The number of events in the current experiment
    UInt_t evtsPerExpt_;

    //! The status of the current fit
    LauAbsFitter::FitStatus fitStatus_;

    //! The worst log likelihood value found so far
    Double_t worstLogLike_;

    //! The fit covariance matrix
    TMatrixD covMatrix_;

    //! The number of successful fits
    UInt_t numberOKFits_;

    //! The number of fit failures
    UInt_t numberBadFits_;

    ClassDef( LauFitObject, 0 )
};

#endif
