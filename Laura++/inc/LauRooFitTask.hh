
/*
Copyright 2017 University of Warwick

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

/*! \file LauRooFitTask.hh
    \brief File containing declaration of LauRooFitTask class.
*/

/*! \class LauRooFitTask
    \brief A class for creating a RooFit-based task process for simultaneous/combined fits

    Implementation of the JFit method described in arXiv:1409.5080 [physics.data-an].
*/

#ifndef LAU_ROO_FIT_TASK
#define LAU_ROO_FIT_TASK

#include "LauSimFitTask.hh"

#include "RooAbsData.h"
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "TFile.h"
#include "TMatrixDfwd.h"
#include "TString.h"

#include <memory>
#include <set>
#include <utility>
#include <vector>

class LauParameter;
class RooAbsPdf;
class RooFormulaVar;
class RooRealVar;
class TTree;

class LauRooFitTask : public LauSimFitTask {

  public:
    //! Constructor
    LauRooFitTask( RooAbsPdf& model,
                   const Bool_t extended,
                   const RooArgSet& vars,
                   const TString& weightVarName = "" );

    //! Generate events from the model
    /*!
        \param [in] dataFileName the name of the file where the generated events are stored
        \param [in] dataTreeName the name of the tree used to store the variables
        \param [in] nExpt the number of experiments to generate
        \param [in] firstExpt the index of the first experiment number to generate
        \param [in] nEventsPerExpt the number of events to generate per experiment (only used if the PDF is not extended)
    */
    virtual void generate( const TString& dataFileName,
                           const TString& dataTreeName,
                           const UInt_t nExpt,
                           const UInt_t firstExpt     = 0,
                           const Int_t nEventsPerExpt = 0 );

    //! Initialise the fit model
    virtual void initialise();

    //! This function sets the parameter values from Minuit
    /*!
        \param [in] par an array storing the various parameter values
        \param [in] npar the number of free parameters
    */
    virtual void setParsFromMinuit( Double_t* par, Int_t npar );

    //! Calculates the total negative log-likelihood
    virtual Double_t getTotNegLogLikelihood();

  protected:
    //! Package the initial fit parameters for transmission to the coordinator
    /*!
        \param [out] array the array to be filled with the LauParameter objects
    */
    virtual void prepareInitialParArray( TObjArray& array );

    //! Convert a RooRealVar into a LauParameter
    /*!
        \param [in] rooParameter the RooRealVar to be converted
        \return the newly created LauParameter
    */
    LauParameter* convertToLauParameter( const RooRealVar* rooParameter ) const;

    //! Convert a RooFormulaVar into LauParameters
    /*!
        \param [in] rooFormula the RooFormulaVar to be converted
        \return the vector of pointers to the RooRealVars and newly created LauParameters
    */
    std::vector<std::pair<RooRealVar*, LauParameter*>> convertToLauParameters(
        const RooFormulaVar* rooFormula ) const;

    //! Perform all finalisation actions
    /*!
        - Receive the results of the fit from the coordinator
        - Perform any finalisation routines
        - Package the finalised fit parameters for transmission back to the coordinator

        \param [in] fitStat the status of the fit, e.g. status code, EDM, NLL
        \param [in] parsFromCoordinator the parameters at the fit minimum
        \param [in] covMat the fit covariance matrix
        \param [out] parsToCoordinator the array to be filled with the finalised LauParameter objects
    */
    virtual void finaliseExperiment( const LauAbsFitter::FitStatus& fitStat,
                                     const TObjArray* parsFromCoordinator,
                                     const TMatrixD* covMat,
                                     TObjArray& parsToCoordinator );

    //! Open the input file and verify that all required variables are present
    /*!
        \param [in] dataFileName the name of the input file
        \param [in] dataTreeName the name of the input tree
    */
    virtual Bool_t verifyFitData( const TString& dataFileName, const TString& dataTreeName );

    //! Read in the data for the current experiment
    /*!
        \return the number of events read in
    */
    virtual UInt_t readExperimentData();

    //! Cache the input data values to calculate the likelihood during the fit
    virtual void cacheInputFitVars();

  private:
    //! Cleanup the data
    void cleanData();

    //! The model
    RooAbsPdf& model_;

    //! The dataset variables
    RooArgSet dataVars_;

    //! The name of the (optional) weight variable in the dataset
    TString weightVarName_;

    //! The data file
    std::unique_ptr<TFile> dataFile_;

    //! The data tree
    TTree* dataTree_ { nullptr };

    //! The data for the current experiment
    std::unique_ptr<RooAbsData> exptData_;

    //! Is the PDF extended?
    const Bool_t extended_ { false };

    //! The experiment category variable
    std::set<UInt_t> iExptSet_;

    //! The NLL variable
    std::unique_ptr<RooAbsReal> nllVar_;

    //! The fit parameters (as RooRealVar's)
    std::vector<RooRealVar*> fitVars_;

    //! The fit parameters (as LauParameter's)
    std::vector<LauParameter*> fitPars_;

    ClassDef( LauRooFitTask, 0 );
};

#endif
