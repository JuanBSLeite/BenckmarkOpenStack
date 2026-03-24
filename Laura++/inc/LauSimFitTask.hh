
/*
Copyright 2015 University of Warwick

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

/*! \file LauSimFitTask.hh
    \brief File containing declaration of LauSimFitTask class.
*/

/*! \class LauSimFitTask
    \brief The base class for any task process for simultaneous/combined fits

    Implementation of the JFit method described in arXiv:1409.5080 [physics.data-an].

    This class acts as the base class from which tasks should inherit.
    This allows any fitting framework to plug in to the JFit method.
*/

#ifndef LAU_SIM_FIT_TASK
#define LAU_SIM_FIT_TASK

#include "LauFitObject.hh"

#include "TMatrixDfwd.h"

class TMessage;
class TSocket;
class TString;
class LauFitNtuple;

class LauSimFitTask : public LauFitObject {

  public:
    //! Constructor
    LauSimFitTask();

    //! Destructor
    virtual ~LauSimFitTask();

    //! Obtain the number of tasks
    UInt_t nTasks() const { return nTasks_; }

    //! Obtain the ID number of this task
    UInt_t taskId() const { return taskId_; }

    //! Start the task process for simultaneous fitting
    /*!
        \param [in] dataFileName the name of the input data file
        \param [in] dataTreeName the name of the tree containing the data
        \param [in] histFileName the file name for the output histograms
        \param [in] tableFileName the file name for the latex output file
        \param [in] addressCoordinator the hostname of the machine running the coordinator process
        \param [in] portCoordinator the port number on which the coordinator process is listening
    */
    virtual void runTask( const TString& dataFileName,
                          const TString& dataTreeName,
                          const TString& histFileName,
                          const TString& tableFileName      = "",
                          const TString& addressCoordinator = "localhost",
                          const UInt_t portCoordinator      = 9090 );

    //! Initialise the fit model
    /*!
        Each class that inherits from this one must implement
        this to do what is appropriate
    */
    virtual void initialise() = 0;

    //! This function sets the parameter values from Minuit
    /*!
        \param [in] par an array storing the various parameter values
        \param [in] npar the number of free parameters
    */
    virtual void setParsFromMinuit( Double_t* par, Int_t npar ) = 0;

    //! Calculates the total negative log-likelihood
    virtual Double_t getTotNegLogLikelihood() = 0;

  protected:
    //! Const access to the fit ntuple
    const LauFitNtuple* fitNtuple() const { return fitNtuple_; }

    //! Access to the fit ntuple
    LauFitNtuple* fitNtuple() { return fitNtuple_; }

    //! Establish the connection to the coordinator process
    /*!
        \param [in] addressCoordinator the hostname of the machine running the coordinator process
        \param [in] portCoordinator the port number on which the coordinator process is listening
    */
    void connectToCoordinator( const TString& addressCoordinator, const UInt_t portCoordinator );

    //! Listen for requests from the coordinator and act accordingly
    void processCoordinatorRequests();

    //! Setup saving of fit results to ntuple/LaTeX table etc.
    /*!
        Provide here a default implementation that produces an ntuple only.
        Derived classes can override as they wish.

        \param [in] histFileName the file name for the output histograms
        \param [in] tableFileName the file name for the latex output file
    */
    virtual void setupResultsOutputs( const TString& histFileName, const TString& tableFileName );

    //! Package the initial fit parameters for transmission to the coordinator
    /*!
        \param [out] array the array to be filled with the LauParameter objects
    */
    virtual void prepareInitialParArray( TObjArray& array ) = 0;

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
                                     TObjArray& parsToCoordinator ) = 0;

    //! Open the input file and verify that all required variables are present
    /*!
        \param [in] dataFileName the name of the input file
        \param [in] dataTreeName the name of the input tree
    */
    virtual Bool_t verifyFitData( const TString& dataFileName, const TString& dataTreeName ) = 0;

    //! Read in the data for the current experiment
    /*!
        \return the number of events read in
    */
    virtual UInt_t readExperimentData() = 0;

    //! Cache the input data values to calculate the likelihood during the fit
    virtual void cacheInputFitVars() = 0;

    //! Write out any fit results
    virtual void writeOutAllFitResults();

  private:
    //! Copy constructor (not implemented)
    LauSimFitTask( const LauSimFitTask& rhs );

    //! Copy assignment operator (not implemented)
    LauSimFitTask& operator=( const LauSimFitTask& rhs );

    //! A socket to enable parallel setup
    TSocket* socketCoordinator_;

    //! Message from coordinator to the tasks
    TMessage* messageFromCoordinator_;

    //! Task id number
    UInt_t taskId_;

    //! The total number of tasks
    UInt_t nTasks_;

    //! Parameter values array (for reading from the coordinator)
    Double_t* parValues_;

    //! The fit ntuple
    LauFitNtuple* fitNtuple_;

    ClassDef( LauSimFitTask, 0 );
};

#endif
