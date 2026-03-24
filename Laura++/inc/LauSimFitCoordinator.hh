
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

/*! \file LauSimFitCoordinator.hh
    \brief File containing declaration of LauSimFitCoordinator class.
*/

/*! \class LauSimFitCoordinator
    \brief The coordinator process for simultaneous/combined fits

    Implementation of the JFit method described in arXiv:1409.5080 [physics.data-an].

    This class acts as the interface between the task processes and the minimiser.
*/

#ifndef LAU_SIM_FIT_COORDINATOR
#define LAU_SIM_FIT_COORDINATOR

#include "LauFitObject.hh"

#include "TMatrixD.h"
#include "TStopwatch.h"
#include "TString.h"

#include <map>
#include <vector>

class TMessage;
class TMonitor;
class TSocket;
class LauAbsRValue;
class LauParameter;
class LauFitNtuple;

class LauSimFitCoordinator : public LauFitObject {

  public:
    //! Constructor
    /*!
        \param [in] numTasks the number of tasks processes to expect connections from
        \param [in] port the port on which to listen for connections from the tasks
    */
    LauSimFitCoordinator( UInt_t numTasks, UInt_t port = 9090 );

    //! Destructor
    virtual ~LauSimFitCoordinator();

    //! Run the fit
    /*!
        \param [in] fitNtupleFileName the file to which the ntuple containing the fit results should be written
        \param [in] useAsymmErrors should asymmetric errors be calculated or not
        \param [in] doTwoStageFit should the fit be performed in two stages or not
    */
    void runSimFit( const TString& fitNtupleFileName,
                    const Bool_t useAsymmErrors = kFALSE,
                    const Bool_t doTwoStageFit  = kFALSE );

    //! Mark that the fit is calculating asymmetric errors
    /*!
        This function has to be public since it is called by
        the fitter interface to mark when entering and exiting
        the asymmetric error calculation.
        It should not be called otherwise!

        \param [in] inAsymErrCalc boolean marking that the fit is calculating the asymmetric errors
    */
    virtual void withinAsymErrorCalc( const Bool_t inAsymErrCalc );

    // Need to unshadow the query method defined in the base class
    using LauFitObject::withinAsymErrorCalc;

    //! This function sets the parameter values from Minuit
    /*!
        This function has to be public since it is called from the global FCN.
        It should not be called otherwise!

        \param [in] par an array storing the various parameter values
        \param [in] npar the number of free parameters
    */
    virtual void setParsFromMinuit( Double_t* par, Int_t npar );

    //! Calculate the new value of the negative log likelihood
    /*!
        This function has to be public since it is called from the global FCN.
        It should not be called otherwise!
    */
    virtual Double_t getTotNegLogLikelihood();

  protected:
    //! Print information on the parameters
    void printParInfo() const;

    //! Initialise
    void initialise();

    //! Initialise socket connections for the tasks
    void initSockets();

    //! Determine/update the parameter initial values from all tasks
    void getParametersFromTasks();

    //! Determine the parameter names and initial values from all tasks
    void getParametersFromTasksFirstTime();

    //! Update and verify the parameter initial values from all tasks
    void updateParametersFromTasks();

    //! Check for compatibility between two same-named parameters, which should therefore be identical
    void checkParameter( const LauParameter* param, UInt_t index ) const;

    //! Add parameters to the list of Gaussian constrained parameters
    void addConParameters();

    //! Calculate the penalty terms to the log likelihood from Gaussian constraints
    Double_t getLogLikelihoodPenalty();

    //! Instruct the tasks to read the input data for the given experiment
    /*!
        \return success/failure of the reading operations
    */
    Bool_t readData();

    //! Instruct the tasks to perform the caching
    /*!
        \return success/failure of the caching operations
    */
    Bool_t cacheInputData();

    //! Perform the fit for the current experiment
    void fitExpt();

    //! Instruct the tasks to update the initial fit parameter values, if required
    void checkInitFitParams();

    //! Return the final parameters to the tasks and instruct them to perform their finalisation
    Bool_t finalise();

    //! Instruct the tasks to write out the fit results
    Bool_t writeOutResults();

  private:
    //! Copy constructor (not implemented)
    LauSimFitCoordinator( const LauSimFitCoordinator& rhs );

    //! Copy assignment operator (not implemented)
    LauSimFitCoordinator& operator=( const LauSimFitCoordinator& rhs );

    //! The number of tasks
    const UInt_t nTasks_;

    //! The requested port
    const UInt_t reqPort_;

    //! The covariance sub-matrices for each task
    std::vector<TMatrixD> covMatrices_;

    //! Parallel setup monitor
    TMonitor* socketMonitor_;

    //! Sockets for each of the tasks
    std::vector<TSocket*> socketTasks_;

    //! Messages to tasks
    std::vector<TMessage*> messagesToTasks_;

    //! Message from tasks to the coordinator
    TMessage* messageFromTask_;

    //! Map of parameter names to index in the values vector
    std::map<TString, UInt_t> parIndices_;

    //! Reverse map of index in the values vector to parameter names
    std::map<UInt_t, TString> parNames_;

    //! Parameters
    std::vector<LauParameter*> params_;

    //! Gaussian constraints
    std::vector<LauAbsRValue*> conVars_;

    //! Parameter values
    std::vector<Double_t> parValues_;

    //! Lists of indices for each task
    std::vector<std::vector<UInt_t>> taskIndices_;

    //! Lists of indices of free parameters for each task
    std::vector<std::vector<UInt_t>> taskFreeIndices_;

    //! Parameter values to send to the tasks
    std::vector<Double_t*> vectorPar_;

    //! Likelihood values returned from the tasks
    std::vector<Double_t> vectorRes_;

    //! The fit timer
    TStopwatch timer_;

    //! The total fit timer
    TStopwatch cumulTimer_;

    //! The fit results ntuple
    LauFitNtuple* fitNtuple_;

    ClassDef( LauSimFitCoordinator, 0 );
};

#endif
