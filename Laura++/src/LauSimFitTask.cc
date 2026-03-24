
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

/*! \file LauSimFitTask.cc
    \brief File containing implementation of LauSimFitTask class.
*/

#include "LauSimFitTask.hh"

#include "LauFitNtuple.hh"

#include "TMatrixD.h"
#include "TMessage.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TSocket.h"
#include "TSystem.h"

#include <cstdlib>
#include <iostream>

LauSimFitTask::LauSimFitTask() :
    socketCoordinator_( 0 ),
    messageFromCoordinator_( 0 ),
    taskId_( 0 ),
    nTasks_( 0 ),
    parValues_( 0 ),
    fitNtuple_( 0 )
{
}

LauSimFitTask::~LauSimFitTask()
{
    delete socketCoordinator_;
    delete messageFromCoordinator_;
    delete[] parValues_;
    delete fitNtuple_;
}

void LauSimFitTask::runTask( const TString& dataFileName,
                             const TString& dataTreeName,
                             const TString& histFileName,
                             const TString& tableFileName,
                             const TString& addressCoordinator,
                             const UInt_t portCoordinator )
{
    // Establish the connection to the coordinator process
    this->connectToCoordinator( addressCoordinator, portCoordinator );

    // Initialise the fit model
    this->initialise();

    // NB call to addConParameters() is intentionally not included here cf.
    // LauAbsFitModel::run() since this has to be dealt with by the coordinator
    // to avoid multiple inclusions of each penalty term
    // Print a warning if constraints on combinations of parameters have been specified
    const auto& storeCon = this->formulaConstraints();
    if ( ! storeCon.empty() ) {
        std::cerr << "WARNING in LauSimFitTask::runTask : Constraints have been added but these will be ignored - they should have been added to the coordinator process"
                  << std::endl;
    }

    // Setup saving of fit results to ntuple/LaTeX table etc.
    this->setupResultsOutputs( histFileName, tableFileName );

    // This reads in the given dataFile and creates an input
    // fit data tree that stores them for all events and experiments.
    Bool_t dataOK = this->verifyFitData( dataFileName, dataTreeName );
    if ( ! dataOK ) {
        std::cerr << "ERROR in LauSimFitTask::runTask : Problem caching the fit data." << std::endl;
        return;
    }

    // Now process the various requests from the coordinator
    this->processCoordinatorRequests();

    std::cout << "INFO in LauSimFitTask::runTask : Fit task " << this->taskId()
              << " has finished successfully" << std::endl;
}

void LauSimFitTask::setupResultsOutputs( const TString& histFileName,
                                         const TString& /*tableFileName*/ )
{
    // Create and setup the fit results ntuple
    std::cout << "INFO in LauSimFitTask::setupResultsOutputs : Creating fit ntuple." << std::endl;
    if ( fitNtuple_ != 0 ) {
        delete fitNtuple_;
        fitNtuple_ = 0;
    }
    fitNtuple_ = new LauFitNtuple( histFileName, this->useAsymmFitErrors() );
}

void LauSimFitTask::connectToCoordinator( const TString& addressCoordinator,
                                          const UInt_t portCoordinator )
{
    if ( socketCoordinator_ != 0 ) {
        std::cerr << "ERROR in LauSimFitTask::connectToCoordinator : coordinator socket already present"
                  << std::endl;
        return;
    }

    // Open connection to coordinator
    socketCoordinator_ = new TSocket( addressCoordinator, portCoordinator );
    socketCoordinator_->Recv( messageFromCoordinator_ );

    messageFromCoordinator_->ReadUInt( taskId_ );
    messageFromCoordinator_->ReadUInt( nTasks_ );

    Bool_t useAsymErrs( kFALSE );
    messageFromCoordinator_->ReadBool( useAsymErrs );
    this->useAsymmFitErrors( useAsymErrs );

    delete messageFromCoordinator_;
    messageFromCoordinator_ = 0;

    std::cout << "INFO in LauSimFitTask::connectToCoordinator : Established connection to coordinator on port "
              << portCoordinator << std::endl;
    std::cout << "                                        : We are task " << taskId_ << " of "
              << nTasks_ << std::endl;
    if ( useAsymErrs ) {
        std::cout << "                                        : The fit will determine asymmetric errors"
                  << std::endl;
    }
}

void LauSimFitTask::processCoordinatorRequests()
{
    // Listen for requests from the coordinator and act accordingly

    TMessage messageToCoordinator( kMESS_ANY );

    while ( kTRUE ) {

        socketCoordinator_->Recv( messageFromCoordinator_ );

        if ( messageFromCoordinator_->What() == kMESS_STRING ) {

            TString msgStr;
            messageFromCoordinator_->ReadTString( msgStr );

            std::cout << "INFO in LauSimFitTask::processCoordinatorRequests : Received message from coordinator: "
                      << msgStr << std::endl;

            if ( msgStr == "Send Parameters" ) {

                // Send the fit parameters

                TObjArray array;
                this->prepareInitialParArray( array );

                // Create array to efficiently exchange parameter values with coordinator
                if ( parValues_ != 0 ) {
                    delete[] parValues_;
                    parValues_ = 0;
                }
                UInt_t nPar = array.GetEntries();
                parValues_  = new Double_t[nPar];

                messageToCoordinator.Reset( kMESS_OBJECT );
                messageToCoordinator.WriteObject( &array );
                socketCoordinator_->Send( messageToCoordinator );

            } else if ( msgStr == "Read Expt" ) {

                // Read the data for this experiment
                UInt_t iExp( 0 );
                messageFromCoordinator_->ReadUInt( iExp );

                this->setCurrentExperiment( iExp );

                UInt_t nEvents = this->readExperimentData();
                if ( nEvents < 1 ) {
                    std::cerr << "WARNING in LauSimFitTask::processCoordinatorRequests : Zero events in experiment "
                              << iExp << ", the coordinator should skip this experiment..."
                              << std::endl;
                }

                messageToCoordinator.Reset( kMESS_ANY );
                messageToCoordinator.WriteUInt( taskId_ );
                messageToCoordinator.WriteUInt( nEvents );
                socketCoordinator_->Send( messageToCoordinator );

            } else if ( msgStr == "Cache" ) {

                // Perform the caching

                this->cacheInputFitVars();

                messageToCoordinator.Reset( kMESS_ANY );
                messageToCoordinator.WriteUInt( taskId_ );
                messageToCoordinator.WriteBool( kTRUE );
                socketCoordinator_->Send( messageToCoordinator );

            } else if ( msgStr == "Asym Error Calc" ) {

                Bool_t asymErrorCalc( kFALSE );
                messageFromCoordinator_->ReadBool( asymErrorCalc );
                this->withinAsymErrorCalc( asymErrorCalc );

                messageToCoordinator.Reset( kMESS_ANY );
                messageToCoordinator.WriteUInt( taskId_ );
                messageToCoordinator.WriteBool( asymErrorCalc );
                socketCoordinator_->Send( messageToCoordinator );

            } else if ( msgStr == "Write Results" ) {

                this->writeOutAllFitResults();

                messageToCoordinator.Reset( kMESS_ANY );
                messageToCoordinator.WriteUInt( taskId_ );
                messageToCoordinator.WriteBool( kTRUE );
                socketCoordinator_->Send( messageToCoordinator );

            } else if ( msgStr == "Finish" ) {

                std::cout << "INFO in LauSimFitTask::processCoordinatorRequests : Message from coordinator to finish"
                          << std::endl;

                delete messageFromCoordinator_;
                messageFromCoordinator_ = 0;

                break;

            } else {

                std::cerr << "ERROR in LauSimFitTask::processCoordinatorRequests : Unexpected message from coordinator"
                          << std::endl;
                gSystem->Exit( EXIT_FAILURE );
            }

        } else if ( messageFromCoordinator_->What() == kMESS_OBJECT ) {

            std::cout << "INFO in LauSimFitTask::processCoordinatorRequests : Received message from coordinator: Finalise"
                      << std::endl;

            Int_t status( 0 );
            Double_t NLL( 0.0 );
            Double_t EDM( 0.0 );
            messageFromCoordinator_->ReadInt( status );
            messageFromCoordinator_->ReadDouble( NLL );
            messageFromCoordinator_->ReadDouble( EDM );

            TObjArray* objarray = dynamic_cast<TObjArray*>(
                messageFromCoordinator_->ReadObject( messageFromCoordinator_->GetClass() ) );
            if ( ! objarray ) {
                std::cerr << "ERROR in LauSimFitTask::processCoordinatorRequests : Error reading parameters from coordinator"
                          << std::endl;
                gSystem->Exit( EXIT_FAILURE );
            }

            TMatrixD* covMat = dynamic_cast<TMatrixD*>(
                messageFromCoordinator_->ReadObject( messageFromCoordinator_->GetClass() ) );
            if ( ! covMat ) {
                std::cerr << "ERROR in LauSimFitTask::processCoordinatorRequests : Error reading covariance matrix from coordinator"
                          << std::endl;
                gSystem->Exit( EXIT_FAILURE );
            }

            TObjArray array;
            LauAbsFitter::FitStatus fitStat { status, NLL, EDM };
            this->finaliseExperiment( fitStat, objarray, covMat, array );

            delete objarray;
            objarray = 0;
            delete covMat;
            covMat = 0;

            // Send the finalised parameters back to the coordinator
            messageToCoordinator.Reset( kMESS_ANY );
            messageToCoordinator.WriteUInt( taskId_ );
            messageToCoordinator.WriteBool( kTRUE );
            messageToCoordinator.WriteObject( &array );
            socketCoordinator_->Send( messageToCoordinator );

        } else if ( messageFromCoordinator_->What() == kMESS_ANY ) {

            UInt_t nPars( 0 );
            UInt_t nFreePars( 0 );
            messageFromCoordinator_->ReadUInt( nPars );
            messageFromCoordinator_->ReadUInt( nFreePars );

            if ( nPars != this->nTotParams() ) {
                std::cerr << "ERROR in LauSimFitTask::processCoordinatorRequests : Unexpected number of parameters received from coordinator"
                          << std::endl;
                std::cerr << "                                               : Received " << nPars
                          << " when expecting " << this->nTotParams() << std::endl;
                gSystem->Exit( EXIT_FAILURE );
            }

            messageFromCoordinator_->ReadFastArray( parValues_, nPars );

            this->setParsFromMinuit( parValues_, nFreePars );

            Double_t negLogLike = this->getTotNegLogLikelihood();

            messageToCoordinator.Reset( kMESS_ANY );
            messageToCoordinator.WriteDouble( negLogLike );
            socketCoordinator_->Send( messageToCoordinator );

        } else {
            std::cerr << "ERROR in LauSimFitTask::processCoordinatorRequests : Unexpected message type"
                      << std::endl;
            gSystem->Exit( EXIT_FAILURE );
        }

        delete messageFromCoordinator_;
        messageFromCoordinator_ = 0;
    }
}

void LauSimFitTask::writeOutAllFitResults()
{
    // Write out histograms at end
    if ( fitNtuple_ != 0 ) {
        fitNtuple_->writeOutFitResults();
    }
}
