
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

/*! \file LauSimFitCoordinator.cc
    \brief File containing implementation of LauSimFitCoordinator class.
*/

#include "LauSimFitCoordinator.hh"

#include "LauAbsFitter.hh"
#include "LauFitNtuple.hh"
#include "LauFitter.hh"
#include "LauFormulaPar.hh"
#include "LauParamFixed.hh"
#include "LauParameter.hh"

#include "TMath.h"
#include "TMatrixD.h"
#include "TMessage.h"
#include "TMonitor.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TServerSocket.h"
#include "TSocket.h"
#include "TSystem.h"

#include <cstdlib>
#include <iostream>
#include <limits>

LauSimFitCoordinator::LauSimFitCoordinator( UInt_t numTasks, UInt_t port ) :
    nTasks_( numTasks ),
    reqPort_( port ),
    socketMonitor_( 0 ),
    messageFromTask_( 0 ),
    fitNtuple_( 0 )
{
    messagesToTasks_.resize( nTasks_ );
    for ( UInt_t iTask( 0 ); iTask < nTasks_; ++iTask ) {
        messagesToTasks_[iTask] = new TMessage();
    }
}

LauSimFitCoordinator::~LauSimFitCoordinator()
{
    delete socketMonitor_;
    socketMonitor_ = 0;

    // Tell all tasks that they are finished and delete corresponding socket
    TString msgStr( "Finish" );
    TMessage message( kMESS_STRING );
    message.WriteTString( msgStr );
    for ( std::vector<TSocket*>::iterator iter = socketTasks_.begin(); iter != socketTasks_.end();
          ++iter ) {
        ( *iter )->Send( message );
        ( *iter )->Close();
        delete ( *iter );
    }
    socketTasks_.clear();

    // Remove all fit parameters
    for ( std::vector<LauParameter*>::iterator iter = params_.begin(); iter != params_.end();
          ++iter ) {
        delete *iter;
    }
    params_.clear();

    for ( std::vector<Double_t*>::iterator iter = vectorPar_.begin(); iter != vectorPar_.end();
          ++iter ) {
        delete[] ( *iter );
    }
    vectorPar_.clear();

    delete messageFromTask_;
    messageFromTask_ = 0;

    for ( std::vector<TMessage*>::iterator iter = messagesToTasks_.begin();
          iter != messagesToTasks_.end();
          ++iter ) {
        delete ( *iter );
    }
    messagesToTasks_.clear();

    delete fitNtuple_;
}

void LauSimFitCoordinator::initSockets()
{
    if ( socketMonitor_ != 0 ) {
        std::cerr << "ERROR in LauSimFitCoordinator::initSockets : Sockets already initialised."
                  << std::endl;
        return;
    }

    //initialise socket connection, then accept a connection and return a full-duplex communication socket.
    socketMonitor_ = new TMonitor();

    TServerSocket* ss  = new TServerSocket( reqPort_, kFALSE );
    UInt_t actual_port = ss->GetLocalPort();

    std::cout << "INFO in LauSimFitCoordinator::initSockets : Waiting for connection with "
              << nTasks_ << " tasks on port " << actual_port << std::endl;

    socketTasks_.resize( nTasks_ );
    for ( UInt_t iTask( 0 ); iTask < nTasks_; ++iTask ) {
        socketTasks_[iTask] = ss->Accept();
        std::cout << "                                     : Added task " << iTask << std::endl;
    }

    // tell the clients to start
    std::cout << "INFO in LauSimFitCoordinator::initSockets : Initialising tasks" << std::endl;
    for ( UInt_t iTask( 0 ); iTask < nTasks_; ++iTask ) {

        TMessage message( kMESS_ANY );
        message.WriteUInt( iTask );
        message.WriteUInt( nTasks_ );
        message.WriteBool( this->useAsymmFitErrors() );

        socketTasks_[iTask]->Send( message );

        socketMonitor_->Add( socketTasks_[iTask] );
    }
    std::cout << "                                     : Now start fit\n" << std::endl;

    ss->Close();
    delete ss;
}

/*
 * OLD VERSION THAT JUST GETS THE NAMES - COULD HAVE A SERIES OF EXCHANGES TO GET THE NAMES, INIT VALUES, RANGES, ETC. INSTEAD OF PASSING PARAMETERS
 * THIS INCREASES THE GENERALITY OF THE CODE, I.E. THERE IS NO NEED FOR THE TASKS TO KNOW ANY LAURA++ CLASS BUT THIS ONE, BUT MAKES IT RATHER MORE DENSE
 * FOR THE MOMENT I WILL STICK WITH THE METHOD OF PASSING LAUPARAMETER OBJECTS AROUND AND CONSIDER GOING BACK TO THIS GENERAL METHOD ONCE EVERYTHING IS WORKING
 *
void LauSimFitCoordinator::getParametersFromTasksFirstTime()
{
    taskIndices_.resize( nTasks_ );

    TSocket* sActive(0);
    for ( UInt_t iTask(0); iTask<nTasks_; ++iTask ) {
        // Send a message to the task, requesting the list of parameter names
        TString msgStr = "Parameter Names";
        TMessage message( kMESS_STRING );
        message.WriteTString( msgStr );
        socketTasks_[iTask]->Send(message);

        // Wait to receive the response and check that it has come from the task we just requested from
        sActive = socketMonitor_->Select();
        if ( sActive != socketTasks_[iTask] ) {
            std::cerr << "ERROR in LauSimFitCoordinator::getParametersFromTasksFirstTime : Received message from a different task than expected!" << std::endl;
            gSystem->Exit(1);
        }

        // Read the object and extract the parameter names
        socketTasks_[iTask]->Recv( messageFromTask_ );
        TObjArray * objarray = dynamic_cast<TObjArray*>( messageFromTask_->ReadObject( messageFromTask_->GetClass() ) );
        if ( ! objarray ) {
            std::cerr << "ERROR in LauSimFitCoordinator::getParametersFromTasksFirstTime : Error reading parameter names from task" << std::endl;
            gSystem->Exit(1);
        }

        Int_t nPars = objarray->GetEntries();
        for ( Int_t iPar(0); iPar < nPars; ++iPar ) {
            TObjString* objstring = dynamic_cast<TObjString*>( (*objarray)[iPar] );
            if ( ! objstring ) {
                std::cerr << "ERROR in LauSimFitCoordinator::getParametersFromTasksFirstTime : Error reading parameter names from task" << std::endl;
                gSystem->Exit(1);
            }
            TString parname = objstring->GetString();

            std::map< TString, UInt_t >::iterator iter = parIndices_.find( parname );
            if ( iter != parIndices_.end() ) {
                UInt_t index = iter->second;
                taskIndices_[iTask].push_back( index );
            } else {
                UInt_t index = parIndices_.size();
                parIndices_.insert( std::make_pair( parname, index ) );
                parNames_.insert( std::make_pair( index, parname ) );
                taskIndices_[iTask].push_back( index );
            }
        }

        delete objarray; objarray = 0;
        delete messageFromTask_; messageFromTask_ = 0;
    }

    UInt_t nPars = parNames_.size();
    parValues_.resize( nPars );
}
*/

void LauSimFitCoordinator::getParametersFromTasks()
{
    if ( socketMonitor_ == 0 ) {
        std::cerr << "ERROR in LauSimFitCoordinator::getParametersFromTasks : Sockets not initialised."
                  << std::endl;
        return;
    }

    if ( params_.empty() ) {
        this->getParametersFromTasksFirstTime();

        // Add variables to Gaussian constrain to a list
        this->addConParameters();
    } else {
        this->updateParametersFromTasks();
    }
}

void LauSimFitCoordinator::updateParametersFromTasks()
{
    TSocket* sActive( 0 );

    // Construct a message, requesting the list of parameter names
    TString msgStr = "Send Parameters";
    TMessage message( kMESS_STRING );
    message.WriteTString( msgStr );

    for ( UInt_t iTask( 0 ); iTask < nTasks_; ++iTask ) {
        // Send the message to the task
        socketTasks_[iTask]->Send( message );

        // Wait to receive the response and check that it has come from the task we just requested from
        sActive = socketMonitor_->Select();
        if ( sActive != socketTasks_[iTask] ) {
            std::cerr << "ERROR in LauSimFitCoordinator::updateParametersFromTasks : Received message from a different task than expected!"
                      << std::endl;
            gSystem->Exit( 1 );
        }

        // Read the object and extract the parameter names
        socketTasks_[iTask]->Recv( messageFromTask_ );
        TObjArray* objarray = dynamic_cast<TObjArray*>(
            messageFromTask_->ReadObject( messageFromTask_->GetClass() ) );
        if ( ! objarray ) {
            std::cerr << "ERROR in LauSimFitCoordinator::updateParametersFromTasks : Error reading parameter names from task"
                      << std::endl;
            gSystem->Exit( 1 );
        }

        // We want to auto-delete the supplied parameters since we only copy their values in this case
        objarray->SetOwner( kTRUE );

        const UInt_t nPars = objarray->GetEntries();
        if ( nPars != taskIndices_[iTask].size() ) {
            std::cerr << "ERROR in LauSimFitCoordinator::updateParametersFromTasks : Unexpected number of parameters received from task"
                      << std::endl;
            gSystem->Exit( 1 );
        }

        for ( UInt_t iPar( 0 ); iPar < nPars; ++iPar ) {
            LauParameter* parameter = dynamic_cast<LauParameter*>( ( *objarray )[iPar] );
            if ( ! parameter ) {
                std::cerr << "ERROR in LauSimFitCoordinator::updateParametersFromTasks : Error reading parameter from task"
                          << std::endl;
                gSystem->Exit( 1 );
            }

            TString parname   = parameter->name();
            Double_t parvalue = parameter->initValue();

            std::map<TString, UInt_t>::iterator iter = parIndices_.find( parname );
            if ( iter == parIndices_.end() ) {
                std::cerr << "ERROR in LauSimFitCoordinator::updateParametersFromTasks : Unexpected parameter name received from task"
                          << std::endl;
                gSystem->Exit( 1 );
            }

            const UInt_t index = iter->second;
            if ( taskIndices_[iTask][iPar] != index ) {
                std::cerr << "ERROR in LauSimFitCoordinator::updateParametersFromTasks : Unexpected parameter received from task"
                          << std::endl;
                gSystem->Exit( 1 );
            }

            params_[index]->initValue( parvalue );
            parValues_[index]       = parvalue;
            vectorPar_[iTask][iPar] = parvalue;
            this->checkParameter( parameter, index );
        }

        delete objarray;
        objarray = 0;
        delete messageFromTask_;
        messageFromTask_ = 0;
    }
}

void LauSimFitCoordinator::getParametersFromTasksFirstTime()
{
    taskIndices_.resize( nTasks_ );
    taskFreeIndices_.resize( nTasks_ );
    vectorPar_.resize( nTasks_ );
    vectorRes_.resize( nTasks_ );

    TSocket* sActive( 0 );

    // Construct a message, requesting the list of parameter names
    TString msgStr = "Send Parameters";
    TMessage message( kMESS_STRING );
    message.WriteTString( msgStr );

    for ( UInt_t iTask( 0 ); iTask < nTasks_; ++iTask ) {
        // Send the message to the task
        socketTasks_[iTask]->Send( message );

        // Wait to receive the response and check that it has come from the task we just requested from
        sActive = socketMonitor_->Select();
        if ( sActive != socketTasks_[iTask] ) {
            std::cerr << "ERROR in LauSimFitCoordinator::getParametersFromTasksFirstTime : Received message from a different task than expected!"
                      << std::endl;
            gSystem->Exit( 1 );
        }

        // Read the object and extract the parameter names
        socketTasks_[iTask]->Recv( messageFromTask_ );
        TObjArray* objarray = dynamic_cast<TObjArray*>(
            messageFromTask_->ReadObject( messageFromTask_->GetClass() ) );
        if ( ! objarray ) {
            std::cerr << "ERROR in LauSimFitCoordinator::getParametersFromTasksFirstTime : Error reading parameters from task"
                      << std::endl;
            gSystem->Exit( 1 );
        }

        const UInt_t nPars = objarray->GetEntries();

        vectorPar_[iTask] = new Double_t[nPars];

        for ( UInt_t iPar( 0 ); iPar < nPars; ++iPar ) {
            LauParameter* parameter = dynamic_cast<LauParameter*>( ( *objarray )[iPar] );
            if ( ! parameter ) {
                std::cerr << "ERROR in LauSimFitCoordinator::getParametersFromTasksFirstTime : Error reading parameter from task"
                          << std::endl;
                gSystem->Exit( 1 );
            }

            TString parname   = parameter->name();
            Double_t parvalue = parameter->initValue();
            Bool_t parfixed   = parameter->fixed();

            std::map<TString, UInt_t>::iterator iter = parIndices_.find( parname );
            if ( iter != parIndices_.end() ) {
                UInt_t index = iter->second;
                taskIndices_[iTask].push_back( index );
                if ( ! parfixed ) {
                    taskFreeIndices_[iTask].push_back( index );
                }
                this->checkParameter( parameter, index );
            } else {
                UInt_t index = parIndices_.size();
                parIndices_.insert( std::make_pair( parname, index ) );
                parNames_.insert( std::make_pair( index, parname ) );
                taskIndices_[iTask].push_back( index );
                if ( ! parfixed ) {
                    taskFreeIndices_[iTask].push_back( index );
                }
                params_.push_back( parameter );
                parValues_.push_back( parvalue );
            }
            vectorPar_[iTask][iPar] = parvalue;
        }

        delete objarray;
        objarray = 0;
        delete messageFromTask_;
        messageFromTask_ = 0;
    }
}

void LauSimFitCoordinator::printParInfo() const
{
    for ( UInt_t iTask( 0 ); iTask < nTasks_; ++iTask ) {
        const std::vector<UInt_t>& indices = taskIndices_[iTask];

        std::cout << "INFO in LauSimFitCoordinator::printParInfo : Task " << iTask
                  << " has the following parameters:\n";
        for ( std::vector<UInt_t>::const_iterator iter = indices.begin(); iter != indices.end();
              ++iter ) {
            const TString& parName  = parNames_.find( *iter )->second;
            Double_t parValue       = parValues_[*iter];
            const LauParameter* par = params_[*iter];
            if ( par->name() != parName || par->initValue() != parValue ) {
                std::cerr << "ERROR in LauSimFitCoordinator::printParInfo : Discrepancy in parameter name and value records, this is very strange!!"
                          << std::endl;
            }

            std::cout << "                                      : " << parName << " = " << parValue
                      << " and has index " << *iter << "\n";
        }

        std::cout << std::endl;
    }

    std::cout << "INFO in LauSimFitCoordinator::printParInfo : "
              << "There are " << params_.size() << " parameters in total" << std::endl;
}

void LauSimFitCoordinator::checkParameter( const LauParameter* param, UInt_t index ) const
{
    const LauParameter* storedPar = params_[index];

    TString parName = storedPar->name();
    if ( param->name() != parName ) {
        std::cerr << "ERROR in LauSimFitCoordinator::checkParameter : Parameter name is different!!  This shouldn't happen!!"
                  << std::endl;
    }
    if ( param->initValue() != storedPar->initValue() ) {
        std::cerr << "WARNING in LauSimFitCoordinator::checkParameter : Initial value for parameter "
                  << parName
                  << " is different, will use the value first set: " << storedPar->initValue()
                  << std::endl;
    }
    if ( param->minValue() != storedPar->minValue() ) {
        std::cerr << "WARNING in LauSimFitCoordinator::checkParameter : Minimum allowed value for parameter "
                  << parName
                  << " is different, will use the value first set: " << storedPar->minValue()
                  << std::endl;
    }
    if ( param->maxValue() != storedPar->maxValue() ) {
        std::cerr << "WARNING in LauSimFitCoordinator::checkParameter : Maximum allowed value for parameter "
                  << parName
                  << " is different, will use the value first set: " << storedPar->maxValue()
                  << std::endl;
    }
    if ( param->fixed() != storedPar->fixed() ) {
        std::cerr << "WARNING in LauSimFitCoordinator::checkParameter : Fixed/floating property of parameter "
                  << parName << " is different, will use the value first set: "
                  << ( storedPar->fixed() ? "fixed" : "floating" ) << std::endl;
    }
    if ( param->secondStage() != storedPar->secondStage() ) {
        std::cerr << "WARNING in LauSimFitCoordinator::checkParameter : Second stage property of parameter "
                  << parName << " is different, will use the value first set: "
                  << ( storedPar->secondStage() ? "true" : "false" ) << std::endl;
    }
}

void LauSimFitCoordinator::initialise()
{
    this->initSockets();
}

void LauSimFitCoordinator::runSimFit( const TString& fitNtupleFileName,
                                      const Bool_t useAsymmErrors,
                                      const Bool_t doTwoStageFit )
{
    // Routine to perform the total fit.

    // First, initialise
    this->useAsymmFitErrors( useAsymmErrors );
    this->twoStageFit( doTwoStageFit );
    this->initialise();

    const UInt_t nExp { this->nExpt() };
    const UInt_t firstExp { this->firstExpt() };
    const Bool_t isToy { this->toyExpts() };

    std::cout << "INFO in LauSimFitCoordinator::runSimFit : First experiment = " << firstExp
              << std::endl;
    std::cout << "INFO in LauSimFitCoordinator::runSimFit : Number of experiments = " << nExp
              << std::endl;
    std::cout << "INFO in LauSimFitCoordinator::runSimFit : Is toy MC = "
              << ( isToy ? "True" : "False" ) << std::endl;

    // Start the cumulative timer
    cumulTimer_.Start();

    this->resetFitCounters();

    // Create and setup the fit results ntuple
    std::cout << "INFO in LauSimFitCoordinator::runSimFit : Creating fit ntuple." << std::endl;
    if ( fitNtuple_ != 0 ) {
        delete fitNtuple_;
        fitNtuple_ = 0;
    }
    fitNtuple_ = new LauFitNtuple( fitNtupleFileName, useAsymmErrors );

    // Loop over the number of experiments
    for ( UInt_t iExp = firstExp; iExp < ( firstExp + nExp ); ++iExp ) {

        // Start the timer to see how long each fit takes
        timer_.Start();

        this->setCurrentExperiment( iExp );

        // Instruct the tasks to read the data for this experiment
        Bool_t readOK = this->readData();
        if ( ! readOK ) {
            std::cerr << "ERROR in LauSimFitCoordinator::runSimFit : One or more tasks reported problems with reading data for experiment "
                      << iExp << ", skipping..." << std::endl;
            timer_.Stop();
            continue;
        }

        // Instruct the tasks to perform the caching
        this->cacheInputData();

        // Do the fit
        this->fitExpt();

        // Stop the timer and see how long the program took so far
        timer_.Stop();
        timer_.Print();

        // Instruct the tasks to finalise the results
        this->finalise();
    }

    // Print out total timing info.
    std::cout << "INFO in LauSimFitCoordinator::runSimFit : Cumulative timing:" << std::endl;
    cumulTimer_.Stop();
    cumulTimer_.Print();

    // Print out stats on OK fits.
    const UInt_t nOKFits  = this->numberOKFits();
    const UInt_t nBadFits = this->numberBadFits();
    std::cout << "INFO in LauSimFitCoordinator::runSimFit : Number of OK Fits = " << nOKFits
              << std::endl;
    std::cout << "INFO in LauSimFitCoordinator::runSimFit : Number of Failed Fits = " << nBadFits
              << std::endl;
    Double_t fitEff( 0.0 );
    if ( nExp != 0 ) {
        fitEff = nOKFits / ( 1.0 * nExp );
    }
    std::cout << "INFO in LauSimFitCoordinator::runSimFit : Fit efficiency = " << fitEff * 100.0
              << "%." << std::endl;

    // Instruct the tasks to write out any fit results (ntuples etc...).
    this->writeOutResults();
}

void LauSimFitCoordinator::withinAsymErrorCalc( const Bool_t inAsymErrCalc )
{
    this->LauFitObject::withinAsymErrorCalc( inAsymErrCalc );

    if ( socketMonitor_ == 0 ) {
        std::cerr << "ERROR in LauSimFitCoordinator::withinAsymErrorCalc : Sockets not initialised."
                  << std::endl;
        return;
    }

    // Construct a message, informing the tasks whether or not we are now within the asymmetric error calculation
    TString msgStr( "Asym Error Calc" );
    const Bool_t asymErrorCalc( this->withinAsymErrorCalc() );
    TMessage message( kMESS_STRING );
    message.WriteTString( msgStr );
    message.WriteBool( asymErrorCalc );

    // Send the message to the tasks
    for ( UInt_t iTask( 0 ); iTask < nTasks_; ++iTask ) {
        socketTasks_[iTask]->Send( message );
    }

    TSocket* sActive( 0 );
    UInt_t responsesReceived( 0 );
    while ( responsesReceived != nTasks_ ) {

        // Get the next queued response
        sActive = socketMonitor_->Select();

        // Extract from the message the ID of the task and the number of events read
        Bool_t response( kTRUE );
        UInt_t iTask( 0 );
        sActive->Recv( messageFromTask_ );
        messageFromTask_->ReadUInt( iTask );
        messageFromTask_->ReadBool( response );

        if ( response != asymErrorCalc ) {
            std::cerr << "WARNING in LauSimFitCoordinator::withinAsymErrorCalc : Problem informing task "
                      << iTask << std::endl;
        }

        ++responsesReceived;
    }
}

Bool_t LauSimFitCoordinator::readData()
{
    if ( socketMonitor_ == 0 ) {
        std::cerr << "ERROR in LauSimFitCoordinator::readData : Sockets not initialised."
                  << std::endl;
        return kFALSE;
    }

    // Construct a message, requesting to read the data for the given experiment
    TString msgStr( "Read Expt" );
    const UInt_t iExp( this->iExpt() );
    TMessage message( kMESS_STRING );
    message.WriteTString( msgStr );
    message.WriteUInt( iExp );

    // Send the message to the tasks
    for ( UInt_t iTask( 0 ); iTask < nTasks_; ++iTask ) {
        socketTasks_[iTask]->Send( message );
    }

    TSocket* sActive( 0 );
    UInt_t responsesReceived( 0 );
    Bool_t ok( kTRUE );
    while ( responsesReceived != nTasks_ ) {

        // Get the next queued response
        sActive = socketMonitor_->Select();

        // Extract from the message the ID of the task and the number of events read
        sActive->Recv( messageFromTask_ );
        UInt_t iTask( 0 );
        UInt_t nEvents( 0 );
        messageFromTask_->ReadUInt( iTask );
        messageFromTask_->ReadUInt( nEvents );

        if ( nEvents <= 0 ) {
            std::cerr << "ERROR in LauSimFitCoordinator::readData : Task " << iTask
                      << " reports no events found for experiment " << iExp << std::endl;
            ok = kFALSE;
        } else {
            std::cerr << "INFO in LauSimFitCoordinator::readData : Task " << iTask << " reports "
                      << nEvents << " events found for experiment " << iExp << std::endl;
        }

        ++responsesReceived;
    }

    return ok;
}

Bool_t LauSimFitCoordinator::cacheInputData()
{
    if ( socketMonitor_ == 0 ) {
        std::cerr << "ERROR in LauSimFitCoordinator::cacheInputData : Sockets not initialised."
                  << std::endl;
        return kFALSE;
    }

    // Construct a message, requesting it to read the data for the given experiment
    TString msgStr( "Cache" );
    TMessage message( kMESS_STRING );
    message.WriteTString( msgStr );

    for ( UInt_t iTask( 0 ); iTask < nTasks_; ++iTask ) {
        // Send the message to the task
        socketTasks_[iTask]->Send( message );
    }

    TSocket* sActive( 0 );
    UInt_t responsesReceived( 0 );
    Bool_t allOK( kTRUE );
    while ( responsesReceived != nTasks_ ) {

        // Get the next queued response
        sActive = socketMonitor_->Select();

        // Extract from the message the ID of the task and the success/failure flag
        sActive->Recv( messageFromTask_ );
        UInt_t iTask( 0 );
        Bool_t ok( kTRUE );
        messageFromTask_->ReadUInt( iTask );
        messageFromTask_->ReadBool( ok );

        if ( ! ok ) {
            std::cerr << "ERROR in LauSimFitCoordinator::cacheInputData : Task " << iTask
                      << " reports an error performing caching" << std::endl;
            allOK = kFALSE;
        }

        ++responsesReceived;
    }

    return allOK;
}

void LauSimFitCoordinator::checkInitFitParams()
{
    this->getParametersFromTasks();
    this->printParInfo();
}

void LauSimFitCoordinator::fitExpt()
{
    // Routine to perform the actual fit for the given experiment

    // Instruct the tasks to update initial fit parameters if required (e.g. if using random numbers).
    this->checkInitFitParams();

    // If we're fitting toy experiments then re-generate the means of any constraints
    this->generateConstraintMeans( conVars_ );

    // Initialise the fitter
    LauFitter::fitter().useAsymmFitErrors( this->useAsymmFitErrors() );
    LauFitter::fitter().twoStageFit( this->twoStageFit() );
    LauFitter::fitter().initialise( this, params_ );

    this->startNewFit( LauFitter::fitter().nParameters(), LauFitter::fitter().nFreeParameters() );

    // Now ready for minimisation step
    std::cout << "\nINFO in LauSimFitCoordinator::fitExpt : Start minimisation...\n";
    LauAbsFitter::FitStatus fitResult = LauFitter::fitter().minimise();

    // If we're doing a two stage fit we can now release (i.e. float)
    // the 2nd stage parameters and re-fit
    if ( this->twoStageFit() ) {
        if ( fitResult.status != 3 ) {
            std::cerr << "ERROR in LauSimFitCoordinator:fitExpt : Not running second stage fit since first stage failed."
                      << std::endl;
            LauFitter::fitter().releaseSecondStageParameters();
        } else {
            LauFitter::fitter().releaseSecondStageParameters();
            this->startNewFit( LauFitter::fitter().nParameters(),
                               LauFitter::fitter().nFreeParameters() );
            fitResult = LauFitter::fitter().minimise();
        }
    }

    const TMatrixD& covMat = LauFitter::fitter().covarianceMatrix();
    this->storeFitStatus( fitResult, covMat );

    // Store the final fit results and errors into protected internal vectors that
    // all sub-classes can use within their own finalFitResults implementation
    // used below (e.g. putting them into an ntuple in a root file)
    LauFitter::fitter().updateParameters();
}

void LauSimFitCoordinator::setParsFromMinuit( Double_t* par, Int_t npar )
{
    // This function sets the internal parameters based on the values
    // that Minuit is using when trying to minimise the total likelihood function.

    // MINOS reports different numbers of free parameters depending on the
    // situation, so disable this check
    if ( ! this->withinAsymErrorCalc() ) {
        const UInt_t nFreePars = this->nFreeParams();
        if ( static_cast<UInt_t>( npar ) != nFreePars ) {
            std::cerr << "ERROR in LauSimFitCoordinator::setParsFromMinuit : Unexpected number of free parameters: "
                      << npar << ".\n";
            std::cerr << "                                              Expected: " << nFreePars
                      << ".\n"
                      << std::endl;
            gSystem->Exit( EXIT_FAILURE );
        }
    }

    // Despite npar being the number of free parameters
    // the par array actually contains all the parameters,
    // free and floating...
    // Update all the parameters with their new values.
    // Change the value in the array to be sent out to the tasks and the
    // parameters themselves (so that constraints are correctly calculated)
    for ( UInt_t i( 0 ); i < this->nTotParams(); ++i ) {
        if ( ! params_[i]->fixed() ) {
            parValues_[i] = par[i];
            params_[i]->value( par[i] );
        }
    }
}

Double_t LauSimFitCoordinator::getTotNegLogLikelihood()
{
    if ( socketMonitor_ == 0 ) {
        std::cerr << "ERROR in LauSimFitCoordinator::getTotNegLogLikelihood : Sockets not initialised."
                  << std::endl;
        return 0.0;
    }

    // Send current values of the parameters to the tasks.
    for ( UInt_t iTask( 0 ); iTask < nTasks_; ++iTask ) {

        std::vector<UInt_t>& indices     = taskIndices_[iTask];
        std::vector<UInt_t>& freeIndices = taskFreeIndices_[iTask];
        UInt_t nPars                     = indices.size();
        UInt_t nFreePars                 = freeIndices.size();
        for ( UInt_t iPar( 0 ); iPar < nPars; ++iPar ) {
            vectorPar_[iTask][iPar] = parValues_[indices[iPar]];
        }

        TMessage* message = messagesToTasks_[iTask];
        message->Reset( kMESS_ANY );
        message->WriteUInt( nPars );
        message->WriteUInt( nFreePars );
        message->WriteFastArray( vectorPar_[iTask], nPars );

        socketTasks_[iTask]->Send( *message );
    }

    Double_t negLogLike( 0.0 );
    TSocket* sActive( 0 );
    UInt_t responsesReceived( 0 );
    Bool_t allOK( kTRUE );
    while ( responsesReceived != nTasks_ ) {

        sActive = socketMonitor_->Select();
        sActive->Recv( messageFromTask_ );

        messageFromTask_->ReadDouble( vectorRes_[responsesReceived] );

        Double_t& nLL = vectorRes_[responsesReceived];
        if ( nLL == 0.0 || TMath::IsNaN( nLL ) || ! TMath::Finite( nLL ) ) {
            allOK = kFALSE;
        }

        negLogLike += vectorRes_[responsesReceived];

        ++responsesReceived;
    }

    // Calculate any penalty terms from Gaussian constrained variables
    const auto& multiDimCons = this->multiDimConstraints();
    if ( ! conVars_.empty() || ! multiDimCons.empty() ) {
        negLogLike += this->getLogLikelihoodPenalty();
    }

    const Double_t worstNegLogLike = -1.0 * this->worstLogLike();
    if ( ! allOK ) {
        std::cerr << "WARNING in LauSimFitCoordinator::getTotNegLogLikelihood : Strange NLL value returned by one or more tasks\n";
        std::cerr << "                                                   : Returning worst NLL found so far to force MINUIT out of this region."
                  << std::endl;
        negLogLike = worstNegLogLike;
    } else if ( negLogLike > worstNegLogLike ) {
        this->worstLogLike( -negLogLike );
    }

    return negLogLike;
}

Double_t LauSimFitCoordinator::getLogLikelihoodPenalty()
{
    Double_t penalty { 0.0 };

    for ( LauAbsRValue* par : conVars_ ) {
        penalty += par->constraintPenalty();
    }

    auto& multiDimCons = this->multiDimConstraints();
    for ( auto& constraint : multiDimCons ) {
        penalty += constraint.constraintPenalty();
    }

    return penalty;
}

void LauSimFitCoordinator::addConParameters()
{
    // Add penalties from the constraints to fit parameters

    // First, constraints on the fit parameters themselves
    for ( LauParameter* param : params_ ) {
        if ( param->gaussConstraint() ) {
            conVars_.push_back( param );
            std::cout << "INFO in LauSimFitCoordinator::addConParameters : Added Gaussian constraint to parameter "
                      << param->name() << std::endl;
        }
    }

    // Second, constraints on arbitrary combinations
    auto& conStore = this->formulaConstraints();
    for ( auto& constraint : conStore ) {
        std::vector<LauParameter*> params;
        for ( const auto& name : constraint.conPars_ ) {
            for ( LauParameter* par : params_ ) {
                if ( par->name() == name ) {
                    params.push_back( par );
                }
            }
        }

        // If the parameters are not found, skip it
        if ( params.size() != constraint.conPars_.size() ) {
            std::cerr << "WARNING in LauSimFitCoordinator::addConParameters: Could not find parameters to constrain in the formula... skipping"
                      << std::endl;
            continue;
        }

        constraint.formulaPar_ =
            std::make_unique<LauFormulaPar>( constraint.formula_, constraint.formula_, params );
        constraint.formulaPar_->addGaussianConstraint( constraint.mean_, constraint.width_ );

        conVars_.push_back( constraint.formulaPar_.get() );

        std::cout << "INFO in LauSimFitCoordinator::addConParameters : Added Gaussian constraint to formula\n";
        std::cout << "                                               : Formula: "
                  << constraint.formula_ << std::endl;
        for ( LauParameter* param : params ) {
            std::cout << "                                               : Parameter: "
                      << param->name() << std::endl;
        }
    }

    // Thirdly, add n-dimensional constraints
    auto& multiDimCons = this->multiDimConstraints();
    for ( auto& constraint : multiDimCons ) {
        for ( auto& parname : constraint.conPars_ ) {
            for ( auto& fitPar : params_ ) {
                if ( parname == fitPar->name() ) {
                    // Check parameters do not have a 1D Gaussian constraint applied
                    if ( fitPar->gaussConstraint() ) {
                        std::cerr << "ERROR in LauSimFitCoordinator::addConParameters: parameter in n-dimensional constraint already has a 1d constraint applied"
                                  << std::endl;
                        gSystem->Exit( EXIT_FAILURE );
                    }
                    constraint.conLauPars_.push_back( fitPar );
                }
            }
        }
        if ( constraint.conLauPars_.size() != constraint.conPars_.size() ) {
            std::cerr << "Error in LauSimFitCoordinator::addConParameters : Could not match parameter names for n-dimensional constraint"
                      << std::endl;
            gSystem->Exit( EXIT_FAILURE );
        }
    }
}

Bool_t LauSimFitCoordinator::finalise()
{
    if ( socketMonitor_ == 0 ) {
        std::cerr << "ERROR in LauSimFitCoordinator::finalise : Sockets not initialised."
                  << std::endl;
        return kFALSE;
    }

    // Prepare the covariance matrices
    const TMatrixD& covMatrix = this->covarianceMatrix();
    covMatrices_.resize( nTasks_ );

    LauParamFixed pred;

    std::map<UInt_t, UInt_t> freeParIndices;

    UInt_t counter( 0 );
    for ( UInt_t iPar( 0 ); iPar < this->nTotParams(); ++iPar ) {
        const LauParameter* par = params_[iPar];
        if ( ! pred( par ) ) {
            freeParIndices.insert( std::make_pair( iPar, counter ) );
            ++counter;
        }
    }

    for ( UInt_t iTask( 0 ); iTask < nTasks_; ++iTask ) {
        const UInt_t nPar = taskIndices_[iTask].size();

        std::vector<UInt_t> freeIndices;
        freeIndices.reserve( nPar );

        for ( UInt_t iPar( 0 ); iPar < nPar; ++iPar ) {
            UInt_t index                                = taskIndices_[iTask][iPar];
            std::map<UInt_t, UInt_t>::iterator freeIter = freeParIndices.find( index );
            if ( freeIter == freeParIndices.end() ) {
                continue;
            }
            UInt_t freeIndex = freeIter->second;
            freeIndices.push_back( freeIndex );
        }

        const UInt_t nFreePars = freeIndices.size();
        TMatrixD& covMat       = covMatrices_[iTask];
        covMat.ResizeTo( nFreePars, nFreePars );

        for ( UInt_t iPar( 0 ); iPar < nFreePars; ++iPar ) {
            for ( UInt_t jPar( 0 ); jPar < nFreePars; ++jPar ) {
                UInt_t i             = freeIndices[iPar];
                UInt_t j             = freeIndices[jPar];
                covMat( iPar, jPar ) = covMatrix( i, j );
            }
        }
    }

    // The array to hold the parameters
    TObjArray array;

    // Send messages to all tasks containing the final parameters and fit status, NLL
    for ( UInt_t iTask( 0 ); iTask < nTasks_; ++iTask ) {

        array.Clear();

        std::vector<UInt_t>& indices = taskIndices_[iTask];
        UInt_t nPars                 = indices.size();
        for ( UInt_t iPar( 0 ); iPar < nPars; ++iPar ) {
            array.Add( params_[indices[iPar]] );
        }

        const Int_t status = this->statusCode();
        const Double_t NLL = this->nll();
        const Double_t EDM = this->edm();
        TMatrixD& covMat   = covMatrices_[iTask];

        TMessage* message = messagesToTasks_[iTask];
        message->Reset( kMESS_OBJECT );
        message->WriteInt( status );
        message->WriteDouble( NLL );
        message->WriteDouble( EDM );
        message->WriteObject( &array );
        message->WriteObject( &covMat );

        socketTasks_[iTask]->Send( *message );
    }

    TSocket* sActive( 0 );
    UInt_t responsesReceived( 0 );
    Bool_t allOK( kTRUE );
    while ( responsesReceived != nTasks_ ) {

        // Get the next queued response
        sActive = socketMonitor_->Select();

        // Extract from the message the ID of the task and the number of events read
        sActive->Recv( messageFromTask_ );
        UInt_t iTask( 0 );
        Bool_t ok( kTRUE );
        messageFromTask_->ReadUInt( iTask );
        messageFromTask_->ReadBool( ok );

        if ( ok ) {
            TObjArray* objarray = dynamic_cast<TObjArray*>(
                messageFromTask_->ReadObject( messageFromTask_->GetClass() ) );
            if ( ! objarray ) {
                std::cerr << "ERROR in LauSimFitCoordinator::finalise : Error reading finalised parameters from task"
                          << std::endl;
                allOK = kFALSE;
            } else {
                // We want to auto-delete the supplied parameters since we only copy their values in this case
                objarray->SetOwner( kTRUE );

                const UInt_t nPars = objarray->GetEntries();
                if ( nPars != taskIndices_[iTask].size() ) {
                    std::cerr << "ERROR in LauSimFitCoordinator::finalise : Unexpected number of finalised parameters received from task"
                              << std::endl;
                    allOK = kFALSE;
                } else {
                    for ( UInt_t iPar( 0 ); iPar < nPars; ++iPar ) {
                        LauParameter* parameter = dynamic_cast<LauParameter*>( ( *objarray )[iPar] );
                        if ( ! parameter ) {
                            std::cerr
                                << "ERROR in LauSimFitCoordinator::finalise : Error reading parameter from task"
                                << std::endl;
                            allOK = kFALSE;
                            continue;
                        }

                        TString parname = parameter->name();

                        std::map<TString, UInt_t>::iterator iter = parIndices_.find( parname );
                        if ( iter == parIndices_.end() ) {
                            std::cerr
                                << "ERROR in LauSimFitCoordinator::finalise : Unexpected parameter name received from task"
                                << std::endl;
                            allOK = kFALSE;
                            continue;
                        }

                        const UInt_t index = iter->second;
                        if ( taskIndices_[iTask][iPar] != index ) {
                            std::cerr
                                << "ERROR in LauSimFitCoordinator::finalise : Unexpected parameter received from task"
                                << std::endl;
                            allOK = kFALSE;
                            continue;
                        }

                        Double_t parvalue = parameter->value();
                        params_[index]->value( parvalue );
                        parValues_[index]       = parvalue;
                        vectorPar_[iTask][iPar] = parvalue;
                    }
                }
                delete objarray;
            }
        } else {
            std::cerr << "ERROR in LauSimFitCoordinator::finalise : Task " << iTask
                      << " reports an error performing finalisation" << std::endl;
            allOK = kFALSE;
        }

        ++responsesReceived;
    }

    // Fill our ntuple as well
    if ( fitNtuple_ != 0 ) {
        for ( std::vector<LauParameter*>::iterator iter = params_.begin(); iter != params_.end();
              ++iter ) {
            if ( ! ( *iter )->fixed() ) {
                ( *iter )->updatePull();
            }
        }
        fitNtuple_->storeParsAndErrors( params_, this->multiDimConstrainedPars(), {} );
        fitNtuple_->storeCorrMatrix( this->iExpt(), this->fitStatus(), this->covarianceMatrix() );
        fitNtuple_->updateFitNtuple();
    }

    return allOK;
}

Bool_t LauSimFitCoordinator::writeOutResults()
{
    if ( socketMonitor_ == 0 ) {
        std::cerr << "ERROR in LauSimFitCoordinator::writeOutResults : Sockets not initialised."
                  << std::endl;
        return kFALSE;
    }

    // Construct a message, requesting to write out the fit results
    TString msgStr( "Write Results" );
    TMessage message( kMESS_STRING );
    message.WriteTString( msgStr );

    // Send the message to the tasks
    for ( UInt_t iTask( 0 ); iTask < nTasks_; ++iTask ) {
        socketTasks_[iTask]->Send( message );
    }

    TSocket* sActive( 0 );
    UInt_t responsesReceived( 0 );
    Bool_t allOK( kTRUE );
    while ( responsesReceived != nTasks_ ) {

        // Get the next queued response
        sActive = socketMonitor_->Select();

        // Extract from the message the ID of the task and the number of events read
        sActive->Recv( messageFromTask_ );
        UInt_t iTask( 0 );
        Bool_t ok( kTRUE );
        messageFromTask_->ReadUInt( iTask );
        messageFromTask_->ReadBool( ok );

        if ( ! ok ) {
            std::cerr << "ERROR in LauSimFitCoordinator::writeOutResults : Task " << iTask
                      << " reports an error performing finalisation" << std::endl;
            allOK = kFALSE;
        }

        ++responsesReceived;
    }

    // Write out our ntuple as well
    if ( fitNtuple_ != 0 ) {
        fitNtuple_->writeOutFitResults();
    }

    return allOK;
}
