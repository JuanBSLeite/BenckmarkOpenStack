
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

/*! \file LauRooFitTask.cc
    \brief File containing implementation of LauRooFitTask class.
*/

#include "LauRooFitTask.hh"

#include "LauFitNtuple.hh"
#include "LauParameter.hh"
#include "LauSimFitTask.hh"

#include "RooAbsPdf.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooFormulaVar.h"
#include "RooRealVar.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

#include <iostream>
#include <vector>

LauRooFitTask::LauRooFitTask( RooAbsPdf& model,
                              const Bool_t extended,
                              const RooArgSet& vars,
                              const TString& weightVarName ) :
    LauSimFitTask {},
    model_ { model },
    dataVars_ { vars },
    weightVarName_ { weightVarName },
    extended_ { extended }
{
}

void LauRooFitTask::cleanData()
{
    if ( dataFile_ ) {
        dataFile_->Close();
    }

    dataFile_.reset();
    dataTree_ = nullptr;

    exptData_.reset();
}

void LauRooFitTask::initialise()
{
    if ( weightVarName_ == "" ) {
        return;
    }

    Bool_t weightVarFound { kFALSE };
    for ( RooAbsArg* param : dataVars_ ) {
        TString name = param->GetName();
        if ( name == weightVarName_ ) {
            weightVarFound = kTRUE;
            break;
        }
    }
    if ( ! weightVarFound ) {
        std::cerr << "ERROR in LauRooFitTask::initialise : The set of data variables does not contain the weighting variable \""
                  << weightVarName_ << std::endl;
        std::cerr << "                                   : Weighting will be disabled." << std::endl;
        weightVarName_ = "";
    }
}

Bool_t LauRooFitTask::verifyFitData( const TString& dataFileName, const TString& dataTreeName )
{
    // Clean-up from any previous runs
    if ( dataFile_ ) {
        this->cleanData();
    }

    // Open the data file
    dataFile_.reset( TFile::Open( dataFileName ) );
    if ( ! dataFile_ ) {
        std::cerr << "ERROR in LauRooFitTask::verifyFitData : Problem opening data file \""
                  << dataFileName << "\"" << std::endl;
        return kFALSE;
    }

    // Retrieve the tree
    dataTree_ = dynamic_cast<TTree*>( dataFile_->Get( dataTreeName ) );
    if ( ! dataTree_ ) {
        std::cerr << "ERROR in LauRooFitTask::verifyFitData : Problem retrieving tree \""
                  << dataTreeName << "\" from data file \"" << dataFileName << "\"" << std::endl;
        this->cleanData();
        return kFALSE;
    }

    // Check that the tree contains branches for all the fit variables
    Bool_t allOK { kTRUE };
    for ( RooAbsArg* param : dataVars_ ) {
        TString name    = param->GetName();
        TBranch* branch = dataTree_->GetBranch( name );
        if ( ! branch ) {
            std::cerr << "ERROR in LauRooFitTask::verifyFitData : The data tree does not contain a branch for fit variable \""
                      << name << std::endl;
            allOK = kFALSE;
        }
    }
    if ( ! allOK ) {
        return kFALSE;
    }

    // Check whether the tree has the branch iExpt
    TBranch* branch = dataTree_->GetBranch( "iExpt" );
    if ( ! branch ) {
        std::cout << "WARNING in LauRooFitTask::verifyFitData : Cannot find branch \"iExpt\" in the tree, will treat all data as being from a single experiment"
                  << std::endl;
        iExptSet_.clear();
    } else {
        // Define the valid values for iExpt
        iExptSet_.clear();
        const UInt_t firstExp = dataTree_->GetMinimum( "iExpt" );
        const UInt_t lastExp  = dataTree_->GetMaximum( "iExpt" );
        for ( UInt_t iExp = firstExp; iExp <= lastExp; ++iExp ) {
            iExptSet_.insert( iExp );
        }
    }

    return kTRUE;
}

void LauRooFitTask::prepareInitialParArray( TObjArray& array )
{
    // Check that the NLL variable has been initialised
    if ( ! nllVar_ ) {
        std::cerr << "ERROR in LauRooFitTask::prepareInitialParArray : NLL var not initialised"
                  << std::endl;
        return;
    }

    // If we already prepared the entries in the fitPars_ vector then we only need to add the contents to the array
    if ( ! fitPars_.empty() ) {
        for ( LauParameter* par : fitPars_ ) {
            array.Add( par );
        }
        return;
    }

    // Store the set of parameters and the total number of parameters
    std::unique_ptr<RooArgSet> varSet { nllVar_->getParameters( *exptData_ ) };
    UInt_t nFreePars { 0 };

    // Loop through the fit parameters
    for ( RooAbsArg* param : *varSet ) {
        // Only consider the free parameters
        if ( ! param->isConstant() ) {
            // Add the parameter
            RooRealVar* rrvar = dynamic_cast<RooRealVar*>( param );
            if ( rrvar ) {
                // Count the number of free parameters
                ++nFreePars;
                // Do the conversion and add it to the array
                LauParameter* lpar = this->convertToLauParameter( rrvar );
                fitVars_.push_back( rrvar );
                fitPars_.push_back( lpar );
                array.Add( lpar );
            } else {
                RooFormulaVar* rfvar = dynamic_cast<RooFormulaVar*>( param );
                if ( ! rfvar ) {
                    std::cerr << "ERROR in LauRooFitTask::prepareInitialParArray : The parameter is neither a RooRealVar nor a RooFormulaVar, don't know what to do"
                              << std::endl;
                    continue;
                }
                auto lpars = this->convertToLauParameters( rfvar );
                for ( auto [rrv, lpar] : lpars ) {
                    if ( ! rrv->isConstant() ) {
                        continue;
                    }

                    // Count the number of free parameters
                    ++nFreePars;
                    // Add the parameter to the array
                    fitVars_.push_back( rrvar );
                    fitPars_.push_back( lpar );
                    array.Add( lpar );
                }
            }
        }
    }

    this->startNewFit( nFreePars, nFreePars );
}

LauParameter* LauRooFitTask::convertToLauParameter( const RooRealVar* rooParameter ) const
{
    // TODO - can we make these unique_ptr's
    return new LauParameter( rooParameter->GetName(),
                             rooParameter->getVal(),
                             rooParameter->getMin(),
                             rooParameter->getMax(),
                             rooParameter->isConstant() );
}

std::vector<std::pair<RooRealVar*, LauParameter*>> LauRooFitTask::convertToLauParameters(
    const RooFormulaVar* rooFormula ) const
{
    // Create the empty vector
    std::vector<std::pair<RooRealVar*, LauParameter*>> lauParameters;

    Int_t parIndex { 0 };
    RooAbsArg* rabsarg { nullptr };
    RooRealVar* rrvar { nullptr };
    RooFormulaVar* rfvar { nullptr };
    // Loop through all the parameters of the formula
    while ( ( rabsarg = rooFormula->getParameter( parIndex ) ) ) {
        // First try converting to a RooRealVar
        rrvar = dynamic_cast<RooRealVar*>( rabsarg );
        if ( rrvar ) {
            // Do the conversion and add it to the array
            LauParameter* lpar = this->convertToLauParameter( rrvar );
            lauParameters.push_back( std::make_pair( rrvar, lpar ) );
            continue;
        }

        // If that didn't work, try converting to a RooFormulaVar
        rfvar = dynamic_cast<RooFormulaVar*>( rabsarg );
        if ( rfvar ) {
            // Do the conversion and add these to the array
            auto lpars = this->convertToLauParameters( rfvar );
            lauParameters.insert( lauParameters.end(), lpars.begin(), lpars.end() );
            continue;
        }

        // If neither of those worked we don't know what to do, so print an error message and continue
        std::cerr << "ERROR in LauRooFitTask::convertToLauParameters : One of the parameters is not a RooRealVar nor a RooFormulaVar, it is a: "
                  << rabsarg->ClassName() << std::endl;
        std::cerr << "                                                : Do not know how to process that - it will be skipped."
                  << std::endl;
    }

    return lauParameters;
}

Double_t LauRooFitTask::getTotNegLogLikelihood()
{
    Double_t nLL = nllVar_ ? nllVar_->getVal() : 0.0;
    return nLL;
}

void LauRooFitTask::setParsFromMinuit( Double_t* par, Int_t npar )
{
    // This function sets the internal parameters based on the values
    // that Minuit is using when trying to minimise the total likelihood function.

    // MINOS reports different numbers of free parameters depending on the
    // situation, so disable this check
    const UInt_t nFreePars = this->nFreeParams();
    if ( ! this->withinAsymErrorCalc() ) {
        if ( static_cast<UInt_t>( npar ) != nFreePars ) {
            std::cerr << "ERROR in LauRooFitTask::setParsFromMinuit : Unexpected number of free parameters: "
                      << npar << ".\n";
            std::cerr << "                                             Expected: " << nFreePars
                      << ".\n"
                      << std::endl;
            gSystem->Exit( EXIT_FAILURE );
        }
    }

    // Despite npar being the number of free parameters
    // the par array actually contains all the parameters,
    // free and floating...

    // Update all the floating ones with their new values
    for ( UInt_t i { 0 }; i < nFreePars; ++i ) {
        if ( ! fitPars_[i]->fixed() ) {
            // Set both the RooRealVars and the LauParameters
            fitPars_[i]->value( par[i] );
            fitVars_[i]->setVal( par[i] );
        }
    }
}

UInt_t LauRooFitTask::readExperimentData()
{
    // check that we're being asked to read a valid index
    const UInt_t exptIndex = this->iExpt();
    if ( iExptSet_.empty() && exptIndex != 0 ) {
        std::cerr << "ERROR in LauRooFitTask::readExperimentData : Invalid experiment number "
                  << exptIndex << ", data contains only one experiment" << std::endl;
        return 0;
    } else if ( ! iExptSet_.empty() && iExptSet_.find( exptIndex ) == iExptSet_.end() ) {
        std::cerr << "ERROR in LauRooFitTask::readExperimentData : Invalid experiment number "
                  << exptIndex << std::endl;
        return 0;
    }

    // cleanup the data from any previous experiment
    exptData_.reset();

    // retrieve the data and find out how many events have been read
    if ( iExptSet_.empty() ) {
        if ( weightVarName_ != "" ) {
            exptData_ = std::make_unique<RooDataSet>( TString::Format( "expt%dData", exptIndex ),
                                                      "",
                                                      dataVars_,
                                                      RooFit::Import( *dataTree_ ),
                                                      RooFit::WeightVar( weightVarName_.Data() ) );
        } else {
            exptData_ = std::make_unique<RooDataSet>( TString::Format( "expt%dData", exptIndex ),
                                                      "",
                                                      dataVars_,
                                                      RooFit::Import( *dataTree_ ) );
        }
    } else {
        const TString selectionString { TString::Format( "iExpt==%d", exptIndex ) };
        TTree* exptTree { dataTree_->CopyTree( selectionString ) };
        if ( weightVarName_ != "" ) {
            exptData_ = std::make_unique<RooDataSet>( TString::Format( "expt%dData", exptIndex ),
                                                      "",
                                                      dataVars_,
                                                      RooFit::Import( *exptTree ),
                                                      RooFit::WeightVar( weightVarName_.Data() ) );
        } else {
            exptData_ = std::make_unique<RooDataSet>( TString::Format( "expt%dData", exptIndex ),
                                                      "",
                                                      dataVars_,
                                                      RooFit::Import( *exptTree ) );
        }
        delete exptTree;
    }

    const UInt_t nEvent = exptData_->numEntries();
    this->eventsPerExpt( nEvent );
    return nEvent;
}

void LauRooFitTask::cacheInputFitVars()
{
    // construct the new NLL variable for this dataset
    nllVar_.reset( model_.createNLL( *exptData_, RooFit::Extended( extended_ ) ) );
}

void LauRooFitTask::finaliseExperiment( const LauAbsFitter::FitStatus& fitStat,
                                        const TObjArray* parsFromCoordinator,
                                        const TMatrixD* covMat,
                                        TObjArray& parsToCoordinator )
{
    // Copy the fit status information
    this->storeFitStatus( fitStat, *covMat );

    // Now process the parameters
    const UInt_t nFreePars = this->nFreeParams();
    UInt_t nPars           = parsFromCoordinator->GetEntries();
    if ( nPars != nFreePars ) {
        std::cerr << "ERROR in LauRooFitTask::finaliseExperiment : Unexpected number of parameters received from coordinator"
                  << std::endl;
        std::cerr << "                                            : Received " << nPars
                  << " when expecting " << nFreePars << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    for ( UInt_t iPar { 0 }; iPar < nPars; ++iPar ) {
        LauParameter* parameter = dynamic_cast<LauParameter*>( ( *parsFromCoordinator )[iPar] );
        if ( ! parameter ) {
            std::cerr << "ERROR in LauRooFitTask::finaliseExperiment : Error reading parameter from coordinator"
                      << std::endl;
            gSystem->Exit( EXIT_FAILURE );
        }

        if ( parameter->name() != fitPars_[iPar]->name() ) {
            std::cerr << "ERROR in LauRooFitTask::finaliseExperiment : Error reading parameter from coordinator"
                      << std::endl;
            gSystem->Exit( EXIT_FAILURE );
        }

        *( fitPars_[iPar] ) = *parameter;

        RooRealVar* rrv = fitVars_[iPar];
        rrv->setVal( parameter->value() );
        rrv->setError( parameter->error() );
        rrv->setAsymError( parameter->negError(), parameter->posError() );
    }

    // Update the pulls and add each finalised fit parameter to the list to
    // send back to the coordinator
    for ( LauParameter* par : fitPars_ ) {
        par->updatePull();
        parsToCoordinator.Add( par );
    }

    // Write the results into the ntuple
    LauFitNtuple* ntuple = this->fitNtuple();
    ntuple->storeParsAndErrors( fitPars_, this->multiDimConstrainedPars(), {} );

    // find out the correlation matrix for the parameters
    ntuple->storeCorrMatrix( this->iExpt(), this->fitStatus(), this->covarianceMatrix() );

    // Fill the data into ntuple
    ntuple->updateFitNtuple();
}

void LauRooFitTask::generate( const TString& dataFileName,
                              const TString& dataTreeName,
                              const UInt_t nExpt,
                              const UInt_t firstExpt,
                              const Int_t nEventsPerExpt )
{
    if ( nExpt == 0 ) {
        return;
    }

    // Create the file in which to store the generated data
    std::unique_ptr<TFile> file { TFile::Open( dataFileName, "recreate" ) };
    if ( ! file ) {
        std::cerr << "ERROR in LauRooFitTask::generate : problem opening file \"" << dataFileName
                  << "\" for writing" << std::endl;
        return;
    }

    // Create the variable for holding the experiment number index
    // and define all experiments as being valid index values
    RooCategory exptIndex( "iExpt", "" );
    for ( UInt_t iExpt { firstExpt }; iExpt < ( firstExpt + nExpt ); iExpt++ ) {
        exptIndex.defineType( Form( "Experiment %u", iExpt ), iExpt );
    }

    RooCmdArg cmdArg = extended_ ? RooFit::Extended() : RooFit::NumEvents( nEventsPerExpt );

    // Generate data for the first experiment
    std::unique_ptr<RooDataSet> data { model_.generate( dataVars_, cmdArg ) };
    exptIndex.setIndex( firstExpt );
    data->addColumn( exptIndex );

    // Generate data for all subsequent experiments and append to the main dataset
    for ( UInt_t iExpt { firstExpt + 1 }; iExpt < ( firstExpt + nExpt ); iExpt++ ) {

        std::unique_ptr<RooDataSet> tmpdata { model_.generate( dataVars_, cmdArg ) };
        exptIndex.setIndex( iExpt );
        tmpdata->addColumn( exptIndex );

        data->append( *tmpdata );
    }

    // Create a TTree from the dataset
    TTree* tree = data->GetClonedTree();

    // Rename the tree to its requested name
    tree->SetName( dataTreeName );

    // Rename the iExpt branch/leaf since it gets "_idx" appended for some reason
    tree->GetBranch( "iExpt_idx" )->SetNameTitle( "iExpt", "iExpt/I" );
    tree->GetLeaf( "iExpt_idx" )->SetNameTitle( "iExpt", "iExpt" );

    // Attach the tree to the file and write it
    tree->SetDirectory( file.get() );
    file->Write();
    file->Close();
}
