
/*
Copyright 2004 University of Warwick

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

/*! \file LauAbsFitModel.cc
    \brief File containing implementation of LauAbsFitModel class.
*/

#include "LauAbsFitModel.hh"

#include "LauAbsFitter.hh"
#include "LauAbsPdf.hh"
#include "LauComplex.hh"
#include "LauFitDataTree.hh"
#include "LauFitter.hh"
#include "LauGenNtuple.hh"
#include "LauParamFixed.hh"
#include "LauParameter.hh"
#include "LauPrint.hh"
#include "LauSPlot.hh"

#include "TMath.h"
#include "TMessage.h"
#include "TMonitor.h"
#include "TSystem.h"
#include "TVirtualFitter.h"

#include <iostream>
#include <limits>
#include <vector>

LauAbsFitModel::LauAbsFitModel() :
    fixParams_( kFALSE ),
    compareFitData_( kFALSE ),
    savePDF_( kFALSE ),
    writeLatexTable_( kFALSE ),
    writeSPlotData_( kFALSE ),
    storeDPEff_( kFALSE ),
    randomFit_( kFALSE ),
    emlFit_( kFALSE ),
    poissonSmear_( kFALSE ),
    enableEmbedding_( kFALSE ),
    usingDP_( kTRUE ),
    pdfsDependOnDP_( kFALSE ),
    inputFitData_( 0 ),
    genNtuple_( 0 ),
    sPlotNtuple_( 0 ),
    nullString_( "" ),
    doSFit_( kFALSE ),
    sWeightBranchName_( "" ),
    sWeightScaleFactor_( 1.0 ),
    outputTableName_( "" ),
    fitToyMCFileName_( "fitToyMC.root" ),
    fitToyMCTableName_( "fitToyMCTable" ),
    fitToyMCScale_( 10 ),
    fitToyMCPoissonSmear_( kFALSE ),
    sPlotFileName_( "" ),
    sPlotTreeName_( "" ),
    sPlotVerbosity_( "" )
{
}

LauAbsFitModel::~LauAbsFitModel()
{
    delete inputFitData_;
    inputFitData_ = 0;
    delete genNtuple_;
    genNtuple_ = 0;
    delete sPlotNtuple_;
    sPlotNtuple_ = 0;
}

void LauAbsFitModel::run( const TString& applicationCode,
                          const TString& dataFileName,
                          const TString& dataTreeName,
                          const TString& histFileName,
                          const TString& tableFileName )
{
    // Chose whether you want to generate or fit events in the Dalitz plot.
    // To generate events choose applicationCode = "gen", to fit events choose
    // applicationCode = "fit".

    TString runCode( applicationCode );
    runCode.ToLower();

    TString histFileNameCopy( histFileName );
    TString tableFileNameCopy( tableFileName );
    TString dataFileNameCopy( dataFileName );
    TString dataTreeNameCopy( dataTreeName );

    // Initialise the fit par vectors. Each class that inherits from this one
    // must implement this sensibly for all vectors specified in clearFitParVectors,
    // i.e. specify parameter names, initial, min, max and fixed values
    this->initialise();

    // Add variables to Gaussian constrain to a list
    this->addConParameters();

    if ( dataFileNameCopy == "" ) {
        dataFileNameCopy = "data.root";
    }
    if ( dataTreeNameCopy == "" ) {
        dataTreeNameCopy = "genResults";
    }

    if ( runCode.Contains( "gen" ) ) {

        if ( histFileNameCopy == "" ) {
            histFileNameCopy = "parInfo.root";
        }
        if ( tableFileNameCopy == "" ) {
            tableFileNameCopy = "genResults";
        }

        this->setGenValues();
        this->generate( dataFileNameCopy, dataTreeNameCopy, histFileNameCopy, tableFileNameCopy );

    } else if ( runCode.Contains( "fit" ) ) {

        if ( histFileNameCopy == "" ) {
            histFileNameCopy = "parInfo.root";
        }
        if ( tableFileNameCopy == "" ) {
            tableFileNameCopy = "fitResults";
        }

        this->fit( dataFileNameCopy, dataTreeNameCopy, histFileNameCopy, tableFileNameCopy );

    } else if ( runCode.Contains( "plot" ) ) {

        this->savePDFPlots( "plot" );

    } else if ( runCode.Contains( "weight" ) ) {

        this->weightEvents( dataFileNameCopy, dataTreeNameCopy );
    }
}

void LauAbsFitModel::doSFit( const TString& sWeightBranchName, Double_t scaleFactor )
{
    if ( sWeightBranchName == "" ) {
        std::cerr << "WARNING in LauAbsFitModel::doSFit : sWeight branch name is empty string, not setting-up sFit."
                  << std::endl;
        return;
    }

    doSFit_             = kTRUE;
    sWeightBranchName_  = sWeightBranchName;
    sWeightScaleFactor_ = scaleFactor;
}

void LauAbsFitModel::setBkgndClassNames( const std::vector<TString>& names )
{
    if ( ! bkgndClassNames_.empty() ) {
        std::cerr << "WARNING in LauAbsFitModel::setBkgndClassNames : Names already stored, not changing them."
                  << std::endl;
        return;
    }

    UInt_t nBkgnds = names.size();
    for ( UInt_t i( 0 ); i < nBkgnds; ++i ) {
        bkgndClassNames_.insert( std::make_pair( i, names[i] ) );
    }

    this->setupBkgndVectors();
}

Bool_t LauAbsFitModel::validBkgndClass( const TString& className ) const
{
    if ( bkgndClassNames_.empty() ) {
        return kFALSE;
    }

    Bool_t found( kFALSE );
    for ( LauBkgndClassMap::const_iterator iter = bkgndClassNames_.begin();
          iter != bkgndClassNames_.end();
          ++iter ) {
        if ( iter->second == className ) {
            found = kTRUE;
            break;
        }
    }

    return found;
}

UInt_t LauAbsFitModel::bkgndClassID( const TString& className ) const
{
    if ( ! this->validBkgndClass( className ) ) {
        std::cerr << "ERROR in LauAbsFitModel::bkgndClassID : Request for ID for invalid background class \""
                  << className << "\"." << std::endl;
        return ( bkgndClassNames_.size() + 1 );
    }

    UInt_t bgID( 0 );
    for ( LauBkgndClassMap::const_iterator iter = bkgndClassNames_.begin();
          iter != bkgndClassNames_.end();
          ++iter ) {
        if ( iter->second == className ) {
            bgID = iter->first;
            break;
        }
    }

    return bgID;
}

const TString& LauAbsFitModel::bkgndClassName( UInt_t classID ) const
{
    LauBkgndClassMap::const_iterator iter = bkgndClassNames_.find( classID );

    if ( iter == bkgndClassNames_.end() ) {
        std::cerr << "ERROR in LauAbsFitModel::bkgndClassName : Request for name of invalid background class ID "
                  << classID << "." << std::endl;
        return nullString_;
    }

    return iter->second;
}

void LauAbsFitModel::clearFitParVectors()
{
    std::cout << "INFO in LauAbsFitModel::clearFitParVectors : Clearing fit variable vectors"
              << std::endl;

    conVars_.clear();
    resVars_.clear();
    fitVars_.clear();
}

void LauAbsFitModel::clearExtraVarVectors()
{
    std::cout << "INFO in LauAbsFitModel::clearExtraVarVectors : Clearing extra ntuple variable vectors"
              << std::endl;
    extraVars_.clear();
}

void LauAbsFitModel::setGenValues()
{
    // makes sure each parameter holds its genValue as its current value
    for ( LauParameterPList::iterator iter = fitVars_.begin(); iter != fitVars_.end(); ++iter ) {
        ( *iter )->value( ( *iter )->genValue() );
    }
    this->propagateParUpdates();
}

void LauAbsFitModel::writeSPlotData( const TString& fileName,
                                     const TString& treeName,
                                     Bool_t storeDPEfficiency,
                                     const TString& verbosity )
{
    if ( this->writeSPlotData() ) {
        std::cerr << "ERROR in LauAbsFitModel::writeSPlotData : Already have an sPlot ntuple setup, not doing it again."
                  << std::endl;
        return;
    }
    writeSPlotData_ = kTRUE;
    sPlotFileName_  = fileName;
    sPlotTreeName_  = treeName;
    sPlotVerbosity_ = verbosity;
    storeDPEff_     = storeDPEfficiency;
}

// TODO : histFileName isn't used here at the moment but could be used for
//        storing the values of the parameters used in the generation.
//        These could then be read and used for setting the "true" values
//        in a subsequent fit.
void LauAbsFitModel::generate( const TString& dataFileName,
                               const TString& dataTreeName,
                               const TString& /*histFileName*/,
                               const TString& tableFileNameBase )
{
    // Create the ntuple for storing the results
    std::cout << "INFO in LauAbsFitModel::generate : Creating generation ntuple." << std::endl;
    if ( genNtuple_ != 0 ) {
        delete genNtuple_;
        genNtuple_ = 0;
    }
    genNtuple_ = new LauGenNtuple( dataFileName, dataTreeName );

    // add branches for storing the experiment number and the number of
    // the event within the current experiment
    this->addGenNtupleIntegerBranch( "iExpt" );
    this->addGenNtupleIntegerBranch( "iEvtWithinExpt" );
    this->setupGenNtupleBranches();

    // Start the cumulative timer
    cumulTimer_.Start();

    const UInt_t firstExp = this->firstExpt();
    const UInt_t nExp     = this->nExpt();

    Bool_t genOK( kTRUE );
    do {
        // Loop over the number of experiments
        for ( UInt_t iExp = firstExp; iExp < ( firstExp + nExp ); ++iExp ) {

            // Start the timer to see how long each experiment takes to generate
            timer_.Start();

            // Store the experiment number in the ntuple
            this->setGenNtupleIntegerBranchValue( "iExpt", iExp );

            // Do the generation for this experiment
            std::cout << "INFO in LauAbsFitModel::generate : Generating experiment number " << iExp
                      << std::endl;
            genOK = this->genExpt();

            // Stop the timer and see how long the program took so far
            timer_.Stop();
            timer_.Print();

            if ( ! genOK ) {
                // delete and recreate an empty tree
                genNtuple_->deleteAndRecreateTree();

                // then break out of the experiment loop
                std::cerr << "WARNING in LauAbsFitModel::generate : Problem in toy MC generation.  Starting again with updated parameters..."
                          << std::endl;
                break;
            }

            if ( this->writeLatexTable() ) {
                TString tableFileName( tableFileNameBase );
                tableFileName += "_";
                tableFileName += iExp;
                tableFileName += ".tex";
                this->writeOutTable( tableFileName );
            }

        }    // Loop over number of experiments
    } while ( ! genOK );

    // Print out total timing info.
    cumulTimer_.Stop();
    std::cout << "INFO in LauAbsFitModel::generate : Finished generating all experiments."
              << std::endl;
    std::cout << "INFO in LauAbsFitModel::generate : Cumulative timing:" << std::endl;
    cumulTimer_.Print();

    // Build the event index
    std::cout << "INFO in LauAbsFitModel::generate : Building experiment:event index." << std::endl;
    // TODO - can test this return value?
    //Int_t nIndexEntries =
    genNtuple_->buildIndex( "iExpt", "iEvtWithinExpt" );

    // Write out toy MC ntuple
    std::cout << "INFO in LauAbsFitModel::generate : Writing data to file " << dataFileName << "."
              << std::endl;
    genNtuple_->writeOutGenResults();
}

void LauAbsFitModel::addGenNtupleIntegerBranch( const TString& name )
{
    genNtuple_->addIntegerBranch( name );
}

void LauAbsFitModel::addGenNtupleDoubleBranch( const TString& name )
{
    genNtuple_->addDoubleBranch( name );
}

void LauAbsFitModel::setGenNtupleIntegerBranchValue( const TString& name, Int_t value )
{
    genNtuple_->setIntegerBranchValue( name, value );
}

void LauAbsFitModel::setGenNtupleDoubleBranchValue( const TString& name, Double_t value )
{
    genNtuple_->setDoubleBranchValue( name, value );
}

Int_t LauAbsFitModel::getGenNtupleIntegerBranchValue( const TString& name ) const
{
    return genNtuple_->getIntegerBranchValue( name );
}

Double_t LauAbsFitModel::getGenNtupleDoubleBranchValue( const TString& name ) const
{
    return genNtuple_->getDoubleBranchValue( name );
}

void LauAbsFitModel::fillGenNtupleBranches()
{
    genNtuple_->fillBranches();
}

void LauAbsFitModel::addSPlotNtupleIntegerBranch( const TString& name )
{
    sPlotNtuple_->addIntegerBranch( name );
}

void LauAbsFitModel::addSPlotNtupleDoubleBranch( const TString& name )
{
    sPlotNtuple_->addDoubleBranch( name );
}

void LauAbsFitModel::setSPlotNtupleIntegerBranchValue( const TString& name, Int_t value )
{
    sPlotNtuple_->setIntegerBranchValue( name, value );
}

void LauAbsFitModel::setSPlotNtupleDoubleBranchValue( const TString& name, Double_t value )
{
    sPlotNtuple_->setDoubleBranchValue( name, value );
}

void LauAbsFitModel::fillSPlotNtupleBranches()
{
    sPlotNtuple_->fillBranches();
}

void LauAbsFitModel::fit( const TString& dataFileName,
                          const TString& dataTreeName,
                          const TString& histFileName,
                          const TString& tableFileNameBase )
{
    // Routine to perform the total fit.

    const UInt_t firstExp = this->firstExpt();
    const UInt_t nExp     = this->nExpt();

    std::cout << "INFO in LauAbsFitModel::fit : First experiment = " << firstExp << std::endl;
    std::cout << "INFO in LauAbsFitModel::fit : Number of experiments = " << nExp << std::endl;

    // Start the cumulative timer
    cumulTimer_.Start();

    this->resetFitCounters();

    // Create and setup the fit results ntuple
    this->setupResultsOutputs( histFileName, tableFileNameBase );

    // Create and setup the sPlot ntuple
    if ( this->writeSPlotData() ) {
        std::cout << "INFO in LauAbsFitModel::fit : Creating sPlot ntuple." << std::endl;
        if ( sPlotNtuple_ != 0 ) {
            delete sPlotNtuple_;
            sPlotNtuple_ = 0;
        }
        sPlotNtuple_ = new LauGenNtuple( sPlotFileName_, sPlotTreeName_ );
        this->setupSPlotNtupleBranches();
    }

    // This reads in the given dataFile and creates an input
    // fit data tree that stores them for all events and experiments.
    Bool_t dataOK = this->verifyFitData( dataFileName, dataTreeName );
    if ( ! dataOK ) {
        std::cerr << "ERROR in LauAbsFitModel::fit : Problem caching the fit data." << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Loop over the number of experiments
    for ( UInt_t iExp = firstExp; iExp < ( firstExp + nExp ); ++iExp ) {

        // Start the timer to see how long each fit takes
        timer_.Start();

        this->setCurrentExperiment( iExp );

        UInt_t nEvents = this->readExperimentData();
        if ( nEvents < 1 ) {
            std::cerr << "WARNING in LauAbsFitModel::fit : Zero events in experiment " << iExp
                      << ", skipping..." << std::endl;
            timer_.Stop();
            continue;
        }

        // Now the sub-classes must implement whatever they need to do
        // to cache any more input fit data they need in order to
        // calculate the likelihoods during the fit.
        // They need to use the inputFitData_ tree as input. For example,
        // inputFitData_ contains m13Sq and m23Sq. The appropriate fit model
        // then caches the resonance dynamics for the signal model, as
        // well as the background likelihood values in the Dalitz plot
        this->cacheInputFitVars();
        if ( this->doSFit() ) {
            this->cacheInputSWeights();
        }

        // If we're fitting toy experiments then re-generate the means of any constraints
        this->generateConstraintMeans( conVars_ );

        // Do the fit for this experiment
        this->fitExpt();

        // Write the results into the ntuple
        this->finaliseFitResults( outputTableName_ );

        // Stop the timer and see how long the program took so far
        timer_.Stop();
        timer_.Print();

        // Store the per-event likelihood values
        if ( this->writeSPlotData() ) {
            this->storePerEvtLlhds();
        }

        // Create a toy MC sample using the fitted parameters so that
        // the user can compare the fit to the data.
        if ( compareFitData_ == kTRUE && this->statusCode() == 3 ) {
            this->createFitToyMC( fitToyMCFileName_, fitToyMCTableName_ );
        }

    }    // Loop over number of experiments

    // Print out total timing info.
    cumulTimer_.Stop();
    std::cout << "INFO in LauAbsFitModel::fit : Cumulative timing:" << std::endl;
    cumulTimer_.Print();

    // Print out stats on OK fits.
    const UInt_t nOKFits  = this->numberOKFits();
    const UInt_t nBadFits = this->numberBadFits();
    std::cout << "INFO in LauAbsFitModel::fit : Number of OK Fits = " << nOKFits << std::endl;
    std::cout << "INFO in LauAbsFitModel::fit : Number of Failed Fits = " << nBadFits << std::endl;
    Double_t fitEff( 0.0 );
    if ( nExp != 0 ) {
        fitEff = nOKFits / ( 1.0 * nExp );
    }
    std::cout << "INFO in LauAbsFitModel::fit : Fit efficiency = " << fitEff * 100.0 << "%."
              << std::endl;

    // Write out any fit results (ntuples etc...).
    this->writeOutAllFitResults();
    if ( this->writeSPlotData() ) {
        this->calculateSPlotData();
    }
}

void LauAbsFitModel::setupResultsOutputs( const TString& histFileName, const TString& tableFileName )
{
    this->LauSimFitTask::setupResultsOutputs( histFileName, tableFileName );

    outputTableName_ = tableFileName;
}

Bool_t LauAbsFitModel::verifyFitData( const TString& dataFileName, const TString& dataTreeName )
{
    // From the input data stream, store the variables into the
    // internal tree inputFitData_ that can be used by the sub-classes
    // in calculating their likelihood functions for the fit
    delete inputFitData_;
    inputFitData_ = new LauFitDataTree( dataFileName, dataTreeName );
    Bool_t dataOK = inputFitData_->findBranches();

    if ( ! dataOK ) {
        delete inputFitData_;
        inputFitData_ = 0;
    }

    return dataOK;
}

void LauAbsFitModel::cacheInputSWeights()
{
    Bool_t hasBranch = inputFitData_->haveBranch( sWeightBranchName_ );
    if ( ! hasBranch ) {
        std::cerr << "ERROR in LauAbsFitModel::cacheInputSWeights : Input data does not contain variable \""
                  << sWeightBranchName_ << "\".\n";
        std::cerr << "                                            : Turning off sFit!" << std::endl;
        doSFit_ = kFALSE;
        return;
    }

    UInt_t nEvents = this->eventsPerExpt();
    sWeights_.clear();
    sWeights_.reserve( nEvents );

    for ( UInt_t iEvt = 0; iEvt < nEvents; ++iEvt ) {

        const LauFitData& dataValues = inputFitData_->getData( iEvt );

        LauFitData::const_iterator iter = dataValues.find( sWeightBranchName_ );

        sWeights_.push_back( iter->second * sWeightScaleFactor_ );
    }
}

void LauAbsFitModel::fitExpt()
{
    // Routine to perform the actual fit for the given experiment

    // Update initial fit parameters if required (e.g. if using random numbers).
    this->checkInitFitParams();

    // Initialise the fitter
    LauFitter::fitter().useAsymmFitErrors( this->useAsymmFitErrors() );
    LauFitter::fitter().twoStageFit( this->twoStageFit() );
    LauFitter::fitter().initialise( this, fitVars_ );

    this->startNewFit( LauFitter::fitter().nParameters(), LauFitter::fitter().nFreeParameters() );

    // Now ready for minimisation step
    std::cout << "\nINFO in LauAbsFitModel::fitExpt : Start minimisation...\n";
    LauAbsFitter::FitStatus fitResult = LauFitter::fitter().minimise();

    // If we're doing a two stage fit we can now release (i.e. float)
    // the 2nd stage parameters and re-fit
    if ( this->twoStageFit() ) {
        if ( fitResult.status != 3 ) {
            std::cerr << "WARNING in LauAbsFitModel:fitExpt : Not running second stage fit since first stage failed."
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

void LauAbsFitModel::calculateSPlotData()
{
    if ( sPlotNtuple_ != 0 ) {
        sPlotNtuple_->addFriendTree( inputFitData_->fileName(), inputFitData_->treeName() );
        sPlotNtuple_->writeOutGenResults();
        LauSPlot splot( sPlotNtuple_->fileName(),
                        sPlotNtuple_->treeName(),
                        this->firstExpt(),
                        this->nExpt(),
                        this->variableNames(),
                        this->freeSpeciesNames(),
                        this->fixdSpeciesNames(),
                        this->twodimPDFs(),
                        this->splitSignal(),
                        this->scfDPSmear() );
        splot.runCalculations( sPlotVerbosity_ );
        splot.writeOutResults();
    }
}

void LauAbsFitModel::compareFitData( UInt_t toyMCScale,
                                     const TString& mcFileName,
                                     const TString& tableFileName,
                                     Bool_t poissonSmearing )
{
    compareFitData_       = kTRUE;
    fitToyMCScale_        = toyMCScale;
    fitToyMCFileName_     = mcFileName;
    fitToyMCTableName_    = tableFileName;
    fitToyMCPoissonSmear_ = poissonSmearing;
}

void LauAbsFitModel::createFitToyMC( const TString& mcFileName, const TString& tableFileName )
{
    // Create a toy MC sample so that the user can compare the fitted
    // result with the data.
    // Generate more toy MC to reduce statistical fluctuations:
    // - use the rescaling value fitToyMCScale_

    // Store the info on the number of experiments, first expt and current expt
    const UInt_t oldNExpt( this->nExpt() );
    const UInt_t oldFirstExpt( this->firstExpt() );
    const UInt_t oldIExpt( this->iExpt() );
    const Bool_t oldToyExpt( this->toyExpts() );

    // Turn off Poisson smearing if required
    const Bool_t poissonSmearing( this->doPoissonSmearing() );
    this->doPoissonSmearing( fitToyMCPoissonSmear_ );

    // Turn off embedding, since we need toy MC, not reco'ed events
    const Bool_t enableEmbeddingOrig( this->enableEmbedding() );
    this->enableEmbedding( kFALSE );

    // Need to make sure that the generation of the DP co-ordinates is
    // switched on if any of our PDFs depend on it
    const Bool_t origUseDP = this->useDP();
    if ( ! origUseDP && this->pdfsDependOnDP() ) {
        this->useDP( kTRUE );
        this->initialiseDPModels();
    }

    // Construct a unique filename for this experiment
    TString exptString( "_expt" );
    exptString += oldIExpt;
    TString fileName( mcFileName );
    fileName.Insert( fileName.Last( '.' ), exptString );

    // Generate the toy MC
    std::cout << "INFO in LauAbsFitModel::createFitToyMC : Generating toy MC in " << fileName
              << " to compare fit with data..." << std::endl;
    std::cout << "                                       : Number of experiments to generate = "
              << fitToyMCScale_ << "." << std::endl;
    std::cout << "                                       : This is to allow the toy MC to be made with reduced statistical fluctuations."
              << std::endl;

    // Set the genValue of each parameter to its current (fitted) value
    // but first store the original genValues for restoring later
    std::vector<Double_t> origGenValues;
    origGenValues.reserve( this->nTotParams() );
    Bool_t blind( kFALSE );
    for ( LauParameterPList::iterator iter = fitVars_.begin(); iter != fitVars_.end(); ++iter ) {
        origGenValues.push_back( ( *iter )->genValue() );
        ( *iter )->genValue( ( *iter )->unblindValue() );
        if ( ( *iter )->blind() ) {
            blind = kTRUE;
        }
    }
    if ( blind ) {
        std::cerr << "WARNING in LauAbsFitModel::createFitToyMC : One or more parameters are blind but the toy will be created using the unblind values - use with caution!!"
                  << std::endl;
    }

    // If we're asked to generate more than 100 experiments then split it
    // up into multiple files since otherwise can run into memory issues
    // when building the index
    // TODO - this obviously depends on the number of events per experiment as well, so should do this properly

    UInt_t totalExpts = fitToyMCScale_;
    if ( totalExpts > 100 ) {
        UInt_t nFiles = totalExpts / 100;
        if ( totalExpts % 100 ) {
            nFiles += 1;
        }
        TString fileNameBase { fileName };
        for ( UInt_t iFile( 0 ); iFile < nFiles; ++iFile ) {

            UInt_t firstExp( iFile * 100 );

            // Set number of experiments and first experiment to generate
            UInt_t nExp = ( ( firstExp + 100 ) > totalExpts ) ? totalExpts - firstExp : 100;
            this->setNExpts( nExp, firstExp, kTRUE );

            // Create a unique filename and generate the events
            fileName           = fileNameBase;
            TString extraname  = "_file";
            extraname         += iFile;
            fileName.Insert( fileName.Last( '.' ), extraname );
            this->generate( fileName, "genResults", "dummy.root", tableFileName );
        }
    } else {
        // Set number of experiments to new value
        this->setNExpts( fitToyMCScale_, 0, kTRUE );
        // Generate the toy
        this->generate( fileName, "genResults", "dummy.root", tableFileName );
    }

    // Reset number of experiments to original value
    this->setNExpts( oldNExpt, oldFirstExpt, oldToyExpt );
    this->setCurrentExperiment( oldIExpt );

    // Restore the Poisson smearing to its former value
    this->doPoissonSmearing( poissonSmearing );

    // Restore the embedding status
    this->enableEmbedding( enableEmbeddingOrig );

    // Restore "useDP" to its former status
    this->useDP( origUseDP );

    // Restore the original genValue to each parameter
    for ( UInt_t i( 0 ); i < this->nTotParams(); ++i ) {
        fitVars_[i]->genValue( origGenValues[i] );
    }

    std::cout << "INFO in LauAbsFitModel::createFitToyMC : Finished in createFitToyMC." << std::endl;
}

Double_t LauAbsFitModel::getTotNegLogLikelihood()
{
    // Calculate the total negative log-likelihood over all events.
    // This function assumes that the fit parameters and data tree have
    // already been set-up correctly.

    // Loop over the data points to calculate the log likelihood
    Double_t logLike = this->getLogLikelihood( 0, this->eventsPerExpt() );

    // Include the Poisson term in the extended likelihood if required
    if ( this->doEMLFit() ) {
        logLike -= this->getEventSum();
    }

    // Calculate any penalty terms from Gaussian constrained variables
    const auto& multiDimCons = this->multiDimConstraints();
    if ( ! conVars_.empty() || ! multiDimCons.empty() ) {
        logLike -= this->getLogLikelihoodPenalty();
    }

    Double_t totNegLogLike = -logLike;
    return totNegLogLike;
}

Double_t LauAbsFitModel::getLogLikelihoodPenalty()
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

Double_t LauAbsFitModel::getLogLikelihood( UInt_t iStart, UInt_t iEnd )
{
    // Calculate the total negative log-likelihood over all events.
    // This function assumes that the fit parameters and data tree have
    // already been set-up correctly.

    // Loop over the data points to calculate the log likelihood
    Double_t logLike( 0.0 );
    const Double_t worstLL = this->worstLogLike();

    // Loop over the number of events in this experiment
    Bool_t ok( kTRUE );
    for ( UInt_t iEvt = iStart; iEvt < iEnd; ++iEvt ) {

        Double_t likelihood = this->getTotEvtLikelihood( iEvt );

        if ( likelihood > std::numeric_limits<Double_t>::min() ) {    // Is the likelihood zero?
            Double_t evtLogLike = TMath::Log( likelihood );
            if ( doSFit_ ) {
                evtLogLike *= sWeights_[iEvt];
            }
            logLike += evtLogLike;
        } else {
            ok = kFALSE;
            std::cerr << "WARNING in LauAbsFitModel::getLogLikelihood : Strange likelihood value for event "
                      << iEvt << ": " << likelihood << "\n";
            this->printEventInfo( iEvt );
            this->printVarsInfo();    //Write the values of the floated variables for which the likelihood is zero
            break;
        }
    }

    if ( ! ok ) {
        std::cerr << "                                                  : Returning worst NLL found so far to force MINUIT out of this region."
                  << std::endl;
        logLike = worstLL;
    } else if ( logLike < worstLL ) {
        this->worstLogLike( logLike );
    }

    return logLike;
}

void LauAbsFitModel::setParsFromMinuit( Double_t* par, Int_t npar )
{
    // This function sets the internal parameters based on the values
    // that Minuit is using when trying to minimise the total likelihood function.

    // MINOS reports different numbers of free parameters depending on the
    // situation, so disable this check
    if ( ! this->withinAsymErrorCalc() ) {
        const UInt_t nFreePars = this->nFreeParams();
        if ( static_cast<UInt_t>( npar ) != nFreePars ) {
            std::cerr << "ERROR in LauAbsFitModel::setParsFromMinuit : Unexpected number of free parameters: "
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
    // Also check if we have any parameters on which the DP integrals depend
    // and whether they have changed since the last iteration
    Bool_t recalcNorm( kFALSE );
    const LauParameterPSet::const_iterator resVarsEnd = resVars_.end();
    for ( UInt_t i( 0 ); i < this->nTotParams(); ++i ) {
        if ( ! fitVars_[i]->fixed() ) {
            if ( resVars_.find( fitVars_[i] ) != resVarsEnd ) {
                if ( fitVars_[i]->value() != par[i] ) {
                    recalcNorm = kTRUE;
                }
            }
            fitVars_[i]->value( par[i] );
        }
    }

    // If so, then recalculate the normalisation
    if ( recalcNorm ) {
        this->recalculateNormalisation();
    }

    this->propagateParUpdates();
}

UInt_t LauAbsFitModel::addFitParameters( LauParameter* param, const Bool_t addFixed )
{
    UInt_t nParsAdded { 0 };

    // if we're dealing with a clone,
    // we should instead deal with its parent
    if ( param->clone() ) {
        param = param->parent();
    }

    // we need to include the parameter if it is either:
    // - floating
    // - currently fixed but will float in the second stage of a two-stage fit
    if ( addFixed || ! param->fixed() || ( this->twoStageFit() && param->secondStage() ) ) {
        // check whether we already have this parameter stored
        auto [_, addedOK] { fitVarsSet_.insert( param ) };
        if ( addedOK ) {
            // if not, add it to the list
            fitVars_.push_back( param );
            ++nParsAdded;
        }
    }

    return nParsAdded;
}

UInt_t LauAbsFitModel::addFitParameters( LauParameterPList& paramList, const Bool_t addFixed )
{
    UInt_t nParsAdded { 0 };

    for ( auto param : paramList ) {
        nParsAdded += this->addFitParameters( param, addFixed );
    }

    return nParsAdded;
}

UInt_t LauAbsFitModel::addFitParameters( LauAbsRValue* param, const Bool_t addFixed )
{
    UInt_t nParsAdded { 0 };

    LauParameterPList pars { param->getPars() };
    for ( auto par : pars ) {
        nParsAdded += this->addFitParameters( par, addFixed );
    }

    return nParsAdded;
}

UInt_t LauAbsFitModel::addFitParameters( LauAbsRValuePList& paramList, const Bool_t addFixed )
{
    UInt_t nParsAdded { 0 };

    for ( auto param : paramList ) {
        nParsAdded += this->addFitParameters( param, addFixed );
    }

    return nParsAdded;
}

UInt_t LauAbsFitModel::addFitParameters( LauPdfPList& pdfList, const Bool_t addFixed )
{
    UInt_t nParsAdded { 0 };

    for ( auto pdf : pdfList ) {
        if ( pdf->isDPDependent() ) {
            this->pdfsDependOnDP( kTRUE );
        }

        LauAbsRValuePList& pars  = pdf->getParameters();
        nParsAdded              += this->addFitParameters( pars, addFixed );
    }

    return nParsAdded;
}

UInt_t LauAbsFitModel::addResonanceParameters( LauParameter* param )
{
    UInt_t nParsAdded { 0 };

    // first check if we have this parameter in the set of resonance parameters
    auto [_, addedOK] { resVars_.insert( param ) };

    if ( addedOK ) {
        nParsAdded += this->addFitParameters( param );
    }

    return nParsAdded;
}

UInt_t LauAbsFitModel::addResonanceParameters( LauParameterPList& paramList )
{
    UInt_t nParsAdded { 0 };

    for ( auto param : paramList ) {
        nParsAdded += this->addResonanceParameters( param );
    }

    return nParsAdded;
}

void LauAbsFitModel::addConParameters()
{
    // Add penalties from the constraints to fit parameters

    // First, constraints on the fit parameters themselves
    for ( LauParameter* param : fitVars_ ) {
        if ( param->gaussConstraint() ) {
            conVars_.push_back( param );
            std::cout << "INFO in LauAbsFitModel::addConParameters : Added Gaussian constraint to parameter "
                      << param->name() << std::endl;
        }
    }

    // Second, constraints on arbitrary combinations
    auto& conStore = this->formulaConstraints();
    for ( auto& constraint : conStore ) {
        std::vector<LauParameter*> params;
        for ( const auto& name : constraint.conPars_ ) {
            for ( LauParameter* par : fitVars_ ) {
                if ( par->name() == name ) {
                    params.push_back( par );
                }
            }
        }

        // If the parameters are not found, skip it
        if ( params.size() != constraint.conPars_.size() ) {
            std::cerr << "WARNING in LauAbsFitModel::addConParameters: Could not find parameters to constrain in the formula... skipping"
                      << std::endl;
            continue;
        }

        constraint.formulaPar_ =
            std::make_unique<LauFormulaPar>( constraint.formula_, constraint.formula_, params );
        constraint.formulaPar_->addGaussianConstraint( constraint.mean_, constraint.width_ );

        conVars_.push_back( constraint.formulaPar_.get() );

        std::cout << "INFO in LauAbsFitModel::addConParameters : Added Gaussian constraint to formula\n";
        std::cout << "                                         : Formula: " << constraint.formula_
                  << std::endl;
        for ( LauParameter* param : params ) {
            std::cout << "                                         : Parameter: " << param->name()
                      << std::endl;
        }
    }

    // Add n-dimensional constraints
    auto& multiDimCons = this->multiDimConstraints();
    for ( auto& constraint : multiDimCons ) {
        for ( auto& parname : constraint.conPars_ ) {
            for ( auto& fitPar : fitVars_ ) {
                if ( parname == fitPar->name() ) {
                    // Check parameters do not have a 1D Gaussian constraint applied
                    if ( fitPar->gaussConstraint() ) {
                        std::cerr << "ERROR in LauAbsFitModel::addConParameters: parameter in n-dimensional constraint already has a 1d constraint applied"
                                  << std::endl;
                        gSystem->Exit( EXIT_FAILURE );
                    }
                    constraint.conLauPars_.push_back( fitPar );
                }
            }
        }
        if ( constraint.conLauPars_.size() != constraint.conPars_.size() ) {
            std::cerr << "Error in LauAbsFitModel::addConParameters : Could not match parameter names for n-dimensional constraint"
                      << std::endl;
            gSystem->Exit( EXIT_FAILURE );
        }
    }
}

void LauAbsFitModel::updateFitParameters( LauPdfPList& pdfList )
{
    for ( LauPdfPList::iterator pdf_iter = pdfList.begin(); pdf_iter != pdfList.end(); ++pdf_iter ) {
        ( *pdf_iter )->updatePulls();
    }
}

void LauAbsFitModel::printFitParameters( const LauPdfPList& pdfList, std::ostream& fout ) const
{
    LauPrint print;
    for ( LauPdfPList::const_iterator pdf_iter = pdfList.begin(); pdf_iter != pdfList.end();
          ++pdf_iter ) {
        const LauAbsRValuePList& pars = ( *pdf_iter )->getParameters();
        for ( LauAbsRValuePList::const_iterator pars_iter = pars.begin(); pars_iter != pars.end();
              ++pars_iter ) {
            LauParameterPList params = ( *pars_iter )->getPars();
            for ( LauParameterPList::iterator params_iter = params.begin();
                  params_iter != params.end();
                  ++params_iter ) {
                if ( ! ( *params_iter )->clone() ) {
                    fout << ( *params_iter )->name() << "  &  $";
                    print.printFormat( fout, ( *params_iter )->value() );
                    if ( ( *params_iter )->fixed() == kTRUE ) {
                        fout << "$ (fixed) \\\\";
                    } else {
                        fout << " \\pm ";
                        print.printFormat( fout, ( *params_iter )->error() );
                        fout << "$ \\\\" << std::endl;
                    }
                }
            }
        }
    }
}

void LauAbsFitModel::cacheInfo( LauPdfPList& pdfList, const LauFitDataTree& theData )
{
    for ( LauPdfPList::iterator pdf_iter = pdfList.begin(); pdf_iter != pdfList.end(); ++pdf_iter ) {
        ( *pdf_iter )->cacheInfo( theData );
    }
}

Double_t LauAbsFitModel::prodPdfValue( LauPdfPList& pdfList, UInt_t iEvt )
{
    Double_t pdfVal = 1.0;
    for ( LauPdfPList::iterator pdf_iter = pdfList.begin(); pdf_iter != pdfList.end(); ++pdf_iter ) {
        ( *pdf_iter )->calcLikelihoodInfo( iEvt );
        pdfVal *= ( *pdf_iter )->getLikelihood();
    }
    return pdfVal;
}

void LauAbsFitModel::printEventInfo( UInt_t iEvt ) const
{
    const LauFitData& data = inputFitData_->getData( iEvt );
    std::cerr << "                                                  : Input data values for this event:"
              << std::endl;
    for ( LauFitData::const_iterator iter = data.begin(); iter != data.end(); ++iter ) {
        std::cerr << " " << iter->first << " = " << iter->second << std::endl;
    }
}

void LauAbsFitModel::printVarsInfo() const
{
    std::cerr << "                                                  : Current values of fit parameters:"
              << std::endl;
    for ( UInt_t i( 0 ); i < this->nTotParams(); ++i ) {
        std::cerr << " " << ( fitVars_[i]->name() ).Data() << " = " << fitVars_[i]->value()
                  << std::endl;
    }
}

void LauAbsFitModel::prepareInitialParArray( TObjArray& array )
{
    // Update initial fit parameters if required (e.g. if using random numbers).
    this->checkInitFitParams();

    // Store the total number of parameters and the number of free parameters
    UInt_t nPars     = fitVars_.size();
    UInt_t nFreePars = 0;

    // Send the fit parameters
    for ( LauParameterPList::iterator iter = fitVars_.begin(); iter != fitVars_.end(); ++iter ) {
        if ( ! ( *iter )->fixed() ) {
            ++nFreePars;
        }
        array.Add( *iter );
    }

    this->startNewFit( nPars, nFreePars );
}

void LauAbsFitModel::finaliseExperiment( const LauAbsFitter::FitStatus& fitStat,
                                         const TObjArray* parsFromCoordinator,
                                         const TMatrixD* covMat,
                                         TObjArray& parsToCoordinator )
{
    // Copy the fit status information
    this->storeFitStatus( fitStat, *covMat );

    // Now process the parameters
    const UInt_t nPars          = this->nTotParams();
    UInt_t nParsFromCoordinator = parsFromCoordinator->GetEntries();
    if ( nParsFromCoordinator != nPars ) {
        std::cerr << "ERROR in LauAbsFitModel::finaliseExperiment : Unexpected number of parameters received from coordinator"
                  << std::endl;
        std::cerr << "                                            : Received "
                  << nParsFromCoordinator << " when expecting " << nPars << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    for ( UInt_t iPar( 0 ); iPar < nParsFromCoordinator; ++iPar ) {
        LauParameter* parameter = dynamic_cast<LauParameter*>( ( *parsFromCoordinator )[iPar] );
        if ( ! parameter ) {
            std::cerr << "ERROR in LauAbsFitModel::finaliseExperiment : Error reading parameter from coordinator"
                      << std::endl;
            gSystem->Exit( EXIT_FAILURE );
        }

        if ( parameter->name() != fitVars_[iPar]->name() ) {
            std::cerr << "ERROR in LauAbsFitModel::finaliseExperiment : Error reading parameter from coordinator"
                      << std::endl;
            gSystem->Exit( EXIT_FAILURE );
        }

        *( fitVars_[iPar] ) = *parameter;
    }

    // Write the results into the ntuple
    this->finaliseFitResults( outputTableName_ );

    // Store the per-event likelihood values
    if ( this->writeSPlotData() ) {
        this->storePerEvtLlhds();
    }

    // Create a toy MC sample using the fitted parameters so that
    // the user can compare the fit to the data.
    if ( compareFitData_ == kTRUE && fitStat.status == 3 ) {
        this->createFitToyMC( fitToyMCFileName_, fitToyMCTableName_ );
    }

    // Send the finalised fit parameters
    for ( LauParameterPList::iterator iter = fitVars_.begin(); iter != fitVars_.end(); ++iter ) {
        parsToCoordinator.Add( *iter );
    }
}

UInt_t LauAbsFitModel::readExperimentData()
{
    // retrieve the data and find out how many events have been read
    const UInt_t exptIndex = this->iExpt();
    inputFitData_->readExperimentData( exptIndex );
    const UInt_t nEvent = inputFitData_->nEvents();
    this->eventsPerExpt( nEvent );
    return nEvent;
}

void LauAbsFitModel::setParametersFromFile( const TString& fileName,
                                            const TString& treeName,
                                            const Bool_t fix )
{
    fixParamFileName_ = fileName;
    fixParamTreeName_ = treeName;
    fixParams_        = fix;
}

void LauAbsFitModel::setParametersFromMap( const std::map<TString, Double_t>& parameters,
                                           const Bool_t fix )
{
    fixParamMap_ = parameters;
    fixParams_   = fix;
}

void LauAbsFitModel::setNamedParameters( const TString& fileName,
                                         const TString& treeName,
                                         const std::set<TString>& parameters,
                                         const Bool_t fix )
{
    fixParamFileName_ = fileName;
    fixParamTreeName_ = treeName;
    fixParamNames_    = parameters;
    fixParams_        = fix;
}

void LauAbsFitModel::setParametersFileFallback( const TString& fileName,
                                                const TString& treeName,
                                                const std::map<TString, Double_t>& parameters,
                                                const Bool_t fix )
{
    fixParamFileName_ = fileName;
    fixParamTreeName_ = treeName;
    fixParamMap_      = parameters;
    fixParams_        = fix;
}
