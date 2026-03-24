
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

/*! \file LauFitDataTree.cc
    \brief File containing implementation of LauFitDataTree class.
*/

#include "LauFitDataTree.hh"

#include "TFile.h"
#include "TLeaf.h"
#include "TString.h"
#include "TSystem.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

LauFitDataTree::LauFitDataTree( const TString& rootFileName, const TString& rootTreeName ) :
    rootFileName_( rootFileName ),
    rootTreeName_( rootTreeName ),
    rootFile_( 0 ),
    rootTree_( 0 ),
    eventList_( 0 )
{
    if ( rootFileName_ != "" && rootTreeName_ != "" ) {
        this->openFileAndTree();
    }
}

LauFitDataTree::~LauFitDataTree()
{
    if ( rootFile_ && rootFile_->IsOpen() ) {
        delete rootTree_;
        rootTree_ = 0;
        delete eventList_;
        eventList_ = 0;
    }
    delete rootFile_;
    rootFile_ = 0;
}

void LauFitDataTree::openFileAndTree()
{
    // first check whether we've already opened up the file or not
    if ( ! rootFile_ ) {
        // if not, first check the filename and if all ok create the file
        if ( rootFileName_ == "" ) {
            std::cerr << "ERROR in LauFitDataTree::openFileAndTree : Bad filename supplied, not creating file or tree."
                      << std::endl;
            return;
        }
        rootFile_ = TFile::Open( rootFileName_, "read" );
        if ( ! rootFile_ || rootFile_->IsZombie() || ! rootFile_->IsOpen() ) {
            std::cerr << "ERROR in LauFitDataTree::openFileAndTree : Problem opening file \""
                      << rootFileName_ << "\" for reading." << std::endl;
            return;
        }
    }
    // check whether we've already created the tree
    if ( ! rootTree_ ) {
        // if not change to the file's directory and grab the tree
        rootFile_->cd();
        rootTree_ = dynamic_cast<TTree*>( rootFile_->Get( rootTreeName_ ) );
        if ( ! rootTree_ ) {
            std::cerr << "ERROR in LauFitDataTree::openFileAndTree : Problem reading tree \""
                      << rootTreeName_ << "\" from file." << std::endl;
            return;
        }
        rootTree_->SetDirectory( rootFile_ );
    }
}

Bool_t LauFitDataTree::findBranches()
{
    if ( ! rootTree_ ) {
        std::cerr << "ERROR in LauFitDataTree::findBranches : Invalid pointer to data tree."
                  << std::endl;
        return kFALSE;
    }

    // this method should only be called once
    if ( ! leafNames_.empty() ) {
        std::cerr << "ERROR in LauFitDataTree::findBranches : Branches already found, not running again."
                  << std::endl;
        return kFALSE;
    }

    leafNames_.clear();
    eventData_.clear();
    eventDataOut_.clear();
    leaves_.clear();
    treeEvents_.clear();
    fakeEvents_.clear();

    const UInt_t numBranches( this->nBranches() );
    eventData_.reserve( numBranches );

    TObjArray* pLeaves = rootTree_->GetListOfLeaves();
    if ( ! pLeaves ) {
        std::cerr << "ERROR in LauFitDataTree::findBranches : Problem retrieving leaves from the tree."
                  << std::endl;
        return kFALSE;
    }
    TObjArray& leaves = *pLeaves;

    const UInt_t numLeaves = static_cast<UInt_t>( leaves.GetEntries() );
    if ( numBranches != numLeaves ) {
        std::cerr << "ERROR in LauFitDataTree::findBranches : List of leaves is not the same size as the number of branches - this is strange!"
                  << std::endl;
        std::cerr << "                                      : There are " << numLeaves
                  << " leaves and " << numBranches << " branches" << std::endl;
        return kFALSE;
    }

    for ( UInt_t iLeaf = 0; iLeaf < numBranches; ++iLeaf ) {

        TLeaf* leaf = dynamic_cast<TLeaf*>( leaves[iLeaf] );

        // we can't deal with arrays
        Int_t size = leaf->GetNdata();
        if ( size != 1 ) {
            std::cerr << "ERROR in LauFitDataTree::findBranches : Tree has array branches, can't deal with those."
                      << std::endl;
            return kFALSE;
        }

        // find the name and type of the leaf
        TString name = leaf->GetName();

        // initialise an entry in the maps to hold the value
        leafNames_[name] = iLeaf;
        leaves_.push_back( leaf );
        eventData_.push_back( 0.0 );
        eventDataOut_[name] = 0.0;
    }

    const UInt_t numEvents( this->nTotalEvents() );

    std::cout << "INFO in LauFitDataTree::findBranches : Finished finding data tree branches."
              << std::endl;
    std::cout << "                                     : Tree contains " << numBranches
              << " branches and a total of " << numEvents << " events." << std::endl;

    return kTRUE;
}

void LauFitDataTree::readExperimentData( UInt_t iExpt )
{
    // Check that we have a valid tree to read from
    if ( ! rootTree_ ) {
        std::cerr << "ERROR in LauFitDataTree::readExperimentData : Invalid pointer to data tree."
                  << std::endl;
        return;
    }

    rootFile_->cd();

    // If the event list doesn't exist, create it
    if ( ! eventList_ ) {
        eventList_ = new TEventList( "elist", "elist" );
        eventList_->SetDirectory( rootFile_ );
    }

    // Reset the event list so it contains no entries
    eventList_->Reset();

    // Check that the tree contains the iExpt branch
    if ( ! this->haveBranch( "iExpt" ) ) {
        // If it doesn't, check which experiment we've been asked to read
        if ( iExpt == 0 ) {
            // If it's experiment 0, assume that we should read all events in the tree (but print a warning since this is an assumption)
            std::cerr << "WARNING in LauFitDataTree::readExperimentData : Tree does not contain \"iExpt\" branch, will read all data in the tree"
                      << std::endl;
            this->readAllData();
            return;
        } else {
            // Otherwise we can't really do anything, so print an error and return
            std::cerr << "ERROR in LauFitDataTree::readExperimentData : Tree does not contain \"iExpt\" branch and experiment requested is > 0, will not read anything"
                      << std::endl;
            return;
        }
    }

    // Fill the event list with the entries for the requested experiment
    TString listName( eventList_->GetName() );
    listName.Prepend( ">>" );
    TString selection( "iExpt==" );
    selection += iExpt;

    std::cout << "INFO in LauFitDataTree::readExperimentData : Setting tree to experiment number "
              << iExpt << "." << std::endl;
    rootTree_->Draw( listName, selection );
    std::cout << "                                           : Found " << this->nEvents()
              << " events." << std::endl;

    // Load the data
    this->loadData();
}

UInt_t LauFitDataTree::nBranches() const
{
    if ( rootTree_ ) {
        return static_cast<UInt_t>( rootTree_->GetNbranches() );
    } else if ( ! fakeEvents_.empty() ) {
        std::vector<LauEventData>::const_iterator fakeIter = fakeEvents_.begin();
        return ( *fakeIter ).size();
    } else {
        return 0;
    }
}

void LauFitDataTree::appendFakePoints( const std::vector<Double_t>& xCoords,
                                       const std::vector<Double_t>& yCoords )
{
    if ( xCoords.size() != yCoords.size() ) {
        std::cerr << "ERROR in LauFitDataTree::appendFakePoints : The co-ordinate vectors are not the same size!"
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // clear the list of fake events and reserve enough space
    fakeEvents_.clear();
    fakeEvents_.reserve( xCoords.size() );

    // loop over the event data zeroing everything
    for ( LauEventData::iterator iter = eventData_.begin(); iter != eventData_.end(); ++iter ) {
        *iter = 0.0;
    }

    // find the vector indices for the DP co-ordinates
    LauNameIndexMap::const_iterator iter = leafNames_.find( "m13Sq" );
    if ( iter == leafNames_.end() ) {
        std::cerr << "ERROR in LauFitDataTree::appendFakePoints : Can't find entry \"m13Sq\" in event data map."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }
    UInt_t m13SqIdx = iter->second;

    iter = leafNames_.find( "m23Sq" );
    if ( iter == leafNames_.end() ) {
        std::cerr << "ERROR in LauFitDataTree::appendFakePoints : Can't find entry \"m23Sq\" in event data map."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }
    UInt_t m23SqIdx = iter->second;

    // TODO - is this really needed?  Comment out for now.
    /*
    // and a fake tagCat which we don't care about, but the program
    // needs a number here for comptibility with conventional events.
    iter = leafNames_.find( "tagCat" );
    if ( iter == leafNames_.end() ) {
        std::cerr << "ERROR in LauFitDataTree::appendFakePoints : Can't find entry \"tagCat\" in event data map." << std::endl;
        gSystem->Exit(EXIT_FAILURE);
    }
    UInt_t tagCatIdx = iter->second;
     */

    for ( UInt_t iEvt = 0; iEvt < xCoords.size(); ++iEvt ) {

        // set the values
        eventData_[m13SqIdx] = xCoords[iEvt];
        eventData_[m23SqIdx] = yCoords[iEvt];
        //eventData_[ tagCatIdx ] = 0.0;

        // append to the list of fake events
        fakeEvents_.push_back( eventData_ );
    }
}

Bool_t LauFitDataTree::haveBranch( const TString& name ) const
{
    LauNameIndexMap::const_iterator iter = leafNames_.find( name );
    return ( iter != leafNames_.end() );
}

void LauFitDataTree::disableAllBranches() const
{
    rootTree_->SetBranchStatus( "*", 0 );
}

void LauFitDataTree::enableAllBranches() const
{
    rootTree_->SetBranchStatus( "*", 1 );
}

void LauFitDataTree::enableBranch( const TString& name ) const
{
    if ( this->haveBranch( name ) ) {
        rootTree_->SetBranchStatus( name, 1 );
    }
}

void LauFitDataTree::disableBranch( const TString& name ) const
{
    if ( this->haveBranch( name ) ) {
        rootTree_->SetBranchStatus( name, 0 );
    }
}

void LauFitDataTree::loadData()
{
    // Clear any data from previous expts
    treeEvents_.clear();

    // Reserve enough space in the vector
    Long64_t nEvts = this->nEvents();
    treeEvents_.reserve( nEvts );

    const UInt_t nLeaves = this->nBranches();

    // Loop through the tree and store the data
    for ( Long64_t iEvt( 0 ); iEvt < nEvts; ++iEvt ) {

        // Find which entry from the full tree contains the requested event
        Long64_t iEntry = eventList_ ? eventList_->GetEntry( iEvt ) : iEvt;
        if ( iEntry < 0 ) {    // this shouldn't happen, but just in case...
            std::cerr << "ERROR in LauFitDataTree::loadData : Requested event not found."
                      << std::endl;
            gSystem->Exit( EXIT_FAILURE );
        }

        // Then retrieve that entry from the tree
        rootTree_->GetEntry( iEntry );

        // Copy the leaf values into eventData_
        for ( UInt_t iLeaf( 0 ); iLeaf < nLeaves; ++iLeaf ) {

            const TLeaf* leaf = leaves_[iLeaf];

            Double_t value = leaf->GetValue();
            if ( std::isnan( value ) || std::isinf( value ) ) {
                std::cerr << "ERROR in LauFitDataTree::loadData : Event " << iEvt
                          << " has infinite or NaN entry for variable " << leaf->GetName()
                          << std::endl;
                gSystem->Exit( EXIT_FAILURE );
            }

            eventData_[iLeaf] = value;
        }

        // Store the event in the vector
        treeEvents_.push_back( eventData_ );
    }
}

void LauFitDataTree::readAllData()
{
    delete eventList_;
    eventList_ = 0;
    this->loadData();
}

const LauFitData& LauFitDataTree::getData( UInt_t iEvt ) const
{
    // For a given event, within the current experiment,
    // retrieve the input fit variables

    const UInt_t numTreeEvents = this->nEvents();
    const UInt_t numFakeEvents = this->nFakeEvents();

    // Does the requested event come from the tree or from the fake events list?
    if ( iEvt < numTreeEvents ) {
        if ( iEvt > treeEvents_.size() ) {    // this shouldn't happen, but just in case...
            std::cerr << "ERROR in LauFitDataTree::getData : Requested event, " << iEvt
                      << ", not found." << std::endl;
            gSystem->Exit( EXIT_FAILURE );
        }
        // Retrieve that event from the vector
        eventData_ = treeEvents_[iEvt];
    } else if ( iEvt < ( numTreeEvents + numFakeEvents ) ) {
        // Find the entry from the fake events list contains the requested event
        UInt_t iEntry = iEvt - numTreeEvents;
        // Copy the entry into eventData_
        eventData_ = fakeEvents_[iEntry];
    } else {
        std::cerr << "ERROR in LauFitDataTree::getData : Requested event " << iEvt
                  << " not found for " << rootTreeName_ << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    for ( LauNameIndexMap::const_iterator iter = leafNames_.begin(); iter != leafNames_.end();
          ++iter ) {
        const TString& name = iter->first;
        UInt_t index        = iter->second;

        LauFitData::iterator outIter = eventDataOut_.find( name );
        if ( outIter == eventDataOut_.end() ) {
            std::cerr << "ERROR in LauFitDataTree::getData : Something very strange happening."
                      << std::endl;
            gSystem->Exit( EXIT_FAILURE );
        }

        outIter->second = eventData_[index];
    }
    return eventDataOut_;
}
