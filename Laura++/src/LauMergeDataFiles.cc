
/*
Copyright 2008 University of Warwick

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

/*! \file LauMergeDataFiles.cc
    \brief File containing implementation of LauMergeDataFiles class.
*/

#include "LauMergeDataFiles.hh"

#include "TLeaf.h"
#include "TObjArray.h"
#include "TSystem.h"

#include <iostream>
#include <map>

LauMergeDataFiles::LauMergeDataFiles( const TString& fileName1,
                                      const TString& fileName2,
                                      const TString& treeName ) :
    fileName1_( fileName1 ),
    fileName2_( fileName2 ),
    treeName_( treeName ),
    inputFile1_( 0 ),
    inputFile2_( 0 ),
    inputTree1_( 0 ),
    inputTree2_( 0 ),
    outputFile_( 0 ),
    outputTree_( 0 )
{
}

LauMergeDataFiles::~LauMergeDataFiles()
{
    if ( inputFile1_ && inputFile1_->IsOpen() ) {
        inputFile1_->Close();
    }
    delete inputFile1_;

    if ( inputFile2_ && inputFile2_->IsOpen() ) {
        inputFile2_->Close();
    }
    delete inputFile2_;

    if ( outputFile_ && outputFile_->IsOpen() ) {
        outputFile_->Close();
    }
    delete outputFile_;
}

void LauMergeDataFiles::openInputFiles()
{
    // open the two input ROOT files
    inputFile1_ = TFile::Open( fileName1_ );
    if ( ! inputFile1_ || inputFile1_->IsZombie() ) {
        std::cerr << "Problem opening file: \"" << fileName1_ << "\", exiting..." << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }
    inputTree1_ = dynamic_cast<TTree*>( inputFile1_->Get( treeName_ ) );
    if ( ! inputTree1_ ) {
        std::cerr << "Problem getting tree called " << treeName_ << "from file: \"" << fileName1_
                  << "\", exiting..." << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    inputFile2_ = TFile::Open( fileName2_ );
    if ( ! inputFile2_ || inputFile2_->IsZombie() ) {
        std::cerr << "Problem opening file: \"" << fileName2_ << "\", exiting..." << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }
    inputTree2_ = dynamic_cast<TTree*>( inputFile2_->Get( treeName_ ) );
    if ( ! inputTree2_ ) {
        std::cerr << "Problem getting tree called " << treeName_ << "from file: \"" << fileName2_
                  << "\", exiting..." << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }
}

void LauMergeDataFiles::setupInputTrees()
{
    TObjArray* leaves1 = inputTree1_->GetListOfLeaves();
    TObjArray* leaves2 = inputTree2_->GetListOfLeaves();
    Int_t nLeaves1     = leaves1->GetEntries();
    Int_t nLeaves2     = leaves2->GetEntries();
    if ( nLeaves1 != nLeaves2 ) {
        std::cerr << "Different number of leaves in the two input trees, not continuing."
                  << std::endl;
        return;
    }

    std::cout << "Setting branches for input trees \"" << treeName_ << "\"..." << std::endl;
    inputTree1_->SetBranchAddress( "iExpt", &iExpt_ );
    inputTree1_->SetBranchAddress( "iEvtWithinExpt", &iEvtWithinExpt_ );
    inputTree2_->SetBranchAddress( "iExpt", &iExpt_ );
    inputTree2_->SetBranchAddress( "iEvtWithinExpt", &iEvtWithinExpt_ );

    for ( Int_t iLeaf( 0 ); iLeaf < nLeaves1; ++iLeaf ) {

        TLeaf* leaf  = dynamic_cast<TLeaf*>( ( *leaves1 )[iLeaf] );
        TString type = leaf->GetTypeName();
        TString name = leaf->GetName();
        Int_t size   = leaf->GetNdata();

        if ( ( name == "iExpt" ) || ( name == "iEvtWithinExpt" ) || ( size != 1 ) ) {
            continue;
        }

        if ( type == "Double_t" ) {
            std::pair<LeafDoubleMap::iterator, bool> result = doubleVars_.insert(
                std::make_pair( name, 0.0 ) );
            LeafDoubleMap::iterator iter = result.first;
            bool ok                      = result.second;
            if ( ok ) {
                inputTree1_->SetBranchAddress( name, &( iter->second ) );
                inputTree2_->SetBranchAddress( name, &( iter->second ) );
            }
        } else if ( type == "Int_t" ) {
            std::pair<LeafIntegerMap::iterator, bool> result = integerVars_.insert(
                std::make_pair( name, 0 ) );
            LeafIntegerMap::iterator iter = result.first;
            bool ok                       = result.second;
            if ( ok ) {
                inputTree1_->SetBranchAddress( name, &( iter->second ) );
                inputTree2_->SetBranchAddress( name, &( iter->second ) );
            }
        }
    }
    std::cout << "Set branch addresses for " << doubleVars_.size() << " Double_t branches.\n";
    std::cout << "Set branch addresses for " << integerVars_.size() << " Int_t branches.\n"
              << std::endl;
}

void LauMergeDataFiles::setupOutputTree()
{
    std::cout << "Creating branches for output tree \"" << outputTree_->GetName() << "\"..."
              << std::endl;

    outputTree_->Branch( "iExpt", &iExpt_, "iExpt/I" );
    outputTree_->Branch( "iEvtWithinExpt", &iEvtWithinExpt_, "iEvtWithinExpt/I" );

    for ( LeafDoubleMap::iterator iter = doubleVars_.begin(); iter != doubleVars_.end(); ++iter ) {
        TString name       = iter->first;
        Double_t* address  = &( iter->second );
        TString thirdBit   = name;
        thirdBit          += "/D";

        outputTree_->Branch( name, address, thirdBit );
    }
    for ( LeafIntegerMap::iterator iter = integerVars_.begin(); iter != integerVars_.end(); ++iter ) {
        TString name      = iter->first;
        Int_t* address    = &( iter->second );
        TString thirdBit  = name;
        thirdBit         += "/I";

        outputTree_->Branch( name, address, thirdBit );
    }
    std::cout << "Created " << doubleVars_.size() + integerVars_.size() << " branches.\n"
              << std::endl;
}

void LauMergeDataFiles::process( const TString& fileName )
{
    this->openInputFiles();

    this->setupInputTrees();

    outputFile_ = TFile::Open( fileName, "recreate" );
    outputTree_ = new TTree( treeName_, "" );
    this->setupOutputTree();

    // loop over the trees and combine the corresponding experiments
    std::cout << "Starting to combine the trees..." << std::endl;

    // Find the first and last entries for each experiment in each tree
    this->findExperiments( inputTree1_, tree1Expts_ );
    this->findExperiments( inputTree2_, tree2Expts_ );

    // Check that the experiments in the two trees match
    if ( ! this->checkExperimentMaps() ) {
        return;
    }

    // Loop through the experiments
    for ( ExptsMap::const_iterator iter1 = tree1Expts_.begin(); iter1 != tree1Expts_.end(); ++iter1 ) {

        // get the map element for tree2
        Int_t expt                     = iter1->first;
        ExptsMap::const_iterator iter2 = tree2Expts_.find( expt );

        // determine the number of entries in tree1
        Int_t nEntriesInTree1 = iter1->second.second - iter1->second.first + 1;

        // read the entries from the trees, filling the output tree
        this->readExperiment( inputTree1_, iter1, 0 );
        this->readExperiment( inputTree2_, iter2, nEntriesInTree1 );
    }

    // Write the output file
    this->writeFile();
}

void LauMergeDataFiles::findExperiments( TTree* tree, ExptsMap& exptsMap )
{
    const Int_t nEntries = tree->GetEntries();

    // loop through the tree
    for ( Int_t iEntry( 0 ); iEntry < nEntries; ++iEntry ) {
        // grab the entry
        tree->GetEntry( iEntry );

        // see if we already have an element in the map for the
        // current experiment
        ExptsMap::iterator iter = exptsMap.find( iExpt_ );
        if ( iter == exptsMap.end() ) {
            // if not, we need to add an element that points to
            // this entry in the tree as the start entry
            exptsMap.insert( std::make_pair( iExpt_, std::make_pair( iEntry, -99 ) ) );

            // also we need to complete the map element for the
            // previous experiment with the previous tree entry
            // as the last entry
            ExptsMap::iterator previter = exptsMap.find( iExpt_ - 1 );
            if ( previter != exptsMap.end() ) {
                previter->second.second = iEntry - 1;
            }
        }
    }
    // need to complete the map element for the final experiment
    exptsMap[iExpt_].second = nEntries - 1;
}

Bool_t LauMergeDataFiles::checkExperimentMaps() const
{
    // first check that the two maps are the same size
    UInt_t size1 = tree1Expts_.size();
    UInt_t size2 = tree2Expts_.size();
    if ( size1 != size2 ) {
        std::cerr << "ERROR in LauMergeDataFiles::checkExperimentMaps : Experiment maps are not the same size.\n";
        std::cerr << "                                                : Tree from " << fileName1_
                  << " has " << size1 << " experiments.";
        std::cerr << "                                                : Tree from " << fileName2_
                  << " has " << size2 << " experiments.";
        return kFALSE;
    }

    for ( ExptsMap::const_iterator iter1 = tree1Expts_.begin(); iter1 != tree1Expts_.end(); ++iter1 ) {
        Int_t expt                     = iter1->first;
        ExptsMap::const_iterator iter2 = tree2Expts_.find( expt );
        if ( iter2 == tree2Expts_.end() ) {
            std::cerr << "ERROR in LauMergeDataFiles::checkExperimentMaps : Cannot find experiment "
                      << expt << " in tree from " << fileName2_ << std::endl;
            return kFALSE;
        }
    }

    return kTRUE;
}

void LauMergeDataFiles::readExperiment( TTree* tree, const ExptsMap::const_iterator& expt, Int_t offset )
{
    // find the first and last entry for this experiment from the map element
    const Int_t firstEntry = expt->second.first;
    const Int_t lastEntry  = expt->second.second;

    // loop through all the entries
    for ( Int_t iEntry( firstEntry ); iEntry <= lastEntry; ++iEntry ) {
        // get the entry from the tree
        tree->GetEntry( iEntry );
        // apply the offset to the "event within experiment" variable
        iEvtWithinExpt_ += offset;
        // fill the output tree
        outputTree_->Fill();
    }
}

void LauMergeDataFiles::writeFile()
{
    std::cout << "Building experiment:event index" << std::endl;
    outputTree_->BuildIndex( "iExpt", "iEvtWithinExpt" );

    std::cout << "Writing data to outputfile " << outputFile_->GetName() << std::endl;
    outputTree_->SetDirectory( outputFile_ );
    outputFile_->Write();

    // clean-up
    outputFile_->Close();
    delete outputFile_;
    outputFile_ = 0;
    outputTree_ = 0;

    inputFile1_->Close();
    delete inputFile1_;
    inputFile1_ = 0;
    inputTree1_ = 0;

    inputFile2_->Close();
    delete inputFile2_;
    inputFile2_ = 0;
    inputTree2_ = 0;

    doubleVars_.clear();
    integerVars_.clear();
    tree1Expts_.clear();
    tree2Expts_.clear();
}
