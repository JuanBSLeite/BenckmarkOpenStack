
/*
Copyright 2005 University of Warwick

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

/*! \file LauResultsExtractor.cc
    \brief File containing implementation of LauResultsExtractor class.
*/

#include "LauResultsExtractor.hh"

#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TLeaf.h"
#include "TObjArray.h"
#include "TSystem.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

LauResultsExtractor::LauResultsExtractor( const TString& inputFileName,
                                          const TString& outputFileName,
                                          const TString& treeName ) :
    inputFileName_( inputFileName ),
    outputFileName_( outputFileName ),
    treeName_( treeName ),
    inputTree_( 0 ),
    outputFile_( 0 ),
    outputTree_( 0 ),
    nEntries_( 0 )
{
}

LauResultsExtractor::~LauResultsExtractor()
{
    this->clearMaps();
    delete inputTree_;
    inputTree_ = 0;
    if ( outputFile_ && outputFile_->IsOpen() ) {
        delete outputTree_;
        outputTree_ = 0;
    }
    delete outputFile_;
    outputFile_ = 0;
}

void LauResultsExtractor::setupInputTree()
{
    TObjArray* leaves = inputTree_->GetListOfLeaves();
    Int_t nLeaves     = leaves->GetEntries();

    std::cout << "Setting branches for input tree \"" << inputTree_->GetName() << "\" with "
              << nEntries_ << " entries..." << std::endl;
    inputTree_->SetBranchAddress( "iExpt", &iExpt_ );
    inputTree_->SetBranchAddress( "fitStatus", &fitStatus_ );
    inputTree_->SetBranchAddress( "NLL", &NLL_ );
    inputTree_->SetBranchAddress( "EDM", &EDM_ );

    for ( Int_t iLeaf( 3 ); iLeaf < nLeaves; ++iLeaf ) {

        TLeaf* leaf  = dynamic_cast<TLeaf*>( ( *leaves )[iLeaf] );
        TString type = leaf->GetTypeName();
        TString name = leaf->GetName();
        Int_t size   = leaf->GetNdata();

        if ( ( type != "Double_t" ) || ( size != 1 ) ) {
            continue;
        }

        std::pair<std::map<TString, Double_t>::iterator, bool> result = otherVars_.insert(
            std::make_pair( name, 0.0 ) );
        std::map<TString, Double_t>::iterator iter = result.first;
        bool ok                                    = result.second;
        if ( ok ) {
            inputTree_->SetBranchAddress( name, &( iter->second ) );
        }
    }

    std::cout << "Set branch addresses for " << otherVars_.size() + 3 << " branches.\n"
              << std::endl;
}

void LauResultsExtractor::setupOutputTree( TTree* tree )
{
    std::cout << "Creating branches for output tree \"" << tree->GetName() << "\"..." << std::endl;

    tree->Branch( "iExpt", &iExpt_, "iExpt/I" );
    tree->Branch( "fitStatus", &fitStatus_, "fitStatus/I" );
    tree->Branch( "NLL", &NLL_, "NLL/D" );
    tree->Branch( "EDM", &EDM_, "EDM/D" );

    for ( std::map<TString, Double_t>::iterator iter = otherVars_.begin(); iter != otherVars_.end();
          ++iter ) {
        TString name       = iter->first;
        Double_t* address  = &( iter->second );
        TString thirdBit   = name;
        thirdBit          += "/D";

        tree->Branch( name, address, thirdBit );
    }
    std::cout << "Created " << otherVars_.size() << " branches.\n" << std::endl;
}

void LauResultsExtractor::setInputTreeBranchStatus( const Bool_t status )
{
    inputTree_->SetBranchStatus( "iExpt", kTRUE );
    inputTree_->SetBranchStatus( "fitStatus", kTRUE );
    inputTree_->SetBranchStatus( "NLL", kTRUE );
    inputTree_->SetBranchStatus( "EDM", kTRUE );

    for ( std::map<TString, Double_t>::iterator iter = otherVars_.begin(); iter != otherVars_.end();
          ++iter ) {
        TString name = iter->first;
        inputTree_->SetBranchStatus( name, status );
    }
}

void LauResultsExtractor::clearMaps()
{
    for ( std::map<Int_t, TH1*>::iterator iter = nllHistos_.begin(); iter != nllHistos_.end();
          ++iter ) {
        delete ( iter->second );
    }
    bestNLL_.clear();
    worstNLL_.clear();
    allNLLs_.clear();
    nllHistos_.clear();
}

void LauResultsExtractor::process( Int_t numExpts )
{
    // open the text file
    std::cout << "\n"
              << "Chaining...\n"
              << std::endl;
    std::ifstream textFile( inputFileName_, std::ios::in );
    if ( ! textFile.good() ) {
        std::cerr << "Problem opening file: \"" << inputFileName_ << "\", exiting..." << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    if ( inputTree_ ) {
        delete inputTree_;
        inputTree_ = 0;
    }
    inputTree_ = new TChain( treeName_ );

    // Read the text file and add each valid entry to the chain
    TString inputFileName = "";
    while ( inputFileName.ReadLine( textFile ) && ( ! inputFileName.IsNull() ) ) {
        if ( inputFileName.EndsWith( ".root" ) && ! inputFileName.BeginsWith( "#" ) ) {
            std::cout << inputFileName << std::endl;
            inputTree_->Add( inputFileName );
        } else {
            std::cout << inputFileName << "\t *** Skipped ***" << std::endl;
        }
    }

    textFile.close();
    std::cout << "\n"
              << "... finished.\n"
              << std::endl;

    nEntries_ = inputTree_->GetEntries();
    this->setupInputTree();

    outputTree_ = new TTree( treeName_, "" );
    this->setupOutputTree( outputTree_ );

    // setup the map:
    // for each experiment there is a pair object holding
    // the best NLL and the tree entry for that NLL value
    // each expt starts out with NLL = 0.0 and entry = -1
    std::cout << "Setting up the map..." << std::flush;
    this->clearMaps();
    for ( Int_t i( 0 ); i < numExpts; ++i ) {
        bestNLL_.insert( std::make_pair( i, std::make_pair( 0.0, -1 ) ) );
        worstNLL_.insert( std::make_pair( i, std::make_pair( 0.0, -1 ) ) );
        allNLLs_.insert( std::make_pair( i, std::vector<Double_t>() ) );
        allNLLs_[i].reserve( nEntries_ );
    }
    std::cout << " done.\n" << std::endl;

    // only read the 3 info branches
    //this->setInputTreeBranchStatus(kFALSE);

    // loop over the tree and store the best entries for each expt
    std::cout << "Starting to store best entry info..." << std::endl;
    for ( Int_t j( 0 ); j < nEntries_; ++j ) {

        if ( ( nEntries_ < 100 ) || ( j % ( nEntries_ / 100 ) == 0 ) ) {
            std::cout << "Examining entry " << j << std::endl;
        }

        inputTree_->GetEntry( j );

        if ( ( fitStatus_ == 3 ) && ( NLL_ > -DBL_MAX / 10.0 ) ) {
            allNLLs_[iExpt_].push_back( NLL_ );

            Double_t curBestNLL = bestNLL_[iExpt_].first;
            Int_t curBestEntry  = bestNLL_[iExpt_].second;
            if ( ( NLL_ < curBestNLL ) || ( curBestEntry == -1 ) ) {
                bestNLL_[iExpt_] = std::make_pair( NLL_, j );
            }

            Double_t curWorstNLL = worstNLL_[iExpt_].first;
            Int_t curWorstEntry  = worstNLL_[iExpt_].second;
            if ( ( NLL_ > curWorstNLL ) || ( curWorstEntry == -1 ) ) {
                worstNLL_[iExpt_] = std::make_pair( NLL_, j );
            }
        }
    }
    std::cout << "Finished storing best entry info.\n" << std::endl;

    std::cout << "Creating NLL histograms..." << std::flush;
    TH1* histo( 0 );
    for ( Int_t i( 0 ); i < numExpts; ++i ) {
        Double_t min   = bestNLL_[i].first;
        Double_t max   = worstNLL_[i].first;
        Double_t range = max - min;
        if ( range < 1e-3 ) {
            min -= 0.005;
            max += 0.005;
        } else {
            min -= range * 0.2;
            max += range * 0.2;
        }
        TString name( "expt" );
        name  += i;
        name  += "NLL";
        histo  = new TH1F( name, "", 100, min, max );
        for ( std::vector<Double_t>::const_iterator iter = allNLLs_[i].begin();
              iter != allNLLs_[i].end();
              ++iter ) {
            histo->Fill( *iter );
        }
        nllHistos_.insert( std::make_pair( i, histo ) );
    }
    std::cout << " done.\n" << std::endl;

    // now need to read all branches
    //this->setInputTreeBranchStatus(kTRUE);

    std::ofstream fout( "best-fit.txt" );

    // loop over the experiments, grab the best entry and store it
    std::cout << "Starting to retrieve best entries and fill output tree." << std::endl;
    for ( Int_t i( 0 ); i < numExpts; ++i ) {
        Int_t bestEntry = bestNLL_[i].second;
        if ( bestEntry != -1 ) {
            inputTree_->GetEntry( bestEntry );
            outputTree_->Fill();
        }
        if ( ( numExpts < 100 ) || ( i % ( numExpts / 100 ) == 0 ) ) {
            std::cout << "Writing out experiment " << i << std::endl;
        }
        TString bestFit( inputTree_->GetCurrentFile()->GetName() );
        bestFit.Remove( 0, 3 );
        Int_t index = bestFit.Index( "_" );
        if ( index < 1 ) {
            index = bestFit.Index( "." );
        }
        bestFit.Remove( index );
        fout << "Experiment " << i << " BestFit " << bestFit << std::endl;
    }
    std::cout << "Finished filling best entries in output tree.\n" << std::endl;
    fout.close();

    std::cout << "Writing output file." << std::endl;

    this->writeFile();
}

void LauResultsExtractor::writeFile()
{
    if ( ! outputFile_ ) {
        outputFile_ = new TFile( outputFileName_, "recreate" );
    }
    for ( std::map<Int_t, TH1*>::iterator iter = nllHistos_.begin(); iter != nllHistos_.end();
          ++iter ) {
        ( iter->second )->SetDirectory( outputFile_ );
    }
    outputTree_->SetDirectory( outputFile_ );
    outputFile_->Write();
    outputFile_->Close();
    delete outputFile_;
    outputFile_ = 0;
    nllHistos_.clear();
}
