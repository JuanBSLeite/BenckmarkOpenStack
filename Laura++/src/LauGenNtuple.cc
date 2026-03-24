
/*
Copyright 2006 University of Warwick

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

/*! \file LauGenNtuple.cc
    \brief File containing implementation of LauGenNtuple class.
*/

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

#include "LauGenNtuple.hh"

#include "TFile.h"
#include "TTree.h"

LauGenNtuple::LauGenNtuple( const TString& rootFileName, const TString& rootTreeName ) :
    rootFileName_( rootFileName ),
    rootTreeName_( rootTreeName ),
    rootFile_( 0 ),
    rootTree_( 0 ),
    definedBranches_( kFALSE )
{
    this->createFileAndTree();
}

LauGenNtuple::~LauGenNtuple()
{
    // seems that closing the file deletes the tree
    // so only delete if the file is still open for some reason
    if ( rootFile_ && rootFile_->IsOpen() ) {
        delete rootTree_;
        rootTree_ = 0;
    }
    delete rootFile_;
    rootFile_ = 0;
}

void LauGenNtuple::createFileAndTree()
{
    // first check whether we've already opened up the file or not
    if ( ! rootFile_ ) {
        // if not, first check the filename and if all ok create the file
        if ( rootFileName_ == "" ) {
            cerr << "ERROR in LauGenNtuple::createFileAndTree : Bad filename supplied, not creating file or tree."
                 << endl;
            return;
        }
        rootFile_ = TFile::Open( rootFileName_, "recreate" );
        if ( ! rootFile_ || rootFile_->IsZombie() || ! rootFile_->IsWritable() ) {
            cerr << "ERROR in LauGenNtuple::createFileAndTree : Problem opening file \""
                 << rootFileName_ << "\" for writing, not creating tree." << endl;
            return;
        }
    }
    // check whether we've already created the tree
    if ( ! rootTree_ ) {
        // if not change to the file's directory and create the tree
        rootFile_->cd();
        rootTree_ = new TTree( rootTreeName_, rootTreeName_ );
        rootTree_->SetDirectory( rootFile_ );
        this->definedBranches( kFALSE );
    }
}

void LauGenNtuple::addIntegerBranch( const TString& name )
{
    if ( this->definedBranches() ) {
        cerr << "ERROR in LauGenNtuple::addIntegerBranch : Already defined branches, can't add further ones."
             << endl;
        return;
    }
    this->setIntegerBranchValue( name, 0 );
}

void LauGenNtuple::addDoubleBranch( const TString& name )
{
    if ( this->definedBranches() ) {
        cerr << "ERROR in LauGenNtuple::addDoubleBranch : Already defined branches, can't add further ones."
             << endl;
        return;
    }
    this->setDoubleBranchValue( name, 0.0 );
}

void LauGenNtuple::setIntegerBranchValue( const TString& name, Int_t value )
{
    intVars_[name] = value;
}

void LauGenNtuple::setDoubleBranchValue( const TString& name, Double_t value )
{
    doubleVars_[name] = value;
}

Int_t LauGenNtuple::getIntegerBranchValue( const TString& name ) const
{
    IntVarMap::const_iterator iter = intVars_.find( name );
    if ( iter == intVars_.end() ) {
        cerr << "ERROR in LauGenNtuple::getIntegerBranchValue : no such branch \"" << name << "\"."
             << endl;
        return -99;
    } else {
        return iter->second;
    }
}

Double_t LauGenNtuple::getDoubleBranchValue( const TString& name ) const
{
    DoubleVarMap::const_iterator iter = doubleVars_.find( name );
    if ( iter == doubleVars_.end() ) {
        cerr << "ERROR in LauGenNtuple::getDoubleBranchValue : no such branch \"" << name << "\"."
             << endl;
        return -99.0;
    } else {
        return iter->second;
    }
}

void LauGenNtuple::defineBranches()
{
    if ( this->definedBranches() ) {
        cerr << "ERROR in LauGenNtuple::defineBranches : Already defined branches, not doing it again."
             << endl;
        return;
    }
    for ( IntVarMap::iterator iter = intVars_.begin(); iter != intVars_.end(); ++iter ) {
        TString name   = iter->first;
        Int_t* pointer = &( iter->second );
        TString thirdPart( name );
        thirdPart += "/I";
        rootTree_->Branch( name, pointer, thirdPart );
    }
    for ( DoubleVarMap::iterator iter = doubleVars_.begin(); iter != doubleVars_.end(); ++iter ) {
        TString name      = iter->first;
        Double_t* pointer = &( iter->second );
        TString thirdPart( name );
        thirdPart += "/D";
        rootTree_->Branch( name, pointer, thirdPart );
    }
    this->definedBranches( kTRUE );
}

void LauGenNtuple::fillBranches()
{
    if ( ! rootTree_ ) {
        cerr << "ERROR in LauGenNtuple::fillBranches : Tree not created, cannot fill branches."
             << endl;
        return;
    } else if ( ! this->definedBranches() ) {
        this->defineBranches();
    }
    rootTree_->Fill();
}

void LauGenNtuple::deleteAndRecreateTree()
{
    if ( rootTree_ ) {
        delete rootTree_;
        rootTree_ = 0;
    }
    this->createFileAndTree();
}

Int_t LauGenNtuple::buildIndex( const TString& majorName, const TString& minorName )
{
    if ( ! rootTree_ ) {
        cerr << "ERROR in LauGenNtuple::buildIndex : Tree not created, cannot build index." << endl;
        return -1;
    }
    return rootTree_->BuildIndex( majorName, minorName );
}

void LauGenNtuple::writeOutGenResults()
{
    // Write out the generated ntuple

    // Check that the file is open
    if ( rootFile_ == 0 ) {
        cerr << "ERROR in LauGenNtuple::writeOutGenResults : File not opened, can't write anything."
             << endl;
        return;
    }

    // Check that the tree exists and if so then make sure we have the
    // up to date pointer to the file since if it splits via the
    // TTree::ChangeFile mechanism we're left with a dangling pointer
    if ( rootTree_ != 0 ) {
        rootFile_ = rootTree_->GetCurrentFile();
    }
    rootFile_->cd();
    rootTree_->Write( "", TObject::kOverwrite );
    rootFile_->Close();
    delete rootFile_;
    rootFile_ = 0;
}

void LauGenNtuple::addFriendTree( const TString& rootFileName, const TString& rootTreeName )
{
    if ( ! rootTree_ ) {
        cerr << "ERROR in LauGenNtuple::addFriendTree : Tree not created, cannot add friend." << endl;
        return;
    }
    rootTree_->AddFriend( rootTreeName, rootFileName );
}
