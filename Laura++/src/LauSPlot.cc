
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

/*! \file LauSPlot.cc
    \brief File containing implementation of LauSPlot class.

    Class for defining the SPlot technique based on TSplot from ROOT by the following authors:
    Muriel Pivk, Anna Kreshuk (10/2005).
    (Original copyright notice below)
*/

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 ROOT Foundation,  CERN/PH-SFT                   *
 *                                                                    *
 **********************************************************************/

#include <cfloat>
#include <iostream>
#include <vector>
using std::cerr;
using std::cout;
using std::endl;

#include "LauSPlot.hh"

#include "TEventList.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVirtualFitter.h"

extern void Yields( Int_t&, Double_t*, Double_t& f, Double_t* x, Int_t iflag );

LauSPlot::LauSPlot( const TString& fileName,
                    const TString& treeName,
                    Int_t firstExpt,
                    Int_t nExpt,
                    const NameSet& variableNames,
                    const NumbMap& freeSpecies,
                    const NumbMap& fixdSpecies,
                    const TwoDMap& twodimPDFs,
                    Bool_t sigSplit,
                    Bool_t scfDPSmeared ) :
    fileName_( fileName ),
    inputTreeName_( treeName ),
    cnTreeName_( "" ),
    sweightTreeName_( "" ),
    file_( 0 ),
    inputTree_( 0 ),
    cnTree_( 0 ),
    sweightTree_( 0 ),
    eventList_( 0 ),
    variableNames_( variableNames ),
    freeSpecies_( freeSpecies ),
    fixdSpecies_( fixdSpecies ),
    origFreeSpecies_( freeSpecies ),
    origFixdSpecies_( fixdSpecies ),
    twodimPDFs_( twodimPDFs ),
    signalSplit_( sigSplit ),
    scfDPSmear_( scfDPSmeared ),
    readInput_( kFALSE ),
    definedCNBranches_( kFALSE ),
    definedSWeightBranches_( kFALSE ),
    firstExpt_( firstExpt ),
    nExpt_( nExpt ),
    iExpt_( 0 ),
    nEvents_( 0 ),
    nDiscVars_( variableNames.size() ),
    nFreeSpecies_( freeSpecies.size() ),
    nFixdSpecies_( fixdSpecies.size() ),
    nSpecies_( freeSpecies.size() + fixdSpecies.size() )
{
    this->openInputFileAndTree();
    this->readInputInfo();
    this->createSWeightTree();
    if ( nFixdSpecies_ > 0 ) {
        this->createCNTree();
    }
}

LauSPlot::~LauSPlot()
{
    // seems that closing the file deletes the tree
    // so only delete if the file is still open
    if ( file_ && file_->IsOpen() ) {
        delete inputTree_;
        inputTree_ = 0;
        delete sweightTree_;
        sweightTree_ = 0;
        delete cnTree_;
        cnTree_ = 0;
    }
    delete file_;
    file_ = 0;
}

void LauSPlot::openInputFileAndTree()
{
    // first check whether we've already opened up the file or not
    if ( ! file_ ) {
        // if not, first check the filename and if all ok create the file
        if ( fileName_ == "" ) {
            cerr << "ERROR in LauSPlot::createFileAndTree : Bad filename supplied, not creating file or tree."
                 << endl;
            return;
        }
        file_ = TFile::Open( fileName_, "update" );
        if ( ! file_ || file_->IsZombie() || ! file_->IsWritable() ) {
            cerr << "ERROR in LauSPlot::createFileAndTree : Problem opening file \"" << fileName_
                 << "\" for updating, can't do anything." << endl;
            return;
        }
    }
    // next open the input tree for reading
    if ( ! inputTree_ ) {
        file_->cd();
        inputTree_ = dynamic_cast<TTree*>( file_->Get( inputTreeName_ ) );
        inputTree_->SetDirectory( file_ );
    }
}

void LauSPlot::readInputInfo()
{
    // Read the tree and then setup the maps so we know which leaves to
    // read from the tree to get the various PDF values

    Bool_t inputOK = this->readInputLeaves();
    if ( ! inputOK ) {
        this->readInput( inputOK );
        return;
    }

    inputOK = this->checkLeaves();
    this->readInput( inputOK );
    return;
}

Bool_t LauSPlot::readInputLeaves()
{
    // Read all the leaves in the tree and store them in the leaves map

    if ( ! inputTree_ ) {
        cerr << "ERROR in LauSPlot::readInputInfo : Invalid pointer to data tree." << endl;
        return kFALSE;
    }

    Int_t nBranches    = inputTree_ ? static_cast<Int_t>( inputTree_->GetNbranches() ) : 0;
    TObjArray* pLeaves = inputTree_->GetListOfLeaves();
    if ( ! pLeaves ) {
        cerr << "ERROR in LauSPlot::readInputInfo : Problem retrieving leaves from the tree." << endl;
        return kFALSE;
    }
    TObjArray& leaves = *pLeaves;

    if ( nBranches > leaves.GetSize() ) {
        cerr << "ERROR in LauSPlot::readInputInfo : List of leaves is smaller than number of branches - this is strange!"
             << endl;
        return kFALSE;
    }

    for ( Int_t iLeaf = 0; iLeaf < nBranches; ++iLeaf ) {

        TLeaf* leaf = dynamic_cast<TLeaf*>( leaves[iLeaf] );

        // we can't deal with arrays
        Int_t size = leaf->GetNdata();
        if ( size != 1 ) {
            cerr << "ERROR in LauSPlot::readInputInfo : Tree has array branches, can't deal with those."
                 << endl;
            return kFALSE;
        }

        // find the name of the leaf
        TString name = leaf->GetName();

        // initialise an entry in the maps to hold the value
        leaves_[name] = leaf;
    }

    return kTRUE;
}

Bool_t LauSPlot::checkLeaves() const
{
    // Analyse the names of the leaves to check that we have a leaf for
    // all bits of information we require, i.e. a likelihood value for
    // each combination of variable and species

    // If we have 2D PDFs then we have to look for some extra leaves
    for ( TwoDMap::const_iterator twodim_iter = twodimPDFs_.begin(); twodim_iter != twodimPDFs_.end();
          ++twodim_iter ) {

        const TString& specName      = twodim_iter->first;
        const TString& firstVarName  = twodim_iter->second.first;
        const TString& secondVarName = twodim_iter->second.second;

        TString expectedName( specName );
        expectedName += firstVarName;
        expectedName += secondVarName;
        expectedName += "Like";

        Bool_t found( kFALSE );
        for ( LeafMap::const_iterator leaf_iter = leaves_.begin(); leaf_iter != leaves_.end();
              ++leaf_iter ) {
            const TString& leafName = leaf_iter->first;
            const TLeaf* leaf       = leaf_iter->second;

            if ( ( leafName == expectedName ) && ( leaf != 0 ) ) {
                found = kTRUE;
                break;
            }
        }

        if ( ! found ) {
            cerr << "ERROR in LauSPlot::checkLeaves : Could not find expected leaf \""
                 << expectedName << "\"." << endl;
            return kFALSE;
        }
    }

    // Search for all other "standard" leaves, i.e. <species><var>Like

    for ( NumbMap::const_iterator fixd_iter = fixdSpecies_.begin(); fixd_iter != fixdSpecies_.end();
          ++fixd_iter ) {

        const TString& specName = fixd_iter->first;

        // if the signal is split we need to do a dedicated search
        // for sigTM and sigSCF, so skip the signal here
        if ( specName == "sig" && this->signalSplit() ) {
            continue;
        }

        for ( NameSet::const_iterator vars_iter = variableNames_.begin();
              vars_iter != variableNames_.end();
              ++vars_iter ) {
            const TString& varName = ( *vars_iter );
            TString expectedName( specName );
            expectedName += varName;
            expectedName += "Like";
            Bool_t found( kFALSE );
            for ( LeafMap::const_iterator leaf_iter = leaves_.begin(); leaf_iter != leaves_.end();
                  ++leaf_iter ) {
                const TString& leafName = leaf_iter->first;
                const TLeaf* leaf       = leaf_iter->second;

                if ( ( leafName == expectedName ) && ( leaf != 0 ) ) {
                    found = kTRUE;
                    break;
                }
            }

            if ( ! found ) {
                cerr << "ERROR in LauSPlot::checkLeaves : Could not find expected leaf \""
                     << expectedName << "\"." << endl;
                return kFALSE;
            }
        }
    }

    for ( NumbMap::const_iterator free_iter = freeSpecies_.begin(); free_iter != freeSpecies_.end();
          ++free_iter ) {

        const TString& specName = free_iter->first;

        // if the signal is split we need to do a dedicated search
        // for sigTM and sigSCF, so skip the signal here
        if ( specName == "sig" && this->signalSplit() ) {
            continue;
        }

        for ( NameSet::const_iterator vars_iter = variableNames_.begin();
              vars_iter != variableNames_.end();
              ++vars_iter ) {
            const TString& varName = ( *vars_iter );
            TString expectedName( specName );
            expectedName += varName;
            expectedName += "Like";
            Bool_t found( kFALSE );
            for ( LeafMap::const_iterator leaf_iter = leaves_.begin(); leaf_iter != leaves_.end();
                  ++leaf_iter ) {
                const TString& leafName = leaf_iter->first;
                const TLeaf* leaf       = leaf_iter->second;

                if ( ( leafName == expectedName ) && ( leaf != 0 ) ) {
                    found = kTRUE;
                    break;
                }
            }

            if ( ! found ) {
                cerr << "ERROR in LauSPlot::checkLeaves : Could not find expected leaf \""
                     << expectedName << "\"." << endl;
                return kFALSE;
            }
        }
    }

    if ( this->signalSplit() ) {

        // now need to search for the sigTM and sigSCF leaves

        for ( NameSet::const_iterator vars_iter = variableNames_.begin();
              vars_iter != variableNames_.end();
              ++vars_iter ) {
            const TString& varName = ( *vars_iter );
            TString expectedName( "sigTM" );
            expectedName += varName;
            expectedName += "Like";
            Bool_t found( kFALSE );
            for ( LeafMap::const_iterator leaf_iter = leaves_.begin(); leaf_iter != leaves_.end();
                  ++leaf_iter ) {
                const TString& leafName = leaf_iter->first;
                const TLeaf* leaf       = leaf_iter->second;

                if ( ( leafName == expectedName ) && ( leaf != 0 ) ) {
                    found = kTRUE;
                    break;
                }
            }

            if ( ! found ) {
                cerr << "ERROR in LauSPlot::checkLeaves : Could not find expected leaf \""
                     << expectedName << "\"." << endl;
                return kFALSE;
            }

            expectedName  = "sigSCF";
            expectedName += varName;
            expectedName += "Like";
            found         = kFALSE;
            for ( LeafMap::const_iterator leaf_iter = leaves_.begin(); leaf_iter != leaves_.end();
                  ++leaf_iter ) {
                const TString& leafName = leaf_iter->first;
                const TLeaf* leaf       = leaf_iter->second;

                if ( ( leafName == expectedName ) && ( leaf != 0 ) ) {
                    found = kTRUE;
                    break;
                }
            }

            if ( ! found ) {
                cerr << "ERROR in LauSPlot::checkLeaves : Could not find expected leaf \""
                     << expectedName << "\"." << endl;
                return kFALSE;
            }

            expectedName = "sigSCFFrac";
            found        = kFALSE;
            for ( LeafMap::const_iterator leaf_iter = leaves_.begin(); leaf_iter != leaves_.end();
                  ++leaf_iter ) {
                const TString& leafName = leaf_iter->first;
                const TLeaf* leaf       = leaf_iter->second;

                if ( ( leafName == expectedName ) && ( leaf != 0 ) ) {
                    found = kTRUE;
                    break;
                }
            }

            if ( ! found ) {
                cerr << "ERROR in LauSPlot::checkLeaves : Could not find expected leaf \""
                     << expectedName << "\"." << endl;
                return kFALSE;
            }
        }
    }

    return kTRUE;
}

void LauSPlot::createCNTree()
{
    // check whether we've already created the tree
    if ( ! cnTree_ ) {
        // if not change to the file's directory and create the tree
        cnTreeName_  = inputTreeName_;
        cnTreeName_ += "_cNCoeffs";
        file_->cd();
        cnTree_ = new TTree( cnTreeName_, cnTreeName_ );
        cnTree_->SetDirectory( file_ );
        this->definedCNBranches( kFALSE );
    }
}

void LauSPlot::createSWeightTree()
{
    // check whether we've already created the tree
    if ( ! sweightTree_ ) {
        // if not change to the file's directory and create the tree
        sweightTreeName_  = inputTreeName_;
        sweightTreeName_ += "_sWeights";
        file_->cd();
        sweightTree_ = new TTree( sweightTreeName_, sweightTreeName_ );
        sweightTree_->SetDirectory( file_ );
        this->definedSWeightBranches( kFALSE );
    }
}

void LauSPlot::defineCNBranches()
{
    if ( this->definedCNBranches() ) {
        cerr << "ERROR in LauSPlot::defineCNBranches : Already defined branches, not doing it again."
             << endl;
        return;
    }
    if ( cN_.empty() ) {
        cerr << "ERROR in LauSPlot::defineCNBranches : No entries in the cN container, can't define branches."
             << endl;
        return;
    }
    if ( ! cnTree_ ) {
        cerr << "ERROR in LauSPlot::defineCNBranches : Invalid pointer to the tree, can't define branches."
             << endl;
        return;
    }

    // In the cN tree there is one entry per experiment, so need to know the experiment number.
    cnTree_->Branch( "iExpt", &iExpt_, "iExpt/I" );

    for ( std::map<TString, NumbMap>::iterator excl_iter = cN_.begin(); excl_iter != cN_.end();
          ++excl_iter ) {
        const TString& exclName = excl_iter->first;
        NumbMap& species        = excl_iter->second;
        for ( NumbMap::iterator spec_iter = species.begin(); spec_iter != species.end();
              ++spec_iter ) {
            const TString& specName = spec_iter->first;
            Double_t* pointer       = &( spec_iter->second );
            TString name( specName );
            name += "_cN";
            if ( exclName == "none" ) {
                name += "_all";
            } else {
                name += "_no";
                name += exclName;
            }
            TString thirdPart( name );
            thirdPart += "/D";
            cnTree_->Branch( name, pointer, thirdPart );
        }
    }
    this->definedCNBranches( kTRUE );
}

void LauSPlot::defineSWeightBranches()
{
    if ( this->definedSWeightBranches() ) {
        cerr << "ERROR in LauSPlot::defineSWeightBranches : Already defined branches, not doing it again."
             << endl;
        return;
    }
    if ( sWeights_.empty() ) {
        cerr << "ERROR in LauSPlot::defineSWeightBranches : No entries in the sWeights container, can't define branches."
             << endl;
        return;
    }
    if ( ! sweightTree_ ) {
        cerr << "ERROR in LauSPlot::defineSWeightBranches : Invalid pointer to the tree, can't define branches."
             << endl;
        return;
    }

    for ( std::map<TString, NumbMap>::iterator excl_iter = sWeightsCurrent_.begin();
          excl_iter != sWeightsCurrent_.end();
          ++excl_iter ) {
        const TString& exclName = excl_iter->first;
        NumbMap& species        = excl_iter->second;
        for ( NumbMap::iterator spec_iter = species.begin(); spec_iter != species.end();
              ++spec_iter ) {
            const TString& specName = spec_iter->first;
            Double_t* pointer       = &( spec_iter->second );
            TString name( specName );
            name += "_sWeight";
            if ( exclName == "none" ) {
                name += "_all";
            } else {
                name += "_no";
                name += exclName;
            }
            TString thirdPart( name );
            thirdPart += "/D";
            sweightTree_->Branch( name, pointer, thirdPart );
        }
    }
    this->definedSWeightBranches( kTRUE );
}

void LauSPlot::setExperiment( Int_t iExpt )
{
    if ( ! inputTree_ ) {
        cerr << "ERROR in LauSPlot::setExperiment : Invalid pointer to data tree." << endl;
        return;
    }

    // create the event list if we haven't already done so
    if ( ! eventList_ ) {
        eventList_ = new TEventList( "splotelist", "splotelist" );
        eventList_->SetDirectory( file_ );
    }

    // fill the event list with this experiment's events
    TString listName( eventList_->GetName() );
    listName.Prepend( ">>" );
    TString selection( "iExpt==" );
    selection += iExpt;
    cout << "LauSPlot::setExperiment : Setting tree to experiment number " << iExpt << "." << endl;
    inputTree_->Draw( listName, selection );

    // find how many events there are
    nEvents_ = eventList_->GetN();
    cout << "                          Found " << nEvents_ << " events." << endl;

    // make sure we have enough space in the per-event value vectors
    pdfTot_.clear();
    pdfTot_.resize( nEvents_ );
    discPdf_.clear();
    discPdf_.resize( nEvents_ );
    scfFrac_.clear();
    scfFrac_.resize( nEvents_ );
    sWeights_.clear();
    sWeights_.resize( nEvents_ );

    // read the info for this experiment
    this->readExpt();
}

void LauSPlot::readExpt()
{
    for ( Int_t iEvt( 0 ); iEvt < nEvents_; ++iEvt ) {
        // Find which entry from the full tree contains the requested event
        Long64_t iEntry = eventList_ ? eventList_->GetEntry( iEvt ) : iEvt;
        if ( iEntry < 0 ) {    // this shouldn't happen, but just in case...
            cerr << "ERROR in LauSPlot::readExpt : Problem retrieving event." << endl;
            gSystem->Exit( EXIT_FAILURE );
        }

        // Then retrieve that entry from the tree
        inputTree_->GetEntry( iEntry );

        // If needed retrieve the SCF fraction values
        if ( this->signalSplit() ) {
            for ( LeafMap::const_iterator leaf_iter = leaves_.begin(); leaf_iter != leaves_.end();
                  ++leaf_iter ) {
                const TString& leafName = leaf_iter->first;
                const TLeaf* leaf       = leaf_iter->second;
                if ( ( leafName == "sigSCFFrac" ) && ( leaf != 0 ) ) {
                    scfFrac_[iEvt] = leaf->GetValue();
                    break;
                }
            }
        }

        // Copy the leaf values into discPdf_
        for ( TwoDMap::const_iterator twodim_iter = twodimPDFs_.begin();
              twodim_iter != twodimPDFs_.end();
              ++twodim_iter ) {
            const TString& specName      = twodim_iter->first;
            const TString& firstVarName  = twodim_iter->second.first;
            const TString& secondVarName = twodim_iter->second.second;

            TString varName = firstVarName + secondVarName;

            TString expectedName( specName );
            expectedName += varName;
            expectedName += "Like";

            for ( LeafMap::const_iterator leaf_iter = leaves_.begin(); leaf_iter != leaves_.end();
                  ++leaf_iter ) {
                const TString& leafName = leaf_iter->first;
                const TLeaf* leaf       = leaf_iter->second;
                if ( ( leafName == expectedName ) && ( leaf != 0 ) ) {
                    discPdf_[iEvt][specName][varName] = leaf->GetValue();
                    break;
                }
            }
        }

        Bool_t needSignalSearch( kFALSE );
        for ( NumbMap::const_iterator fixd_iter = fixdSpecies_.begin();
              fixd_iter != fixdSpecies_.end();
              ++fixd_iter ) {
            const TString& specName = fixd_iter->first;

            // if the signal is split we need to do a dedicated search
            // for sigTM and sigSCF, so skip the signal here
            if ( specName == "sig" && this->signalSplit() ) {
                needSignalSearch = kTRUE;
                continue;
            }

            for ( NameSet::const_iterator vars_iter = variableNames_.begin();
                  vars_iter != variableNames_.end();
                  ++vars_iter ) {
                const TString& varName = ( *vars_iter );
                TString expectedName( specName );
                expectedName += varName;
                expectedName += "Like";
                for ( LeafMap::const_iterator leaf_iter = leaves_.begin(); leaf_iter != leaves_.end();
                      ++leaf_iter ) {
                    const TString& leafName = leaf_iter->first;
                    const TLeaf* leaf       = leaf_iter->second;
                    if ( ( leafName == expectedName ) && ( leaf != 0 ) ) {
                        discPdf_[iEvt][specName][varName] = leaf->GetValue();
                        break;
                    }
                }
            }
        }

        for ( NumbMap::const_iterator free_iter = freeSpecies_.begin();
              free_iter != freeSpecies_.end();
              ++free_iter ) {
            const TString& specName = free_iter->first;

            // if the signal is split we need to do a dedicated search
            // for sigTM and sigSCF, so skip the signal here
            if ( specName == "sig" && this->signalSplit() ) {
                needSignalSearch = kTRUE;
                continue;
            }

            for ( NameSet::const_iterator vars_iter = variableNames_.begin();
                  vars_iter != variableNames_.end();
                  ++vars_iter ) {
                const TString& varName = ( *vars_iter );
                TString expectedName( specName );
                expectedName += varName;
                expectedName += "Like";
                for ( LeafMap::const_iterator leaf_iter = leaves_.begin(); leaf_iter != leaves_.end();
                      ++leaf_iter ) {
                    const TString& leafName = leaf_iter->first;
                    const TLeaf* leaf       = leaf_iter->second;
                    if ( ( leafName == expectedName ) && ( leaf != 0 ) ) {
                        discPdf_[iEvt][specName][varName] = leaf->GetValue();
                        break;
                    }
                }
            }
        }

        if ( needSignalSearch ) {
            for ( NameSet::const_iterator vars_iter = variableNames_.begin();
                  vars_iter != variableNames_.end();
                  ++vars_iter ) {
                const TString& varName = ( *vars_iter );
                TString specName( "sigTM" );
                TString expectedName( specName );
                expectedName += varName;
                expectedName += "Like";
                for ( LeafMap::const_iterator leaf_iter = leaves_.begin(); leaf_iter != leaves_.end();
                      ++leaf_iter ) {
                    const TString& leafName = leaf_iter->first;
                    const TLeaf* leaf       = leaf_iter->second;
                    if ( ( leafName == expectedName ) && ( leaf != 0 ) ) {
                        discPdf_[iEvt][specName][varName] = leaf->GetValue();
                        break;
                    }
                }
                specName      = "sigSCF";
                expectedName  = specName;
                expectedName += varName;
                expectedName += "Like";
                for ( LeafMap::const_iterator leaf_iter = leaves_.begin(); leaf_iter != leaves_.end();
                      ++leaf_iter ) {
                    const TString& leafName = leaf_iter->first;
                    const TLeaf* leaf       = leaf_iter->second;
                    if ( ( leafName == expectedName ) && ( leaf != 0 ) ) {
                        discPdf_[iEvt][specName][varName] = leaf->GetValue();
                        break;
                    }
                }
            }
        }
    }
}

void LauSPlot::runCalculations( const TString& option )
{
    // Calculates the sWeights and cN coeffs
    // The option controls the print level:
    // "Q" - no print out (default)
    // "V" - prints the estimated # of events in species
    // "VV" - as "V" + the MINUIT printing + sums of weights for control

    if ( ! this->readInput() ) {
        cerr << "ERROR in LauSPlot::runCalculations : The input ntuple has not been successfully read, can't calculate anything."
             << endl;
        return;
    }

    if ( freeSpecies_.empty() ) {
        cerr << "ERROR in LauSPlot::runCalculations : Numbers of events in species have not been set."
             << endl;
        return;
    }

    TString opt( option );
    opt.ToUpper();
    opt.ReplaceAll( "VV", "W" );

    // Make sure that global fitter is MINUIT
    this->checkFitter();

    // Loop over experiments
    for ( iExpt_ = firstExpt_; iExpt_ < ( firstExpt_ + nExpt_ ); ++iExpt_ ) {

        this->setExperiment( iExpt_ );

        if ( nEvents_ < 1 ) {
            cerr << "ERROR in LauSPlot::runCalculations : Zero events in experiment " << iExpt_
                 << ", skipping..." << endl;
            continue;
        }

        // Now loop over the PDFs to exclude, including the case where none are to be excluded.
        NameSet excludePdf;
        if ( variableNames_.size() < 2 ) {
            excludePdf.insert( "none" );
        } else {
            excludePdf = variableNames_;
            excludePdf.insert( "none" );
        }

        for ( NameSet::const_iterator excl_iter = excludePdf.begin(); excl_iter != excludePdf.end();
              ++excl_iter ) {

            const TString& exclName = ( *excl_iter );

            cout << "LauSPlot::runCalculations : Calculating sWeights, excluding PDF: " << exclName
                 << "." << endl;

            // Calculate the per-event total PDF values for each species.
            this->calcTotPDFValues( exclName );

            // Reset the fitter
            this->initialiseFitter( opt );

            // Set the parameters to their initial values
            this->setFitParameters();

            // Run the fit
            this->runFit();

            // Get the fitted parameter values back
            this->retrieveFittedParameters( opt );

            // Get the covariance matrix
            Bool_t covMatOK = this->calcCovMatrix();
            Double_t* covmat( 0 );
            if ( ! covMatOK ) {
                TVirtualFitter* fitter = TVirtualFitter::GetFitter();
                covmat                 = fitter->GetCovarianceMatrix();
            }
            if ( opt.Contains( "W" ) ) {
                this->printCovMatrixElements( covmat );
            }

            // calculate the cN and sWeights from the covariance matrix
            if ( nFixdSpecies_ > 0 ) {
                this->calcCNCoeffs( exclName, covmat );
            }
            this->calcSWeights( exclName, covmat );

            // print verbose info if required
            if ( opt.Contains( "W" ) ) {
                this->printSumOfWeights( exclName );
            }
        }

        // Finally fill all the branches for this experiment
        if ( nFixdSpecies_ > 0 ) {
            this->fillCNBranches();
        }
        this->fillSWeightBranches();
    }
}

void LauSPlot::checkFitter() const
{
    TString minuitName( "TFitter" );
    TVirtualFitter* fitter = TVirtualFitter::GetFitter();
    if ( fitter ) {
        TString fitterName( fitter->IsA()->GetName() );
        if ( fitterName != minuitName ) {
            delete fitter;
            fitter = 0;
        }
    }
    fitter = TVirtualFitter::Fitter( 0 );
}

void LauSPlot::initialiseFitter( const TString& opt )
{
    TVirtualFitter* fitter = TVirtualFitter::GetFitter();
    fitter->Clear();
    fitter->SetFCN( Yields );
    fitter->SetObjectFit( this );

    // Set the print level
    Double_t arglist[10];
    if ( opt.Contains( "Q" ) ) {
        arglist[0] = -1;
    }
    if ( opt.Contains( "V" ) ) {
        arglist[0] = 0;
    }
    if ( opt.Contains( "W" ) ) {
        arglist[0] = 1;
    }
    fitter->ExecuteCommand( "SET PRINT", arglist, 1 );

    // Need to set the "SET ERR" command to +0.5 for +/-1 sigma errors
    // for maximum likelihood fit.
    arglist[0] = 0.5;
    fitter->ExecuteCommand( "SET ERR", arglist, 1 );
}

void LauSPlot::setFitParameters() const
{
    TVirtualFitter* fitter = TVirtualFitter::GetFitter();

    // must add the parameters in the same order as they are stored in pdfTot_
    Int_t ispecies( 0 );
    const NumbMap& species = pdfTot_.front();

    for ( NumbMap::const_iterator spec_iter = species.begin(); spec_iter != species.end();
          ++spec_iter ) {
        const TString& name( spec_iter->first );

        // starting parameters should be the original values,
        // not those that came out of the last fit
        NumbMap::const_iterator free_iter = origFreeSpecies_.find( name );
        NumbMap::const_iterator fixd_iter = origFixdSpecies_.find( name );
        Bool_t fixed                      = fixd_iter != origFixdSpecies_.end();

        Double_t value( 0.0 );
        if ( fixed ) {
            value = fixd_iter->second;
        } else {
            value = free_iter->second;
        }
        fitter->SetParameter( ispecies, name, value, 1, 0, 0 );
        if ( fixed ) {
            fitter->FixParameter( ispecies );
        }
        ++ispecies;
    }
}

void LauSPlot::runFit()
{
    TVirtualFitter* fitter = TVirtualFitter::GetFitter();

    Double_t arglist[10];
    arglist[0] = 1000 * nFreeSpecies_;    // maximum iterations
    arglist[1] = 0.05;                    // tolerance : min EDM = 0.001*tolerance

    // Execute MIGRAD
    Int_t fitStatus = fitter->ExecuteCommand( "MIGRAD", arglist, 2 );

    if ( fitStatus != 0 ) {
        cerr << "ERROR in LauSPlot::runFit : Error during MIGRAD minimisation." << endl;
    } else {
        // Execute HESSE
        fitStatus = fitter->ExecuteCommand( "HESSE", arglist, 1 );
        if ( fitStatus != 0 ) {
            cerr << "ERROR in LauSPlot::runFit : Error during HESSE error calculation." << endl;
        }
    }
}

void LauSPlot::printCovMatrixElements( const Double_t* covmat ) const
{
    Int_t a( 0 );
    cout << endl;
    for ( NumbMap::const_iterator iter_a = freeSpecies_.begin(); iter_a != freeSpecies_.end();
          ++iter_a ) {
        Int_t b( 0 );
        for ( NumbMap::const_iterator iter_b = freeSpecies_.begin(); iter_b != freeSpecies_.end();
              ++iter_b ) {
            if ( covmat ) {
                cout << "CovMat Elem: [" << iter_a->first << "," << iter_b->first
                     << "] = " << covmat[a * nSpecies_ + b] << endl;
            } else {
                cout << "CovMat Elem: [" << iter_a->first << "," << iter_b->first
                     << "] = " << covMat_( a, b ) << endl;
            }
            ++b;
        }
        ++a;
    }
    cout << endl;

    if ( ! covmat ) {
        covMat_.Print();
    }
}

void LauSPlot::retrieveFittedParameters( const TString& opt )
{
    TVirtualFitter* fitter = TVirtualFitter::GetFitter();

    // remember fit parameters are in same order as in pdfTot_
    Int_t ispecies( 0 );
    const NumbMap& species = pdfTot_.front();

    for ( NumbMap::const_iterator spec_iter = species.begin(); spec_iter != species.end();
          ++spec_iter ) {
        const TString& name( spec_iter->first );
        NumbMap::iterator free_iter = freeSpecies_.find( name );
        if ( free_iter != freeSpecies_.end() ) {
            free_iter->second = fitter->GetParameter( ispecies );
            if ( ! opt.Contains( "Q" ) ) {
                cout << "Estimated # of events in species " << name << " = " << free_iter->second
                     << endl;
            }
        }
        ++ispecies;
    }
    if ( ! opt.Contains( "Q" ) ) {
        cout << endl;
    }
}

void LauSPlot::printSumOfWeights( const TString& exclName ) const
{
    for ( NumbMap::const_iterator spec_iter = freeSpecies_.begin(); spec_iter != freeSpecies_.end();
          ++spec_iter ) {
        const TString& specName = spec_iter->first;
        Double_t sumweight( 0.0 );
        for ( Int_t iEvt( 0 ); iEvt < nEvents_; ++iEvt ) {
            const NumbMap& specWeightMap  = sWeights_[iEvt].find( exclName )->second;
            Double_t weight               = specWeightMap.find( specName )->second;
            sumweight                    += weight;
        }
        cout << "Sum of sWeights for species \"" << specName << "\" = " << sumweight << endl;
    }
    cout << endl;
}

Bool_t LauSPlot::calcCovMatrix()
{
    // Calculate our inverse covariance matrix from the various PDFs

    TMatrixD invMatrix( nFreeSpecies_, nFreeSpecies_ );

    // First calculate the denominator, which is common to all elements
    std::vector<Double_t> denominator( nEvents_ );
    for ( Int_t iEvt( 0 ); iEvt < nEvents_; ++iEvt ) {
        denominator[iEvt] = 0.0;
        for ( NumbMap::const_iterator spec_iter = freeSpecies_.begin();
              spec_iter != freeSpecies_.end();
              ++spec_iter ) {
            const TString& specName  = spec_iter->first;
            Double_t specNumEvents   = spec_iter->second;
            denominator[iEvt]       += specNumEvents * pdfTot_[iEvt][specName];
        }
        for ( NumbMap::const_iterator spec_iter = fixdSpecies_.begin();
              spec_iter != fixdSpecies_.end();
              ++spec_iter ) {
            const TString& specName  = spec_iter->first;
            Double_t specNumEvents   = spec_iter->second;
            denominator[iEvt]       += specNumEvents * pdfTot_[iEvt][specName];
        }
        // Square to get the final denominator
        denominator[iEvt] *= denominator[iEvt];
    }

    // Then calculate each element
    Int_t i( 0 );
    for ( NumbMap::const_iterator spec_iter_i = freeSpecies_.begin();
          spec_iter_i != freeSpecies_.end();
          ++spec_iter_i ) {
        const TString& specName_i = spec_iter_i->first;
        Int_t j( 0 );
        for ( NumbMap::const_iterator spec_iter_j = freeSpecies_.begin();
              spec_iter_j != freeSpecies_.end();
              ++spec_iter_j ) {
            const TString& specName_j = spec_iter_j->first;
            invMatrix( i, j )         = 0.0;
            for ( Int_t iEvt( 0 ); iEvt < nEvents_; ++iEvt ) {
                Double_t numerator  = pdfTot_[iEvt][specName_i] * pdfTot_[iEvt][specName_j];
                invMatrix( i, j )  += numerator / denominator[iEvt];
            }
            ++j;
        }
        ++i;
    }

    // Check for a singular matrix
    if ( invMatrix.Determinant() < 1e-15 ) {
        cerr << "ERROR in LauSPlot::calcCovMatrix : Calculated inverse covariance matrix is singular, can't invert it."
             << endl;
        return kFALSE;
    }

    // Invert and store in the covariance matrix
    covMat_.ResizeTo( nFreeSpecies_, nFreeSpecies_ );
    covMat_ = TMatrixD( TMatrixD::kInverted, invMatrix );

    return kTRUE;
}

void LauSPlot::calcTotPDFValues( const TString& exclName )
{
    for ( Int_t iEvt( 0 ); iEvt < nEvents_; iEvt++ ) {

        Bool_t needSignalSearch( kFALSE );

        for ( NumbMap::const_iterator spec_iter = fixdSpecies_.begin();
              spec_iter != fixdSpecies_.end();
              ++spec_iter ) {
            const TString& specName = spec_iter->first;

            // if the signal is split we need to do a dedicated search
            // for sigTM and sigSCF, so skip the signal here
            if ( specName == "sig" && this->signalSplit() ) {
                needSignalSearch = kTRUE;
                continue;
            }

            pdfTot_[iEvt][specName] = 1.0;

            // loop through the 2D histo list
            NameSet skipList;
            for ( TwoDMap::const_iterator twodim_iter = twodimPDFs_.begin();
                  twodim_iter != twodimPDFs_.end();
                  ++twodim_iter ) {
                // if the entry doesn't refer to this
                // species then skip on
                if ( specName != twodim_iter->first ) {
                    continue;
                }

                // retrieve the two variable names
                const TString& firstVarName  = twodim_iter->second.first;
                const TString& secondVarName = twodim_iter->second.second;
                if ( firstVarName != exclName && secondVarName != exclName ) {
                    // if neither name is the one being excluded then...
                    // add them both to the skip list
                    skipList.insert( firstVarName );
                    skipList.insert( secondVarName );
                    // and multiply the total by the combined PDF value
                    TString varName          = firstVarName + secondVarName;
                    pdfTot_[iEvt][specName] *= discPdf_[iEvt][specName][varName];
                }
            }

            // loop through all the variables
            for ( NameSet::const_iterator var_iter = variableNames_.begin();
                  var_iter != variableNames_.end();
                  ++var_iter ) {
                const TString& varName = ( *var_iter );
                // if the variable isn't the one being excluded
                if ( exclName != varName ) {
                    // and it's not involved in a 2D PDF
                    if ( skipList.find( varName ) == skipList.end() ) {
                        // multiply the total by its PDF value
                        pdfTot_[iEvt][specName] *= discPdf_[iEvt][specName][varName];
                    }
                }
            }
        }

        for ( NumbMap::const_iterator spec_iter = freeSpecies_.begin();
              spec_iter != freeSpecies_.end();
              ++spec_iter ) {
            const TString& specName = spec_iter->first;

            // if the signal is split we need to do a dedicated search
            // for sigTM and sigSCF, so skip the signal here
            if ( specName == "sig" && this->signalSplit() ) {
                needSignalSearch = kTRUE;
                continue;
            }

            pdfTot_[iEvt][specName] = 1.0;

            // loop through the 2D histo list
            NameSet skipList;
            for ( TwoDMap::const_iterator twodim_iter = twodimPDFs_.begin();
                  twodim_iter != twodimPDFs_.end();
                  ++twodim_iter ) {
                // if the entry doesn't refer to this
                // species then skip on
                if ( specName != twodim_iter->first ) {
                    continue;
                }

                // retrieve the two variable names
                const TString& firstVarName  = twodim_iter->second.first;
                const TString& secondVarName = twodim_iter->second.second;
                if ( firstVarName != exclName && secondVarName != exclName ) {
                    // if neither name is the one being excluded then...
                    // add them both to the skip list
                    skipList.insert( firstVarName );
                    skipList.insert( secondVarName );
                    // and multiply the total by the combined PDF value
                    TString varName          = firstVarName + secondVarName;
                    pdfTot_[iEvt][specName] *= discPdf_[iEvt][specName][varName];
                }
            }

            // loop through all the variables
            for ( NameSet::const_iterator var_iter = variableNames_.begin();
                  var_iter != variableNames_.end();
                  ++var_iter ) {
                const TString& varName = ( *var_iter );
                // if the variable isn't the one being excluded
                if ( exclName != varName ) {
                    // and it's not involved in a 2D PDF
                    if ( skipList.find( varName ) == skipList.end() ) {
                        // multiply the total by its PDF value
                        pdfTot_[iEvt][specName] *= discPdf_[iEvt][specName][varName];
                    }
                }
            }
        }

        if ( needSignalSearch ) {

            // loop through the 2D histo list
            TString specName( "sigTM" );
            Double_t tmPDFVal( 1.0 );
            NameSet skipList;
            for ( TwoDMap::const_iterator twodim_iter = twodimPDFs_.begin();
                  twodim_iter != twodimPDFs_.end();
                  ++twodim_iter ) {
                // if the entry doesn't refer to this
                // species then skip on
                if ( specName != twodim_iter->first ) {
                    continue;
                }

                // retrieve the two variable names
                const TString& firstVarName  = twodim_iter->second.first;
                const TString& secondVarName = twodim_iter->second.second;
                if ( firstVarName != exclName && secondVarName != exclName ) {
                    // if neither name is the one being excluded then...
                    // add them both to the skip list
                    skipList.insert( firstVarName );
                    skipList.insert( secondVarName );
                    // and multiply the total by the combined PDF value
                    TString varName  = firstVarName + secondVarName;
                    tmPDFVal        *= discPdf_[iEvt][specName][varName];
                }
            }

            // loop through all the variables
            for ( NameSet::const_iterator var_iter = variableNames_.begin();
                  var_iter != variableNames_.end();
                  ++var_iter ) {
                const TString& varName = ( *var_iter );
                // if the variable isn't the one being excluded
                if ( exclName != varName ) {
                    // and it's not involved in a 2D PDF
                    if ( skipList.find( varName ) == skipList.end() ) {
                        // multiply the total by its PDF value
                        tmPDFVal *= discPdf_[iEvt][specName][varName];
                    }
                }
            }

            tmPDFVal *= ( 1.0 - scfFrac_[iEvt] );

            // loop through the 2D histo list
            specName = "sigSCF";
            Double_t scfPDFVal( 1.0 );
            skipList.clear();
            for ( TwoDMap::const_iterator twodim_iter = twodimPDFs_.begin();
                  twodim_iter != twodimPDFs_.end();
                  ++twodim_iter ) {
                // if the entry doesn't refer to this
                // species then skip on
                if ( specName != twodim_iter->first ) {
                    continue;
                }

                // retrieve the two variable names
                const TString& firstVarName  = twodim_iter->second.first;
                const TString& secondVarName = twodim_iter->second.second;
                if ( firstVarName != exclName && secondVarName != exclName ) {
                    // if neither name is the one being excluded then...
                    // add them both to the skip list
                    skipList.insert( firstVarName );
                    skipList.insert( secondVarName );
                    // and multiply the total by the combined PDF value
                    TString varName  = firstVarName + secondVarName;
                    scfPDFVal       *= discPdf_[iEvt][specName][varName];
                }
            }

            // loop through all the variables
            for ( NameSet::const_iterator var_iter = variableNames_.begin();
                  var_iter != variableNames_.end();
                  ++var_iter ) {
                const TString& varName = ( *var_iter );
                // if the variable isn't the one being excluded
                if ( exclName != varName ) {
                    // and it's not involved in a 2D PDF
                    if ( skipList.find( varName ) == skipList.end() ) {
                        // multiply the total by its PDF value
                        scfPDFVal *= discPdf_[iEvt][specName][varName];
                    }
                }
            }

            if ( exclName == "DP" || ! this->scfDPSmear() ) {
                scfPDFVal *= scfFrac_[iEvt];
            }

            pdfTot_[iEvt]["sig"] = tmPDFVal + scfPDFVal;
        }
    }
}

void LauSPlot::calcCNCoeffs( const TString& exclName, const Double_t* covmat )
{
    // Computes the cN for the extended sPlots from the covariance matrix

    if ( nFixdSpecies_ <= 0 ) {
        return;
    }

    Int_t species_n( 0 );
    for ( NumbMap::const_iterator iter_n = freeSpecies_.begin(); iter_n != freeSpecies_.end();
          ++iter_n ) {
        Int_t species_j( 0 );
        const TString& specName = iter_n->first;
        Double_t value          = iter_n->second;
        cN_[exclName][specName] = value;
        for ( NumbMap::const_iterator iter_j = freeSpecies_.begin(); iter_j != freeSpecies_.end();
              ++iter_j ) {
            if ( covmat ) {
                cN_[exclName][specName] -= covmat[species_n * nSpecies_ + species_j];
            } else {
                cN_[exclName][specName] -= covMat_( species_n, species_j );
            }
            ++species_j;
        }
        ++species_n;
    }
}

void LauSPlot::calcSWeights( const TString& exclName, Double_t* covmat )
{
    // Computes the sWeights from the covariance matrix
    // NB for the extended sPlot the sum in the denominator is still over all species,
    // while that in the numerator is only over the free species.
    // Similarly the sWeights can only be calculated for the free species.

    Double_t numerator( 0.0 ), denominator( 0.0 );
    for ( Int_t iEvent = 0; iEvent < nEvents_; ++iEvent ) {
        denominator = 0.0;
        for ( NumbMap::const_iterator free_iter = freeSpecies_.begin();
              free_iter != freeSpecies_.end();
              ++free_iter ) {
            denominator += free_iter->second * pdfTot_[iEvent][free_iter->first];
        }
        for ( NumbMap::const_iterator fixd_iter = fixdSpecies_.begin();
              fixd_iter != fixdSpecies_.end();
              ++fixd_iter ) {
            denominator += fixd_iter->second * pdfTot_[iEvent][fixd_iter->first];
        }
        Int_t species_n( 0 );
        for ( NumbMap::const_iterator iter_n = freeSpecies_.begin(); iter_n != freeSpecies_.end();
              ++iter_n ) {
            numerator = 0.0;
            Int_t species_j( 0 );
            for ( NumbMap::const_iterator iter_j = freeSpecies_.begin(); iter_j != freeSpecies_.end();
                  ++iter_j ) {
                if ( covmat ) {
                    numerator += covmat[species_n * nSpecies_ + species_j] *
                                 pdfTot_[iEvent][iter_j->first];
                } else {
                    numerator += covMat_( species_n, species_j ) * pdfTot_[iEvent][iter_j->first];
                }
                ++species_j;
            }
            sWeights_[iEvent][exclName][iter_n->first] = numerator / denominator;
            ++species_n;
        }
    }
}

void LauSPlot::fillCNBranches()
{
    if ( ! cnTree_ ) {
        cerr << "ERROR in LauSPlot::fillCNBranches : Tree not created, cannot fill branches." << endl;
        return;
    } else if ( ! this->definedCNBranches() ) {
        this->defineCNBranches();
    }
    cnTree_->Fill();
}

void LauSPlot::fillSWeightBranches()
{
    if ( ! sweightTree_ ) {
        cerr << "ERROR in LauSPlot::fillSWeightBranches : Tree not created, cannot fill branches."
             << endl;
        return;
    } else if ( sWeights_.empty() ) {
        cerr << "ERROR in LauSPlot::fillSWeightBranches : No sWeights calculated, can't fill branches."
             << endl;
        return;
    } else if ( ! this->definedSWeightBranches() ) {
        this->copyEventWeights( 0 );
        this->defineSWeightBranches();
    }

    for ( Int_t iEvent = 0; iEvent < nEvents_; ++iEvent ) {
        this->copyEventWeights( iEvent );
        sweightTree_->Fill();
    }
}

void LauSPlot::copyEventWeights( Int_t iEvent )
{
    for ( std::map<TString, NumbMap>::const_iterator excl_iter = sWeights_[iEvent].begin();
          excl_iter != sWeights_[iEvent].end();
          ++excl_iter ) {
        const TString& exclName = excl_iter->first;
        const NumbMap& species  = excl_iter->second;
        for ( NumbMap::const_iterator spec_iter = species.begin(); spec_iter != species.end();
              ++spec_iter ) {
            const TString& specName              = spec_iter->first;
            sWeightsCurrent_[exclName][specName] = spec_iter->second;
        }
    }
}

void LauSPlot::writeOutResults()
{
    // write out the results

    // remove the transient objects that we don't want (re-)written to the file
    if ( eventList_ ) {
        delete eventList_;
        eventList_ = 0;
    }
    if ( inputTree_ ) {
        delete inputTree_;
        inputTree_ = 0;
    }

    // first add the input tree as a friend of the output tree
    this->addFriendTree();

    // then write everything to the file and clean up
    file_->cd();
    file_->Write();
    file_->Close();
    delete file_;
    file_ = 0;
}

void LauSPlot::addFriendTree()
{
    if ( ! sweightTree_ ) {
        cerr << "ERROR in LauSPlot::addFriendTree : Tree not created, cannot add friend." << endl;
        return;
    }
    sweightTree_->AddFriend( inputTreeName_, fileName_ );
}

void Yields( Int_t&, Double_t*, Double_t& f, Double_t* x, Int_t /*iflag*/ )
{
    // FCN-function for Minuit

    f = 0.0;

    TVirtualFitter* fitter = TVirtualFitter::GetFitter();
    LauSPlot* theModel     = dynamic_cast<LauSPlot*>( fitter->GetObjectFit() );
    const std::vector<LauSPlot::NumbMap>& pdfTot = theModel->totalPdf();

    Double_t ntot( 0.0 );
    for ( std::vector<LauSPlot::NumbMap>::const_iterator evt_iter = pdfTot.begin();
          evt_iter != pdfTot.end();
          ++evt_iter ) {    // loop over events
        const LauSPlot::NumbMap& species = ( *evt_iter );
        Int_t ispecies( 0 );
        Double_t lik( 0.0 );
        ntot = 0.0;
        for ( LauSPlot::NumbMap::const_iterator spec_iter = species.begin();
              spec_iter != species.end();
              ++spec_iter ) {    // loop over species
            lik  += x[ispecies] * spec_iter->second;
            ntot += x[ispecies];
            ++ispecies;
        }
        if ( lik < 0.0 ) {
            // make f the worst possible value to force MINUIT
            // out of this region of parameter space
            f = -DBL_MAX;
            break;
        } else {
            f += TMath::Log( lik );
        }
    }

    // extended likelihood
    f = ( ntot - f );
}
