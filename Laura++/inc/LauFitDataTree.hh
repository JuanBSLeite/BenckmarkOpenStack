
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

/*! \file LauFitDataTree.hh
    \brief File containing declaration of LauFitDataTree class.
*/

/*! \class LauFitDataTree
    \brief Class to store the input fit variables.

    Events are loaded from a tree and fake events may be added manually.
*/

#ifndef LAU_FIT_DATA_TREE
#define LAU_FIT_DATA_TREE

#include "TEventList.h"
#include "TTree.h"

#include <map>
#include <vector>

class TFile;
class TLeaf;

//! Type for holding event data
typedef std::map<TString, Double_t> LauFitData;

class LauFitDataTree {

  public:
    //! Constructor
    /*!
        \param [in] rootFileName the name of the file containing the data
        \param [in] rootTreeName the name of the tree containing the data
    */
    LauFitDataTree( const TString& rootFileName, const TString& rootTreeName );

    //! Destructor
    virtual ~LauFitDataTree();

    //! Retrieve the file name
    /*!
        \return the file name
    */
    const TString& fileName() const { return rootFileName_; }

    //! Retrieve the tree name
    /*!
        \return the tree name
    */
    const TString& treeName() const { return rootTreeName_; }

    //! Find all of the branches in the tree
    /*!
        \returns true if sucessful, otherwise returns false and prints an error message
    */
    Bool_t findBranches();

    //! Read all events from the tree
    void readAllData();

    //! Read events only for the given experiment
    /*!
        \param [in] iExpt the experiment to read
    */
    void readExperimentData( UInt_t iExpt );

    //! Add fake events to the data
    /*!
        Used, for example, for storing the information on the histogram bin centres for calculating the DP smearing of mis-reconstructed events, see, for example, LauSimpleFitModel::splitSignalComponent and LauScfMap

        \param [in] xCoords the values of m13 squared for the fake events
        \param [in] yCoords the values of m23 squared for the fake events
    */
    void appendFakePoints( const std::vector<Double_t>& xCoords,
                           const std::vector<Double_t>& yCoords );

    //! Obtain the number of branches in the tree
    /*!
        \return the number of branches
    */
    UInt_t nBranches() const;

    //! Retrieve the number of events in the tree
    /*!
        \return the number of events in the tree
    */
    UInt_t nTreeEvents() const
    {
        return rootTree_ ? static_cast<UInt_t>( rootTree_->GetEntries() ) : 0;
    }

    //! Retrieve the number of events
    /*!
        \return the number of events in the current experiment if one is selected, otherwise the total number of events in the tree
    */
    UInt_t nEvents() const
    {
        return eventList_ ? static_cast<UInt_t>( eventList_->GetN() ) : this->nTreeEvents();
    }

    //! Retrieve the number of fake events
    /*!
        \return the number of fake events
    */
    UInt_t nFakeEvents() const { return fakeEvents_.size(); }

    //! Retrieve the total number of events
    /*!
        \return the total number of events
    */
    UInt_t nTotalEvents() const { return this->nTreeEvents() + this->nFakeEvents(); }

    //! Check if the named branch is stored
    /*!
        \param [in] name the name of the branch
        \return true if the branch is stored, otherwise return false
    */
    Bool_t haveBranch( const TString& name ) const;

    //! Retrieve the data for a given event
    /*!
        \param [in] iEvt the index of the event
        \return the event data
    */
    const LauFitData& getData( UInt_t iEvt ) const;

    //! Disable all branches
    void disableAllBranches() const;

    //! Enable all branches
    void enableAllBranches() const;

    //! Enable the named branch
    /*!
        \param [in] name the name of the branch
    */
    void enableBranch( const TString& name ) const;

    //! Disable the named branch
    /*!
        \param [in] name the name of the branch
    */
    void disableBranch( const TString& name ) const;

  protected:
    //! Open the file and tree
    void openFileAndTree();

    //! Load events from the tree
    /*!
        If an event list is present, load only those events, otherwise load all events
    */
    void loadData();

  private:
    //! Copy constructor (not implemented)
    LauFitDataTree( const LauFitDataTree& rhs );

    //! Copy assignment operator (not implemented)
    LauFitDataTree& operator=( const LauFitDataTree& rhs );

    //! The type used to map the leaf names to the vector indices
    typedef std::map<TString, UInt_t> LauNameIndexMap;

    //! The type used to contain the data for each event
    typedef std::vector<Double_t> LauEventData;

    //! The type used to hold the leaves
    typedef std::vector<TLeaf*> LauLeafList;

    //! The name of the file containing the data
    TString rootFileName_;

    //! The name of the tree ocntaining the data
    TString rootTreeName_;

    //! The file containing the data
    TFile* rootFile_;

    //! The tree containing the data
    TTree* rootTree_;

    //! A list of the events in the current experiment
    TEventList* eventList_;

    //! Stores the mapping from the leaf names to the vector indices
    LauNameIndexMap leafNames_;

    //! Stores the current event
    mutable LauEventData eventData_;

    //! Stores the current event (for external use)
    mutable LauFitData eventDataOut_;

    //! The leaf objects
    LauLeafList leaves_;

    //! The events read from the tree
    std::vector<LauEventData> treeEvents_;

    //! The fake events, which are not from the tree
    std::vector<LauEventData> fakeEvents_;

    ClassDef( LauFitDataTree, 0 )
};

#endif
