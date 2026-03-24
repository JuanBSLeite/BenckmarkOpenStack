
/*
Copyright 2007 University of Warwick

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

/*! \file LauEmbeddedData.hh
    \brief File containing declaration of LauEmbeddedData class.
*/

/*! \class LauEmbeddedData
    \brief Class to store the data for embedding in toy experiments

    Class to store the data from external data sources, such as full simulation, for embedding in toy experiments
*/

#ifndef LAU_EMBEDDED_DATA
#define LAU_EMBEDDED_DATA

#include "LauFitDataTree.hh"

#include <set>
#include <vector>

class LauKinematics;
class LauIsobarDynamics;

class LauEmbeddedData {

  public:
    //! Constructor
    /*!
        \param [in] fileName the ntuple file name
        \param [in] treeName the tree name
        \param [in] allowReuseOfEvents flag to indicated whether events are allowed to be sampled more than once
    */
    LauEmbeddedData( const TString& fileName, const TString& treeName, Bool_t allowReuseOfEvents );

    //! Destructor
    virtual ~LauEmbeddedData();

    //! Find and read the branches in data tree
    /*!
        \return success/failure flag
    */
    Bool_t findBranches();

    //! Retrieve the number of branches
    /*!
        \return the number of branches
    */
    UInt_t nBranches() const { return theDataTree_ ? theDataTree_->nBranches() : 0; }

    //! Retrieve the number of events
    /*!
        \return the number of events
    */
    UInt_t nEvents() const { return theDataTree_ ? theDataTree_->nEvents() : 0; }

    //! Retrieve the number of events that have already been sampled
    /*!
        \return the number of used events
    */
    UInt_t nUsedEvents() const { return usedEvents_.size(); }

    //! Boolean to determine whether branch exists
    /*!
        \param [in] name the branch name
    */
    Bool_t haveBranch( const TString& name ) const
    {
        return theDataTree_ ? theDataTree_->haveBranch( name ) : kFALSE;
    }

    //! Retrieve an event from the data sample
    /*!
        \param [in] kinematics object in which to store the kinematics of the Dalitz plot
    */
    void getEmbeddedEvent( LauKinematics* kinematics );

    //! Retrieve an event from the data sample, applying an accept/reject based on the given DP model
    /*!
        \param [in] dynamics the amplitude model
    */
    Bool_t getReweightedEvent( LauIsobarDynamics* dynamics );

    //! Get the value of a specified branch
    /*!
        \param [in] name the branch name
        \return the branch value
    */
    Double_t getValue( const TString& name ) const;

    //! Get values of specified branches
    /*!
        \param [in] names the branch names
        \return the branch values
    */
    LauFitData getValues( const std::vector<TString>& names ) const;

    //! Clear the list of used events
    void clearUsedList() { usedEvents_.clear(); }

  protected:
    //! Boolean determining whether events should be reused
    /*!
        \return flag determining whether events should be reused
    */
    Bool_t reuseEvents() const { return reuseEvents_; }

  private:
    //! Copy constructor (not implemented)
    LauEmbeddedData( const LauEmbeddedData& rhs );

    //! Copy assignment operator (not implemented)
    LauEmbeddedData& operator=( const LauEmbeddedData& rhs );

    //! The structure containing the data
    LauFitDataTree* theDataTree_;
    //! The data for the currently retrieved event
    LauFitData theData_;
    //! Flag whether events can be reused
    Bool_t reuseEvents_;
    //! Used events
    std::set<UInt_t> usedEvents_;

    ClassDef( LauEmbeddedData, 0 )    // a non-persistent bare-bones complex class
};

#endif
