
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

/*! \file LauGenNtuple.hh
    \brief File containing declaration of LauGenNtuple class.
*/

/*! \class LauGenNtuple
    \brief Class to store the results from the toy MC generation into an ntuple

    Class to store the results from the toy MC generation into an ntuple
*/

#ifndef LAU_GEN_NTUPLE
#define LAU_GEN_NTUPLE

#include "TString.h"

#include <map>

class TFile;
class TTree;

class LauGenNtuple {

  public:
    //! Constructor
    /*!
        \param [in] rootFileName the name for the ntuple
        \param [in] rootTreeName the name for the tree in the ntuple
    */
    LauGenNtuple( const TString& rootFileName, const TString& rootTreeName );

    //! Destructor
    virtual ~LauGenNtuple();

    //! Ntuple file name
    /*!
        \return the name of the ntuple
    */
    const TString& fileName() const { return rootFileName_; }

    //! Ntuple tree name
    /*!
        \return the name of the tree
    */
    const TString& treeName() const { return rootTreeName_; }

    //! Add integer branch to tree
    /*!
        \param [in] name the name of the branch
    */
    void addIntegerBranch( const TString& name );

    //! Add double branch to tree
    /*!
        \param [in] name the name of the branch
    */
    void addDoubleBranch( const TString& name );

    //! Set value of an integer branch
    /*!
        \param [in] name the name of the branch
        \param [in] value the value to set the branch
    */
    void setIntegerBranchValue( const TString& name, Int_t value );

    //! Set value of a double branch
    /*!
        \param [in] name the name of the branch
        \param [in] value the value to set the branch
    */
    void setDoubleBranchValue( const TString& name, Double_t value );

    //! Get value of an integer branch
    /*!
        \param [in] name the name of the branch
        \return the value of the branch
    */
    Int_t getIntegerBranchValue( const TString& name ) const;

    //! Get value of a double branch
    /*!
        \param [in] name the name of the branch
        \return the value of the branch
    */
    Double_t getDoubleBranchValue( const TString& name ) const;

    //! Fill branches in the ntuple
    void fillBranches();

    //! Delete and recreate tree
    void deleteAndRecreateTree();

    //! Create an index table using leaves of the tree
    /*!
        \param [in] majorName the majorname
        \param [in] minorName the minorname
        \return the number of entries in the index (< 0 indicates failure)
    */
    Int_t buildIndex( const TString& majorName, const TString& minorName = "0" );

    //! Write out the results from the generation
    void writeOutGenResults();

    //! Add a friend tree
    /*!
        \param [in] rootFileName the name of the root file
        \param [in] rootTreeName the name of the root tree
    */
    void addFriendTree( const TString& rootFileName, const TString& rootTreeName );

  protected:
    //! Create ntuple file and the tree
    void createFileAndTree();

    //! Define branches of the tree
    void defineBranches();

    //! Flags whether branches have been defined
    /*!
        \return boolean to determine whether branches have been defined
    */
    Bool_t definedBranches() const { return definedBranches_; }

    //! Flags whether branches have been defined
    /*!
        \param [in] defined boolean to determine whether branches have been defined
    */
    void definedBranches( Bool_t defined ) { definedBranches_ = defined; }

  private:
    //! Copy constructor (not implemented)
    LauGenNtuple( const LauGenNtuple& rhs );

    //! Copy assignment operator (not implemented)
    LauGenNtuple& operator=( const LauGenNtuple& rhs );

    //! Name of root file
    TString rootFileName_;
    //! Name of root tree
    TString rootTreeName_;
    //! Root file
    TFile* rootFile_;
    //! Root tree
    TTree* rootTree_;

    //! Flags whether branches are defined
    Bool_t definedBranches_;

    //! Type to hold integer variables
    typedef std::map<TString, Int_t> IntVarMap;
    //! Type to hold double precision floating point variables
    typedef std::map<TString, Double_t> DoubleVarMap;

    //! Integer variables
    IntVarMap intVars_;
    //! Double variables
    DoubleVarMap doubleVars_;

    ClassDef( LauGenNtuple, 0 )    // Generated toyMC ntuple
};

#endif
