
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

#include "TFile.h"
#include "TString.h"
#include "TTree.h"

#include <cstdlib>
#include <map>

/*! \file LauMergeDataFiles.hh
    \brief File containing declaration of LauMergeDataFiles class.
*/

/*! \class LauMergeDataFiles
    \brief Utility class to allow the merging of data files on a expt-by-expt basis

    The files are merged such that events for expt 0 from tree 1 will be followed
    by events for expt 0 from tree 2, then expt 1 from tree1, expt 1 from tree 2, etc.
*/

class LauMergeDataFiles {
  public:
    //! Constructor
    /*!
        \param [in] fileName1 name of first file to be merged
        \param [in] fileName2 name of second file to be merged
        \param [in] treeName name of the tree to read from the input files
    */
    LauMergeDataFiles( const TString& fileName1, const TString& fileName2, const TString& treeName );

    //! Destructor
    virtual ~LauMergeDataFiles();

    //! Do the merge
    /*!
        \param [in] fileName name of file to which to write the merged tree
    */
    void process( const TString& fileName );

  protected:
    //! Type to relate leaf names with their double-precision value
    typedef std::map<TString, Double_t> LeafDoubleMap;
    //! Type to relate leaf names with their integer value
    typedef std::map<TString, Int_t> LeafIntegerMap;
    //! Type to hold for each experiment the first and last entry numbers in a tree
    typedef std::map<Int_t, std::pair<Int_t, Int_t>> ExptsMap;

    //! Open the specified input files and check that the trees can be read
    void openInputFiles();
    //! Read the structure of the input trees, create appropriate storage and set the branch addresses
    void setupInputTrees();
    //! Create the structure of the output tree
    void setupOutputTree();
    //! Determine the experiments stored a given tree
    void findExperiments( TTree* tree, ExptsMap& exptsMap );
    //! Check that the experiments in each tree match
    Bool_t checkExperimentMaps() const;
    //! Read the entries for a given experiment from the given tree and store in the output tree
    void readExperiment( TTree* tree, const ExptsMap::const_iterator& exptsMap, Int_t offset );
    //! Write the output file
    void writeFile();

  private:
    //! Name of file 1
    TString fileName1_;
    //! Name of file 2
    TString fileName2_;
    //! Name of the tree
    TString treeName_;

    //! Input file 1
    TFile* inputFile1_;
    //! Input file 2
    TFile* inputFile2_;

    //! Input tree 1
    TTree* inputTree1_;
    //! Input tree 2
    TTree* inputTree2_;

    //! Output file
    TFile* outputFile_;
    //! Output tree
    TTree* outputTree_;

    // Tree variables
    //! Storage for the experiment index variable
    Int_t iExpt_;
    //! Storage for the event-within-experiment index variable
    Int_t iEvtWithinExpt_;

    //! Storage for double-precision leaves
    LeafDoubleMap doubleVars_;
    //! Storage for integer leaves
    LeafIntegerMap integerVars_;

    //! Experiment -> first and last tree entry for tree 1
    ExptsMap tree1Expts_;
    //! Experiment -> first and last tree entry for tree 2
    ExptsMap tree2Expts_;

    ClassDef( LauMergeDataFiles, 0 )
};
