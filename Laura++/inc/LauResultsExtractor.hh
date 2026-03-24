
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

#include "TString.h"

#include <map>
#include <vector>

class TChain;
class TFile;
class TH1;
class TTree;

/*! \file LauResultsExtractor.hh
    \brief File containing declaration of LauResultsExtractor class.
*/

/*! \class LauResultsExtractor
    \brief Utility class to allow the extraction of the best fit from a series of fits to a given data sample

    A utility class to allow the extraction of the best fit from a series of fits to a given data sample.

    When fitting amplitude models, the likelihood parameter space is highly complex.
    Hence the fitter can often wander into local minima.
    To mitigate this effect a data sample can be fitted many times with randomised starting values of the isobar parameters.
    It is then necessary to determine which of these fits gives the best solution, i.e. the minimum NLL.
    This class performs this task, reading in a series of fits for each data sample and writing out a single file that contains the results of the best fit for each data sample.
*/

class LauResultsExtractor {
  public:
    //! Constructor
    /*!
        \param [in] inputFileName name of text file containing the input ROOT files
        \param [in] outputFileName name of the file to which the best fit info should be written
        \param [in] treeName name of the tree to read from the input files
    */
    LauResultsExtractor( const TString& inputFileName,
                         const TString& outputFileName,
                         const TString& treeName );

    //! Destructor
    virtual ~LauResultsExtractor();

    //! Run the calculations
    /*!
        \param [in] numExpts the number of experiments to process
    */
    void process( const Int_t numExpts );

  protected:
    //! Create storage for leaves and call SetBranchAddress for each
    void setupInputTree();
    //! Create branches in the output tree
    void setupOutputTree( TTree* tree );
    //! Toggle branch status of input tree branches (except iExpt, fitStatus and NLL)
    void setInputTreeBranchStatus( const Bool_t status );
    //! Clear all information
    void clearMaps();
    //! Write the output file
    void writeFile();

  private:
    //! Name of text file containing list of input files
    TString inputFileName_;
    //! Name of output ROOT file
    TString outputFileName_;
    //! Name of tree in input ROOT files
    TString treeName_;

    //! Chain of inputs
    TChain* inputTree_;

    //! Output file
    TFile* outputFile_;
    //! Output tree
    TTree* outputTree_;

    //! Number of entries in the input chain
    Int_t nEntries_;

    // Tree variables
    //! Storage for experiment ID variable
    Int_t iExpt_;
    //! Storage for fit status variable
    Int_t fitStatus_;
    //! Storage for NLL variable
    Double_t NLL_;
    //! Storage for EDM variable
    Double_t EDM_;
    //! Storage for other input variables
    std::map<TString, Double_t> otherVars_;

    //! Best NLL and corresponding tree entries for each experiment
    std::map<Int_t, std::pair<Double_t, Int_t>> bestNLL_;
    //! Worst NLL and corresponding tree entries for each experiment
    std::map<Int_t, std::pair<Double_t, Int_t>> worstNLL_;
    //! All NLL values for each experiment
    std::map<Int_t, std::vector<Double_t>> allNLLs_;

    //! Histograms of the NLL values for each experiment
    std::map<Int_t, TH1*> nllHistos_;

    ClassDef( LauResultsExtractor, 0 )
};
