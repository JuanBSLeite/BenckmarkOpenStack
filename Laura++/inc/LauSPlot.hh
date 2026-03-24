
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

/*! \file LauSPlot.hh
    \brief File containing declaration of LauSPlot class.
*/

/*! \class LauSPlot
    \brief Class for defining the SPlot technique

    Class for defining the SPlot technique based on TSplot from ROOT by the following authors:
    Muriel Pivk, Anna Kreshuk (10/2005).
    (Original copyright notice below)

    Code extended to deal with the following two extra scenarios:
    - Extended sPlots (see appendix B of sPlot paper)
    - 2D PDFs, i.e. one PDF for 2 variables

    When performing a multidimensional fit, the sWeights are calculated
    excluding each dimension in turn and excluding none of the variables.  This
    allows sPlots of the fit variables to be made by using the weight
    calculated when that variable is excluded, while variables not in the fit
    can be plotted from the complete information.
*/

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 ROOT Foundation,  CERN/PH-SFT                   *
 *                                                                    *
 **********************************************************************/

#ifndef LAU_SPLOT
#define LAU_SPLOT

#include "TMatrixD.h"
#include "TString.h"

#include <map>
#include <set>
#include <vector>

class TEventList;
class TLeaf;
class TFile;
class TTree;

class LauSPlot : public TObject {

  public:
    //! Type to store names, e.g. of the discriminating/control variables
    typedef std::set<TString> NameSet;

    //! Type to associate a category name with a double precision number, e.g. a yield or PDF value for a given species
    typedef std::map<TString, Double_t> NumbMap;

    //! Type to associate a variable name with the leaf of the input tree
    typedef std::map<TString, TLeaf*> LeafMap;

    //! Type to associate the name of the species that have 2D PDFs with the names of the two variables involved in each such PDF
    typedef std::multimap<TString, std::pair<TString, TString>> TwoDMap;

    //! Constructor
    /*!
        \param [in] fileName the name of the data file (containing the input tree)
        \param [in] treeName the name of the input tree (containing the per-event llhds)
        \param [in] firstExpt the first experiment
        \param [in] nExpt the number of experiments
        \param [in] variableNames the names of the discriminating variables
        \param [in] freeSpecies the species that are free to float in the fit
        \param [in] fixdSpecies the species that are fixed in the fit
        \param [in] twodimPDFs the species that have 2D PDFs and the name of the variables involved in each such PDF
        \param [in] sigSplit boolean flag to check whether the signal is split into Truth Matched and Self Cross Feed
        \param [in] scfDPSmeared boolean flag to check whether the SCF is smeared in the DP
    */
    LauSPlot( const TString& fileName,
              const TString& treeName,
              Int_t firstExpt,
              Int_t nExpt,
              const NameSet& variableNames,
              const LauSPlot::NumbMap& freeSpecies,
              const LauSPlot::NumbMap& fixdSpecies,
              const TwoDMap& twodimPDFs,
              Bool_t sigSplit     = kFALSE,
              Bool_t scfDPSmeared = kFALSE );

    //! Destructor
    virtual ~LauSPlot();

    //! Method to calculate the sWeights and cN coeffs.
    /*!
        \param [in] option control print level (default set to no print out)\n
        "Q" - no print out (default)\n
        "V" - prints the estimated # of events in species\n
        "VV" - as "V" + the MINUIT printing + sums of weights for control
    */
    void runCalculations( const TString& option = "q" );

    //! Save the sWeight results as a friend tree to the input tree (in the same file)
    void writeOutResults();

    //! Access the per-event total PDF values for each species
    /*!
        \return the per-event total PDF values
    */
    const std::vector<LauSPlot::NumbMap>& totalPdf() const { return pdfTot_; }

  protected:
    //! Check whether the input tree has been successfully read
    /*!
        \return true/false whether the task have been successfully performed
    */
    Bool_t readInput() const { return readInput_; }

    //! Set that the input tree has been successfully read
    /*!
        \param [in] ok set that the input tree has been succeessfully read
    */
    void readInput( Bool_t ok ) { readInput_ = ok; }

    //! Check whether the signal is split into Truth Matched and Self Cross Feed
    /*!
        \return true/false whether the signal is split into TM and SCF
    */
    Bool_t signalSplit() const { return signalSplit_; }

    //! Check whether the Self Cross Feed is smeared in the DP
    /*!
        \return true/false whether the Self Cross Feed is smeared in the DP
    */
    Bool_t scfDPSmear() const { return scfDPSmear_; }

    //! Check whether the cN branches have been already created
    /*!
        \return true/false whether cN branches have been already created
    */
    Bool_t definedCNBranches() const { return definedCNBranches_; }

    //! Set that the cN branches have been already defined
    /*!
        \param [in] defined set that the cN branches have been already defined
    */
    void definedCNBranches( Bool_t defined ) { definedCNBranches_ = defined; }

    //! Check whether the sWeights branches have been already created
    /*!
        \return true/false whether sWeights branches have been already created
    */
    Bool_t definedSWeightBranches() const { return definedSWeightBranches_; }

    //! Set that the sWeights branches have been already defined
    /*!
        \param [in] defined set that the sWeights branches have been already defined
    */
    void definedSWeightBranches( Bool_t defined ) { definedSWeightBranches_ = defined; }

    //! Method to open the file in "update" mode and grab the input tree for reading
    void openInputFileAndTree();

    //! Read the leaf structure from the tree and check the status of the read (calls LauSPlot::readInputLeaves and LauSPlot::checkLeaves)
    void readInputInfo();

    //! Read the leaf structure from the tree and setup the leaf map
    Bool_t readInputLeaves();

    //! Check whether the leaf structure makes sense given the PDFs we are expecting
    Bool_t checkLeaves() const;

    //! Create (if not already done) the tree for storing the cN coeffs
    void createCNTree();
    //! Create the branches for each cN coefficient
    void defineCNBranches();

    //! Create (if not already done) the tree for storing the sWeights
    void createSWeightTree();
    //! Create the branches to store the sWeights
    void defineSWeightBranches();

    //! Set the event list to contain only events from the given experiment
    /*!
        \param [in] iExpt the required experiment number
    */
    void setExperiment( Int_t iExpt );

    //! Reads the values of each PDF likelihood for every event in the experiment
    void readExpt();

    //! Make sure that we're using Minuit
    void checkFitter() const;

    //! Initialise Minuit, set the verbosity
    /*!
        \param [in] opt option to set the print level (opt = "Q", "V" or "W" == "VV")
    */
    void initialiseFitter( const TString& opt );

    //! Add the species yields as fit parameters and fix them as appropriate
    void setFitParameters() const;

    //! Perform the minimisation wrt the yields alone
    void runFit();

    //! Update the yields with the newly fitted values and print them (unless print option is "Q").
    /*!
        \param [in] opt option to set the print level (opt = "Q", "V" or "W")
    */
    void retrieveFittedParameters( const TString& opt );

    //! Print the supplied covariance matrix or, if pointer is null, the one previously calculated
    /*!
        \param [in] covmat 2D array containing the covariace matrix elements
    */
    void printCovMatrixElements( const Double_t* covmat = 0 ) const;

    //! Print the sum of sWeights for all species
    /*!
        \param [in] exclName the name of variable excluded (or "none")
    */
    void printSumOfWeights( const TString& exclName ) const;

    // Calculate the covariance matrix from the various PDFs
    /*!
        \return true if calculation succeeded or false in case the covariance matrix can not be inverted
    */
    Bool_t calcCovMatrix();

    //! Calculate the total likelihood for each species by multiply together all the PDFs for that species
    /*!
        \param [in] exclName the name of excluded variable (or "none")
    */
    void calcTotPDFValues( const TString& exclName );

    //! Computes the cN for the extended sPlots from the covariance matrix
    /*!
        \param [in] exclName the name of excluded variable (or "none")
        \param [in] covmat 2D array containing the covariace matrix elements
    */
    void calcCNCoeffs( const TString& exclName, const Double_t* covmat = 0 );

    //! Computes the sWeights from the PDFs and covariance matrix
    /*!
        \param [in] exclName the name of excluded variable (or "none")
        \param [in] covmat 2D array containing the covariace matrix elements
    */
    void calcSWeights( const TString& exclName, Double_t* covmat = 0 );

    //! Copy the sWeight of a given event into LauSPlot::sWeightsCurrent_, from which they can be stored in the output tree
    /*!
        \param [in] iEvent the requested event
    */
    void copyEventWeights( Int_t iEvent );

    //! Fill the cN branches
    void fillCNBranches();

    //! Fill the sWeights branches
    void fillSWeightBranches();

    //! Add the sWeightTree as a friend tree of the input tree
    void addFriendTree();

  private:
    //! Copy constructor (not implemented)
    LauSPlot( const LauSPlot& rhs );

    //! Copy assignment operator (not implemented)
    LauSPlot& operator=( const LauSPlot& rhs );

    //! The name of the data file
    TString fileName_;
    //! The name of the input tree (containing the per-event llhds)
    TString inputTreeName_;
    //! The name of the cn tree (containing the cN coefficients)
    TString cnTreeName_;
    //! The name of the sweight tree (containing the sWeights)
    TString sweightTreeName_;
    //! Pointer to the data file object
    TFile* file_;
    //! Pointer to the input tree
    TTree* inputTree_;
    //! Pointer to the output tree containing the cN coefficients
    TTree* cnTree_;
    //! Pointer to the output tree containing the sWeights
    TTree* sweightTree_;

    //! Pointer to an event list, that is used to loop through the experiments
    TEventList* eventList_;

    //! Collection to hold pointers to the leaves of the input tree
    LeafMap leaves_;

    //! The names of the discriminating variables
    NameSet variableNames_;
    //! The names and estimated yields of the free species
    NumbMap freeSpecies_;
    //! The names and estimated yields of the fixed species
    NumbMap fixdSpecies_;
    //! The names and estimated yields of the free species - need to keep the original values
    NumbMap origFreeSpecies_;
    //! The names and estimated yields of the fixed species - need to keep the original values
    NumbMap origFixdSpecies_;
    //! The names of the species that have 2D PDFs and the names of the variables involved
    TwoDMap twodimPDFs_;

    //! Is the signal split into TM and SCF?
    const Bool_t signalSplit_;
    //! If so then is the SCF smeared in the DP?
    const Bool_t scfDPSmear_;

    //! Flag whether the input tree has been successfully read
    Bool_t readInput_;
    //! Flag whether the cN branches have already been created
    Bool_t definedCNBranches_;
    //! Flag whether the sWeights branches have already been created
    Bool_t definedSWeightBranches_;

    //! First experiment
    Int_t firstExpt_;
    //! Number of experiments
    Int_t nExpt_;
    //! The current experiment
    Int_t iExpt_;
    //! Number of events in current experiment
    Int_t nEvents_;
    //! Number of discriminating variables
    Int_t nDiscVars_;
    //! Number of species free to float in the fit
    Int_t nFreeSpecies_;
    //! Number of species fixed in the fit
    Int_t nFixdSpecies_;
    //! Total number of species (free + fixed)
    Int_t nSpecies_;

    //! The per-event values of the total PDF for each species
    std::vector<NumbMap> pdfTot_;
    //! The per-event values of the PDFs for each species for each disc variable
    std::vector<std::map<TString, NumbMap>> discPdf_;
    //! The per-event values of the SCF fraction
    std::vector<Double_t> scfFrac_;

    //! The calculated covariance matrix
    TMatrixD covMat_;

    //! The per-event values of the computed sWeights (for each species and for each combination of excluded vars)
    std::vector<std::map<TString, NumbMap>> sWeights_;
    //! The current-event values of the computed sWeights
    std::map<TString, NumbMap> sWeightsCurrent_;
    //! The extended sPlot coefficients (for each species and for each combination of excluded vars)
    std::map<TString, NumbMap> cN_;

    ClassDef( LauSPlot, 0 )
};

#endif
