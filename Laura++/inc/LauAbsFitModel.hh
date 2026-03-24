
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

/*! \file LauAbsFitModel.hh
    \brief File containing declaration of LauAbsFitModel class.
*/

/*! \class LauAbsFitModel
    \brief Abstract interface to the fitting and toy MC model

    Abstract interface to the fitting and toy MC model

    Any class inheriting from this must implement the following functions:
    - cacheInputFitVars()
    - checkInitFitParams()
    - finaliseFitResults()
    - fixdSpeciesNames()
    - freeSpeciesNames()
    - genExpt()
    - getEventSum()
    - getTotEvtLikelihood()
    - initialise()
    - initialiseDPModels()
    - propagateParUpdates()
    - recalculateNormalisation()
    - scfDPSmear()
    - setAmpCoeffSet()
    - setNBkgndEvents()
    - setNSigEvents()
    - setupBkgndVectors()
    - setupGenNtupleBranches()
    - setupSPlotNtupleBranches()
    - splitSignal()
    - storePerEvtLlhds()
    - twodimPDFs()
    - updateCoeffs()
    - variableNames()
    - weightEvents()
    - writeOutTable()
*/

#ifndef LAU_ABS_FIT_MODEL
#define LAU_ABS_FIT_MODEL

#include "LauComplex.hh"
#include "LauFitObject.hh"
#include "LauFormulaPar.hh"
#include "LauSimFitTask.hh"

#include "TMatrixDfwd.h"
#include "TStopwatch.h"
#include "TString.h"

#include <iosfwd>
#include <set>
#include <vector>
// LauSPlot included to get LauSPlot::NameSet typedef
#include "LauSPlot.hh"

class LauAbsCoeffSet;
class LauAbsPdf;
class LauFitDataTree;
class LauGenNtuple;
class LauAbsRValue;
class LauParameter;

class LauAbsFitModel : public LauSimFitTask {

  public:
    //! Constructor
    LauAbsFitModel();

    //! Destructor
    virtual ~LauAbsFitModel();

    //! Is the Dalitz plot term in the likelihood
    Bool_t useDP() const { return usingDP_; }

    //! Switch on/off the Dalitz plot term in the Likelihood (allows fits to other quantities, e.g. B mass)
    /*!
        \param [in] usingDP the boolean flag
    */
    void useDP( Bool_t usingDP ) { usingDP_ = usingDP; }

    //! Return the flag to store the status of using an sFit or not
    Bool_t doSFit() const { return doSFit_; }

    //! Do an sFit (use sWeights to isolate signal decays rather than using background histograms)
    /*!
        \param [in] sWeightBranchName name of the branch of the tree containing the sWeights
        \param [in] scaleFactor scaling factor to get the uncertainties correct
    */
    void doSFit( const TString& sWeightBranchName, Double_t scaleFactor = 1.0 );

    //! Determine whether an extended maximum likelihood fit it being performed
    Bool_t doEMLFit() const { return emlFit_; }

    //! Choice to perform an extended maximum likelihood fit
    /*!
        \param [in] emlFit boolean specifying whether or not to perform the EML
    */
    void doEMLFit( Bool_t emlFit ) { emlFit_ = emlFit; }

    //! Determine whether Poisson smearing is enabled for the toy MC generation
    Bool_t doPoissonSmearing() const { return poissonSmear_; }

    //! Turn Poisson smearing (for the toy MC generation) on or off
    /*!
        \param [in] poissonSmear boolean specifying whether or not to do Poisson smearing
    */
    void doPoissonSmearing( Bool_t poissonSmear ) { poissonSmear_ = poissonSmear; }

    //! Determine whether embedding of events is enabled in the generation
    Bool_t enableEmbedding() const { return enableEmbedding_; }

    //! Turn on or off embedding of events in the generation
    /*!
        \param [in] enable boolean specifying whether to embed events
    */
    void enableEmbedding( Bool_t enable ) { enableEmbedding_ = enable; }

    //! Determine whether writing out of the latex table is enabled
    Bool_t writeLatexTable() const { return writeLatexTable_; }

    //! Turn on or off the writing out of the latex table
    /*!
        \param [in] writeTable boolean specifying whether or not the latex table should be written out
    */
    void writeLatexTable( Bool_t writeTable ) { writeLatexTable_ = writeTable; }

    //! save files containing graphs of the resonance's PDFs
    Bool_t saveFilePDF() const { return savePDF_; }

    //! Turn on or off the save of files containing graphs of the resonance's PDFs
    /*!
        \param [in] savePDF boolean specifying whether or not the save of files containing graphs of the resonance's PDFs
    */
    void saveFilePDF( Bool_t savePDF ) { savePDF_ = savePDF; }

    //! Set up the sPlot ntuple
    /*!
        \param [in] fileName the sPlot file name
        \param [in] treeName the sPlot tree name
        \param [in] storeDPEfficiency whether or not to store the efficiency information too
        \param [in] verbosity define the level of output
    */
    void writeSPlotData( const TString& fileName,
                         const TString& treeName,
                         Bool_t storeDPEfficiency,
                         const TString& verbosity = "q" );

    //! Determine whether the sPlot data is to be written out
    Bool_t writeSPlotData() const { return writeSPlotData_; }

    //! Determine whether the efficiency information should be stored in the sPlot ntuple
    Bool_t storeDPEff() const { return storeDPEff_; }

    //! Determine whether the initial values of the fit parameters, in particular the isobar coefficient parameters, are to be randomised
    Bool_t useRandomInitFitPars() const { return randomFit_; }

    //! Randomise the initial values of the fit parameters, in particular the isobar coefficient parameters
    void useRandomInitFitPars( Bool_t boolean ) { randomFit_ = boolean; }

    //! Setup the background class names
    /*!
        \param [in] names a vector of all the background names
    */
    virtual void setBkgndClassNames( const std::vector<TString>& names );

    //! Returns the number of background classes
    inline UInt_t nBkgndClasses() const { return bkgndClassNames_.size(); }

    //! Set the number of signal events
    /*!
        \param [in] nSigEvents contains the signal yield and option to fix it
    */
    virtual void setNSigEvents( LauParameter* nSigEvents ) = 0;

    //! Set the number of background events
    /*!
        The name of the parameter must be that of the corresponding background category (so that it can be correctly assigned)

        \param [in] nBkgndEvents contains the name, yield and option to fix the yield of the background

    */
    virtual void setNBkgndEvents( LauAbsRValue* nBkgndEvents ) = 0;

    //! Set the DP amplitude coefficients
    /*!
        The name of the coeffSet must match the name of one of the resonances in the DP model.

        The supplied order of coefficients will be rearranged to match the order in which the resonances are stored in the dynamics, see LauIsobarDynamics::addResonance.

        \param [in] coeffSet the set of coefficients
    */
    virtual void setAmpCoeffSet( LauAbsCoeffSet* coeffSet ) = 0;

    //! Specify that a toy MC sample should be created for a successful fit to an experiment
    /*!
        Generation uses the fitted parameters so that the user can compare the fit to the data

        \param [in] toyMCScale the scale factor to get the number of events to generate
        \param [in] mcFileName the file name where the toy sample will be stored
        \param [in] tableFileName name of the output tex file
        \param [in] poissonSmearing turn smearing on or off
    */
    void compareFitData( UInt_t toyMCScale            = 10,
                         const TString& mcFileName    = "fitToyMC.root",
                         const TString& tableFileName = "fitToyMCTable.tex",
                         Bool_t poissonSmearing       = kTRUE );

    //! Start the toy generation / fitting
    /*!
        \param [in] applicationCode specifies what to do, perform a fit ("fit") or generate toy MC ("gen")
        \param [in] dataFileName the name of the input data file
        \param [in] dataTreeName the name of the tree containing the data
        \param [in] histFileName the file name for the output histograms
        \param [in] tableFileName the file name for the latex output file
    */
    void run( const TString& applicationCode,
              const TString& dataFileName,
              const TString& dataTreeName,
              const TString& histFileName,
              const TString& tableFileName = "" );

    //! This function sets the parameter values from Minuit
    /*!
        This function has to be public since it is called from the global FCN.
        It should not be called otherwise!

        \param [in] par an array storing the various parameter values
        \param [in] npar the number of free parameters
    */
    virtual void setParsFromMinuit( Double_t* par, Int_t npar );

    //! Calculates the total negative log-likelihood
    /*!
        This function has to be public since it is called from the global FCN.
        It should not be called otherwise!
    */
    virtual Double_t getTotNegLogLikelihood();

    //! Set model parameters from a file
    /*!
        \param [in] fileName the name of the file with parameters to set
        \param [in] treeName the name of the tree in the file corresponding to the parameters to set
        \param [in] fix whether to fix (set constant) the loaded parameters, or leave them floating
    */
    void setParametersFromFile( const TString& fileName, const TString& treeName, const Bool_t fix );

    //! Set model parameters from a given std::map
    /*!
        Only parameters named in the map are imported, all others are set to their values specified
        in the model configuration.

        \param [in] parameters map from parameter name to imported value
        \param [in] fix whether to fix (set constant) the loaded parameters, or leave them floating
    */
    void setParametersFromMap( const std::map<TString, Double_t>& parameters, const Bool_t fix );

    //! Set named model parameters from a file
    /*!
        Identical to setParametersFromFile, but only import parameters named from parameters set.
        All others are set to their values specified in the model configuration.

        \param [in] fileName the name of the file with parameters to set
        \param [in] treeName the name of the tree in the file corresponding to the parameters to set
        \param [in] parameters the set of parameters to import from the file
        \param [in] fix whether to fix (set constant) the loaded parameters, or leave them floating
    */
    void setNamedParameters( const TString& fileName,
                             const TString& treeName,
                             const std::set<TString>& parameters,
                             const Bool_t fix );

    //! Set named model parameters from a given std::map, with fallback to those from a file
    /*!
        Parameters named in the map are imported with their specified values. All other parameters are set
        to the values corresponding to the value in the given file.

        \param [in] fileName the name of the file with parameters to set
        \param [in] treeName the name of the tree in the file corresponding to the parameters to set
        \param [in] parameters map from parameter name to imported value (override parameters form the file)
        \param [in] fix whether to fix (set constant) the loaded parameters, or leave them floating
    */
    void setParametersFileFallback( const TString& fileName,
                                    const TString& treeName,
                                    const std::map<TString, Double_t>& parameters,
                                    const Bool_t fix );

    //! Calculate the DP amplitude(s) for a given DP position
    /*!
        If not already done, this function will initialise the DP models

        \param [in] m13Sq the invariant mass squared of children 1 and 3
        \param [in] m23Sq the invariant mass squared of children 2 and 3
        \return a container of complex amplitudes, labelled to indicate to which component they belong
    */
    virtual std::map<TString, LauComplex> getDPAmps( const Double_t m13Sq, const Double_t m23Sq ) = 0;

    //! Calculate the DP likelihood(s) for a given DP position
    /*!
        If not already done, this function will initialise the DP models

        \param [in] m13Sq the invariant mass squared of children 1 and 3
        \param [in] m23Sq the invariant mass squared of children 2 and 3
        \return a container of likelihood values, labelled to indicate to which component they belong
    */
    virtual std::map<TString, Double_t> getDPLikelihoods( const Double_t m13Sq,
                                                          const Double_t m23Sq ) = 0;

  protected:
    // Some typedefs

    //! List of Pdfs
    typedef std::vector<LauAbsPdf*> LauPdfPList;
    //! List of parameter pointers
    typedef std::vector<LauParameter*> LauParameterPList;
    //! List of parameter pointers
    typedef std::vector<LauAbsRValue*> LauAbsRValuePList;
    //! Set of parameter pointers
    typedef std::set<LauParameter*> LauParameterPSet;
    //! List of parameters
    typedef std::vector<LauParameter> LauParameterList;
    //! A type to store background classes
    typedef std::map<UInt_t, TString> LauBkgndClassMap;

    //! Clear the vectors containing fit parameters
    void clearFitParVectors();

    //! Clear the vectors containing extra ntuple variables
    void clearExtraVarVectors();

    //! Weighting - allows e.g. MC events to be weighted by the DP model
    /*!
        \param [in] dataFileName the name of the data file
        \param [in] dataTreeName the name of the tree containing the data
    */
    virtual void weightEvents( const TString& dataFileName, const TString& dataTreeName ) = 0;

    //! Generate toy MC
    /*!
        \param [in] dataFileName the name of the file where the generated events are stored
        \param [in] dataTreeName the name of the tree used to store the variables
        \param [in] histFileName the name of the histogram output file (currently not used)
        \param [in] tableFileNameBase the name the latex output file
    */
    virtual void generate( const TString& dataFileName,
                           const TString& dataTreeName,
                           const TString& histFileName,
                           const TString& tableFileNameBase );

    //! The method that actually generates the toy MC events for the given experiment
    /*!
        \return the success/failure flag of the generation procedure
    */
    virtual Bool_t genExpt() = 0;

    //! Perform the total fit
    /*!
        \param [in] dataFileName the name of the data file
        \param [in] dataTreeName the name of the tree containing the data
        \param [in] histFileName the name of the histogram output file
        \param [in] tableFileNameBase the name the of latex output file
    */
    void fit( const TString& dataFileName,
              const TString& dataTreeName,
              const TString& histFileName,
              const TString& tableFileNameBase );

    //! Routine to perform the actual fit for a given experiment
    void fitExpt();

    //! Routine to perform the minimisation
    /*!
        \return the success/failure flag of the fit
    */
    Bool_t runMinimisation();

    //! Create a toy MC sample from the fitted parameters
    /*!
        \param [in] mcFileName the file name where the toy sample will be stored
        \param [in] tableFileName name of the output tex file
    */
    void createFitToyMC( const TString& mcFileName, const TString& tableFileName );

    //! Read in the data for the current experiment
    /*!
        \return the number of events read in
    */
    virtual UInt_t readExperimentData();

    //! Open the input file and verify that all required variables are present
    /*!
        \param [in] dataFileName the name of the input file
        \param [in] dataTreeName the name of the input tree
    */
    virtual Bool_t verifyFitData( const TString& dataFileName, const TString& dataTreeName );

    //! Cache the input data values to calculate the likelihood during the fit
    virtual void cacheInputFitVars() = 0;

    //! Cache the value of the sWeights to be used in the sFit
    virtual void cacheInputSWeights();

    //! Initialise the fit par vectors
    /*!
        Each class that inherits from this one must implement
        this sensibly for all vectors specified in
        clearFitParVectors, i.e. specify parameter names,
        initial, min, max and fixed values
    */
    virtual void initialise() = 0;

    //! Recalculate normalisation the signal DP model(s)
    virtual void recalculateNormalisation() = 0;

    //! Initialise the DP models
    virtual void initialiseDPModels() = 0;

    /*!
        For each amp in the fit this function takes its
        particular parameters and from them calculates the
        single complex number that is its coefficient.
        The vector of these coeffs can then be passed to the
        signal dynamics.
    */
    virtual void updateCoeffs() = 0;

    //! This function (specific to each model) calculates anything that depends on the fit parameter values
    virtual void propagateParUpdates() = 0;

    //! Calculate the sum of the log-likelihood over the specified events
    /*!
        \param [in] iStart the event number of the first event to be considered
        \param [in] iEnd the event number of the final event to be considered
    */
    Double_t getLogLikelihood( UInt_t iStart, UInt_t iEnd );

    //! Calculate the penalty terms to the log likelihood from Gaussian constraints
    Double_t getLogLikelihoodPenalty();

    //! Calculates the likelihood for a given event
    /*!
        \param [in] iEvt the event number
    */
    virtual Double_t getTotEvtLikelihood( UInt_t iEvt ) = 0;

    //! Returns the sum of the expected events over all hypotheses; used in the EML fit scenario
    virtual Double_t getEventSum() const = 0;

    //! Prints the values of all the fit variables for the specified event - useful for diagnostics
    /*!
        \param [in] iEvt the event number
    */
    virtual void printEventInfo( UInt_t iEvt ) const;

    //! Same as printEventInfo, but printing out the values of the variables in the fit
    virtual void printVarsInfo() const;

    //! Update initial fit parameters if required
    virtual void checkInitFitParams() = 0;

    //! Setup saving of fit results to ntuple/LaTeX table etc.
    /*!
        \param [in] histFileName the file name for the output histograms
        \param [in] tableFileName the file name for the latex output file
    */
    virtual void setupResultsOutputs( const TString& histFileName, const TString& tableFileName );

    //! Package the initial fit parameters for transmission to the coordinator
    /*!
        \param [out] array the array to be filled with the LauParameter objects
    */
    virtual void prepareInitialParArray( TObjArray& array );

    //! Perform all finalisation actions
    /*!
        - Receive the results of the fit from the coordinator
        - Perform any finalisation routines
        - Package the finalised fit parameters for transmission back to the coordinator

        \param [in] fitStat the status of the fit, e.g. status code, EDM, NLL
        \param [in] parsFromCoordinator the parameters at the fit minimum
        \param [in] covMat the fit covariance matrix
        \param [out] parsToCoordinator the array to be filled with the finalised LauParameter objects
    */
    virtual void finaliseExperiment( const LauAbsFitter::FitStatus& fitStat,
                                     const TObjArray* parsFromCoordinator,
                                     const TMatrixD* covMat,
                                     TObjArray& parsToCoordinator );

    //! Write the results of the fit into the ntuple
    /*!
        \param [in] tableFileName the structure containing the results of the fit
    */
    virtual void finaliseFitResults( const TString& tableFileName ) = 0;

    //! Save the pdf Plots for all the resonances of experiment number fitExp
    /*!
        \param [in] label  prefix for the file name to be saved
    */
    virtual void savePDFPlots( const TString& label ) = 0;

    //! Save the pdf Plots for the sum of ressonances correspondint to "sin" of experiment number fitExp
    /*!
        \param [in] label  prefix for the file name to be saved
        \param [in] spin   spin of the wave to be saved
    */
    virtual void savePDFPlotsWave( const TString& label, const Int_t& spin ) = 0;

    //! Write the latex table
    /*!
        \param [in] outputFile the name of the output file
    */
    virtual void writeOutTable( const TString& outputFile ) = 0;

    //! Store the per-event likelihood values
    virtual void storePerEvtLlhds() = 0;

    //! Calculate the sPlot data
    virtual void calculateSPlotData();

    //! Make sure all parameters hold their genValue as the current value
    void setGenValues();

    //! Method to set up the storage for background-related quantities called by setBkgndClassNames
    virtual void setupBkgndVectors() = 0;

    //! Check if the given background class is in the list
    /*!
        \param [in] className the name of the class to check
        \return true or false
    */
    Bool_t validBkgndClass( const TString& className ) const;

    //! The number assigned to a background class
    /*!
        \param [in] className the name of the class to check
        \return the background class ID number
    */
    UInt_t bkgndClassID( const TString& className ) const;

    //! Get the name of a background class from the number
    /*!
        \param [in] classID the ID number of the background class
        \return the class name
    */
    const TString& bkgndClassName( UInt_t classID ) const;

    //! Setup the generation ntuple branches
    virtual void setupGenNtupleBranches() = 0;

    //! Add a branch to the gen tree for storing an integer
    /*!
        \param [in] name the name of the branch
    */
    virtual void addGenNtupleIntegerBranch( const TString& name );

    //! Add a branch to the gen tree for storing a double
    /*!
        \param [in] name the name of the branch
    */
    virtual void addGenNtupleDoubleBranch( const TString& name );

    //! Set the value of an integer branch in the gen tree
    /*!
        \param [in] name the name of the branch
        \param [in] value the value to be stored
    */
    virtual void setGenNtupleIntegerBranchValue( const TString& name, Int_t value );

    //! Set the value of a double branch in the gen tree
    /*!
        \param [in] name the name of the branch
        \param [in] value the value to be stored
    */
    virtual void setGenNtupleDoubleBranchValue( const TString& name, Double_t value );

    //! Get the value of an integer branch in the gen tree
    /*!
        \param [in] name the name of the branch
        \return the value of the parameter
    */
    virtual Int_t getGenNtupleIntegerBranchValue( const TString& name ) const;

    //! Get the value of a double branch in the gen tree
    /*!
        \param [in] name the name of the branch
        \return the value of the parameter
    */
    virtual Double_t getGenNtupleDoubleBranchValue( const TString& name ) const;

    //! Fill the gen tuple branches
    virtual void fillGenNtupleBranches();

    //! Setup the branches of the sPlot tuple
    virtual void setupSPlotNtupleBranches() = 0;

    //! Add a branch to the sPlot tree for storing an integer
    /*!
        \param [in] name the name of the branch
    */
    virtual void addSPlotNtupleIntegerBranch( const TString& name );

    //! Add a branch to the sPlot tree for storing a double
    /*!
        \param [in] name the name of the branch
    */
    virtual void addSPlotNtupleDoubleBranch( const TString& name );

    //! Set the value of an integer branch in the sPlot tree
    /*!
        \param [in] name the name of the branch
        \param [in] value the value to be stored
    */
    virtual void setSPlotNtupleIntegerBranchValue( const TString& name, Int_t value );

    //! Set the value of a double branch in the sPlot tree
    /*!
        \param [in] name the name of the branch
        \param [in] value the value to be stored
    */
    virtual void setSPlotNtupleDoubleBranchValue( const TString& name, Double_t value );

    //! Fill the sPlot tuple
    virtual void fillSPlotNtupleBranches();

    //! Returns the names of all variables in the fit
    virtual LauSPlot::NameSet variableNames() const = 0;

    //! Returns the names and yields of species that are free in the fit
    virtual LauSPlot::NumbMap freeSpeciesNames() const = 0;

    //! Returns the names and yields of species that are fixed in the fit
    virtual LauSPlot::NumbMap fixdSpeciesNames() const = 0;

    //! Returns the species and variables for all 2D PDFs in the fit
    virtual LauSPlot::TwoDMap twodimPDFs() const = 0;

    //! Check if the signal is split into well-reconstructed and mis-reconstructed types
    virtual Bool_t splitSignal() const = 0;

    //! Check if the mis-reconstructed signal is to be smeared in the DP
    virtual Bool_t scfDPSmear() const = 0;

    //! Add the given parameter to the list of all fit parameters
    /*!
        \param [in] param the parameter
        \param [in] addFixed if true add the parameter even if it is not floating
        \return the number of parameters added
    */
    UInt_t addFitParameters( LauParameter* param, const Bool_t addFixed = kFALSE );

    //! Add the given parameter(s) to the list of all fit parameters
    /*!
        \param [in] param the parameter (which may depend on other parameters)
        \param [in] addFixed if true add the parameter even if it is not floating
        \return the number of parameters added
    */
    UInt_t addFitParameters( LauAbsRValue* param, const Bool_t addFixed = kFALSE );

    //! Add the given parameter(s) to the list of all fit parameters
    /*!
        \param [in] paramList a list of parameters
        \param [in] addFixed if true add the parameter even if it is not floating
        \return the number of parameters added
    */
    UInt_t addFitParameters( LauParameterPList& paramList, const Bool_t addFixed = kFALSE );

    //! Add the given parameter(s) to the list of all fit parameters
    /*!
        \param [in] paramList a list of parameters
        \param [in] addFixed if true add the parameter even if it is not floating
        \return the number of parameters added
    */
    UInt_t addFitParameters( LauAbsRValuePList& paramList, const Bool_t addFixed = kFALSE );

    //! Add parameters of the PDFs in the list to the list of all fit parameters
    /*!
        \param [in] pdfList a list of Pdfs
        \param [in] addFixed if true add the parameter even if it is not floating
        \return the number of parameters added
    */
    UInt_t addFitParameters( LauPdfPList& pdfList, const Bool_t addFixed = kFALSE );

    //! Add the given parameter to the list of resonance parameters and the list of all fit parameters
    /*!
        \param [in] param a resonance parameter
        \return the number of parameters added
    */
    UInt_t addResonanceParameters( LauParameter* param );

    //! Add the given parameter(s) to the list of resonance parameters and the list of all fit parameters
    /*!
        \param [in] paramList a list of resonance parameters
        \return the number of parameters added
    */
    UInt_t addResonanceParameters( LauParameterPList& paramList );

    //! Add parameters to the list of Gaussian constrained parameters
    void addConParameters();

    //! Print the fit parameters for all PDFs in the list
    /*!
        \param [in] pdfList a list of Pdfs
        \param [in] fout the output stream to write to
    */
    void printFitParameters( const LauPdfPList& pdfList, std::ostream& fout ) const;

    //! Update the fit parameters for all PDFs in the list
    /*!
        \param [in] pdfList a list of Pdfs
    */
    void updateFitParameters( LauPdfPList& pdfList );

    //! Have all PDFs in the list cache the data
    /*!
        \param [in] pdfList the list of pdfs
        \param [in] theData the data from the fit
    */
    void cacheInfo( LauPdfPList& pdfList, const LauFitDataTree& theData );

    //! Calculate the product of the per-event likelihoods of the PDFs in the list
    /*!
        \param [in] pdfList the list of pdfs
        \param [in] iEvt the event number
    */
    Double_t prodPdfValue( LauPdfPList& pdfList, UInt_t iEvt );

    //! Do any of the PDFs have a dependence on the DP?
    /*!
        \return the flag to indicated if there is a DP dependence
    */
    Bool_t pdfsDependOnDP() const { return pdfsDependOnDP_; }

    //! Do any of the PDFs have a dependence on the DP?
    /*!
        \param [in] dependOnDP the flag to indicated if there is a DP dependence
    */
    void pdfsDependOnDP( Bool_t dependOnDP ) { pdfsDependOnDP_ = dependOnDP; }

    //! Const access to the fit variables
    const LauParameterPList& fitPars() const { return fitVars_; }

    //! Const access the fit variables which affect the DP normalisation
    const LauParameterPSet& resPars() const { return resVars_; }

    //! Const access the extra variables
    const LauParameterList& extraPars() const { return extraVars_; }

    //! Non-const access the extra variables
    LauParameterList& extraPars() { return extraVars_; }

    //! Const access the Gaussian constrained variables
    const LauAbsRValuePList& conPars() const { return conVars_; }

    //! Const access the gen ntuple
    const LauGenNtuple* genNtuple() const { return genNtuple_; }

    //! Access the gen ntuple
    LauGenNtuple* genNtuple() { return genNtuple_; }

    //! Const access the sPlot ntuple
    const LauGenNtuple* sPlotNtuple() const { return sPlotNtuple_; }

    //! Access the sPlot ntuple
    LauGenNtuple* sPlotNtuple() { return sPlotNtuple_; }

    //! Const access the data store
    const LauFitDataTree* fitData() const { return inputFitData_; }

    //! Access the data store
    LauFitDataTree* fitData() { return inputFitData_; }

    //! Imported parameters file name
    TString fixParamFileName_;

    //! Imported parameters tree name
    TString fixParamTreeName_;

    //! Map from imported parameter name to value
    std::map<TString, Double_t> fixParamMap_;

    //! Imported parameter names
    std::set<TString> fixParamNames_;

    //! Whether to fix the loaded parameters (kTRUE) or leave them floating (kFALSE)
    Bool_t fixParams_;

    //! The set of parameters that are imported (either from a file or by value) and not
    //  set to be fixed in the fit. In addition to those from fixParamNames_, these
    //  include those imported from a file.
    std::set<LauParameter*> allImportedFreeParams_;

  private:
    //! Copy constructor (not implemented)
    LauAbsFitModel( const LauAbsFitModel& rhs );

    //! Copy assignment operator (not implemented)
    LauAbsFitModel& operator=( const LauAbsFitModel& rhs );

    // Various control booleans

    //! Option to make toy from 1st successful experiment
    Bool_t compareFitData_;
    //! Option to output a .C file of PDF's
    Bool_t savePDF_;
    //! Option to output a Latex format table
    Bool_t writeLatexTable_;
    //! Option to write sPlot data
    Bool_t writeSPlotData_;
    //! Option to store DP efficiencies in the sPlot ntuple
    Bool_t storeDPEff_;
    //! Option to randomise the initial values of the fit parameters
    Bool_t randomFit_;
    //! Option to perform an extended ML fit
    Bool_t emlFit_;
    //! Option to perform Poisson smearing
    Bool_t poissonSmear_;
    //! Option to enable embedding
    Bool_t enableEmbedding_;
    //! Option to include the DP as part of the fit
    Bool_t usingDP_;
    //! Option to state if pdfs depend on DP position
    Bool_t pdfsDependOnDP_;

    // Info on number of experiments and number of events

    //! Internal vector of fit parameters
    LauParameterPList fitVars_;

    //! Internal set of the same fit parameters (used to check uniqueness)
    LauParameterPSet fitVarsSet_;

    //! Internal set of fit parameters upon which the DP normalisation depends
    LauParameterPSet resVars_;

    //! Extra variables that aren't in the fit but are stored in the ntuple
    LauParameterList extraVars_;

    //! Internal vectors of Gaussian parameters
    LauAbsRValuePList conVars_;

    // Input data and output ntuple

    //! The input data
    LauFitDataTree* inputFitData_;
    //! The generated ntuple
    LauGenNtuple* genNtuple_;
    //! The sPlot ntuple
    LauGenNtuple* sPlotNtuple_;

    // Background class names

    //! The background class names
    LauBkgndClassMap bkgndClassNames_;
    //! An empty string
    const TString nullString_;

    // sFit related variables

    //! Option to perfom the sFit
    Bool_t doSFit_;
    //! The name of the sWeight branch
    TString sWeightBranchName_;
    //! The vector of sWeights
    std::vector<Double_t> sWeights_;
    //! The sWeight scaling factor
    Double_t sWeightScaleFactor_;

    // Fit timers

    //! The fit timer
    TStopwatch timer_;
    //! The total fit timer
    TStopwatch cumulTimer_;

    //! The output table name
    TString outputTableName_;

    // Comparison toy MC related variables

    //! The output file name for Toy MC
    TString fitToyMCFileName_;
    //! The output table name for Toy MC
    TString fitToyMCTableName_;
    //! The scaling factor (toy vs data statistics)
    UInt_t fitToyMCScale_;
    //! Option to perform Poisson smearing
    Bool_t fitToyMCPoissonSmear_;

    // sPlot related variables

    //! The name of the sPlot file
    TString sPlotFileName_;
    //! The name of the sPlot tree
    TString sPlotTreeName_;
    //! Control the verbosity of the sFit
    TString sPlotVerbosity_;

    ClassDef( LauAbsFitModel, 0 )    // Abstract interface to fit/toyMC model
};

#endif
