
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

/*! \file LauSimpleFitModel.hh
    \brief File containing declaration of LauSimpleFitModel class.
*/

/*! \class LauSimpleFitModel
    \brief Define a Dalitz plot according to the isobar model.

    A class that allows the user to define a three-body B meson
    Dalitz plot according to the isobar model, i.e. defining a set of resonances
    that have complex amplitudes that can interfere with each other.
*/

#ifndef LAU_SIMPLE_FIT_MODEL
#define LAU_SIMPLE_FIT_MODEL

#include "LauAbsFitModel.hh"
#include "LauComplex.hh"
#include "LauParameter.hh"

#include "TString.h"

#include <vector>

class TH2;
class LauAbsBkgndDPModel;
class LauAbsCoeffSet;
class LauIsobarDynamics;
class LauAbsPdf;
class LauEffModel;
class LauEmbeddedData;
class LauKinematics;
class LauScfMap;

class LauSimpleFitModel : public LauAbsFitModel {

  public:
    //! Constructor
    /*!
        \param [in] sigDPModel the signal DP model
    */
    explicit LauSimpleFitModel( LauIsobarDynamics* sigDPModel );

    //! Destructor
    virtual ~LauSimpleFitModel();

    //! Set the signal event yield
    /*!
        \param [in] nSigEvents contains the signal yield and boolean to fix it or not
    */
    virtual void setNSigEvents( LauParameter* nSigEvents );

    //! Set the background event yield(s)
    /*!
        The name of the parameter must be that of the corresponding background category (so that it can be correctly assigned)

        \param [in] nBkgndEvents contains the name, yield and option to fix the background yield
    */
    virtual void setNBkgndEvents( LauAbsRValue* nBkgndEvents );

    //! Set the background DP models
    /*!
        \param [in] bkgndClass the name of the background class
        \param [in] bkgndModel the DP model of the background
    */
    void setBkgndDPModel( const TString& bkgndClass, LauAbsBkgndDPModel* bkgndModel );

    //! Split the signal component into well-reconstructed and mis-reconstructed parts
    /*!
        The nomenclature used here is TM (truth-matched) and SCF (self cross feed)
        In this option, the SCF fraction is DP-dependent
        Can also optionally provide a smearing matrix to smear the SCF DP PDF

        \param [in] dpHisto the DP histogram of the SCF fraction value
        \param [in] upperHalf boolean flag to specify that the supplied histogram contains only the upper half of a symmetric DP (or lower half if using square DP coordinates)
        \param [in] fluctuateBins whether the bins on the histogram should be varied in accordance with their uncertainties (for evaluation of systematic uncertainties)
        \param [in] scfMap the (optional) smearing matrix
    */
    void splitSignalComponent( const TH2* dpHisto,
                               const Bool_t upperHalf     = kFALSE,
                               const Bool_t fluctuateBins = kFALSE,
                               LauScfMap* scfMap          = 0 );

    //! Split the signal component into well reconstructed and mis-reconstructed parts
    /*!
        The nomenclature used here is TM (truth-matched) and SCF (self cross feed)
        In this option, the SCF fraction is a single global number

        \param [in] scfFrac the SCF fraction value
        \param [in] fixed whether the SCF fraction is fixed or floated in the fit
    */
    void splitSignalComponent( const Double_t scfFrac, const Bool_t fixed );

    //! Determine whether we are splitting the signal into TM and SCF parts
    Bool_t useSCF() const { return useSCF_; }

    //! Determine whether the SCF fraction is DP-dependent
    Bool_t useSCFHist() const { return useSCFHist_; }

    //! Determine if we are smearing the SCF DP PDF
    Bool_t smearSCFDP() const { return ( scfMap_ != 0 ); }

    //! Set the signal PDF for a given variable
    /*!
        \param [in] pdf the PDF to be added to the signal model
    */
    void setSignalPdf( LauAbsPdf* pdf );

    //! Set the SCF PDF for a given variable
    /*!
        \param [in] pdf the PDF to be added to the SCF model
    */
    void setSCFPdf( LauAbsPdf* pdf );

    //! Set the background PDF
    /*!
        \param [in] bkgndClass the name of the background class
        \param [in] pdf the PDF to be added to the background model
    */
    void setBkgndPdf( const TString& bkgndClass, LauAbsPdf* pdf );

    //! Embed full simulation events for the signal, rather than generating toy from the PDFs
    /*!
        \param [in] fileName the name of the file containing the events
        \param [in] treeName the name of the tree
        \param [in] reuseEventsWithinEnsemble sample with replacement but only replace events once each experiment has been generated
        \param [in] reuseEventsWithinExperiment sample with immediate replacement
        \param [in] useReweighting perform an accept/reject routine using the configured signal amplitude model based on the MC-truth DP coordinate
    */
    void embedSignal( const TString& fileName,
                      const TString& treeName,
                      Bool_t reuseEventsWithinEnsemble,
                      Bool_t reuseEventsWithinExperiment = kFALSE,
                      Bool_t useReweighting              = kFALSE );

    //! Embed full simulation events for the given background class, rather than generating toy from the PDFs
    /*!
        \param [in] bkgndClass the name of the background class
        \param [in] fileName the name of the file containing the events
        \param [in] treeName the name of the tree
        \param [in] reuseEventsWithinEnsemble sample with replacement but only replace events once each experiment has been generated
        \param [in] reuseEventsWithinExperiment sample with immediate replacement
    */
    void embedBkgnd( const TString& bkgndClass,
                     const TString& fileName,
                     const TString& treeName,
                     Bool_t reuseEventsWithinEnsemble,
                     Bool_t reuseEventsWithinExperiment = kFALSE );

    //! Set the DP amplitude coefficients
    /*!
        The name of the coeffSet must match the name of one of the resonances in the DP model.

        The supplied order of coefficients will be rearranged to match the order in which the resonances are stored in the dynamics, see LauIsobarDynamics::addResonance.

        \param [in] coeffSet the set of coefficients
    */
    virtual void setAmpCoeffSet( LauAbsCoeffSet* coeffSet );

    //! Calculate the DP amplitude(s) for a given DP position
    /*!
        If not already done, this function will initialise the DP models

        \param [in] m13Sq the invariant mass squared of children 1 and 3
        \param [in] m23Sq the invariant mass squared of children 2 and 3
        \return a container of complex amplitudes, labelled to indicate to which component they belong
    */
    virtual std::map<TString, LauComplex> getDPAmps( const Double_t m13Sq, const Double_t m23Sq );

    //! Calculate the DP likelihood(s) for a given DP position
    /*!
        If not already done, this function will initialise the DP models

        \param [in] m13Sq the invariant mass squared of children 1 and 3
        \param [in] m23Sq the invariant mass squared of children 2 and 3
        \return a container of likelihood values, labelled to indicate to which component they belong
    */
    virtual std::map<TString, Double_t> getDPLikelihoods( const Double_t m13Sq, const Double_t m23Sq );

  protected:
    //! Define a map to be used to store a category name and numbers
    typedef std::map<TString, std::pair<Int_t, Double_t>> LauGenInfo;

    //! Typedef for a vector of background DP models
    typedef std::vector<LauAbsBkgndDPModel*> LauBkgndDPModelList;

    //! Typedef for a vector of background PDFs
    typedef std::vector<LauPdfPList> LauBkgndPdfsList;

    //! Typedef for a vector of background yields
    typedef std::vector<LauAbsRValue*> LauBkgndYieldList;

    //! Typedef for a vector of embedded data objects
    typedef std::vector<LauEmbeddedData*> LauBkgndEmbDataList;

    //! Typedef for a vector of booleans to flag if events are reused
    typedef std::vector<Bool_t> LauBkgndReuseEventsList;

    //! Weight events based on the DP model
    /*!
        \param [in] dataFileName the name of the data file
        \param [in] dataTreeName the name of the data tree
    */
    virtual void weightEvents( const TString& dataFileName, const TString& dataTreeName );

    //! Initialise the fit
    virtual void initialise();

    //! Initialise the signal DP model
    virtual void initialiseDPModels();

    //! Recalculate Normalization the signal DP models
    virtual void recalculateNormalisation();

    //! Update the coefficients
    virtual void updateCoeffs();

    //! Toy MC generation and fitting overloaded functions
    virtual Bool_t genExpt();

    //! Calculate things that depend on the fit parameters after they have been updated by Minuit
    virtual void propagateParUpdates();

    //! Read in the input fit data variables, e.g. m13Sq and m23Sq
    virtual void cacheInputFitVars();

    //! Check the initial fit parameters
    virtual void checkInitFitParams();

    //! Get the fit results and store them
    /*!
        \param [in] tablePrefixName prefix for the name of the output file
    */
    virtual void finaliseFitResults( const TString& tablePrefixName );

    //! Print the fit fractions, total DP rate and mean efficiency
    /*!
        \param [out] output the stream to which to print
    */
    virtual void printFitFractions( std::ostream& output );

    //! Write the fit results in latex table format
    /*!
        \param [in] outputFile the name of the output file
    */
    virtual void writeOutTable( const TString& outputFile );

    //! Save the pdf Plots for all the resonances of experiment number fitExp
    /*!
        \param [in] label  prefix for the file name to be saved
    */
    virtual void savePDFPlots( const TString& label );

    //! Save the pdf Plots for the sum of ressonances correspondint to "sin" of experiment number fitExp
    /*!
        \param [in] label  prefix for the file name to be saved
        \param [in] spin   spin of the wave to be saved
    */
    virtual void savePDFPlotsWave( const TString& label, const Int_t& spin );

    //! Store the per event likelihood values
    virtual void storePerEvtLlhds();

    // Methods to do with calculating the likelihood functions
    // and manipulating the fitting parameters.

    //! Get the total likelihood for each event
    /*!
        \param [in] iEvt the event number
    */
    virtual Double_t getTotEvtLikelihood( UInt_t iEvt );

    //! Calculate the signal and background likelihoods for the DP for a given event
    /*!
        \param [in] iEvt the event number
    */
    virtual void getEvtDPLikelihood( UInt_t iEvt );

    //! Calculate the SCF likelihood for the DP for a given event
    /*!
        \param [in] iEvt the event number
    */
    virtual Double_t getEvtSCFDPLikelihood( UInt_t iEvt );

    //! Determine the signal and background likelihood for the extra variables for a given event
    /*!
        \param [in] iEvt the event number
    */
    virtual void getEvtExtraLikelihoods( UInt_t iEvt );

    //! Get the total number of events
    /*!
        \return the total number of events
    */
    virtual Double_t getEventSum() const;

    //! Set the fit parameters for the DP model
    void setSignalDPParameters();

    //! Set the fit parameters for the extra PDFs
    void setExtraPdfParameters();

    //! Set the initial yields
    void setFitNEvents();

    //! Set-up other parameters that are derived from the fit results, e.g. fit fractions
    void setExtraNtupleVars();

    //! Randomise the initial fit parameters
    void randomiseInitFitPars();

    //! Define the length of the background vectors
    virtual void setupBkgndVectors();

    //! Determine the number of events to generate for each hypothesis
    std::pair<LauGenInfo, Bool_t> eventsToGenerate();

    //! Generate signal event
    Bool_t generateSignalEvent();

    //! Generate background event
    /*!
        \param [in] bkgndID ID number of the background class
    */
    Bool_t generateBkgndEvent( UInt_t bkgndID );

    //! Setup the required ntuple branches
    void setupGenNtupleBranches();

    //! Store all of the DP information
    void setDPBranchValues();

    //! Generate from the extra PDFs
    /*!
        \param [in] extraPdfs the list of extra PDFs
        \param [in] embeddedData the embedded data sample
    */
    void generateExtraPdfValues( LauPdfPList* extraPdfs, LauEmbeddedData* embeddedData );

    //! Add sPlot branches for the extra PDFs
    /*!
        \param [in] extraPdfs the list of extra PDFs
        \param [in] prefix the list of prefixes for the branch names
    */
    void addSPlotNtupleBranches( const LauPdfPList* extraPdfs, const TString& prefix );

    //! Set the branches for the sPlot ntuple with extra PDFs
    /*!
        \param [in] extraPdfs the list of extra PDFs
        \param [in] prefix the list of prefixes for the branch names
        \param [in] iEvt the event number
    */
    Double_t setSPlotNtupleBranchValues( LauPdfPList* extraPdfs, const TString& prefix, UInt_t iEvt );

    //! Update the signal events after Minuit sets background parameters
    void updateSigEvents();

    //! Add branches to store experiment number and the event number within the experiment
    virtual void setupSPlotNtupleBranches();

    //! Returns the names of all variables in the fit
    virtual LauSPlot::NameSet variableNames() const;

    //! Returns the names and yields of species that are free in the fit
    virtual LauSPlot::NumbMap freeSpeciesNames() const;

    //! Returns the names and yields of species that are fixed in the fit
    virtual LauSPlot::NumbMap fixdSpeciesNames() const;

    //! Returns the species and variables for all 2D PDFs in the fit
    virtual LauSPlot::TwoDMap twodimPDFs() const;

    //! Check if the signal is split into well-reconstructed and mis-reconstructed types
    virtual Bool_t splitSignal() const { return this->useSCF(); }

    //! Check if the mis-reconstructed signal is to be smeared in the DP
    virtual Bool_t scfDPSmear() const { return ( scfMap_ != 0 ); }

    //! Append fake data points to the inputData for each bin in the SCF smearing matrix
    /*!
        We'll be caching the DP amplitudes and efficiencies of the centres of the true bins.
        To do so, we attach some fake points at the end of inputData, the number of the entry
        minus the total number of events corresponding to the number of the histogram for that
        given true bin in the LauScfMap object. (What this means is that when Laura is provided with
        the LauScfMap object by the user, it's the latter who has to make sure that it contains the
        right number of histograms and in exactly the right order!)

        \param [in] inputData the fit data
    */
    void appendBinCentres( LauFitDataTree* inputData );

  private:
    //! Copy constructor (not implemented)
    LauSimpleFitModel( const LauSimpleFitModel& rhs );

    //! Copy assignment operator (not implemented)
    LauSimpleFitModel& operator=( const LauSimpleFitModel& rhs );

    //! The signal Dalitz plot model
    LauIsobarDynamics* sigDPModel_;

    //! The background Dalitz Plot model
    LauBkgndDPModelList bkgndDPModels_;

    //! The Dalitz plot kinematics object
    LauKinematics* kinematics_;

    //! The signal PDFs
    LauPdfPList signalPdfs_;

    //! The SCF PDFs
    LauPdfPList scfPdfs_;

    //! The background PDFs
    LauBkgndPdfsList bkgndPdfs_;

    //! Background boolean
    Bool_t usingBkgnd_;

    //! Number of signal components
    UInt_t nSigComp_;

    //! Number of signal DP parameters
    UInt_t nSigDPPar_;

    //! Number of extra PDF parameters
    UInt_t nExtraPdfPar_;

    //! Number of normalisation parameters (i.e. yields)
    UInt_t nNormPar_;

    //! Magnitudes and Phases
    std::vector<LauAbsCoeffSet*> coeffPars_;

    //! Fit fractions
    LauParArray fitFrac_;

    //! Fit fractions (uncorrected for the efficiency)
    LauParArray fitFracEffUnCorr_;

    //! The mean efficiency
    LauParameter meanEff_;

    //! The average DP rate
    LauParameter dpRate_;

    //! Signal yield
    LauParameter* signalEvents_;

    //! Background yield(s)
    LauBkgndYieldList bkgndEvents_;

    //! Is the signal split into TM and SCF
    Bool_t useSCF_;

    //! Is the SCF fraction DP-dependent
    Bool_t useSCFHist_;

    //! The (global) SCF fraction parameter
    LauParameter scfFrac_;

    //! The histogram giving the DP-dependence of the SCF fraction
    LauEffModel* scfFracHist_;

    //! The smearing matrix for the SCF DP PDF
    LauScfMap* scfMap_;

    //! The cached values of the SCF fraction for each event
    std::vector<Double_t> recoSCFFracs_;

    //! The cached values of the SCF fraction for each bin centre
    std::vector<Double_t> fakeSCFFracs_;

    //! The cached values of the sqDP jacobians for each event
    std::vector<Double_t> recoJacobians_;

    //! The cached values of the sqDP jacobians for each true bin
    std::vector<Double_t> fakeJacobians_;

    //! Run choice variables
    Bool_t compareFitData_;

    //! The complex coefficients
    std::vector<LauComplex> coeffs_;

    // Embedding full simulation events

    //! The signal event tree
    LauEmbeddedData* signalTree_;

    //! The background event tree
    LauBkgndEmbDataList bkgndTree_;

    //! Boolean to reuse signal events
    Bool_t reuseSignal_;

    //! Boolean to use reweighting
    Bool_t useReweighting_;

    //! Vector of booleans to reuse background events
    LauBkgndReuseEventsList reuseBkgnd_;

    //! Signal likelihood value
    Double_t sigDPLike_;

    //! SCF likelihood value
    Double_t scfDPLike_;

    //! Background likelihood value(s)
    std::vector<Double_t> bkgndDPLike_;

    //! Signal likelihood from extra PDFs
    Double_t sigExtraLike_;

    //! SCF likelihood from extra PDFs
    Double_t scfExtraLike_;

    //! Background likelihood value(s) from extra PDFs
    std::vector<Double_t> bkgndExtraLike_;

    //! Total signal likelihood
    Double_t sigTotalLike_;

    //! Total SCF likelihood
    Double_t scfTotalLike_;

    //! Total background likelihood(s)
    std::vector<Double_t> bkgndTotalLike_;

    ClassDef( LauSimpleFitModel, 0 )    // Total fit/ToyMC model
};

#endif
