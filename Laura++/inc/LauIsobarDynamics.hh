
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

/*! \file LauIsobarDynamics.hh
    \brief File containing declaration of LauIsobarDynamics class.
*/

/*! \class LauIsobarDynamics
    \brief Class for defining signal dynamics using the isobar model.
*/

#ifndef LAU_ISOBAR_DYNAMICS
#define LAU_ISOBAR_DYNAMICS

#include "LauAbsResonance.hh"
#include "LauComplex.hh"

#include "TString.h"

#include <set>
#include <vector>

class LauCacheData;
class LauDaughters;
class LauAbsEffModel;
class LauAbsIncohRes;
class LauFitDataTree;
class LauKMatrixPropagator;
class LauDPPartialIntegralInfo;
class LauKinematics;

class LauIsobarDynamics {

  public:
    //! The type used for containing multiple self cross feed fraction models for different categories (e.g. tagging categories)
    typedef std::map<Int_t, LauAbsEffModel*> LauTagCatScfFractionModelMap;

    //! The possible statuses for toy MC generation
    enum ToyMCStatus {
        GenOK, /*!< Generation completed OK */
        MaxIterError, /*!< Maximum allowed number of iterations completed without success (ASqMax is too high) */
        ASqMaxError /*!< An amplitude squared value was returned that was larger than the maximum expected (ASqMax is too low) */
    };

    //! Constructor
    /*!
        \param [in] daughters the daughters of the decay
        \param [in] effModel the model to describe the efficiency across the Dalitz plot
        \param [in] scfFractionModel the model to describe the fraction of poorly constructed events (the self cross feed fraction) across the Dalitz plot
    */
    LauIsobarDynamics( LauDaughters* daughters,
                       LauAbsEffModel* effModel,
                       LauAbsEffModel* scfFractionModel = 0 );

    //! Constructor
    /*!
        \param [in] daughters the daughters of the decay
        \param [in] effModel the model to describe the efficiency across the Dalitz plot
        \param [in] scfFractionModel the models to describe the fraction of poorly constructed events (the self cross feed fraction) across the Dalitz plot for various tagging categories
    */
    LauIsobarDynamics( LauDaughters* daughters,
                       LauAbsEffModel* effModel,
                       LauTagCatScfFractionModelMap scfFractionModel );

    //! Destructor
    virtual ~LauIsobarDynamics();

    //! Initialise the Dalitz plot dynamics
    /*!
        \param [in] coeffs the complex coefficients for the resonances
    */
    void initialise( const std::vector<LauComplex>& coeffs );

    //! recalculate Normalization
    void recalculateNormalisation();

    //! Set the name of the file to which to save the results of the integrals
    /*!
        \param [in] fileName the name of the file
    */
    inline void setIntFileName( const TString& fileName ) { intFileName_ = fileName; }

    // Integration
    //! Set the widths of the bins to use when integrating across the Dalitz plot or square Dalitz plot
    /*!
        Specify the bin widths required when performing the DP integration.
        Note that the integration is not performed in m13^2 vs m23^2 space
        but in either m13 vs m23 space or mPrime vs thetaPrime space,
        with the appropriate Jacobian applied.
        The default bin widths in m13 vs m23 space are 0.005 GeV.
        The default bin widths in mPrime vs thetaPrime space are 0.001.

        \param [in] m13BinWidth the bin width to use when integrating over m13
        \param [in] m23BinWidth the bin width to use when integrating over m23
        \param [in] mPrimeBinWidth the bin width to use when integrating over mPrime
        \param [in] thPrimeBinWidth the bin width to use when integrating over thetaPrime
    */
    void setIntegralBinWidths( const Double_t m13BinWidth,
                               const Double_t m23BinWidth,
                               const Double_t mPrimeBinWidth  = 0.001,
                               const Double_t thPrimeBinWidth = 0.001 );

    //! Set the value below which a resonance width is considered to be narrow
    /*!
        Narrow resonances trigger different integration behaviour - dividing the DP into regions where a finer binning is used.
        This can cause high memory usage, so use this method and LauIsobarDynamics::setIntegralBinningFactor to tune this behaviour, if needed.

        \param [in] narrowWidth the value below which a resonance is considered to be narrow (defaults to 0.02 GeV/c2)
    */
    void setNarrowResonanceThreshold( const Double_t narrowWidth ) { narrowWidth_ = narrowWidth; }

    //! Set the factor relating the width of a narrow resonance and the binning size in its integration region
    /*!
        Narrow resonances trigger different integration behaviour - dividing the DP into regions where a finer binning is used.
        This can cause high memory usage, so use this method and LauIsobarDynamics::setNarrowResonanceThreshold to tune this behaviour, if needed.

        \param [in] binningFactor the factor by which the resonance width is divided to obtain the bin size (defaults to 100)
    */
    void setIntegralBinningFactor( const Double_t binningFactor )
    {
        binningFactor_ = binningFactor;
    }

    //! Force the symmetrisation of the integration in m13 <-> m23 for non-symmetric but flavour-conjugate final states
    /*!
        This can be necessary for time-dependent fits (where interference terms between A and Abar need to be integrated)

        \param [in] force toggle forcing symmetrisation of the integration for apparently flavour-conjugate final states
    */
    void forceSymmetriseIntegration( const Bool_t force ) { forceSymmetriseIntegration_ = force; }

    //! Add a resonance to the Dalitz plot
    /*!
        NB the stored order of resonances is:
        - Firstly, all coherent resonances (i.e. those added using addResonance() or addKMatrixProdPole() or addKMatrixProdSVP()) in order of addition
        - Followed by all incoherent resonances (i.e. those added using addIncoherentResonance()) in order of addition

        \param [in] resName the name of the resonant particle
        \param [in] resPairAmpInt the index of the daughter not produced by the resonance
        \param [in] resType the model for the resonance dynamics
        \param [in] bwCategory the Blatt-Weisskopf barrier factor category
        \return the newly created resonance
    */
    LauAbsResonance* addResonance( const TString& resName,
                                   const Int_t resPairAmpInt,
                                   const LauAbsResonance::LauResonanceModel resType,
                                   const LauBlattWeisskopfFactor::BlattWeisskopfCategory bwCategory =
                                       LauBlattWeisskopfFactor::Default );

    //! Add an incoherent resonance to the Dalitz plot
    /*!
        NB the stored order of resonances is:
        - Firstly, all coherent resonances (i.e. those added using addResonance() or addKMatrixProdPole() or addKMatrixProdSVP()) in order of addition
        - Followed by all incoherent resonances (i.e. those added using addIncoherentResonance()) in order of addition

        \param [in] resName the name of the resonant particle
        \param [in] resPairAmpInt the index of the daughter not produced by the resonance
        \param [in] resType the model for the resonance dynamics
        \return the newly created resonance
    */
    LauAbsResonance* addIncoherentResonance( const TString& resName,
                                             const Int_t resPairAmpInt,
                                             const LauAbsResonance::LauResonanceModel resType );

    //! Define a new K-matrix Propagator
    /*!
        \param [in] propName the name of the propagator
        \param [in] paramFileName the file that defines the propagator
        \param [in] resPairAmpInt the index of the bachelor
        \param [in] nChannels the number of channels
        \param [in] nPoles the number of poles
        \param [in] rowIndex the index of the row to be used when summing over all amplitude channels: S-wave corresponds to rowIndex = 1.
    */
    LauKMatrixPropagator* defineKMatrixPropagator( const TString& propName,
                                                   const TString& paramFileName,
                                                   Int_t resPairAmpInt,
                                                   Int_t nChannels,
                                                   Int_t nPoles,
                                                   Int_t rowIndex = 1 );

    //! Add a K-matrix production pole term to the model
    /*!
        NB the stored order of resonances is:
        - Firstly, all coherent resonances (i.e. those added using addResonance() or addKMatrixProdPole() or addKMatrixProdSVP()) in order of addition
        - Followed by all incoherent resonances (i.e. those added using addIncoherentResonance()) in order of addition

        \param [in] poleName the name of the pole
        \param [in] propName the name of the propagator to use
        \param [in] poleIndex the index of the pole within the propagator
        \param [in] useProdAdler boolean to turn on/off the production Adler zero factor (default = off)
    */
    void addKMatrixProdPole( const TString& poleName,
                             const TString& propName,
                             Int_t poleIndex,
                             Bool_t useProdAdler = kFALSE );

    //! Add a K-matrix slowly-varying part (SVP) term to the model
    /*!
        NB the stored order of resonances is:
        - Firstly, all coherent resonances (i.e. those added using addResonance() or addKMatrixProdPole() or addKMatrixProdSVP()) in order of addition
        - Followed by all incoherent resonances (i.e. those added using addIncoherentResonance()) in order of addition

        \param [in] SVPName the name of the term
        \param [in] propName the name of the propagator to use
        \param [in] channelIndex the index of the channel within the propagator
        \param [in] useProdAdler boolean to turn on/off the production Adler zero factor (default = off)
    */
    void addKMatrixProdSVP( const TString& SVPName,
                            const TString& propName,
                            Int_t channelIndex,
                            Bool_t useProdAdler = kFALSE );

    //! Set the maximum value of A squared to be used in the accept/reject
    /*!
        Disables the automatic determination of ASqMax

        \param [in] value the new value
    */
    inline void setASqMaxValue( Double_t value )
    {
        aSqMaxSet_  = value;
        aSqMaxAuto_ = kFALSE;
    }

    //! Retrieve the maximum value of A squared to be used in the accept/reject
    /*!
        \return the maximum value of A squared
    */
    inline Double_t getASqMaxSetValue() const { return aSqMaxSet_; }

    //! Retrieve the maximum of A squared that has been found while generating
    /*!
        \return the maximum of A squared that has been found
    */
    inline Double_t getASqMaxVarValue() const { return aSqMaxVar_; }

    //! Generate a toy MC signal event
    /*!
        \return kTRUE if the event is successfully generated, kFALSE otherwise
    */
    Bool_t generate();

    //! Check the status of the toy MC generation
    /*!
        \param [in] printErrorMessages whether error messages should be printed
        \param [in] printInfoMessages whether info messages should be printed
        \return the status of the toy MC generation
    */
    ToyMCStatus checkToyMC( Bool_t printErrorMessages = kTRUE, Bool_t printInfoMessages = kFALSE );

    //! Retrieve the maximum number of iterations allowed when generating an event
    /*!
        \return the maximum number of iterations allowed
    */
    inline Int_t maxGenIterations() const { return iterationsMax_; }

    //! Calculate the likelihood (and all associated information) for the given event number
    /*!
        \param [in] iEvt the event number
    */
    void calcLikelihoodInfo( const UInt_t iEvt );

    //! Calculate the likelihood (and all associated information) given values of the Dalitz plot coordinates
    /*!
        \param [in] m13Sq the invariant mass squared of the first and third daughters
        \param [in] m23Sq the invariant mass squared of the second and third daughters
    */
    void calcLikelihoodInfo( const Double_t m13Sq, const Double_t m23Sq );

    //! Calculate the likelihood (and all associated information) given values of the Dalitz plot coordinates and the tagging category
    /*!
        Also obtain the self cross feed fraction to cache with the rest of the Dalitz plot quantities.

        \param [in] m13Sq the invariant mass squared of the first and third daughters
        \param [in] m23Sq the invariant mass squared of the second and third daughters
        \param [in] tagCat the tagging category
    */
    void calcLikelihoodInfo( const Double_t m13Sq, const Double_t m23Sq, const Int_t tagCat );

    //! Calculate the fit fractions, mean efficiency and total DP rate
    /*!
        \param [in] init whether the calculated values should be stored as the initial/generated values or the fitted values
    */
    void calcExtraInfo( const Bool_t init = kFALSE );

    //! Calculates whether an event with the current kinematics should be accepted in order to produce a distribution of events that matches the model e.g. when reweighting embedded data
    /*!
        Uses the accept/reject method.

        \return kTRUE if the event has been accepted, kFALSE otherwise
    */
    Bool_t gotReweightedEvent();

    //! Calculate the acceptance rate, for events with the current kinematics, when generating events according to the model
    /*!
        \return the weight for the current kinematics
    */
    Double_t getEventWeight();

    //! Retrieve the total amplitude for the current event
    /*!
        \return the total amplitude
    */
    inline const LauComplex& getEvtDPAmp() const { return totAmp_; }

    //! Retrieve the invariant mass squared of the first and third daughters in the current event
    /*!
        \return the invariant mass squared of the first and third daughters in the current event
    */
    inline Double_t getEvtm13Sq() const { return m13Sq_; }

    //! Retrieve the invariant mass squared of the second and third daughters in the current event
    /*!
        \return the invariant mass squared of the second and third daughters in the current event
    */
    inline Double_t getEvtm23Sq() const { return m23Sq_; }

    //! Retrieve the square Dalitz plot coordinate, m', for the current event
    /*!
        \return the square Dalitz plot coordinate, m', for the current event
    */
    inline Double_t getEvtmPrime() const { return mPrime_; }

    //! Retrieve the square Dalitz plot coordinate, theta', for the current event
    /*!
        \return the square Dalitz plot coordinate, theta', for the current event
    */
    inline Double_t getEvtthPrime() const { return thPrime_; }

    //! Retrieve the efficiency for the current event
    /*!
        \return the efficiency for the current event
    */
    inline Double_t getEvtEff() const { return eff_; }

    //! Retrieve the fraction of events that are poorly reconstructed (the self cross feed fraction) for the current event
    /*!
        \return the self cross feed fraction for the current event
    */
    inline Double_t getEvtScfFraction() const { return scfFraction_; }

    //! Retrieve the Jacobian, for the transformation into square DP coordinates, for the current event
    /*!
        \return the Jacobian for the current event
    */
    inline Double_t getEvtJacobian() const { return jacobian_; }

    //! Retrieve the total intensity multiplied by the efficiency for the current event
    /*!
        \return the total intensity multiplied by the efficiency for the current event
    */
    inline Double_t getEvtIntensity() const { return ASq_; }

    //! Retrieve the likelihood for the current event
    /*!
        The likelihood is the normalised total intensity:
        evtLike_ = ASq_/DPNorm_

        \return the likelihood for the current event
    */
    inline Double_t getEvtLikelihood() const { return evtLike_; }

    //! Retrieve the normalised dynamic part of the amplitude of the given amplitude component at the current point in the Dalitz plot
    /*!
        \param [in] resID the index of the component within the model
        \return the amplitude of the given component
    */
    inline LauComplex getDynamicAmp( const Int_t resID ) const
    {
        return ff_[resID].scale( fNorm_[resID] );
    }

    //! Retrieve the Amplitude  of resonance resID
    /*!
        \param [in] resID the index of the component within the model
        \return the amplitude of the given component
    */
    inline LauComplex getFullAmplitude( const Int_t resID ) const
    {
        return Amp_[resID] * this->getDynamicAmp( resID );
    }

    //! Retrieve the event-by-event running totals of amplitude cross terms for all pairs of amplitude components
    /*!
        \return the event-by-event running totals of amplitude cross terms
    */
    inline const std::vector<std::vector<LauComplex>>& getFiFjSum() const { return fifjSum_; }

    //! Retrieve the event-by-event running totals of efficiency corrected amplitude cross terms for all pairs of amplitude components
    /*!
        \return the event-by-event running totals of amplitude cross terms with efficiency corrections applied
    */
    inline const std::vector<std::vector<LauComplex>>& getFiFjEffSum() const { return fifjEffSum_; }

    //! Retrieve the normalisation factors for the dynamic parts of the amplitudes for all of the amplitude components
    /*!
        \return the normalisation factors
    */
    inline const std::vector<Double_t>& getFNorm() const { return fNorm_; }

    //! Fill the internal data structure that caches the resonance dynamics
    /*!
        \param [in] fitDataTree the data source
    */
    void fillDataTree( const LauFitDataTree& fitDataTree );

    //! Recache the amplitude values for those that have changed
    void modifyDataTree();

    //! Check whether this model includes a named resonance
    /*!
        \param [in] resName the resonance
        \return true if the resonance is present, false otherwise
    */
    Bool_t hasResonance( const TString& resName ) const;

    //! Retrieve the index for the given resonance
    /*!
        \param [in] resName the resonance
        \return the index of the resonance if it is present, -1 otherwise
    */
    Int_t resonanceIndex( const TString& resName ) const;

    //! Retrieve the name of the charge conjugate of a named resonance
    /*!
        \param [in] resName the resonance
        \return the name of the charge conjugate
    */
    TString getConjResName( const TString& resName ) const;

    //! Retrieve the named resonance
    /*!
        \param [in] resName the name of the resonance to retrieve
        \return the requested resonance
    */
    const LauAbsResonance* findResonance( const TString& resName ) const;

    //! Retrieve a resonance by its index
    /*!
        \param [in] resIndex the index of the resonance to retrieve
        \return the requested resonance
    */
    const LauAbsResonance* getResonance( const UInt_t resIndex ) const;

    //! Update the complex coefficients for the resonances
    /*!
        \param [in] coeffs the new set of coefficients
    */
    void updateCoeffs( const std::vector<LauComplex>& coeffs );

    //! Collate the resonance parameters to initialise (or re-initialise) the model
    /*!
        NB: This has been factored out of the initialise() method to allow for use in the
        importation of parameters in LauAbsFitModel
    */
    void collateResonanceParameters();

    //! Set the helicity flip flag for new amplitude components
    /*!
        \param [in] boolean the helicity flip flag
    */
    inline void flipHelicityForCPEigenstates( const Bool_t boolean ) { flipHelicity_ = boolean; }

    //! Retrieve the mean efficiency across the Dalitz plot
    /*!
        \return the mean efficiency across the Dalitz plot
    */
    inline const LauParameter& getMeanEff() const { return meanDPEff_; }

    //! Retrieve the overall Dalitz plot rate
    /*!
        \return the overall Dalitz plot rate
    */
    inline const LauParameter& getDPRate() const { return DPRate_; }

    //! Retrieve the fit fractions for the amplitude components
    /*!
        \return the fit fractions
    */
    inline const LauParArray& getFitFractions() const { return fitFrac_; }

    //! Retrieve the fit fractions for the amplitude components
    /*!
        \return the fit fractions
    */
    inline const LauParArray& getFitFractionsEfficiencyUncorrected() const
    {
        return fitFracEffUnCorr_;
    }

    //! Retrieve the total number of amplitude components
    /*!
        \return the total number of amplitude components
    */
    inline UInt_t getnTotAmp() const { return nAmp_ + nIncohAmp_; }

    //! Retrieve the number of coherent amplitude components
    /*!
        \return the number of coherent amplitude components
    */
    inline UInt_t getnCohAmp() const { return nAmp_; }

    //! Retrieve the number of incoherent amplitude components
    /*!
        \return the number of incoherent amplitude components
    */
    inline UInt_t getnIncohAmp() const { return nIncohAmp_; }

    //! Retrieve the normalisation factor for the log-likelihood function
    /*!
        \return the normalisation factor
    */
    inline Double_t getDPNorm() const { return DPNorm_; }

    //! Retrieve the daughters
    /*!
        \return the daughters
    */
    inline const LauDaughters* getDaughters() const { return daughters_; }

    //! Retrieve the Dalitz plot kinematics
    /*!
        \return the Dalitz plot kinematics
    */
    inline const LauKinematics* getKinematics() const { return kinematics_; }

    //! Retrieve the Dalitz plot kinematics
    /*!
        \return the Dalitz plot kinematics
    */
    inline LauKinematics* getKinematics() { return kinematics_; }

    //! Retrieve the model for the efficiency across the Dalitz plot
    /*!
        \return the efficiency model
    */
    inline const LauAbsEffModel* getEffModel() const { return effModel_; }

    //! Check whether a self cross feed fraction model is being used
    /*!
        \return true if a self cross feed fraction model is being used, false otherwise
    */
    inline Bool_t usingScfModel() const { return ! scfFractionModel_.empty(); }

    //! Retrieve any extra parameters/quantities (e.g. K-matrix total fit fractions)
    /*!
        \return any extra parameters
    */
    inline const std::vector<LauParameter>& getExtraParameters() const { return extraParameters_; }

    //! Retrieve the floating parameters of the resonance models
    /*!
        \return the list of floating parameters
    */
    inline std::vector<LauParameter*>& getFloatingParameters() { return resonancePars_; }

    //! Whether to calculate separate rho and omega fit-fractions from LauRhoOmegaMix
    inline void calculateRhoOmegaFitFractions( const Bool_t calcFF )
    {
        calculateRhoOmegaFitFractions_ = calcFF;
    }

  protected:
    //! Print a summary of the model to be used
    void initSummary();

    //! Initialise the internal storage for this model
    void initialiseVectors();

    //! Zero the various values used to store integrals info
    void resetNormVectors();

    //! Calculate the Dalitz plot normalisation integrals across the whole Dalitz plot
    void calcDPNormalisation();

    //! Form the regions that are produced by the spaces between narrow resonances
    /*!
        \param [in] regions the regions defined around narrow resonances
        \param [in] min the minimum value of the invariant mass
        \param [in] max the maximum value of the invariant mass
        \return vector of pointers to LauDPPartialIntegralInfo objects that contain the individual coarse regions
    */
    std::vector<std::pair<Double_t, Double_t>> formGapsFromRegions(
        const std::vector<std::pair<Double_t, Double_t>>& regions,
        const Double_t min,
        const Double_t max ) const;

    //! Removes entries in the vector of LauDPPartialIntegralInfo* that are null
    /*!
        \param [in] regions the list of region pointers
    */
    void cullNullRegions( std::vector<LauDPPartialIntegralInfo*>& regions ) const;

    //! Wrapper for LauDPPartialIntegralInfo constructor
    /*!
        \param [in] minm13 the minimum of the m13 range
        \param [in] maxm13 the maximum of the m13 range
        \param [in] minm23 the minimum of the m23 range
        \param [in] maxm23 the maximum of the m23 range
        \param [in] m13BinWidth the m13 bin width
        \param [in] m23BinWidth the m23 bin width
        \param [in] precision the precision required for the Gauss-Legendre weights
        \param [in] nAmp the number of coherent amplitude components
        \param [in] nIncohAmp the number of incoherent amplitude components
        \return 0 if the integration region has no internal points, otherwise returns a pointer to the newly constructed LauDPPartialIntegralInfo object
    */
    LauDPPartialIntegralInfo* newDPIntegrationRegion( const Double_t minm13,
                                                      const Double_t maxm13,
                                                      const Double_t minm23,
                                                      const Double_t maxm23,
                                                      const Double_t m13BinWidth,
                                                      const Double_t m23BinWidth,
                                                      const Double_t precision,
                                                      const UInt_t nAmp,
                                                      const UInt_t nIncohAmp ) const;

    //! Correct regions to ensure that the finest integration grid takes precedence
    /*!
        \param [in] regions the windows in invariant mass
        \param [in] binnings the corresponding binnings for each window
    */
    void correctDPOverlap( std::vector<std::pair<Double_t, Double_t>>& regions,
                           const std::vector<Double_t>& binnings ) const;

    //! Create the integration grid objects for the m23 narrow resonance regions, including the overlap regions with the m13 narrow resonances
    /*!
        The overlap regions will have an m13Binnings x m23Binnings grid.
        The other regions will have a defaultBinning x m23Binnings grid.

        \param [in] m13Regions the limits of each narrow-resonance region in m13
        \param [in] m23Regions the limits of each narrow-resonance region in m23
        \param [in] m13Binnings the binning of each narrow-resonance region in m13
        \param [in] m23Binnings the binning of each narrow-resonance region in m23
        \param [in] precision the precision required for the Gauss-Legendre weights
        \param [in] defaultBinning the binning used in the bulk of the phase space
        \return vector of pointers to LauDPPartialIntegralInfo objects that contain the individual regions
    */
    std::vector<LauDPPartialIntegralInfo*> m23IntegrationRegions(
        const std::vector<std::pair<Double_t, Double_t>>& m13Regions,
        const std::vector<std::pair<Double_t, Double_t>>& m23Regions,
        const std::vector<Double_t>& m13Binnings,
        const std::vector<Double_t>& m23Binnings,
        const Double_t precision,
        const Double_t defaultBinning ) const;

    //! Create the integration grid objects for the m13 narrow resonance regions, excluding the overlap regions with the m23 narrow resonances
    /*!
        The regions will have a m13Binnings x defaultBinning grid.
        The overlap regions are created by the m23IntegrationRegions function.

        \param [in] m13Regions the limits of each narrow-resonance region in m13
        \param [in] m23Regions the limits of each narrow-resonance region in m23
        \param [in] m13Binnings the binning of each narrow-resonance region in m13
        \param [in] precision the precision required for the Gauss-Legendre weights
        \param [in] defaultBinning the binning used in the bulk of the phase space
        \return vector of pointers to LauDPPartialIntegralInfo objects that contain the individual regions
    */
    std::vector<LauDPPartialIntegralInfo*> m13IntegrationRegions(
        const std::vector<std::pair<Double_t, Double_t>>& m13Regions,
        const std::vector<std::pair<Double_t, Double_t>>& m23Regions,
        const std::vector<Double_t>& m13Binnings,
        const Double_t precision,
        const Double_t defaultBinning ) const;

    //! Calculate the Dalitz plot normalisation integrals across the whole Dalitz plot
    void calcDPNormalisationScheme();

    //! Determine which amplitudes and integrals need to be recalculated
    void findIntegralsToBeRecalculated();

    //! Calculate the Dalitz plot normalisation integrals over a given range
    /*!
        \param [in] intInfo the integration information object
    */
    void calcDPPartialIntegral( LauDPPartialIntegralInfo* intInfo );

    //! Write the results of the integrals (and related information) to a file
    void writeIntegralsFile();

    //! Set the dynamic part of the amplitude for a given amplitude component at the current point in the Dalitz plot
    /*!
        \param [in] index the index of the amplitude component
        \param [in] realPart the real part of the amplitude
        \param [in] imagPart the imaginary part of the amplitude
    */
    void setFFTerm( const UInt_t index, const Double_t realPart, const Double_t imagPart );

    //! Set the dynamic part of the intensity for a given incoherent amplitude component at the current point in the Dalitz plot
    /*!
        \param [in] index the index of the incoherent amplitude component
        \param [in] value the intensity
    */
    void setIncohIntenTerm( const UInt_t index, const Double_t value );

    //! Calculate the amplitudes for all resonances for the current kinematics
    void calculateAmplitudes();

    //! Calculate or retrieve the cached value of the amplitudes for all resonances at the specified integration grid point
    /*!
        \param [in,out] intInfo the integration information object
        \param [in] m13Point the grid index in m13
        \param [in] m23Point the grid index in m23
    */
    void calculateAmplitudes( LauDPPartialIntegralInfo* intInfo,
                              const UInt_t m13Point,
                              const UInt_t m23Point );

    //! Add the amplitude values (with the appropriate weight) at the current grid point to the running integral values
    /*!
        \param [in] weight the weight to apply
    */
    void addGridPointToIntegrals( const Double_t weight );

    //! Calculate the total Dalitz plot amplitude at the current point in the Dalitz plot
    /*!
        \param [in] useEff whether to apply efficiency corrections
    */
    void calcTotalAmp( const Bool_t useEff );

    //! Obtain the efficiency of the current event from the model
    /*!
        \return the efficiency
    */
    Double_t retrieveEfficiency();

    //! Obtain the self cross feed fraction of the current event from the model
    /*!
        \param [in] tagCat the tagging category of the current event
        \return the self cross feed fraction
    */
    Double_t retrieveScfFraction( Int_t tagCat );

    //! Set the maximum of A squared that has been found
    /*!
        \param [in] value the new value
    */
    inline void setASqMaxVarValue( Double_t value ) { aSqMaxVar_ = value; }

    //! Calculate the normalisation factor for the log-likelihood function
    /*!
        \return the normalisation factor
    */
    Double_t calcSigDPNorm();

    //! Calculate the dynamic part of the amplitude for a given component at the current point in the Dalitz plot
    /*!
        \param [in] index the index of the amplitude component within the model
    */
    LauComplex resAmp( const UInt_t index );

    //! Calculate the dynamic part of the intensity for a given incoherent component at the current point in the Dalitz plot
    /*!
        \param [in] index the index of the incoherent component within the model
    */
    Double_t incohResAmp( const UInt_t index );

    //! Load the data for a given event
    /*!
        \param [in] iEvt the number of the event
    */
    void setDataEventNo( UInt_t iEvt );

    //! Retrieve the named resonance
    /*!
        \param [in] resName the name of the resonance to retrieve
        \return the requested resonance
    */
    LauAbsResonance* findResonance( const TString& resName );

    //! Retrieve a resonance by its index
    /*!
        \param [in] resIndex the index of the resonance to retrieve
        \return the requested resonance
    */
    LauAbsResonance* getResonance( const UInt_t resIndex );

    //! Remove the charge from the given particle name
    /*!
        \param [in,out] string the particle name
    */
    void removeCharge( TString& string ) const;

    //! Check whether a resonance is a K-matrix component of a given propagator
    /*!
        \param [in] resAmpInt the index of the resonance within the model
        \param [in] propName the name of the K-matrix propagator
        \return true if the resonance is a component of the given propagator, otherwise return false
    */
    Bool_t gotKMatrixMatch( UInt_t resAmpInt, const TString& propName ) const;

  private:
    //! Copy constructor (not implemented)
    LauIsobarDynamics( const LauIsobarDynamics& rhs );

    //! Copy assignment operator (not implemented)
    LauIsobarDynamics& operator=( const LauIsobarDynamics& rhs );

    //! The type used for containing the K-matrix propagators
    typedef std::map<TString, LauKMatrixPropagator*> KMPropMap;

    //! The type used for mapping K-matrix components to their propagators
    typedef std::map<TString, TString> KMStringMap;

    //! The daughters of the decay
    LauDaughters* daughters_;

    //! The kinematics of the decay
    LauKinematics* kinematics_;

    //! The efficiency model across the Dalitz plot
    LauAbsEffModel* effModel_;

    //! The self cross feed fraction models across the Dalitz plot
    /*!
        These model the fraction of signal events that are poorly reconstructed (the self cross feed fraction) as a function of Dalitz plot position.
        If the self cross feed is depependent on the tagging category then seperate models can be defined.
    */
    LauTagCatScfFractionModelMap scfFractionModel_;

    //! The number of amplitude components
    UInt_t nAmp_;

    //! The number of incoherent amplitude components
    UInt_t nIncohAmp_;

    //! The complex coefficients for the amplitude components
    std::vector<LauComplex> Amp_;

    //! The normalisation factor for the log-likelihood function
    Double_t DPNorm_;

    //! The fit fractions for the amplitude components
    LauParArray fitFrac_;

    //! The efficiency-uncorrected fit fractions for the amplitude components
    LauParArray fitFracEffUnCorr_;

    //! The overall Dalitz plot rate
    LauParameter DPRate_;

    //! The mean efficiency across the Dalitz plot
    LauParameter meanDPEff_;

    //! The cached data for all events
    std::vector<LauCacheData*> data_;

    //! The cached data for the current event
    LauCacheData* currentEvent_;

    //! any extra parameters/quantities (e.g. K-matrix total fit fractions)
    std::vector<LauParameter> extraParameters_;

    //! The resonances in the model
    std::vector<LauAbsResonance*> sigResonances_;

    //! The incoherent resonances in the model
    std::vector<LauAbsIncohRes*> sigIncohResonances_;

    //! The K-matrix propagators
    KMPropMap kMatrixPropagators_;

    //! The names of the M-matrix components in the model mapped to their propagators
    KMStringMap kMatrixPropSet_;

    //! The resonance types of all of the amplitude components
    std::vector<TString> resTypAmp_;

    //! The index of the daughter not produced by the resonance for each amplitude component
    std::vector<Int_t> resPairAmp_;

    //! The resonance types of all of the incoherent amplitude components
    std::vector<TString> incohResTypAmp_;

    //! The index of the daughter not produced by the resonance for each incoherent amplitude component
    std::vector<Int_t> incohResPairAmp_;

    //! The PDG codes of the daughters
    std::vector<Int_t> typDaug_;

    //! Whether the Dalitz plot is symmetrical
    Bool_t symmetricalDP_;

    //! Whether the Dalitz plot is fully symmetric
    Bool_t fullySymmetricDP_;

    //! Whether the Dalitz plot is a flavour-conjugate final state
    Bool_t flavConjDP_;

    //! Whether the integrals have been performed
    Bool_t integralsDone_;

    //! Whether the scheme for the integration has been determined
    Bool_t normalizationSchemeDone_;

    //! Force the symmetrisation of the integration in m13 <-> m23 for non-symmetric but flavour-conjugate final states
    Bool_t forceSymmetriseIntegration_;

    //! The storage of the integration scheme
    std::vector<LauDPPartialIntegralInfo*> dpPartialIntegralInfo_;

    //! The name of the file to save integrals to
    TString intFileName_;

    //! The bin width to use when integrating over m13
    Double_t m13BinWidth_;

    //! The bin width to use when integrating over m23
    Double_t m23BinWidth_;

    //! The bin width to use when integrating over mPrime
    Double_t mPrimeBinWidth_;

    //! The bin width to use when integrating over thetaPrime
    Double_t thPrimeBinWidth_;

    //! The value below which a resonance width is considered to be narrow
    Double_t narrowWidth_;

    //! The factor relating the width of the narrowest resonance and the binning size
    Double_t binningFactor_;

    //! The invariant mass squared of the first and third daughters
    Double_t m13Sq_;

    //! The invariant mass squared of the second and third daughters
    Double_t m23Sq_;

    //! The square Dalitz plot coordinate, m'
    Double_t mPrime_;

    //! The square Dalitz plot coordinate theta'
    Double_t thPrime_;

    //! The tagging category
    Int_t tagCat_;

    //! The efficiency at the current point in the Dalitz plot
    Double_t eff_;

    //!The fraction of events that are poorly reconstructed (the self cross feed fraction) at the current point in the Dalitz plot
    Double_t scfFraction_;

    //! The Jacobian, for the transformation into square DP coordinates at the current point in the Dalitz plot
    Double_t jacobian_;

    //! The value of A squared for the current event
    Double_t ASq_;

    //! The normalised likelihood for the current event
    Double_t evtLike_;

    //! The total amplitude for the current event
    LauComplex totAmp_;

    //! The event-by-event running total of efficiency corrected amplitude cross terms for each pair of amplitude components
    /*!
        Calculated as the sum of ff_[i]*ff_[j]*efficiency for all events
    */
    std::vector<std::vector<LauComplex>> fifjEffSum_;

    //! The event-by-event running total of the amplitude cross terms for each pair of amplitude components
    /*!
        Calculated as the sum of ff_[i]*ff_[j] for all events
    */
    std::vector<std::vector<LauComplex>> fifjSum_;

    //! The dynamic part of the amplitude for each amplitude component at the current point in the Dalitz plot
    std::vector<LauComplex> ff_;

    //! The dynamic part of the intensity for each incoherent amplitude component at the current point in the Dalitz plot
    std::vector<Double_t> incohInten_;

    //! The event-by-event running total of the dynamical amplitude squared for each amplitude component
    std::vector<Double_t> fSqSum_;

    //! The event-by-event running total of the dynamical amplitude squared for each amplitude component
    std::vector<Double_t> fSqEffSum_;

    //! The normalisation factors for the dynamic parts of the amplitude for each amplitude component
    std::vector<Double_t> fNorm_;

    //! The maximum allowed number of attempts when generating an event
    Int_t iterationsMax_;

    //! The number of unsucessful attempts to generate an event so far
    Int_t nSigGenLoop_;

    //! The maximum allowed value of A squared
    Double_t aSqMaxSet_;

    //! The maximum value of A squared that has been seen so far while generating
    Double_t aSqMaxVar_;

    //! Flag to generate aSqMaxSet_ once generate is called
    Bool_t aSqMaxAuto_ { kTRUE };

    //! The helicity flip flag for new amplitude components
    Bool_t flipHelicity_;

    //! Flag to recalculate the normalisation
    Bool_t recalcNormalisation_;

    //! List of floating resonance parameters
    std::vector<LauParameter*> resonancePars_;

    //! List of floating resonance parameter values from previous calculation
    std::vector<Double_t> resonanceParValues_;

    //! Indices in sigResonances_ to point to the corresponding signal resonance(s) for each floating parameter
    std::vector<std::vector<UInt_t>> resonanceParResIndex_;

    //! Resonance indices for which the amplitudes and integrals should be recalculated
    std::set<UInt_t> integralsToBeCalculated_;

    //! Whether to calculate separate rho and omega fit fractions from the LauRhoOmegaMix model
    Bool_t calculateRhoOmegaFitFractions_;

    ClassDef( LauIsobarDynamics, 0 )
};

#endif
