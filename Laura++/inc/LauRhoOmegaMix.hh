
/*
Copyright 2016 University of Warwick

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

/*! \file LauRhoOmegaMix.hh
    \brief File containing declaration of LauRhoOmegaMix class.
*/

/*! \class LauRhoOmegaMix
    \brief Class for defining the rho-omega resonance mixing model

    Formulae from Paul Rensing thesis, SLAC Report 421 and Bill Dunwoodie's note
    http://www.slac.stanford.edu/~wmd/omega-rho_mixing/omega-rho_mixing.note
*/

#ifndef LAU_RHO_OMEGA_MIX
#define LAU_RHO_OMEGA_MIX

#include "LauAbsResonance.hh"
#include "LauComplex.hh"

#include "TString.h"

#include <vector>

class LauResonanceInfo;
class LauDaughters;
class LauParameter;

class LauRhoOmegaMix : public LauAbsResonance {

  public:
    //! Constructor
    /*!
        \param [in] resInfo the object containing information on the resonance name, mass, width, spin, charge, etc.
        \param [in] resType the model of the resonance
        \param [in] resPairAmpInt the number of the daughter not produced by the resonance
        \param [in] daughters the daughter particles
    */
    LauRhoOmegaMix( LauResonanceInfo* resInfo,
                    const LauAbsResonance::LauResonanceModel resType,
                    const Int_t resPairAmpInt,
                    const LauDaughters* daughters );

    //! Destructor
    virtual ~LauRhoOmegaMix();

    //! Initialise the model
    virtual void initialise();

    //! Initialise the rho resonance
    void initialiseRho();

    //! Initialise the omega resonance
    void initialiseOmega();

    //! Get the complex dynamical amplitude
    /*!
        \param [in] kinematics the kinematic variables of the current event
        \return the complex amplitude
    */
    virtual LauComplex amplitude( const LauKinematics* kinematics );

    //! Get the resonance model type
    /*!
        \return the resonance model type
    */
    virtual LauAbsResonance::LauResonanceModel getResonanceModel() const { return model_; }

    //! Set value of the various parameters
    /*!
        \param [in] name the name of the parameter to be changed
        \param [in] value the new parameter value
    */
    virtual void setResonanceParameter( const TString& name, const Double_t value );

    //! Allow the various parameters to float in the fit
    /*!
        \param [in] name the name of the parameter to be floated
    */
    virtual void floatResonanceParameter( const TString& name );

    //! Access the given resonance parameter
    /*!
        \param [in] name the name of the parameter
        \return the corresponding parameter
    */
    virtual LauParameter* getResonanceParameter( const TString& name );

    //! Retrieve the resonance parameters, e.g. so that they can be loaded into a fit
    /*!
        \return floating parameters of the resonance
    */
    virtual const std::vector<LauParameter*>& getFloatingParameters();

    //! Set which rho/omega amplitude to calculate for FF
    void setWhichAmpSq( Int_t which ) { whichAmpSq_ = which; }

  protected:
    //! Set the omega pole mass parameter
    /*!
        \param [in] mOmega new value for the omega mass parameter
    */
    void setmOmegaValue( const Double_t mOmega );

    //! Get the omega pole mass parameter value
    /*!
        \return value of the omega pole mass parameter
    */
    Double_t getmOmegaValue() const { return ( mOmega_ != 0 ) ? mOmega_->unblindValue() : 0.0; }

    //! Fix the omega pole mass parameter value
    /*!
        \return kTRUE if the omega pole mass parameter is fixed, kFALSE otherwise
    */
    Bool_t fixmOmegaValue() const { return ( mOmega_ != 0 ) ? mOmega_->fixed() : 0.0; }

    //! Set the omega pole width parameter
    /*!
        \param [in] wOmega new value for the omega width parameter
    */
    void setwOmegaValue( const Double_t wOmega );

    //! Get the omega pole width parameter value
    /*!
        \return value of the omega pole width parameter
    */
    Double_t getwOmegaValue() const { return ( wOmega_ != 0 ) ? wOmega_->unblindValue() : 0.0; }

    //! Fix the omega pole width parameter value
    /*!
        \return kTRUE if the omega pole width parameter is fixed, kFALSE otherwise
    */
    Bool_t fixwOmegaValue() const { return ( wOmega_ != 0 ) ? wOmega_->fixed() : 0.0; }

    //! Set the omega B magnitude mixing parameter
    /*!
        \param [in] magB new value for the omega B magnitude mixing parameter
    */
    void setmagBValue( const Double_t magB );

    //! Get the omega B magnitude mixing parameter
    /*!
        \return value of the omega B magnitude mixing parameter
    */
    Double_t getmagBValue() const { return ( magB_ != 0 ) ? magB_->unblindValue() : 0.0; }

    //! Fix the omega B magnitude mixing parameter value
    /*!
        \return kTRUE if the omega B magnitude mixing parameter is fixed, kFALSE otherwise
    */
    Bool_t fixmagBValue() const { return ( magB_ != 0 ) ? magB_->fixed() : 0.0; }

    //! Set the omega B phase mixing parameter
    /*!
        \param [in] phiB new value for the omega B phase mixing parameter
    */
    void setphiBValue( const Double_t phiB );

    //! Get the omega B phase mixing parameter
    /*!
        \return value of the omega B phase mixing parameter
    */
    Double_t getphiBValue() const { return ( phiB_ != 0 ) ? phiB_->unblindValue() : 0.0; }

    //! Fix the omega B phase mixing parameter value
    /*!
        \return kTRUE if the omega B phase mixing parameter is fixed, kFALSE otherwise
    */
    Bool_t fixphiBValue() const { return ( phiB_ != 0 ) ? phiB_->fixed() : 0.0; }

    //! Set the omega mixing parameter delta
    /*!
        \param [in] delta new value for the omega mixing parameter delta
    */
    void setdeltaValue( const Double_t delta );

    //! Get the omega mixing parameter delta
    /*!
        \return value of the omega mixing parameter delta
    */
    Double_t getdeltaValue() const { return ( delta_ != 0 ) ? delta_->unblindValue() : 0.0; }

    //! Fix the omega mixing parameter delta value
    /*!
        \return kTRUE if the omega mixing parameter delta is fixed, kFALSE otherwise
    */
    Bool_t fixdeltaValue() const { return ( delta_ != 0 ) ? delta_->fixed() : 0.0; }

    //! This is not called, amplitude is used directly instead
    virtual LauComplex resAmp( Double_t mass, Double_t spinTerm );

    //! Check that both daughters are the same type of particle
    void checkDaughterTypes() const;

  private:
    //! Copy constructor (not implemented)
    LauRhoOmegaMix( const LauRhoOmegaMix& rhs );

    //! Copy assignment operator (not implemented)
    LauRhoOmegaMix& operator=( const LauRhoOmegaMix& rhs );

    //! The model to use
    LauAbsResonance::LauResonanceModel model_;

    //! Previous value of the pole mass of the rho resonance
    Double_t rhoMass_;

    //! Previous value of the barrier radius of the rho resonance
    Double_t rhoResRadius_;

    //! Previous value of the parents barrier radius of the rho resonance
    Double_t rhoParRadius_;

    //! Pole mass of the omega contribution
    LauParameter* mOmega_;

    //! Initial default value of the omega pole mass from LauResonanceMaker
    Double_t mOmega0_;

    //! Current value of the omega pole mass (floating or fixed)
    Double_t mOmegaCur_;

    //! Pole width of the omega contribution
    LauParameter* wOmega_;

    //! Initial default value of the omega pole width from LauResonanceMaker
    Double_t wOmega0_;

    //! Current value of the omega pole mass (floating or fixed)
    Double_t wOmegaCur_;

    //! B magnitude parameter of the omega mixing contribution
    LauParameter* magB_;

    //! B phase parameter of the omega mixing contribution
    LauParameter* phiB_;

    //! delta parameter of the omega mixing contribution
    LauParameter* delta_;

    //! Boolean to specify if we want to use the denominator factor
    Bool_t useDenom_;

    //! Boolean to specify if we have performed the first initialisation
    Bool_t doneFirstInit_;

    //! Pointer to the rho (or first) resonance lineshape
    LauAbsResonance* rhoRes_;

    //! Pointer to the omega (second) resonance lineshape
    LauAbsResonance* omegaRes_;

    //! Which amplitude to calculate for rho/omega fit fractions
    Int_t whichAmpSq_;

    ClassDef( LauRhoOmegaMix, 0 )    // Rho-omega mixing model
};

#endif
