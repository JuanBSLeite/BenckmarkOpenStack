
/*
Copyright 2013 University of Warwick

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

/*! \file LauBelleSymNR.hh
    \brief File containing declaration of LauBelleSymNR class.
*/

/*! \class LauBelleSymNR
    \brief Class for defining the symmetric Belle Non Resonant model

    Defines the Non Resonant model from the Belle collaboration for the case
    of a symmetrical DP
*/

#ifndef LAU_BELLE_SYM_NR
#define LAU_BELLE_SYM_NR

#include "LauAbsResonance.hh"
#include "LauComplex.hh"

#include "TString.h"

class LauKinematics;

class LauBelleSymNR : public LauAbsResonance {

  public:
    //! Constructor
    /*!
        \param [in] resInfo the object containing information on the resonance name, mass, width, spin, charge, etc.
        \param [in] resType the model of the resonance
        \param [in] resPairAmpInt the number of the daughter not produced by the resonance
        \param [in] daughters the daughter particles
    */
    LauBelleSymNR( LauResonanceInfo* resInfo,
                   const LauAbsResonance::LauResonanceModel resType,
                   const Int_t resPairAmpInt,
                   const LauDaughters* daughters );

    //! Destructor
    virtual ~LauBelleSymNR();

    //! Initialise
    virtual void initialise();

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

    //! Is the amplitude pre-symmetrised?
    /*!
        \return kTRUE, indicating that the amplitude is already symmetrised
    */
    virtual Bool_t preSymmetrised() const { return kTRUE; }

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

  protected:
    //! Set the parameter alpha, the effective range
    /*!
        \param [in] alpha the new effective range parameter
    */
    void setAlpha( const Double_t alpha );

    //! Get the effective range parameter
    /*!
        \return the effective range parameter
    */
    Double_t getAlpha() const { return ( alpha_ != 0 ) ? alpha_->unblindValue() : 0.0; }

    //! See if the alpha parameter is fixed or floating
    Bool_t fixAlpha() const { return ( alpha_ != 0 ) ? alpha_->fixed() : kTRUE; }

    //! This is not called, amplitude is used directly instead
    virtual LauComplex resAmp( Double_t mass, Double_t spinTerm );

  private:
    //! Copy constructor (not implemented)
    LauBelleSymNR( const LauBelleSymNR& rhs );

    //! Copy assignment operator (not implemented)
    LauBelleSymNR& operator=( const LauBelleSymNR& rhs );

    //! The range parameter
    LauParameter* alpha_;

    //! The model to use
    LauAbsResonance::LauResonanceModel model_;

    ClassDef( LauBelleSymNR, 0 )    // Belle Non-resonant model
};

#endif
