
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

/*! \file LauBelleNR.hh
    \brief File containing declaration of LauBelleNR class.
*/

/*! \class LauBelleNR
    \brief Class for defining the Belle nonresonant model.

    Defines the nonresonant model from the Belle collaboration.
    arXiv:hep-ex/0412066
*/

#ifndef LAU_BELLE_NR
#define LAU_BELLE_NR

#include "LauAbsResonance.hh"
#include "LauComplex.hh"

#include "TString.h"

class LauKinematics;
class LauParameter;

class LauBelleNR : public LauAbsResonance {

  public:
    //! Constructor
    /*!
        \param [in] resInfo the object containing information on the resonance name, mass, width, spin, charge, etc.
        \param [in] resType the model of the resonance
        \param [in] resPairAmpInt the number of the daughter not produced by the resonance
        \param [in] daughters the daughter particles
    */
    LauBelleNR( LauResonanceInfo* resInfo,
                const LauAbsResonance::LauResonanceModel resType,
                const Int_t resPairAmpInt,
                const LauDaughters* daughters );

    //! Destructor
    virtual ~LauBelleNR();

    //! Override the enforcement of pure Legendre polynomial spin factors
    /*!
        By default this model uses pure Legendre polynomial
        spin factors, regardless of the default type set in
        LauResonanceMaker or any specific request from the user.
        This function allows the enforcement to be overridden.

        \param [in] forceLegendre boolean flag (kTRUE, the default, implies enforcement of pure Legendre spin factors, kFALSE overrides this to allow use of other formalisms)
    */
    void enforceLegendreSpinFactors( const Bool_t forceLegendre )
    {
        forceLegendre_ = forceLegendre;
    }

    //! Initialise the model
    virtual void initialise();

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

  protected:
    //! Set the parameter alpha, the effective range
    /*!
        \param [in] alpha the new effective range parameter value
    */
    void setAlpha( const Double_t alpha );

    //! Get the effective range parameter
    /*!
        \return the effective range parameter
    */
    Double_t getAlpha() const { return ( alpha_ != 0 ) ? alpha_->unblindValue() : 0.0; }

    //! See if the alpha parameter is fixed or floating
    /*!
        \return kTRUE if the effective range parameter is fixed, kFALSE otherwise
    */
    Bool_t fixAlpha() const { return ( alpha_ != 0 ) ? alpha_->fixed() : kTRUE; }

    //! Complex resonant amplitude
    /*!
        \param [in] mass appropriate invariant mass for the resonance
        \param [in] spinTerm spin term
    */
    virtual LauComplex resAmp( Double_t mass, Double_t spinTerm );

  private:
    //! Copy constructor (not implemented)
    LauBelleNR( const LauBelleNR& rhs );

    //! Copy assignment operator (not implemented)
    LauBelleNR& operator=( const LauBelleNR& rhs );

    //! The range parameter
    LauParameter* alpha_;

    //! The model to use
    LauAbsResonance::LauResonanceModel model_;

    //! Force use of Legendre spin factors
    Bool_t forceLegendre_;

    ClassDef( LauBelleNR, 0 )
};

#endif
