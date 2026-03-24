
/*
Copyright 2010 University of Warwick

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

/*! \file LauDabbaRes.hh
    \brief File containing declaration of LauDabbaRes class.
*/

/*! \class LauDabbaRes
    \brief Class for defining the Dabba resonance model

    Class for defining the Dabba resonance model
    Formulae and data values from arXiv:0901.2217 - author D.V.Bugg
*/

#ifndef LAU_DABBA_RES
#define LAU_DABBA_RES

#include "LauAbsResonance.hh"
#include "LauComplex.hh"

#include "TString.h"

class LauDabbaRes : public LauAbsResonance {

  public:
    //! Constructor
    /*!
        \param [in] resInfo the object containing information on the resonance name, mass, width, spin, charge, etc.
        \param [in] resPairAmpInt the number of the daughter not produced by the resonance
        \param [in] daughters the daughter particles
    */
    LauDabbaRes( LauResonanceInfo* resInfo, const Int_t resPairAmpInt, const LauDaughters* daughters );

    //! Destructor
    virtual ~LauDabbaRes();

    //! Initialise the model
    virtual void initialise();

    //! Get the resonance model type
    /*!
        \return the resonance model type
    */
    virtual LauAbsResonance::LauResonanceModel getResonanceModel() const
    {
        return LauAbsResonance::Dabba;
    }

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
    //! Set the b parameter
    /*!
        \param [in] b new value for b parameter
    */
    void setBValue( const Double_t b );

    //! Set the alpha parameter
    /*!
        \param [in] alpha new value for alpha parameter
    */
    void setAlphaValue( const Double_t alpha );

    //! Set the beta parameter
    /*!
        \param [in] beta new value for beta parameter
    */
    void setBetaValue( const Double_t beta );

    //! Get the b parameter value
    /*!
        \return value of the b parameter
    */
    Double_t getBValue() const { return ( b_ != 0 ) ? b_->unblindValue() : 0.0; }

    //! Get the alpha parameter value
    /*!
        \return value of the alpha parameter
    */
    Double_t getAlphaValue() const { return ( alpha_ != 0 ) ? alpha_->unblindValue() : 0.0; }

    //! Get the beta parameter value
    /*!
        \return value of the beta parameter
    */
    Double_t getBetaValue() const { return ( beta_ != 0 ) ? beta_->unblindValue() : 0.0; }

    //! Fix the b parameter value
    /*!
        \return kTRUE if the b parameter is fixed, kFALSE otherwise
    */
    Bool_t fixBValue() const { return ( b_ != 0 ) ? b_->fixed() : 0.0; }

    //! Fix the alpha parameter value
    /*!
        \return kTRUE if the alpha parameter is fixed, kFALSE otherwise
    */
    Bool_t fixAlphaValue() const { return ( alpha_ != 0 ) ? alpha_->fixed() : 0.0; }

    //! Fix the beta parameter value
    /*!
        \return kTRUE if the beta parameter is fixed, kFALSE otherwise
    */
    Bool_t fixBetaValue() const { return ( beta_ != 0 ) ? beta_->fixed() : 0.0; }

    //! Complex resonant amplitude
    /*!
        \param [in] mass appropriate invariant mass for the resonance
        \param [in] spinTerm spin term
    */
    virtual LauComplex resAmp( Double_t mass, Double_t spinTerm );

    //! Check that the daughter particles are D and pi
    void checkDaughterTypes() const;

  private:
    //! Copy constructor (not implemented)
    LauDabbaRes( const LauDabbaRes& rhs );

    //! Copy assignment operator (not implemented)
    LauDabbaRes& operator=( const LauDabbaRes& rhs );

    //! Defined as mD + mPi all squared
    Double_t mSumSq_;
    //! Defined as mD*mD - 0.5*mPi*mPi
    Double_t sAdler_;

    //! Constant factor
    LauParameter* b_;
    //! Constant factor
    LauParameter* alpha_;
    //! Constant factor
    LauParameter* beta_;

    ClassDef( LauDabbaRes, 0 )    // Dabba resonance model
};

#endif
