
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

/*! \file LauKappaRes.hh
    \brief File containing declaration of LauKappaRes class.
*/

/*! \class LauKappaRes
    \brief Class for defining the Kappa resonance model

    Class for defining the Kappa resonance model.
    Formulae and data values from
    Phys.Lett.B 572, 1 (2003) - author D.V.Bugg
*/

#ifndef LAU_KAPPA_RES
#define LAU_KAPPA_RES

#include "LauAbsResonance.hh"
#include "LauComplex.hh"

#include "TString.h"

class LauKappaRes : public LauAbsResonance {

  public:
    //! Constructor
    /*!
        \param [in] resInfo the object containing information on the resonance name, mass, width, spin, charge, etc.
        \param [in] resPairAmpInt the number of the daughter not produced by the resonance
        \param [in] daughters the daughter particles
    */
    LauKappaRes( LauResonanceInfo* resInfo, const Int_t resPairAmpInt, const LauDaughters* daughters );

    //! Destructor
    virtual ~LauKappaRes();

    //! Initialise the model
    virtual void initialise();

    //! Get the resonance model type
    /*!
        \return the resonance model type
    */
    virtual LauAbsResonance::LauResonanceModel getResonanceModel() const
    {
        return LauAbsResonance::Kappa;
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
    //! Set the b1 parameter
    /*!
        \param [in] b1 new value for b1 parameter
    */
    void setB1Value( const Double_t b1 );

    //! Set the b2 parameter
    /*!
        \param [in] b2 new value for b2 parameter
    */
    void setB2Value( const Double_t b2 );

    //! Set the A parameter
    /*!
        \param [in] A new value for A parameter
    */
    void setAValue( const Double_t A );

    //! Set the m0 parameter
    /*!
        \param [in] m0 new value for m0 parameter
    */
    void setM0Value( const Double_t m0 );

    //! Get the b1 parameter value
    /*!
        \return value of the b1 parameter
    */
    Double_t getB1Value() const { return ( b1_ != 0 ) ? b1_->unblindValue() : 0.0; }

    //! Get the b2 parameter value
    /*!
        \return value of the b2 parameter
    */
    Double_t getB2Value() const { return ( b2_ != 0 ) ? b2_->unblindValue() : 0.0; }

    //! Get the A parameter value
    /*!
        \return value of the A parameter
    */
    Double_t getAValue() const { return ( a_ != 0 ) ? a_->unblindValue() : 0.0; }

    //! Get the m0 parameter value
    /*!
        \return value of the m0 parameter
    */
    Double_t getM0Value() const { return ( m0_ != 0 ) ? m0_->unblindValue() : 0.0; }

    //! Fix the b1 parameter value
    /*!
        \return kTRUE if the b1 parameter is fixed, kFALSE otherwise
    */
    Bool_t fixB1Value() const { return ( b1_ != 0 ) ? b1_->fixed() : 0.0; }

    //! Fix the b2 parameter value
    /*!
        \return kTRUE if the b2 parameter is fixed, kFALSE otherwise
    */
    Bool_t fixB2Value() const { return ( b2_ != 0 ) ? b2_->fixed() : 0.0; }

    //! Fix the A parameter value
    /*!
        \return kTRUE if the A parameter is fixed, kFALSE otherwise
    */
    Bool_t fixAValue() const { return ( a_ != 0 ) ? a_->fixed() : 0.0; }

    //! Fix the m0 parameter value
    /*!
        \return kTRUE if the m0 parameter is fixed, kFALSE otherwise
    */
    Bool_t fixM0Value() const { return ( m0_ != 0 ) ? m0_->fixed() : 0.0; }

    //! Complex resonant amplitude
    /*!
        \param [in] mass appropriate invariant mass for the resonance
        \param [in] spinTerm spin term
    */
    virtual LauComplex resAmp( Double_t mass, Double_t spinTerm );

    //! Check that Kappa daughters are K and pi
    void checkDaughterTypes() const;

  private:
    //! Copy constructor (not implemented)
    LauKappaRes( const LauKappaRes& rhs );

    //! Copy assignment operator (not implemented)
    LauKappaRes& operator=( const LauKappaRes& rhs );

    //! Square of (mK + mPi)
    const Double_t mSumSq_;
    //! Defined as mK*mK - 0.5*mPi*mPi
    const Double_t sAdler_;

    //! Factor from BES data
    LauParameter* b1_;
    //! Factor from BES data
    LauParameter* b2_;
    //! Factor from BES data
    LauParameter* a_;
    //! Factor from BES data
    LauParameter* m0_;

    ClassDef( LauKappaRes, 0 )    // Kappa resonance model
};

#endif
