
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

/*! \file LauFlatteRes.hh
    \brief File containing declaration of LauFlatteRes class.
*/

/*! \class LauFlatteRes
    \brief Class for defining the Flatte resonance model

    Class for defining the Flatte resonance model.
    For use with the f_0(980) resonance.
*/

#ifndef LAU_FLATTE_RES
#define LAU_FLATTE_RES

#include "LauAbsResonance.hh"
#include "LauComplex.hh"

#include "TString.h"

class LauFlatteRes : public LauAbsResonance {

  public:
    //! Constructor
    /*!
        \param [in] resInfo the object containing information on the resonance name, mass, width, spin, charge, etc.
        \param [in] resPairAmpInt the number of the daughter not produced by the resonance
        \param [in] daughters the daughter particles
    */
    LauFlatteRes( LauResonanceInfo* resInfo, const Int_t resPairAmpInt, const LauDaughters* daughters );

    //! Destructor
    virtual ~LauFlatteRes();

    //! Initialise the model
    virtual void initialise();

    //! Get the resonance model type
    /*!
        \return the resonance model type
    */
    virtual LauAbsResonance::LauResonanceModel getResonanceModel() const
    {
        return LauAbsResonance::Flatte;
    }

    //! Set value of a resonance parameter
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
    //! Set the g1 parameter
    /*!
        \param [in] g1 constant factor
    */
    void setg1Parameter( const Double_t g1 );

    //! Set the g2 parameter
    /*!
        \param [in] g2 constant factor
    */
    void setg2Parameter( const Double_t g2 );

    //! Get the g1 parameter
    /*!
        \return constant factor g1
    */
    Double_t getg1Parameter() const { return ( g1_ != 0 ) ? g1_->unblindValue() : 0.0; }

    //! Get the g2 parameter
    /*!
        \return constant factor g2
    */
    Double_t getg2Parameter() const { return ( g2_ != 0 ) ? g2_->unblindValue() : 0.0; }

    //! See if the g1 parameter is fixed or floating
    /*!
        \return kTRUE if the g1 parameter is fixed, kFALSE otherwise
    */
    Double_t fixg1Parameter() const { return ( g1_ != 0 ) ? g1_->fixed() : kTRUE; }

    //! See if the g2 parameter is fixed or floating
    /*!
        \return kTRUE if the g2 parameter is fixed, kFALSE otherwise
    */
    Double_t fixg2Parameter() const { return ( g2_ != 0 ) ? g2_->fixed() : kTRUE; }

    //! Complex resonant amplitude
    /*!
        \param [in] mass appropriate invariant mass for the resonance
        \param [in] spinTerm spin term
    */
    virtual LauComplex resAmp( Double_t mass, Double_t spinTerm );

  private:
    //! Copy constructor (not implemented)
    LauFlatteRes( const LauFlatteRes& rhs );

    //! Copy assignment operator (not implemented)
    LauFlatteRes& operator=( const LauFlatteRes& rhs );

    //! Channel 1 coupling parameter
    LauParameter* g1_;
    //! Channel 1 coupling parameter
    LauParameter* g2_;

    //! Channel 1, subchannel 1 invariant mass
    Double_t mSumSq0_;
    //! Channel 1, subchannel 2 invariant mass
    Double_t mSumSq1_;
    //! Channel 2, subchannel 1 invariant mass
    Double_t mSumSq2_;
    //! Channel 2, subchannel 2 invariant mass
    Double_t mSumSq3_;

    //! Flag to turn on Adler term in the width
    Bool_t useAdlerTerm_;

    //! The Adler zero
    Double_t sA_;

    //! Flag to specify whether the couplings absorb the m_0 factor
    Bool_t absorbM0_;

    ClassDef( LauFlatteRes, 0 )
};

#endif
