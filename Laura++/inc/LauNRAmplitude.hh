
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

/*! \file LauNRAmplitude.hh
    \brief File containing declaration of LauNRAmplitude class.
*/

/*! \class LauNRAmplitude
    \brief Class for defining the NR amplitude model.

    Class for defining the NR amplitude model.
    Formulae and data values from arXiv:0709.0075v1 [hep-ph].
*/

#ifndef LAU_NR_AMPLITUDE
#define LAU_NR_AMPLITUDE

#include "LauAbsResonance.hh"
#include "LauComplex.hh"

#include "TMath.h"
#include "TString.h"

class LauKinematics;

class LauNRAmplitude : public LauAbsResonance {

  public:
    //! Constructor
    /*!
        \param [in] resInfo the object containing information on the resonance name, mass, width, spin, charge, etc.
        \param [in] resPairAmpInt the number of the daughter not produced by the resonance
        \param [in] daughters the daughter particles
    */
    LauNRAmplitude( LauResonanceInfo* resInfo,
                    const Int_t resPairAmpInt,
                    const LauDaughters* daughters );

    //! Destructor
    virtual ~LauNRAmplitude();

    //! Initialise the model
    virtual void initialise();

    //! Complex resonant amplitude
    /*!
        \param [in] kinematics kinematic parameters of the parent and daughter particles
        \return the complex amplitude
    */
    virtual LauComplex amplitude( const LauKinematics* kinematics );

    //! Get the resonance model type
    /*!
        \return the resonance model type
    */
    virtual LauAbsResonance::LauResonanceModel getResonanceModel() const
    {
        return LauAbsResonance::NRModel;
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
    //! This is not meant to be called
    virtual LauComplex resAmp( Double_t mass, Double_t spinTerm );

    //! Evaluate the expression 1.0 / (1.0 + TMath::Exp( c * (s-p) ))
    /*!
        \param [in] s invariant mass squared
        \param [in] c parameter from arXiv:0709.0075v1 [hep-ph]
        \param [in] p parameter from arXiv:0709.0075v1 [hep-ph]
        \return the value of the expression 1.0 / (1.0 + TMath::Exp( c * (s-p) ))
    */
    Double_t f( const Double_t s, const Double_t c, const Double_t p ) const;

    //! Set the d parameter value
    /*!
        \param [in] d the new d parameter value
    */
    void setdParameter( const Double_t d );

    //! Set the c1 parameter value
    /*!
        \param [in] c1 the new c1 parameter value
    */
    void setc1Parameter( const Double_t c1 );

    //! Set the c2 parameter value
    /*!
        \param [in] c2 the new c2 parameter value
    */
    void setc2Parameter( const Double_t c2 );

    //! Set the p1 parameter value
    /*!
        \param [in] p1 the new p1 parameter value
    */
    void setp1Parameter( const Double_t p1 );

    //! Set the p2 parameter value
    /*!
        \param [in] p2 the new p2 parameter value
    */
    void setp2Parameter( const Double_t p2 );

    //! Get the d parameter value
    /*!
        \return the new d parameter value
    */
    Double_t getdParameter() const { return ( d_ != 0 ) ? d_->unblindValue() : 0.0; }

    //! Get the c1 parameter value
    /*!
        \return the new c1 parameter value
    */
    Double_t getc1Parameter() const { return ( c1_ != 0 ) ? c1_->unblindValue() : 0.0; }

    //! Get the c2 parameter value
    /*!
        \return the new c2 parameter value
    */
    Double_t getc2Parameter() const { return ( c2_ != 0 ) ? c2_->unblindValue() : 0.0; }

    //! Get the p1 parameter value
    /*!
        \return the new p1 parameter value
    */
    Double_t getp1Parameter() const { return ( p1_ != 0 ) ? p1_->unblindValue() : 0.0; }

    //! Get the p2 parameter value
    /*!
        \return the new p2 parameter value
    */
    Double_t getp2Parameter() const { return ( p2_ != 0 ) ? p2_->unblindValue() : 0.0; }

    //! See if the d parameter is fixed or floating
    /*!
        \return kTRUE if the d parameter is fixed, kFALSE otherwise
    */
    Bool_t fixdParameter() const { return ( d_ != 0 ) ? d_->fixed() : kTRUE; }

    //! See if the c1 parameter is fixed or floating
    /*!
        \return kTRUE if the c1 parameter is fixed, kFALSE otherwise
    */
    Bool_t fixc1Parameter() const { return ( c1_ != 0 ) ? c1_->fixed() : kTRUE; }

    //! See if the c2 parameter is fixed or floating
    /*!
        \return kTRUE if the c2 parameter is fixed, kFALSE otherwise
    */
    Bool_t fixc2Parameter() const { return ( c2_ != 0 ) ? c2_->fixed() : kTRUE; }

    //! See if the p1 parameter is fixed or floating
    /*!
        \return kTRUE if the p1 parameter is fixed, kFALSE otherwise
    */
    Bool_t fixp1Parameter() const { return ( p1_ != 0 ) ? p1_->fixed() : kTRUE; }

    //! See if the p2 parameter is fixed or floating
    /*!
        \return kTRUE if the p2 parameter is fixed, kFALSE otherwise
    */
    Bool_t fixp2Parameter() const { return ( p2_ != 0 ) ? p2_->fixed() : kTRUE; }

  private:
    //! Copy constructor (not implemented)
    LauNRAmplitude( const LauNRAmplitude& rhs );

    //! Copy assignment operator (not implemented)
    LauNRAmplitude& operator=( const LauNRAmplitude& rhs );

    //! Parameter from arXiv:0709.0075v1 [hep-ph]
    LauParameter* d_;
    //! Parameter from arXiv:0709.0075v1 [hep-ph]
    LauParameter* c1_;
    //! Parameter from arXiv:0709.0075v1 [hep-ph]
    LauParameter* c2_;
    //! Parameter from arXiv:0709.0075v1 [hep-ph]
    LauParameter* p1_;
    //! Parameter from arXiv:0709.0075v1 [hep-ph]
    LauParameter* p2_;

    ClassDef( LauNRAmplitude, 0 )    // Non-resonant amplitude model
};

#endif
