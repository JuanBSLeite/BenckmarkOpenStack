/*
Copyright 2018 University of Warwick

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

/*! \file LauRescattering2Res.hh
    \brief File containing declaration of LauRescattering2Res class.
*/

/*! \class LauRescattering2Res
    \brief Class for defining an alternative rescattering model.

    Model for pipi SWave proposed by J.Schechter as used by Cleo (PRD76,01200(2007)
    to replace the sigma and f0(980) contributions.
*/

#ifndef LAU_RESCATTERING2_RES
#define LAU_RESCATTERING2_RES

#include "LauAbsResonance.hh"
#include "LauComplex.hh"

#include "TString.h"

class LauKinematics;

class LauRescattering2Res : public LauAbsResonance {

  public:
    //! Constructor
    /*!
        \param [in] resInfo the object containing information on the resonance name, mass, width, spin, charge, etc.
        \param [in] resPairAmpInt the number of the daughter not produced by the resonance
        \param [in] daughters the daughter particles
    */
    LauRescattering2Res( LauResonanceInfo* resInfo,
                         const Int_t resPairAmpInt,
                         const LauDaughters* daughters );

    //! Destructor
    virtual ~LauRescattering2Res();

    //! Initialise the model
    virtual void initialise();

    //! Get the resonance model type
    /*!
        \return the resonance model type
    */
    virtual LauAbsResonance::LauResonanceModel getResonanceModel() const
    {
        return LauAbsResonance::Rescattering2;
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

    Double_t pn( const Double_t x, const Double_t n ) const;
    Double_t x( const Double_t sqr_t, const Int_t i ) const;
    Double_t phi00( const Double_t sqr_t, const Int_t i ) const;
    Double_t g00( const Double_t sqr_t, const Int_t i ) const;

    void setB1Parameter( const Double_t B1 );
    void setB2Parameter( const Double_t B2 );
    void setB3Parameter( const Double_t B3 );
    void setC1Parameter( const Double_t C1 );
    void setC2Parameter( const Double_t C2 );
    void setC3Parameter( const Double_t C3 );
    void setC4Parameter( const Double_t C4 );
    void setC5Parameter( const Double_t C5 );
    void setD0Parameter( const Double_t D0 );
    void setD1Parameter( const Double_t D1 );
    void setD2Parameter( const Double_t D2 );
    void setD3Parameter( const Double_t D3 );
    void setF1Parameter( const Double_t F1 );
    void setF2Parameter( const Double_t F2 );
    void setF3Parameter( const Double_t F3 );
    void setF4Parameter( const Double_t F4 );

    Double_t getB1Parameter() const { return ( B1_ != 0 ) ? B1_->value() : 0.0; }

    Bool_t fixB1Parameter() const { return ( B1_ != 0 ) ? B1_->fixed() : kTRUE; }

    Double_t getB2Parameter() const { return ( B2_ != 0 ) ? B2_->value() : 0.0; }

    Bool_t fixB2Parameter() const { return ( B2_ != 0 ) ? B2_->fixed() : kTRUE; }

    Double_t getB3Parameter() const { return ( B3_ != 0 ) ? B3_->value() : 0.0; }

    Bool_t fixB3Parameter() const { return ( B3_ != 0 ) ? B3_->fixed() : kTRUE; }

    Double_t getC1Parameter() const { return ( C1_ != 0 ) ? C1_->value() : 0.0; }

    Bool_t fixC1Parameter() const { return ( C1_ != 0 ) ? C1_->fixed() : kTRUE; }

    Double_t getC2Parameter() const { return ( C2_ != 0 ) ? C2_->value() : 0.0; }

    Bool_t fixC2Parameter() const { return ( C2_ != 0 ) ? C2_->fixed() : kTRUE; }

    Double_t getC3Parameter() const { return ( C3_ != 0 ) ? C3_->value() : 0.0; }

    Bool_t fixC3Parameter() const { return ( C3_ != 0 ) ? C3_->fixed() : kTRUE; }

    Double_t getC4Parameter() const { return ( C4_ != 0 ) ? C4_->value() : 0.0; }

    Bool_t fixC4Parameter() const { return ( C4_ != 0 ) ? C4_->fixed() : kTRUE; }

    Double_t getC5Parameter() const { return ( C5_ != 0 ) ? C5_->value() : 0.0; }

    Bool_t fixC5Parameter() const { return ( C5_ != 0 ) ? C5_->fixed() : kTRUE; }

    Double_t getD0Parameter() const { return ( D0_ != 0 ) ? D0_->value() : 0.0; }

    Bool_t fixD0Parameter() const { return ( D0_ != 0 ) ? D0_->fixed() : kTRUE; }

    Double_t getD1Parameter() const { return ( D1_ != 0 ) ? D1_->value() : 0.0; }

    Bool_t fixD1Parameter() const { return ( D1_ != 0 ) ? D1_->fixed() : kTRUE; }

    Double_t getD2Parameter() const { return ( D2_ != 0 ) ? D2_->value() : 0.0; }

    Bool_t fixD2Parameter() const { return ( D2_ != 0 ) ? D2_->fixed() : kTRUE; }

    Double_t getD3Parameter() const { return ( D3_ != 0 ) ? D3_->value() : 0.0; }

    Bool_t fixD3Parameter() const { return ( D3_ != 0 ) ? D3_->fixed() : kTRUE; }

    Double_t getF1Parameter() const { return ( F1_ != 0 ) ? F1_->value() : 0.0; }

    Bool_t fixF1Parameter() const { return ( F1_ != 0 ) ? F1_->fixed() : kTRUE; }

    Double_t getF2Parameter() const { return ( F2_ != 0 ) ? F2_->value() : 0.0; }

    Bool_t fixF2Parameter() const { return ( F2_ != 0 ) ? F2_->fixed() : kTRUE; }

    Double_t getF3Parameter() const { return ( F3_ != 0 ) ? F3_->value() : 0.0; }

    Bool_t fixF3Parameter() const { return ( F3_ != 0 ) ? F3_->fixed() : kTRUE; }

    Double_t getF4Parameter() const { return ( F4_ != 0 ) ? F4_->value() : 0.0; }

    Bool_t fixF4Parameter() const { return ( F4_ != 0 ) ? F4_->fixed() : kTRUE; }

  private:
    //! Copy constructor (not implemented)
    LauRescattering2Res( const LauRescattering2Res& rhs );

    //! Copy assignment operator (not implemented)
    LauRescattering2Res& operator=( const LauRescattering2Res& rhs );

    //! Parameter
    LauParameter* B1_;
    LauParameter* B2_;
    LauParameter* B3_;
    LauParameter* C1_;
    LauParameter* C2_;
    LauParameter* C3_;
    LauParameter* C4_;
    LauParameter* C5_;
    LauParameter* D0_;
    LauParameter* D1_;
    LauParameter* D2_;
    LauParameter* D3_;
    LauParameter* F1_;
    LauParameter* F2_;
    LauParameter* F3_;
    LauParameter* F4_;

    Double_t sqr_tmin[3], sqr_tmax[3];
    Double_t B0_, C0_, F0_;

    ClassDef( LauRescattering2Res, 0 )    // pipi S wave model by Schechter amplitude model
};

#endif
