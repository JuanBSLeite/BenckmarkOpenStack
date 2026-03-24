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

/*! \file LauRescatteringRes.hh
    \brief File containing declaration of LauRescatteringRes class.
*/

/*! \class LauRescatteringRes
    \brief Class for defining the rescattering model.

    Defines the Rescatering models from PiPi-KK Inelastic Scatering :
    2005: J.R. Pelaez, F. J. Ynduráin: PHYSICAL REVIEW D 71, 074016 (2005)
    2015: J. H. Alvarenga Nogueira, I. Bediaga, A. B. R. Cavalcante, T. Frederico, and O. Lourenço: PHYSICAL REVIEW D 92, 054010 (2015)
    2018: J.R. Pelaez,  A.Rodas: Unpublished yet PiPi -> KK scattering up to 1.47 GeV with hyperbolic dispersion relations.
*/

#ifndef LAU_RESCATTERING_RES
#define LAU_RESCATTERING_RES

#include "LauAbsResonance.hh"
#include "LauComplex.hh"

#include "TString.h"

class LauKinematics;
class LauParameter;

class LauRescatteringRes : public LauAbsResonance {

  public:
    //! Constructor
    /*!
        \param [in] resInfo the object containing information on the resonance name, mass, width, spin, charge, etc.
        \param [in] resType the model of the resonance
        \param [in] resPairAmpInt the number of the daughter not produced by the resonance
        \param [in] daughters the daughter particles
    */
    LauRescatteringRes( LauResonanceInfo* resInfo,
                        const LauAbsResonance::LauResonanceModel resType,
                        const Int_t resPairAmpInt,
                        const LauDaughters* daughters );

    //! Destructor
    virtual ~LauRescatteringRes();

    //! Initialise the model
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
    //! Set the parameter lambdaPiPi, the term for the PiPi
    /*!
        \param [in] lambda the term for the PiPi
    */
    void setLambdaPiPi( const Double_t lambda );

    //! Get the lambdaPiPi, the term for the PiPi
    /*!
        \return lambdaPiPi, the term for the PiPi
    */
    Double_t getLambdaPiPi() const { return ( lambdaPiPi_ != 0 ) ? lambdaPiPi_->value() : 0.0; }

    //! See if the lambdaPiPi parameter is fixed or floating
    /*!
        \return kTRUE if the lambdaPiPi  parameter is fixed, kFALSE otherwise
    */
    Bool_t fixLambdaPiPi() const { return ( lambdaPiPi_ != 0 ) ? lambdaPiPi_->fixed() : kTRUE; }

    //! Set the parameter lambdaKK, the term for the KK
    /*!
        \param [in] lambda the term for the KK
    */
    void setLambdaKK( const Double_t lambda );

    //! Get the lambdaKK, the term for the KK
    /*!
        \return lambdaKK, the term for the KK
    */
    Double_t getLambdaKK() const { return ( lambdaKK_ != 0 ) ? lambdaKK_->value() : 0.0; }

    //! See if the lambdaKK parameter is fixed or floating
    /*!
        \return kTRUE if the lambdaKK  parameter is fixed, kFALSE otherwise
    */
    Bool_t fixLambdaKK() const { return ( lambdaKK_ != 0 ) ? lambdaKK_->fixed() : kTRUE; }

    //! Set the parameter Ms
    /*!
        \param [in] Ms
    */
    void setMs( const Double_t Ms );

    //! Get the Ms
    /*!
        \return the Ms
    */
    Double_t getMs() const { return ( Ms_ != 0 ) ? Ms_->value() : 0.0; }

    //! See if the Ms parameter is fixed or floating
    /*!
        \return kTRUE if the Ms  parameter is fixed, kFALSE otherwise
    */
    Bool_t fixMs() const { return ( Ms_ != 0 ) ? Ms_->fixed() : kTRUE; }

    //! Set the parameter Mf
    /*!
        \param [in] Mf
    */
    void setMf( const Double_t Mf );

    //! Get the Mf
    /*!
        \return the Mf
    */
    Double_t getMf() const { return ( Mf_ != 0 ) ? Mf_->value() : 0.0; }

    //! See if the Mf parameter is fixed or floating
    /*!
        \return kTRUE if the Mf  parameter is fixed, kFALSE otherwise
    */
    Bool_t fixMf() const { return ( Mf_ != 0 ) ? Mf_->fixed() : kTRUE; }

    //! Set the parameter Mprime
    /*!
        \param [in] Mprime
    */
    void setMprime( const Double_t Mprime );

    //! Get the Mprime
    /*!
        \return the Mprime
    */
    Double_t getMprime() const { return ( Mprime_ != 0 ) ? Mprime_->value() : 0.0; }

    //! See if the Mprime parameter is fixed or floating
    /*!
        \return kTRUE if the Mprime  parameter is fixed, kFALSE otherwise
    */
    Bool_t fixMprime() const { return ( Mprime_ != 0 ) ? Mprime_->fixed() : kTRUE; }

    //! Set the parameter Eps1
    /*!
        \param [in] Eps1
    */
    void setEps1( const Double_t Eps1 );

    //! Get the Eps1
    /*!
        \return the Eps1
    */
    Double_t getEps1() const { return ( Eps1_ != 0 ) ? Eps1_->value() : 0.0; }

    //! See if the Eps1 parameter is fixed or floating
    /*!
        \return kTRUE if the Eps1  parameter is fixed, kFALSE otherwise
    */
    Bool_t fixEps1() const { return ( Eps1_ != 0 ) ? Eps1_->fixed() : kTRUE; }

    //! Set the parameter Eps2
    /*!
        \param [in] Eps2
    */
    void setEps2( const Double_t Eps2 );

    //! Get the Eps2
    /*!
        \return the Eps2
    */
    Double_t getEps2() const { return ( Eps2_ != 0 ) ? Eps2_->value() : 0.0; }

    //! See if the Eps2 parameter is fixed or floating
    /*!
        \return kTRUE if the Eps2  parameter is fixed, kFALSE otherwise
    */
    Bool_t fixEps2() const { return ( Eps2_ != 0 ) ? Eps2_->fixed() : kTRUE; }

    //! Set the parameter C0
    /*!
        \param [in] C0
    */
    void setC0( const Double_t C0 );

    //! Get the C0
    /*!
        \return the C0
    */
    Double_t getC0() const { return ( C0_ != 0 ) ? C0_->value() : 0.0; }

    //! See if the C0 parameter is fixed or floating
    /*!
        \return kTRUE if the C0  parameter is fixed, kFALSE otherwise
    */
    Bool_t fixC0() const { return ( C0_ != 0 ) ? C0_->fixed() : kTRUE; }

    //! Complex resonant amplitude
    /*!
        \param [in] mass appropriate invariant mass for the resonance
        \param [in] spinTerm Zemach spin term
    */
    virtual LauComplex resAmp( Double_t mass, Double_t spinTerm );

  private:
    //! Copy constructor (not implemented)
    LauRescatteringRes( const LauRescatteringRes& rhs );

    //! Copy assignment operator (not implemented)
    LauRescatteringRes& operator=( const LauRescatteringRes& rhs );

    //! the term for the PiPi
    LauParameter* lambdaPiPi_;

    //! the term for the KK
    LauParameter* lambdaKK_;

    //! the term for the Mf_
    LauParameter* Mf_;

    //! the term for the Ms
    LauParameter* Ms_;

    //! the term for the Mprime
    LauParameter* Mprime_;

    //! the term for the Eps1
    LauParameter* Eps1_;

    //! the term for the Eps2
    LauParameter* Eps2_;

    //! the term for the C0
    LauParameter* C0_;

    //! The model to use
    LauAbsResonance::LauResonanceModel model_;

    ClassDef( LauRescatteringRes, 0 )
};

#endif
