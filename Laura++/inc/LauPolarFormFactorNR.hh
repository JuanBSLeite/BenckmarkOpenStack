
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

/*! \file LauPolarFormFactorNR.hh
    \brief File containing declaration of LauPolarFormFactorNR class.
*/

/*! \class LauPolarFormFactorNR
    \brief Class for defining a nonresonant form factor model

    Defines the nonresonant form factor model from:
    Nogueira, Bediaga, Cavalcante, Frederico, Lorenco: Phys. Rev. D92 (2015) 054010, arXiv:1506.08332 [hep-ph]
    Pelaez, Yndurain: Phys. Rev. D71 (2005) 074016, arXiv:hep-ph/0411334
*/

#ifndef LAU_POLAR_FORM_FACTOR_NR
#define LAU_POLAR_FORM_FACTOR_NR

#include "LauAbsResonance.hh"
#include "LauComplex.hh"

#include "TString.h"

class LauKinematics;
class LauParameter;

class LauPolarFormFactorNR : public LauAbsResonance {

  public:
    //! Constructor
    /*!
        \param [in] resInfo the object containing information on the resonance name, mass, width, spin, charge, etc.
        \param [in] resType the model of the resonance
        \param [in] resPairAmpInt the number of the daughter not produced by the resonance
        \param [in] daughters the daughter particles
    */
    LauPolarFormFactorNR( LauResonanceInfo* resInfo,
                          const LauAbsResonance::LauResonanceModel resType,
                          const Int_t resPairAmpInt,
                          const LauDaughters* daughters );

    //! Destructor
    virtual ~LauPolarFormFactorNR();

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
    //! Set the parameter lambda, the NR shape parameter
    /*!
        \param [in] lambda the NR shape parameter
    */
    void setLambda( const Double_t lambda );

    //! Get the parameter lambda, the NR shape parameter
    /*!
        \return lambda, the NR shape parameter
    */
    Double_t getLambda() const { return ( lambda_ != 0 ) ? lambda_->value() : 0.0; }

    //! See if the lambda parameter is fixed or floating
    /*!
        \return kTRUE if the lambda  parameter is fixed, kFALSE otherwise
    */
    Bool_t fixLambda() const { return ( lambda_ != 0 ) ? lambda_->fixed() : kTRUE; }

    //! Complex resonant amplitude
    /*!
        \param [in] mass appropriate invariant mass for the resonance
        \param [in] spinTerm Zemach spin term
    */
    virtual LauComplex resAmp( Double_t mass, Double_t spinTerm );

  private:
    //! Copy constructor (not implemented)
    LauPolarFormFactorNR( const LauPolarFormFactorNR& rhs );

    //! Copy assignment operator (not implemented)
    LauPolarFormFactorNR& operator=( const LauPolarFormFactorNR& rhs );

    //! The NR shape parameter
    LauParameter* lambda_;

    //! The model to use
    LauAbsResonance::LauResonanceModel model_;

    ClassDef( LauPolarFormFactorNR, 0 )
};

#endif
