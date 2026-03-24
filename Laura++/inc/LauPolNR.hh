
/*
Copyright 2014 University of Warwick

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

/*! \file LauPolNR.hh
    \brief File containing declaration of LauPolNR class.
*/

/*! \class LauPolNR
    \brief Class for defining the terms of Babar nonresonant model.

    Defines the nonresonant model from the Babar collaboration.
    arXiv:1201.5897
*/

#ifndef LAU_POL_NR
#define LAU_POL_NR

#include "LauAbsResonance.hh"
#include "LauComplex.hh"

#include "TString.h"

class LauKinematics;

class LauPolNR : public LauAbsResonance {

  public:
    //! Constructor
    /*!
        \param [in] resInfo the object containing information on the resonance name, mass, width, spin, charge, etc.
        \param [in] resPairAmpInt the number of the daughter not produced by the resonance
        \param [in] daughters the daughter particles
    */
    LauPolNR( LauResonanceInfo* resInfo, const Int_t resPairAmpInt, const LauDaughters* daughters );

    //! Destructor
    virtual ~LauPolNR();

    //! Initialise the model
    virtual void initialise();

    //! Get the resonance model type
    /*!
        \return the resonance model type
    */
    virtual LauAbsResonance::LauResonanceModel getResonanceModel() const
    {
        return LauAbsResonance::PolNR;
    }

    //! Set the parameter omega, the offset parameter
    /*!
        \param [in] omega the new offset parameter
    */
    virtual void setOmega( const Double_t omega ) { omega_ = omega; }

    //! Get the offset parameter
    /*!
        \return the offset parameter
    */
    virtual Double_t getOmega() const { return omega_; }

  protected:
    //! Complex resonant amplitude
    /*!
        \param [in] mass appropriate invariant mass for the resonance
        \param [in] spinTerm spin term
    */
    virtual LauComplex resAmp( Double_t mass, Double_t spinTerm );

  private:
    //! Copy constructor (not implemented)
    LauPolNR( const LauPolNR& rhs );

    //! Copy assignment operator (not implemented)
    LauPolNR& operator=( const LauPolNR& rhs );

    //! The offset parameter
    Double_t omega_;

    //! The order in the mass
    UInt_t order_;

    ClassDef( LauPolNR, 0 )    // Polynomial Non-resonant model
};

#endif
