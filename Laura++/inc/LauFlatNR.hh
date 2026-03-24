
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

/*! \file LauFlatNR.hh
    \brief File containing declaration of LauFlatNR class.
*/

/*! \class LauFlatNR
    \brief Class for defining a uniform nonresonant amplitude
*/

#ifndef LAU_FLAT_NR
#define LAU_FLAT_NR

#include "LauAbsResonance.hh"
#include "LauComplex.hh"

#include "TString.h"

class LauKinematics;

class LauFlatNR : public LauAbsResonance {

  public:
    //! Constructor
    /*!
        \param [in] resInfo the object containing information on the resonance name, mass, width, spin, charge, etc.
        \param [in] resPairAmpInt the number of the daughter not produced by the resonance
        \param [in] daughters the daughter particles
    */
    LauFlatNR( LauResonanceInfo* resInfo, const Int_t resPairAmpInt, const LauDaughters* daughters );

    //! Destructor
    virtual ~LauFlatNR();

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
        return LauAbsResonance::FlatNR;
    }

  protected:
    //! This is not meant to be called
    virtual LauComplex resAmp( Double_t mass, Double_t spinTerm );

  private:
    //! Copy constructor (not implemented)
    LauFlatNR( const LauFlatNR& rhs );

    //! Copy assignment operator (not implemented)
    LauFlatNR& operator=( const LauFlatNR& rhs );

    ClassDef( LauFlatNR, 0 )    // Uniform non-resonant model
};

#endif
