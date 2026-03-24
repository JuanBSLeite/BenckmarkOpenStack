
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

/*! \file LauAbsIncohRes.hh
    \brief File containing declaration of LauAbsIncohRes class.
*/

/*! \class LauAbsIncohRes
    \brief Abstract class for defining incoherent resonant amplitude models
*/

#ifndef LAU_ABS_INCOH_RES
#define LAU_ABS_INCOH_RES

#include "LauAbsResonance.hh"

#include "TString.h"

class LauKinematics;

class LauAbsIncohRes : public LauAbsResonance {

  public:
    //! Constructor
    /*!
        \param [in] resInfo the object containing information on the resonance name, mass, width, spin, charge, etc.
        \param [in] resPairAmpInt the number of the daughter not produced by the resonance
        \param [in] daughters the daughter particles
    */
    LauAbsIncohRes( LauResonanceInfo* resInfo,
                    const Int_t resPairAmpInt,
                    const LauDaughters* daughters );

    //! Destructor
    virtual ~LauAbsIncohRes();

    //! Get intensity factor
    /*!
        \param [in] kinematics kinematic parameters of the parent and daughter particles
        \return the intensity factor
    */
    virtual Double_t intensityFactor( const LauKinematics* kinematics ) = 0;

  protected:
    //! Complex resonant amplitude
    /*!
        \param [in] mass appropriate invariant mass for the resonance
        \param [in] spinTerm spin term
    */
    virtual LauComplex resAmp( Double_t mass, Double_t spinTerm );

  private:
    //! Copy constructor (not implemented)
    LauAbsIncohRes( const LauAbsIncohRes& rhs );

    //! Copy assignment operator (not implemented)
    LauAbsIncohRes& operator=( const LauAbsIncohRes& rhs );

    ClassDef( LauAbsIncohRes, 0 )    // Abstract incoherent resonance class
};

#endif
