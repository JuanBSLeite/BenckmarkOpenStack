
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

/*! \file LauGaussIncohRes.hh
    \brief File containing declaration of LauGaussIncohRes class.
*/

/*! \class LauGaussIncohRes
    \brief Class for defining an incoherent resonance with a Gaussian mass dependence

    Useful for modelling resolution-dominated sharp structures such as D*(2010)
*/

#ifndef LAU_GAUSS_INCOH_RES
#define LAU_GAUSS_INCOH_RES

#include "LauAbsIncohRes.hh"
#include "LauComplex.hh"

#include "TString.h"

class LauKinematics;

class LauGaussIncohRes : public LauAbsIncohRes {

  public:
    //! Constructor
    /*!
        \param [in] resInfo the object containing information on the resonance name, mass, width, spin, charge, etc.
        \param [in] resPairAmpInt the number of the daughter not produced by the resonance
        \param [in] daughters the daughter particles
    */
    LauGaussIncohRes( LauResonanceInfo* resInfo,
                      const Int_t resPairAmpInt,
                      const LauDaughters* daughters );

    //! Destructor
    virtual ~LauGaussIncohRes();

    //! Initialise the model
    virtual void initialise();

    //! Get intensity factor
    /*!
        \param [in] kinematics kinematic parameters of the parent and daughter particles
        \return the intensity factor
    */
    virtual Double_t intensityFactor( const LauKinematics* kinematics );

    //! Get the resonance model type
    /*!
        \return the resonance model type
    */
    virtual LauAbsResonance::LauResonanceModel getResonanceModel() const
    {
        return LauAbsResonance::GaussIncoh;
    }

    //! Retrieve the resonance parameters, e.g. so that they can be loaded into a fit
    /*!
        \return floating parameters of the resonance
    */
    virtual const std::vector<LauParameter*>& getFloatingParameters();

  private:
    //! Copy constructor (not implemented)
    LauGaussIncohRes( const LauGaussIncohRes& rhs );

    //! Copy assignment operator (not implemented)
    LauGaussIncohRes& operator=( const LauGaussIncohRes& rhs );

    ClassDef( LauGaussIncohRes, 0 )    // Gaussian incoherent resonance model
};

#endif
