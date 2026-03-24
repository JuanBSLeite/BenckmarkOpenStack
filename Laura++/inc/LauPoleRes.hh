
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

/*! \file LauPoleRes.hh
    \brief File containing declaration of LauPoleRes class.
*/

/*! \class LauPoleRes
    \brief Class for defining a pole-type resonance model

    Class for defining a pole-type resonance model.
    This is a pole type resonance 1/((m_s)^2-m^2); m_s=m0-iw0, good to represent the sigma (f0(500)).
    J. A. Oller, Final state interactions in hadronic D decays, Phys. Rev. D71 (2005) 054030, arXiv:hep-ph/0411105.
*/

#ifndef LAU_POLE_RES
#define LAU_POLE_RES

#include "LauAbsResonance.hh"
#include "LauComplex.hh"

#include "TString.h"

class LauPoleRes : public LauAbsResonance {

  public:
    //! Constructor
    /*!
        \param [in] resInfo the object containing information on the resonance name, mass, width, spin, charge, etc.
        \param [in] resPairAmpInt the number of the daughter not produced by the resonance
        \param [in] daughters the daughter particles
    */
    LauPoleRes( LauResonanceInfo* resInfo, const Int_t resPairAmpInt, const LauDaughters* daughters );

    //! Destructor
    virtual ~LauPoleRes();

    //! Initialise the model
    virtual void initialise();

    //! Get the resonance model type
    /*!
        \return the resonance model type
    */
    virtual LauAbsResonance::LauResonanceModel getResonanceModel() const
    {
        return LauAbsResonance::Pole;
    }

    //! Retrieve the resonance parameters, e.g. so that they can be loaded into a fit
    /*!
        \return floating parameters of the resonance
    */
    virtual const std::vector<LauParameter*>& getFloatingParameters();

  protected:
    //! Complex resonant amplitude
    /*!
        \param [in] mass appropriate invariant mass for the resonance
        \param [in] spinTerm Zemach spin term
    */
    virtual LauComplex resAmp( Double_t mass, Double_t spinTerm );

  private:
    //! Copy constructor (not implemented)
    LauPoleRes( const LauPoleRes& rhs );

    //! Copy assignment operator (not implemented)
    LauPoleRes& operator=( const LauPoleRes& rhs );

    ClassDef( LauPoleRes, 0 )    // Pole resonance model
};

#endif
