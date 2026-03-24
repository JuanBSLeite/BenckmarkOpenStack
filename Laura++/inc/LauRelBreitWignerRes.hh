
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

/*! \file LauRelBreitWignerRes.hh
    \brief File containing declaration of LauRelBreitWignerRes class.
*/

/*! \class LauRelBreitWignerRes
    \brief Class for defining the relativistic Breit-Wigner resonance model

    Class for defining the relativistic Breit-Wigner resonance model, which
    includes the use of Blatt-Weisskopf barrier factors.
*/

#ifndef LAU_REL_BREIT_WIGNER_RES
#define LAU_REL_BREIT_WIGNER_RES

#include "LauAbsResonance.hh"
#include "LauComplex.hh"

#include "TString.h"

class LauRelBreitWignerRes : public LauAbsResonance {

  public:
    //! Constructor
    /*!
        \param [in] resInfo the object containing information on the resonance name, mass, width, spin, charge, etc.
        \param [in] resPairAmpInt the number of the daughter not produced by the resonance
        \param [in] daughters the daughter particles
    */
    LauRelBreitWignerRes( LauResonanceInfo* resInfo,
                          const Int_t resPairAmpInt,
                          const LauDaughters* daughters );

    //! Destructor
    virtual ~LauRelBreitWignerRes();

    //! Initialise the model
    virtual void initialise();

    //! Get the resonance model type
    /*!
        \return the resonance model type
    */
    virtual LauAbsResonance::LauResonanceModel getResonanceModel() const
    {
        return LauAbsResonance::RelBW;
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
        \param [in] spinTerm spin term
    */
    virtual LauComplex resAmp( Double_t mass, Double_t spinTerm );

    //! Calculate the form factor for the resonance
    /*!
        \param [in] z particle momentum multipled by the barrier radius
        \return value of the form factor
    */
    Double_t calcFFactor( Double_t z );

  private:
    //! Copy constructor (not implemented)
    LauRelBreitWignerRes( const LauRelBreitWignerRes& rhs );

    //! Copy assignment operator (not implemented)
    LauRelBreitWignerRes& operator=( const LauRelBreitWignerRes& rhs );

    //! Momentum of the daughters in the resonance rest frame (at pole mass)
    Double_t q0_;
    //! The resonance mass
    Double_t resMass_;
    //! Square of the resonance mass
    Double_t resMassSq_;
    //! The resonance width
    Double_t resWidth_;
    //! The resonance barrier radius
    Double_t resRadius_;
    //! The parent barrier radius
    Double_t parRadius_;
    //! Sum of the two daughter masses
    Double_t mDaugSum_;
    //! Square of the sum of the two daughter masses
    Double_t mDaugSumSq_;
    //! Difference of the two daughter masses
    Double_t mDaugDiff_;
    //! Square of the difference of the two daughter masses
    Double_t mDaugDiffSq_;
    //! Square of the parent mass
    Double_t mParentSq_;
    //! Square of the bachelor mass
    Double_t mBachSq_;
    //! Value of the form factor for resonance decay (at pole mass)
    Double_t FR0_;

    ClassDef( LauRelBreitWignerRes, 0 )    // Relativistic Breit-Wigner resonance model
};

#endif
