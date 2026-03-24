
/*
Copyright 2006 University of Warwick

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

/*! \file LauGounarisSakuraiRes.hh
    \brief File containing declaration of LauGounarisSakuraiRes class.
*/

/*! \class LauGounarisSakuraiRes
    \brief Class for defininf the Gounaris-Sakurai resonance model

    Class for defining the relativistic Gounaris-Sakurai resonance model, which
    includes the use of Blatt-Weisskopf barrier factors.
*/

#ifndef LAU_GOUNARIS_SAKURAI_RES
#define LAU_GOUNARIS_SAKURAI_RES

#include "LauAbsResonance.hh"
#include "LauComplex.hh"

#include "TString.h"

class LauGounarisSakuraiRes : public LauAbsResonance {

  public:
    //! Constructor
    /*!
        \param [in] resInfo the object containing information on the resonance name, mass, width, spin, charge, etc.
        \param [in] resPairAmpInt the number of the daughter not produced by the resonance
        \param [in] daughters the daughter particles
    */
    LauGounarisSakuraiRes( LauResonanceInfo* resInfo,
                           const Int_t resPairAmpInt,
                           const LauDaughters* daughters );

    //! Destructor
    virtual ~LauGounarisSakuraiRes();

    //! Initialise the model
    virtual void initialise();

    //! Get the resonance model type
    /*!
        \return the resonance model type
    */
    virtual LauAbsResonance::LauResonanceModel getResonanceModel() const
    {
        return LauAbsResonance::GS;
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

  private:
    //! Copy constructor (not implemented)
    LauGounarisSakuraiRes( const LauGounarisSakuraiRes& rhs );

    //! Copy assignment operator (not implemented)
    LauGounarisSakuraiRes& operator=( const LauGounarisSakuraiRes& rhs );

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
    //! Difference between the two daughter masses
    Double_t mDaugDiff_;
    //! Square of the difference of the two daughter masses
    Double_t mDaugDiffSq_;
    //! Square of the parent mass
    Double_t mParentSq_;
    //! Square of the bachelor mass
    Double_t mBachSq_;
    //! Extra parameter required by GS shape
    Double_t h0_;
    //! Extra parameter required by GS shape
    Double_t dhdm0_;
    //! Extra parameter required by GS shape
    Double_t d_;
    //! Value of the form factor for resonance decay (at pole mass)
    Double_t FR0_;

    ClassDef( LauGounarisSakuraiRes, 0 )    // Gounaris-Sakurai resonance model
};

#endif
