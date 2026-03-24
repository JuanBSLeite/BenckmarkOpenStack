
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

/*! \file LauDaughters.hh
    \brief File containing declaration of LauDaughters class.
*/

/*! \class LauDaughters
    \brief Class that defines the particular 3-body decay under study.

    The decay has the form P -> h1 h2 h3 (P stands for the parent particle, and h for the daughters).
    The constructor accepts both string name and PDG code types to describe the particles.
*/

#ifndef LAU_DAUGHTERS
#define LAU_DAUGHTERS

#include "LauKinematics.hh"

#include "TString.h"

#include <vector>

class LauParticlePDG;

class LauDaughters {

  public:
    //! Constructor from PDG codes
    /*!
        \param [in] codeParent the parent particle PDG code
        \param [in] code1 the first daughter PDG code
        \param [in] code2 the second daughter PDG code
        \param [in] code3 the third daughter PDG code
        \param [in] useSquareDP the boolean flag decision to use Square Dalitz plot (kTRUE) or standard Dalitz plot (kFALSE). Default value is kFALSE.
    */
    LauDaughters( Int_t codeParent, Int_t code1, Int_t code2, Int_t code3, Bool_t useSquareDP = kFALSE );

    //! Constructor from particle names
    /*!
        \param [in] nameParent the parent particle string name
        \param [in] name1 the first daughter string name
        \param [in] name2 the second daughter string name
        \param [in] name3 the third daughter string name
        \param [in] useSquareDP the boolean flag decision to use Square Dalitz plot (kTRUE) or standard Dalitz plot (kFALSE). Default value is kFALSE.
    */
    LauDaughters( const TString& nameParent,
                  const TString& name1,
                  const TString& name2,
                  const TString& name3,
                  Bool_t useSquareDP = kFALSE );

    //! Destructor
    virtual ~LauDaughters();

    //! Copy constructor
    LauDaughters( const LauDaughters& rhs );

    //! Is Dalitz plot symmetric, i.e. 2 identical particles
    /*!
        \return true/false whether the DP is symmetric
    */
    Bool_t gotSymmetricalDP() const { return symmetricalDP_; }

    //! Is Dalitz plot fully symmetric, i.e. 3 identical particles
    /*!
        \return true/false whether the DP is fully symmetric
    */
    Bool_t gotFullySymmetricDP() const { return fullySymmetricDP_; }

    //! Is Dalitz plot flavour-conjugate, i.e. CP(d1) = d2 and CP(d3) = d3
    /*!
        \return true/false whether the DP is flavour-conjugate
    */
    Bool_t gotFlavourConjugateDP() const { return flavourConjugateDP_; }

    //! Determine to use or not the square Dalitz plot
    /*!
        \return true/false to use the squareDP model
    */
    Bool_t squareDP() const { return kinematics_->squareDP(); }

    //! Get mass of first daughter particle
    Double_t getMassDaug1() const;
    //! Get mass of second daughter particle
    Double_t getMassDaug2() const;
    //! Get mass of third daughter particle
    Double_t getMassDaug3() const;
    //! Get mass of the parent particle
    Double_t getMassParent() const;

    //! Get name of the first daughter particle
    TString getNameDaug1() const;
    //! Get name of the second daughter particle
    TString getNameDaug2() const;
    //! Get name of the third daughter particle
    TString getNameDaug3() const;
    //! Get name of the parent particle
    TString getNameParent() const;

    //! Get sanitised name of the first daughter particle
    TString getSanitisedNameDaug1() const;
    //! Get sanitised name of the second daughter particle
    TString getSanitisedNameDaug2() const;
    //! Get sanitised name of the third daughter particle
    TString getSanitisedNameDaug3() const;
    //! Get sanitised name of the parent particle
    TString getSanitisedNameParent() const;

    //! Get PDG code of the first daughter particle
    Int_t getTypeDaug1() const;
    //! Get PDG code of the second daughter particle
    Int_t getTypeDaug2() const;
    //! Get PDG code of the third daughter particle
    Int_t getTypeDaug3() const;
    //! Get PDG code of the parent particle
    Int_t getTypeParent() const;

    //! Get charge of the first daughter particle
    Int_t getChargeDaug1() const;
    //! Get charge of the second daughter particle
    Int_t getChargeDaug2() const;
    //! Get charge of the third daughter particle
    Int_t getChargeDaug3() const;
    //! Get charge of the parent particle
    Int_t getChargeParent() const;
    //! Get charge of a particular two-daughter combination
    /*!
        \param [in] resPairAmpInt the index of the daughter not in the combination
        \return the charge of the two-daughter combination
    */
    Int_t getCharge( Int_t resPairAmpInt ) const;

    //! Retrieve the Dalitz plot kinematics
    /*!
        \return the Dalitz plot kinematics
    */
    LauKinematics* getKinematics() { return kinematics_; }

    //! Retrieve the Dalitz plot kinematics
    /*!
        \return the Dalitz plot kinematics
    */
    const LauKinematics* getKinematics() const { return kinematics_; }

  protected:
    //! Create list of all the allowed parent/daughter particles
    void createParticleLists();

    //! Set the parent particle type
    /*!
        \param [in] nameParent the name of the parent particle
    */
    void setParentType( const TString& nameParent );

    //! Set the three daughter types
    /*!
        \param [in] name1 the name of the first daughter
        \param [in] name2 the name of the second daughter
        \param [in] name3 the name of the third daughter
    */
    void setDaugType( const TString& name1, const TString& name2, const TString& name3 );

    //! Check whether there is a symmetrical Dalitz plot
    void testDPSymmetry();

    //! Check masses and charges of daughters
    void sanityCheck();

  private:
    //! Copy assignment operator (not implemented)
    LauDaughters& operator=( const LauDaughters& rhs );

    //! Dalitz plot kinematics
    LauKinematics* kinematics_;

    //! All possible daughter types
    std::vector<const LauParticlePDG*> allowedDaughters_;

    //! All possible parent types
    std::vector<const LauParticlePDG*> allowedParents_;

    //! The parent particle
    const LauParticlePDG* parent_;

    //! The daughter particles
    std::vector<const LauParticlePDG*> daughters_;

    //! Boolean flag for symmetrical Dalitz plot
    Bool_t symmetricalDP_;

    //! Boolean flag for fully symmetric Dalitz plot
    Bool_t fullySymmetricDP_;

    //! Flavour-conjugate Dalitz plot
    Bool_t flavourConjugateDP_;

    ClassDef( LauDaughters, 0 )
};

#endif
