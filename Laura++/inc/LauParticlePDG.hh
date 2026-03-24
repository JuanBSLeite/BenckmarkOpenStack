
/*
Copyright 2008 University of Warwick

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

/*! \file LauParticlePDG.hh
    \brief File containing declaration of LauParticlePDG class.
*/

/*! \class LauParticlePDG
    \brief Class for storing particle properties

    Class that provides the properties of particles, including mapping particles and their anti-particles.
*/

#ifndef LAU_PARTICLE_PDG
#define LAU_PARTICLE_PDG

#include "TMath.h"
#include "TString.h"

class TParticlePDG;

class LauParticlePDG {

  public:
    //! Destructor
    virtual ~LauParticlePDG();

    //! Method to access the record for the charge conjugate of this particle
    /*!
        \return pointer to the anti-particle
    */
    const LauParticlePDG* antiParticle() const { return antiParticle_; }

    //! The particle PDG code
    /*!
        \return the PDG code of the particle
    */
    Int_t code() const;

    //! The anti-particle PDG code
    /*!
        \return the PDG code of the charge conjugate of the particle
    */
    Int_t conjugateCode() const;

    //! The PDG code for the particle/antiparticle
    /*!
        \return the modulus of the PDG code of the particle
    */
    Int_t type() const { return TMath::Abs( this->code() ); }

    //! The charge of the particle
    /*!
        \return the charge of the particle
    */
    Int_t charge() const;

    //! The mass of the particle
    /*!
        \return the mass of the particle
    */
    Double_t mass() const { return mass_; }

    //! The mass square of the particle
    /*!
        \return the mass square of the particle
    */
    Double_t massSq() const { return mass_ * mass_; }

    //! Check whether the particle is a kaon (K+, K-, K0, K0bar, KS or KL)
    /*!
        \return true if the particle is a kaon, false otherwise
    */
    Bool_t isKaon() const { return kaon_; }

    //! Check whether the particle is a pion (pi+, pi- or pi0)
    /*!
        \return true if the particle is a pion, false otherwise
    */
    Bool_t isPion() const { return pion_; }

    //! Check whether the particle is a KS
    /*!
        \return true if the particle is a KS, false otherwise
    */
    Bool_t isKs() const { return ks_; }

    //! Check whether the particle is a pi0
    /*!
        \return true if the particle is a pi0, false otherwise
    */
    Bool_t isPi0() const { return pi0_; }

    //! Check whether the particle is a D
    /*!
        \return true if the particle is a D, false otherwise
    */
    Bool_t isD() const { return D_; }

    //! Check whether the particle is a B
    /*!
        \return true if the particle is a B, false otherwise
    */
    Bool_t isB() const { return B_; }

    //! Particle name
    /*!
        \return the name of the particle
    */
    TString string() const;

    //! Particle name, containing only alphanumeric characters
    /*!
        \return the name of the particle, removing all non-alphanumeric characters
    */
    TString stringAlphaNum() const;

  protected:

  private:
    //! LauDataBasePDG is a friend class
    friend class LauDatabasePDG;

    //! Constructor
    /*!
        \param [in] particle the particle
    */
    explicit LauParticlePDG( const TParticlePDG* particle );

    //! Copy constructor (not implemented)
    LauParticlePDG( const LauParticlePDG& rhs );

    //! Copy assignment operator (not implemented)
    LauParticlePDG& operator=( const LauParticlePDG& rhs );

    //! Set the antiparticle
    /*!
        \param [in] antiPart the antiparticle
    */
    void setAntiParticle( const LauParticlePDG* antiPart ) { antiParticle_ = antiPart; }

    //! The particle
    const TParticlePDG* theParticle_;

    //!The antiparticle
    const LauParticlePDG* antiParticle_;

    //! The mass of the particle
    Double_t mass_;

    //! Flag whether the particle is a kaon
    Bool_t kaon_;
    //! Flag whether the particle is a pion
    Bool_t pion_;
    //! Flag whether the particle is a KS
    Bool_t ks_;
    //! Flag whether the particle is a pi0
    Bool_t pi0_;
    //! Flag whether the particle is a D
    Bool_t D_;
    //! Flag whether the particle is a B
    Bool_t B_;

    ClassDef( LauParticlePDG, 0 )
};

#endif
