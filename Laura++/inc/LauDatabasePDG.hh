
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

/*! \file LauDatabasePDG.hh
    \brief File containing declaration of LauDatabasePDG class.
*/

/*! \class LauDatabasePDG
    \brief Singleton class that provides access to particle property records.

    Singleton class that provides access to particle property records from the
    ROOT database and adds missing particles to the default setup.
    Particles may be accessed by both string and PDG code.
*/

#ifndef LAU_DATABASE_PDG
#define LAU_DATABASE_PDG

#include "TString.h"

#include <map>

class TDatabasePDG;
class LauParticlePDG;

class LauDatabasePDG {

  public:
    //! Get particle object based on the PDG code
    /*!
        \param [in] code PDG code
        \return pointer to the particle record
    */
    static const LauParticlePDG* particle( Int_t code );

    //! Get particle object based on the particle name string
    /*!
        \param [in] string particle name string
        \return pointer to the particle record
    */
    static const LauParticlePDG* particle( const TString& string );

    //! Method to convert from a particle name string into a PDG code
    /*!
        \param [in] string the particle name string
        \return the corresponding PDG code
    */
    static Int_t code( const TString& string );

    //! Method to convert from a PDG code to a particle name string
    /*!
        \param [in] code the particle PDG code
        \return the corresponding particle name
    */
    static TString string( Int_t code );

  private:
    // private methods

    //! Constructor
    LauDatabasePDG();

    //! Destructor
    virtual ~LauDatabasePDG();

    //! Copy constructor (not inplemented)
    LauDatabasePDG( const LauDatabasePDG& other );

    //! Copy assignment (not implemented)
    LauDatabasePDG& operator=( const LauDatabasePDG& other );

    //! Get the database instance
    static LauDatabasePDG& get();

    //! Add particles to the ROOT database that are missing from the default setup
    void addMissingParticles();

    //! Create particle in the database for a given PDG code and return the new particle
    /*!
        \param [in] theCode PDG code of the associated particle
        \return pointer to the particle record
    */
    const LauParticlePDG* createParticle( Int_t theCode );

    //! Create particle in the database for a given string and return the new particle
    /*!
        \param [in] theString particle string to be created
        \return pointer to the particle record
    */
    const LauParticlePDG* createParticle( const TString& theString );

    //! Get the associated pdgCode for a given particle name string
    /*!
        \param [in] theString particle string to be associated to the pdgCode
        \return the corresponding PDG code
    */
    Int_t pdgCode( const TString& theString );

    //! Get the associated name string for a given pdgCode
    /*!
        \param [in] theCode particle pdgCode
        \return the corresponding particle name
    */
    TString pdgName( Int_t theCode );

    // private data

    //! The singleton instance
    static LauDatabasePDG* theInstance_;

    //! The ROOT database
    TDatabasePDG* theDatabase_;

    //! The type to hold the list of wrapped particle objects
    typedef std::map<Int_t, const LauParticlePDG*> LauPDGList;

    //! The list of particle objects
    LauPDGList particles_;

    //! Debug flag - set in constructor
    const Bool_t debug_;

    ClassDef( LauDatabasePDG, 0 )
};

#endif
