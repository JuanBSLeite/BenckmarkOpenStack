
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

/*! \file LauDatabasePDG.cc
    \brief File containing implementation of LauDatabasePDG class.
*/

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

#include "LauDatabasePDG.hh"
#include "LauParticlePDG.hh"

#include "TDatabasePDG.h"

LauDatabasePDG* LauDatabasePDG::theInstance_ = 0;

// Constructor

LauDatabasePDG::LauDatabasePDG() :
    theDatabase_( new TDatabasePDG() ),
    particles_(),
    debug_( kFALSE )
{
}

// Destructor

LauDatabasePDG::~LauDatabasePDG()
{
    delete theDatabase_;
    theDatabase_ = 0;
    for ( LauPDGList::iterator iter = particles_.begin(); iter != particles_.end(); ++iter ) {
        delete iter->second;
        iter->second = 0;
    }
}

// Private static method to access the instance
// Will create it and add the missing particles if it doesn't already exist

LauDatabasePDG& LauDatabasePDG::get()
{
    if ( theInstance_ == 0 ) {
        theInstance_ = new LauDatabasePDG();
        theInstance_->addMissingParticles();
    }
    return *theInstance_;
}

// All static methods grab the instance and call the private non-static
// methods - this ensures that the database exists and that the missing
// particles are added to it.

const LauParticlePDG* LauDatabasePDG::particle( Int_t theCode )
{
    return LauDatabasePDG::get().createParticle( theCode );
}

const LauParticlePDG* LauDatabasePDG::particle( const TString& theString )
{
    return LauDatabasePDG::get().createParticle( theString );
}

Int_t LauDatabasePDG::code( const TString& theString )
{
    return LauDatabasePDG::get().pdgCode( theString );
}

TString LauDatabasePDG::string( Int_t theCode )
{
    return LauDatabasePDG::get().pdgName( theCode );
}

// The following method adds particles to the ROOT database that are
// missing from it's default setup.  This should only be called once,
// immediately after the creation of the instance.

void LauDatabasePDG::addMissingParticles()
{
    TParticlePDG* part( 0 );
    TParticlePDG* antipart( 0 );

    // format is:
    // name, title, mass, stable?, width, 3*Q, category, PDG code

    part =
        theDatabase_->AddParticle( "deuteron", "deuteron", 1.8756134, kTRUE, 0.0, 3, "Nucleus", 1011 );
    antipart = theDatabase_->AddAntiParticle( "anti-deuteron", -1011 );
    part->SetAntiParticle( antipart );
    antipart->SetAntiParticle( part );

    part = theDatabase_->AddParticle( "tritium", "tritium", 2.80925, kTRUE, 0.0, 3, "Nucleus", 1021 );
    antipart = theDatabase_->AddAntiParticle( "anti-tritium", -1021 );
    part->SetAntiParticle( antipart );
    antipart->SetAntiParticle( part );

    part     = theDatabase_->AddParticle( "He3", "He3", 2.80923, kTRUE, 0.0, 6, "Nucleus", 1012 );
    antipart = theDatabase_->AddAntiParticle( "anti-He3", -1012 );
    part->SetAntiParticle( antipart );
    antipart->SetAntiParticle( part );

    part = theDatabase_->AddParticle( "alpha", "alpha", 3.727417, kTRUE, 0.0, 3, "Nucleus", 1022 );
    antipart = theDatabase_->AddAntiParticle( "anti-alpha", -1022 );
    part->SetAntiParticle( antipart );
    antipart->SetAntiParticle( part );

    part =
        theDatabase_->AddParticle( "eta_c(2S)", "eta_c(2S)", 3.628, kTRUE, 0.0195, 0, "Meson", 20441 );

    part =
        theDatabase_->AddParticle( "psi(3770)", "psi(3770)", 3.7699, kTRUE, 0.0236, 0, "Meson", 40443 );

    part = theDatabase_->AddParticle( "rho(2S)0", "rho(2S)0", 1.46, kTRUE, 0.31, 0, "Meson", 30113 );
    part = theDatabase_->AddParticle( "rho(2S)+", "rho(2S)+", 1.46, kTRUE, 0.31, 3, "Meson", 30213 );
    antipart = theDatabase_->AddAntiParticle( "rho(2S)-", -30213 );
    part->SetAntiParticle( antipart );
    antipart->SetAntiParticle( part );

    part = theDatabase_->AddParticle( "rho(3S)0", "rho(3S)0", 1.46, kTRUE, 0.31, 0, "Meson", 40113 );
    part = theDatabase_->AddParticle( "rho(3S)+", "rho(3S)+", 1.46, kTRUE, 0.31, 3, "Meson", 40213 );
    antipart = theDatabase_->AddAntiParticle( "rho(3S)-", -40213 );
    part->SetAntiParticle( antipart );
    antipart->SetAntiParticle( part );

    part     = theDatabase_->AddParticle( "Xsu", "Xsu", 1.6, kTRUE, 0.5, 3, "Meson", 30353 );
    antipart = theDatabase_->AddAntiParticle( "anti-Xsu", -30353 );
    part->SetAntiParticle( antipart );
    antipart->SetAntiParticle( part );

    part     = theDatabase_->AddParticle( "Xu0", "Xu0", 1., kTRUE, 0., 0, "Special", 41 );
    part     = theDatabase_->AddParticle( "Xu+", "Xu+", 1., kTRUE, 0., 3, "Special", 42 );
    antipart = theDatabase_->AddAntiParticle( "Xu-", -42 );
    part->SetAntiParticle( antipart );
    antipart->SetAntiParticle( part );

    part     = theDatabase_->AddParticle( "Xsd", "Xsd", 1.6, kTRUE, 0.5, 0, "Meson", 30343 );
    antipart = theDatabase_->AddAntiParticle( "anti-Xsd", -30343 );
    part->SetAntiParticle( antipart );
    antipart->SetAntiParticle( part );

    part     = theDatabase_->AddParticle( "K'*0", "K'*0", 1.414, kTRUE, 0.232, 0, "Meson", 100313 );
    antipart = theDatabase_->AddAntiParticle( "anti-K'*0", -100313 );
    part->SetAntiParticle( antipart );
    antipart->SetAntiParticle( part );

    part     = theDatabase_->AddParticle( "K'*+", "K'*+", 1.414, kTRUE, 0.232, 3, "Meson", 100323 );
    antipart = theDatabase_->AddAntiParticle( "K'*-", -100323 );
    part->SetAntiParticle( antipart );
    antipart->SetAntiParticle( part );

    part = theDatabase_->AddParticle( "K''*0", "K''*0", 1.717, kTRUE, 0.322, 0, "Meson", 30313 );
    antipart = theDatabase_->AddAntiParticle( "anti-K''*0", -30313 );
    part->SetAntiParticle( antipart );
    antipart->SetAntiParticle( part );

    part = theDatabase_->AddParticle( "K''*+", "K''*+", 1.717, kTRUE, 0.322, 3, "Meson", 30323 );
    antipart = theDatabase_->AddAntiParticle( "K''*-", -30323 );
    part->SetAntiParticle( antipart );
    antipart->SetAntiParticle( part );
}

// These next two methods are the ones that actually build-up our database.

const LauParticlePDG* LauDatabasePDG::createParticle( Int_t theCode )
{
    // if we've already created this one simply return it
    LauPDGList::iterator iter = particles_.find( theCode );
    if ( iter != particles_.end() ) {
        return iter->second;
    }

    // otherwise we'll have to make it...

    // first get the particle from the database
    TParticlePDG* theTParticle = theDatabase_->GetParticle( theCode );
    if ( theTParticle == 0 ) {
        cerr << "ERROR in LauDatabasePDG::createParticle : Invalid PDG code: " << theCode << endl;
        return 0;
    }

    // then create a new particle from it and store it in the map
    LauParticlePDG* theParticle                        = new LauParticlePDG( theTParticle );
    std::pair<LauPDGList::iterator, bool> insertResult = particles_.insert(
        std::make_pair( theCode, theParticle ) );
    if ( ! insertResult.second ) {
        cerr << "ERROR in LauDatabasePDG::createParticle : Something went wrong inserting particle "
             << theCode << " into the map." << endl;
        delete theParticle;
        return 0;
    }

    // then get the anti-particle from the database
    TParticlePDG* theTAntiParticle = theTParticle->AntiParticle();

    // check if we have a self-conjugate particle
    if ( theTAntiParticle == 0 || theTAntiParticle == theTParticle ) {
        // if so then set it to be it's own anti-particle and return
        theParticle->setAntiParticle( theParticle );
        if ( debug_ ) {
            cout << "INFO in LauDatabasePDG::createParticle : Created self-conjugate particle "
                 << theTParticle->GetName() << endl;
        }
        return theParticle;
    }

    // otherwise make the particle's antiparticle and store it in the map
    LauParticlePDG* antiParticle = new LauParticlePDG( theTAntiParticle );
    Int_t conjCode               = antiParticle->code();
    insertResult                 = particles_.insert( std::make_pair( conjCode, antiParticle ) );
    if ( ! insertResult.second ) {
        cerr << "ERROR in LauDatabasePDG::createParticle : Something went wrong inserting anti-particle "
             << conjCode << " into the map." << endl;
        particles_.erase( theCode );
        delete theParticle;
        delete antiParticle;
        return 0;
    }

    // set the particle and anti-particle to be each others conjugate
    theParticle->setAntiParticle( antiParticle );
    antiParticle->setAntiParticle( theParticle );

    // and finally return the particle we've created
    if ( debug_ ) {
        cout << "INFO in LauDatabasePDG::createParticle : Created particle "
             << theTParticle->GetName() << " and its conjugate " << theTAntiParticle->GetName()
             << endl;
    }
    return theParticle;
}

const LauParticlePDG* LauDatabasePDG::createParticle( const TString& theString )
{
    // turn the string into a PDG code
    Int_t theCode = this->pdgCode( theString );
    if ( theCode == 0 ) {
        cerr << "ERROR in LauDatabasePDG::createParticle : Invalid particle string: " << theString
             << endl;
        return 0;
    }
    // self-message to get the particle
    return this->createParticle( theCode );
}

// These last two methods have to be independent of our own database since
// they are used by the methods that create the particles and so we'd end
// up with an infinite loop.

Int_t LauDatabasePDG::pdgCode( const TString& theString )
{
    TParticlePDG* theParticle = theDatabase_->GetParticle( theString );
    if ( theParticle == 0 ) {
        cerr << "ERROR in LauDatabasePDG::pdgCode : Invalid particle string: " << theString << endl;
        return 0;
    }
    return theParticle->PdgCode();
}

TString LauDatabasePDG::pdgName( Int_t theCode )
{
    TParticlePDG* theParticle = theDatabase_->GetParticle( theCode );
    if ( theParticle == 0 ) {
        cerr << "ERROR in LauDatabasePDG::pdgName : Invalid particle code: " << theCode << endl;
        return "";
    }
    return theParticle->GetName();
}
