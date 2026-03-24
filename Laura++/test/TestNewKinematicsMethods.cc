
/*
Copyright 2017 University of Warwick

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

#include "LauDaughters.hh"
#include "LauKinematics.hh"
#include "LauRandom.hh"

#include "TFile.h"
#include "TMath.h"
#include "TRandom.h"
#include "TTree.h"

#include <cstdlib>
#include <iostream>
#include <vector>

struct DPpoint {
    DPpoint();
    DPpoint( const LauKinematics* kine );
    void reset();
    void set( const LauKinematics* kine );
    bool operator==( const DPpoint& rhs ) const;
    Double_t m12Sq;
    Double_t m13Sq;
    Double_t m23Sq;
    Double_t m12;
    Double_t m13;
    Double_t m23;
    Double_t c12;
    Double_t c13;
    Double_t c23;
    Double_t mPrime;
    Double_t thPrime;
    Double_t p1_12;
    Double_t p3_12;
    Double_t p2_23;
    Double_t p1_23;
    Double_t p1_13;
    Double_t p2_13;
    Double_t p1_Parent;
    Double_t p2_Parent;
    Double_t p3_Parent;
};

DPpoint::DPpoint() :
    m12Sq( 0.0 ),
    m13Sq( 0.0 ),
    m23Sq( 0.0 ),
    m12( 0.0 ),
    m13( 0.0 ),
    m23( 0.0 ),
    c12( 0.0 ),
    c13( 0.0 ),
    c23( 0.0 ),
    mPrime( 0.0 ),
    thPrime( 0.0 ),
    p1_12( 0.0 ),
    p3_12( 0.0 ),
    p2_23( 0.0 ),
    p1_23( 0.0 ),
    p1_13( 0.0 ),
    p2_13( 0.0 ),
    p1_Parent( 0.0 ),
    p2_Parent( 0.0 ),
    p3_Parent( 0.0 )
{
}

DPpoint::DPpoint( const LauKinematics* kine ) :
    m12Sq( kine->getm12Sq() ),
    m13Sq( kine->getm13Sq() ),
    m23Sq( kine->getm23Sq() ),
    m12( kine->getm12() ),
    m13( kine->getm13() ),
    m23( kine->getm23() ),
    c12( kine->getc12() ),
    c13( kine->getc13() ),
    c23( kine->getc23() ),
    mPrime( kine->getmPrime() ),
    thPrime( kine->getThetaPrime() ),
    p1_12( kine->getp1_12() ),
    p3_12( kine->getp3_12() ),
    p2_23( kine->getp2_23() ),
    p1_23( kine->getp1_23() ),
    p1_13( kine->getp1_13() ),
    p2_13( kine->getp2_13() ),
    p1_Parent( kine->getp1_Parent() ),
    p2_Parent( kine->getp2_Parent() ),
    p3_Parent( kine->getp3_Parent() )
{
}

void DPpoint::reset()
{
    m12Sq     = 0.0;
    m13Sq     = 0.0;
    m23Sq     = 0.0;
    m12       = 0.0;
    m13       = 0.0;
    m23       = 0.0;
    c12       = 0.0;
    c13       = 0.0;
    c23       = 0.0;
    mPrime    = 0.0;
    thPrime   = 0.0;
    p1_12     = 0.0;
    p3_12     = 0.0;
    p2_23     = 0.0;
    p1_23     = 0.0;
    p1_13     = 0.0;
    p2_13     = 0.0;
    p1_Parent = 0.0;
    p2_Parent = 0.0;
    p3_Parent = 0.0;
}

void DPpoint::set( const LauKinematics* kine )
{
    m12Sq     = kine->getm12Sq();
    m13Sq     = kine->getm13Sq();
    m23Sq     = kine->getm23Sq();
    m12       = kine->getm12();
    m13       = kine->getm13();
    m23       = kine->getm23();
    c12       = kine->getc12();
    c13       = kine->getc13();
    c23       = kine->getc23();
    mPrime    = kine->getmPrime();
    thPrime   = kine->getThetaPrime();
    p1_12     = kine->getp1_12();
    p3_12     = kine->getp3_12();
    p2_23     = kine->getp2_23();
    p1_23     = kine->getp1_23();
    p1_13     = kine->getp1_13();
    p2_13     = kine->getp2_13();
    p1_Parent = kine->getp1_Parent();
    p2_Parent = kine->getp2_Parent();
    p3_Parent = kine->getp3_Parent();
}

bool DPpoint::operator==( const DPpoint& rhs ) const
{
    bool ok = true;
    if ( TMath::Abs( m12Sq - rhs.m12Sq ) > 1e-6 ) {
        std::cerr << "Problem with m12Sq: " << m12Sq << " != " << rhs.m12Sq << std::endl;
        ok = false;
    }
    if ( TMath::Abs( m13Sq - rhs.m13Sq ) > 1e-6 ) {
        std::cerr << "Problem with m13Sq: " << m13Sq << " != " << rhs.m13Sq << std::endl;
        ok = false;
    }
    if ( TMath::Abs( m23Sq - rhs.m23Sq ) > 1e-6 ) {
        std::cerr << "Problem with m23Sq: " << m23Sq << " != " << rhs.m23Sq << std::endl;
        ok = false;
    }
    if ( TMath::Abs( m12 - rhs.m12 ) > 1e-6 ) {
        std::cerr << "Problem with m12: " << m12 << " != " << rhs.m12 << std::endl;
        ok = false;
    }
    if ( TMath::Abs( m13 - rhs.m13 ) > 1e-6 ) {
        std::cerr << "Problem with m13: " << m13 << " != " << rhs.m13 << std::endl;
        ok = false;
    }
    if ( TMath::Abs( m23 - rhs.m23 ) > 1e-6 ) {
        std::cerr << "Problem with m23: " << m23 << " != " << rhs.m23 << std::endl;
        ok = false;
    }
    if ( TMath::Abs( c12 - rhs.c12 ) > 1e-6 ) {
        std::cerr << "Problem with c12: " << c12 << " != " << rhs.c12 << std::endl;
        ok = false;
    }
    if ( TMath::Abs( c13 - rhs.c13 ) > 1e-6 ) {
        std::cerr << "Problem with c13: " << c13 << " != " << rhs.c13 << std::endl;
        ok = false;
    }
    if ( TMath::Abs( c23 - rhs.c23 ) > 1e-6 ) {
        std::cerr << "Problem with c23: " << c23 << " != " << rhs.c23 << std::endl;
        ok = false;
    }
    if ( TMath::Abs( mPrime - rhs.mPrime ) > 1e-6 ) {
        std::cerr << "Problem with mPrime: " << mPrime << " != " << rhs.mPrime << std::endl;
        ok = false;
    }
    if ( TMath::Abs( thPrime - rhs.thPrime ) > 1e-6 ) {
        std::cerr << "Problem with thPrime: " << thPrime << " != " << rhs.thPrime << std::endl;
        ok = false;
    }
    if ( TMath::Abs( p1_12 - rhs.p1_12 ) > 1e-6 ) {
        std::cerr << "Problem with p1_12: " << p1_12 << " != " << rhs.p1_12 << std::endl;
        ok = false;
    }
    if ( TMath::Abs( p3_12 - rhs.p3_12 ) > 1e-6 ) {
        std::cerr << "Problem with p3_12: " << p3_12 << " != " << rhs.p3_12 << std::endl;
        ok = false;
    }
    if ( TMath::Abs( p2_23 - rhs.p2_23 ) > 1e-6 ) {
        std::cerr << "Problem with p2_23: " << p2_23 << " != " << rhs.p2_23 << std::endl;
        ok = false;
    }
    if ( TMath::Abs( p1_23 - rhs.p1_23 ) > 1e-6 ) {
        std::cerr << "Problem with p1_23: " << p1_23 << " != " << rhs.p1_23 << std::endl;
        ok = false;
    }
    if ( TMath::Abs( p1_13 - rhs.p1_13 ) > 1e-6 ) {
        std::cerr << "Problem with p1_13: " << p1_13 << " != " << rhs.p1_13 << std::endl;
        ok = false;
    }
    if ( TMath::Abs( p2_13 - rhs.p2_13 ) > 1e-6 ) {
        std::cerr << "Problem with p2_13: " << p2_13 << " != " << rhs.p2_13 << std::endl;
        ok = false;
    }
    if ( TMath::Abs( p1_Parent - rhs.p1_Parent ) > 1e-6 ) {
        std::cerr << "Problem with p1_Parent: " << p1_Parent << " != " << rhs.p1_Parent << std::endl;
        ok = false;
    }
    if ( TMath::Abs( p2_Parent - rhs.p2_Parent ) > 1e-6 ) {
        std::cerr << "Problem with p2_Parent: " << p2_Parent << " != " << rhs.p2_Parent << std::endl;
        ok = false;
    }
    if ( TMath::Abs( p3_Parent - rhs.p3_Parent ) > 1e-6 ) {
        std::cerr << "Problem with p3_Parent: " << p3_Parent << " != " << rhs.p3_Parent << std::endl;
        ok = false;
    }

    return ok;
}

int main( /*int argc, char** argv*/ )
{
    // If you want to use square DP histograms for efficiency,
    // backgrounds or you just want the square DP co-ordinates
    // stored in the toy MC ntuple then set this to kTRUE
    Bool_t squareDP = kTRUE;

    // This defines the DP => decay is B+ -> pi+ pi+ pi-
    // Particle 1 = pi+
    // Particle 2 = pi+
    // Particle 3 = pi-
    // The DP is defined in terms of m13Sq and m23Sq
    LauDaughters* daughters = new LauDaughters( "B+", "pi+", "pi+", "pi-", squareDP );

    LauKinematics* kinematics = daughters->getKinematics();

    Double_t m13Sq( 0.0 ), m23Sq( 0.0 );

    DPpoint pos1;
    DPpoint pos2;

    bool ok( true );

    for ( UInt_t i( 0 ); i < 10; ++i ) {
        pos1.reset();
        pos2.reset();

        kinematics->genFlatPhaseSpace( m13Sq, m23Sq );
        kinematics->updateKinematics( m13Sq, m23Sq );

        pos1.set( kinematics );

        const Double_t m12 = kinematics->getm12();
        const Double_t c12 = kinematics->getc12();
        const Double_t m13 = kinematics->getm13();
        const Double_t c13 = kinematics->getc13();
        const Double_t m23 = kinematics->getm23();
        const Double_t c23 = kinematics->getc23();

        // this is just to make sure that we reset all variables so we can notice if things don't get set properly
        kinematics->genFlatPhaseSpace( m13Sq, m23Sq );
        kinematics->updateKinematics( m13Sq, m23Sq );

        kinematics->updateKinematicsFrom12( m12, c12 );
        pos2.set( kinematics );
        ok = pos1 == pos2;
        if ( ! ok ) {
            std::cerr << "Problem setting from m12, c12" << std::endl;
        }

        // this is just to make sure that we reset all variables so we can notice if things don't get set properly
        kinematics->genFlatPhaseSpace( m13Sq, m23Sq );
        kinematics->updateKinematics( m13Sq, m23Sq );

        kinematics->updateKinematicsFrom13( m13, c13 );
        pos2.set( kinematics );
        ok = pos1 == pos2;
        if ( ! ok ) {
            std::cerr << "Problem setting from m13, c13" << std::endl;
        }

        // this is just to make sure that we reset all variables so we can notice if things don't get set properly
        kinematics->genFlatPhaseSpace( m13Sq, m23Sq );
        kinematics->updateKinematics( m13Sq, m23Sq );

        kinematics->updateKinematicsFrom23( m23, c23 );
        pos2.set( kinematics );
        ok = pos1 == pos2;
        if ( ! ok ) {
            std::cerr << "Problem setting from m23, c23" << std::endl;
        }
    }

    if ( ! ok ) {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
