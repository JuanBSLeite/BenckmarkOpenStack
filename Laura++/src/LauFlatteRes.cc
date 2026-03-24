
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

/*! \file LauFlatteRes.cc
    \brief File containing implementation of LauFlatteRes class.
*/

#include "LauFlatteRes.hh"

#include "LauConstants.hh"
#include "LauResonanceInfo.hh"

#include "TMath.h"
#include "TSystem.h"

#include <iostream>

LauFlatteRes::LauFlatteRes( LauResonanceInfo* resInfo,
                            const Int_t resPairAmpInt,
                            const LauDaughters* daughters ) :
    LauAbsResonance( resInfo, resPairAmpInt, daughters ),
    g1_( 0 ),
    g2_( 0 ),
    mSumSq0_( 0.0 ),
    mSumSq1_( 0.0 ),
    mSumSq2_( 0.0 ),
    mSumSq3_( 0.0 ),
    useAdlerTerm_( kFALSE ),
    sA_( 0.0 ),
    absorbM0_( kFALSE )
{
    Double_t resMass( 0.0 );
    Double_t g1Val( 0.0 );
    Double_t g2Val( 0.0 );

    const TString& resName = this->getResonanceName();

    if ( resName == "f_0(980)" ) {

        // Set the mass threshold values
        mSumSq0_ = ( LauConstants::mPi0 + LauConstants::mPi0 ) *
                   ( LauConstants::mPi0 + LauConstants::mPi0 );
        mSumSq1_ = ( LauConstants::mPi + LauConstants::mPi ) *
                   ( LauConstants::mPi + LauConstants::mPi );
        mSumSq2_ = ( LauConstants::mK + LauConstants::mK ) * ( LauConstants::mK + LauConstants::mK );
        mSumSq3_ = ( LauConstants::mK0 + LauConstants::mK0 ) *
                   ( LauConstants::mK0 + LauConstants::mK0 );

        // Set values of mass and coupling constants from:
        // Phys. Lett. B 607, 243 (2005)
        resMass = 0.965;           // 0.965 +/- 0.008 +/- 0.006 GeV/c^2
        g1Val   = 0.165;           // 0.165 +/- 0.010 +/- 0.015 GeV^2
        g2Val   = 4.21 * g1Val;    // 4.21  +/- 0.25  +/- 0.21

        useAdlerTerm_ = kFALSE;
        sA_           = 0.0;

        absorbM0_ = kTRUE;

    } else if ( resName == "K*0_0(1430)" ) {

        // Set the mass threshold values
        mSumSq0_ = ( LauConstants::mK0 + LauConstants::mPi0 ) *
                   ( LauConstants::mK0 + LauConstants::mPi0 );
        mSumSq1_ = ( LauConstants::mK + LauConstants::mPi ) *
                   ( LauConstants::mK + LauConstants::mPi );
        mSumSq2_ = ( LauConstants::mK0 + LauConstants::mEtaPrime ) *
                   ( LauConstants::mK0 + LauConstants::mEtaPrime );
        mSumSq3_ = ( LauConstants::mK0 + LauConstants::mEtaPrime ) *
                   ( LauConstants::mK0 + LauConstants::mEtaPrime );

        // Set values of mass and coupling constants from:
        // Phys. Lett. B 572, 1 (2003)
        resMass = 1.513;    // GeV/c^2
        g1Val   = 0.304;    // GeV
        g2Val   = 0.380;    // GeV

        useAdlerTerm_ = kTRUE;
        sA_           = 0.234;

        absorbM0_ = kFALSE;

    } else if ( resName == "K*+_0(1430)" || resName == "K*-_0(1430)" ) {

        // Set the mass threshold values
        mSumSq0_ = ( LauConstants::mK + LauConstants::mPi0 ) *
                   ( LauConstants::mK + LauConstants::mPi0 );
        mSumSq1_ = ( LauConstants::mK0 + LauConstants::mPi ) *
                   ( LauConstants::mK0 + LauConstants::mPi );
        mSumSq2_ = ( LauConstants::mK + LauConstants::mEtaPrime ) *
                   ( LauConstants::mK + LauConstants::mEtaPrime );
        mSumSq3_ = ( LauConstants::mK + LauConstants::mEtaPrime ) *
                   ( LauConstants::mK + LauConstants::mEtaPrime );

        // Set values of mass and coupling constants from:
        // Phys. Lett. B 572, 1 (2003)
        resMass = 1.513;    // GeV/c^2
        g1Val   = 0.304;    // GeV
        g2Val   = 0.380;    // GeV

        useAdlerTerm_ = kTRUE;
        sA_           = 0.234;

        absorbM0_ = kFALSE;

    } else if ( resName == "a0_0(980)" ) {

        // Set the mass threshold values
        mSumSq0_ = ( LauConstants::mEta + LauConstants::mPi0 ) *
                   ( LauConstants::mEta + LauConstants::mPi0 );
        mSumSq1_ = ( LauConstants::mEta + LauConstants::mPi0 ) *
                   ( LauConstants::mEta + LauConstants::mPi0 );
        mSumSq2_ = ( LauConstants::mK + LauConstants::mK ) * ( LauConstants::mK + LauConstants::mK );
        mSumSq3_ = ( LauConstants::mK0 + LauConstants::mK0 ) *
                   ( LauConstants::mK0 + LauConstants::mK0 );

        // Set values of mass and coupling constants from:
        // Phys. Rev. D 57, 3860 (1998)
        resMass = 0.982;    // 0.982 +/- 0.003 GeV/c^2
        g1Val   = 0.324 *
                0.324;    // 0.324 +/- 0.015 GeV     (NB this value needs to be squared since the paper uses g_1^2 and g_2^2)
        g2Val = 1.03 *
                g1Val;    // 1.03  +/- 0.14          (NB this is indeed the ratio of what the paper refers to as g_2^2 and g_1^2, so can be used unchanged here)

        useAdlerTerm_ = kFALSE;
        sA_           = 0.0;

        absorbM0_ = kTRUE;

    } else if ( resName == "a+_0(980)" || resName == "a-_0(980)" ) {

        // Set the mass threshold values
        mSumSq0_ = ( LauConstants::mEta + LauConstants::mPi ) *
                   ( LauConstants::mEta + LauConstants::mPi );
        mSumSq1_ = ( LauConstants::mEta + LauConstants::mPi ) *
                   ( LauConstants::mEta + LauConstants::mPi );
        mSumSq2_ = ( LauConstants::mK + LauConstants::mK0 ) *
                   ( LauConstants::mK + LauConstants::mK0 );
        mSumSq3_ = ( LauConstants::mK + LauConstants::mK0 ) *
                   ( LauConstants::mK + LauConstants::mK0 );

        // Set values of mass and coupling constants from:
        // Phys. Rev. D 57, 3860 (1998)
        resMass = 0.982;    // 0.982 +/- 0.003 GeV/c^2
        g1Val   = 0.324 *
                0.324;    // 0.324 +/- 0.015 GeV     (NB this value needs to be squared since the paper uses g_1^2 and g_2^2)
        g2Val = 1.03 *
                g1Val;    // 1.03  +/- 0.14          (NB this is indeed the ratio of what the paper refers to as g_2^2 and g_1^2, so can be used unchanged here)

        useAdlerTerm_ = kFALSE;
        sA_           = 0.0;

        absorbM0_ = kTRUE;
    }

    const TString couplingUnits = ( absorbM0_ ) ? "GeV^2" : "GeV";
    std::cout << "INFO in LauFlatteRes::LauFlatteRes : Setting default parameters for " << resName
              << ":\n";
    std::cout << "                                   : mass = " << resMass << " GeV/c^2\n";
    std::cout << "                                   : g1   = " << g1Val << " " << couplingUnits
              << "\n";
    std::cout << "                                   : g2   = " << g2Val << " " << couplingUnits
              << "\n";
    if ( absorbM0_ ) {
        std::cout << "                                   : Will absorb m0 into couplings\n";
    } else {
        std::cout << "                                   : Will not absorb m0 into couplings\n";
    }
    if ( useAdlerTerm_ ) {
        std::cout << "                                   : Will use Adler zero term\n";
        std::cout << "                                   : sA   = " << sA_ << " GeV^2/c^4\n";
    } else {
        std::cout << "                                   : Will not use Adler zero term\n";
    }

    // Set the mass value
    LauParameter* massPar = this->getMassPar();
    if ( massPar ) {
        massPar->valueAndRange( resMass, 0.0, 3.0 * resMass );
        massPar->initValue( resMass );
        massPar->genValue( resMass );
    } else {
        std::cerr << "ERROR in LauFlatteRes::LauFlatteRes : Unable to retrieve mass parameter"
                  << std::endl;
    }

    // Create the parameters for the coupling constants
    const TString& parNameBase = this->getSanitisedName();

    TString g1Name( parNameBase );
    g1Name += "_g1";
    g1_     = resInfo->getExtraParameter( g1Name );
    if ( g1_ == 0 ) {
        g1_ = new LauParameter( g1Name, g1Val, 0.0, 10.0, kTRUE );
        g1_->secondStage( kTRUE );
        resInfo->addExtraParameter( g1_ );
    }

    TString g2Name( parNameBase );
    g2Name += "_g2";
    g2_     = resInfo->getExtraParameter( g2Name );
    if ( g2_ == 0 ) {
        g2_ = new LauParameter( g2Name, g2Val, 0.0, 10.0, kTRUE );
        g2_->secondStage( kTRUE );
        resInfo->addExtraParameter( g2_ );
    }
}

LauFlatteRes::~LauFlatteRes()
{
}

void LauFlatteRes::initialise()
{
    const TString& resName = this->getResonanceName();
    if ( resName != "f_0(980)" && resName != "K*0_0(1430)" && resName != "K*+_0(1430)" &&
         resName != "K*-_0(1430)" && resName != "a0_0(980)" && resName != "a+_0(980)" &&
         resName != "a-_0(980)" ) {
        std::cerr << "ERROR in LauFlatteRes::initialise : Unexpected resonance name \"" << resName
                  << "\" for Flatte shape." << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    Int_t resSpin = this->getSpin();
    if ( resSpin != 0 ) {
        std::cerr << "WARNING in LauFlatteRes::amplitude : Resonance spin is " << resSpin << "."
                  << std::endl;
        std::cerr << "                                   : Flatte amplitude is only defined for scalers, resetting spin to 0."
                  << std::endl;
        this->changeResonance( -1.0, -1.0, 0 );
    }
}

LauComplex LauFlatteRes::resAmp( Double_t mass, Double_t spinTerm )
{
    // This function returns the complex dynamical amplitude for a Flatte distribution
    // given the invariant mass and cos(helicity) values.

    const Double_t resMass   = this->getMass();
    const Double_t resMassSq = resMass * resMass;
    const Double_t s         = mass * mass;    // Invariant mass squared combination for the system

    const Double_t g1Val = this->getg1Parameter();
    const Double_t g2Val = this->getg2Parameter();

    Double_t dMSq = resMassSq - s;
    Double_t rho1( 0.0 ), rho2( 0.0 );
    if ( s > mSumSq0_ ) {
        rho1 = TMath::Sqrt( 1.0 - mSumSq0_ / s ) / 3.0;
        if ( s > mSumSq1_ ) {
            rho1 += 2.0 * TMath::Sqrt( 1.0 - mSumSq1_ / s ) / 3.0;
            if ( s > mSumSq2_ ) {
                rho2 = 0.5 * TMath::Sqrt( 1.0 - mSumSq2_ / s );
                if ( s > mSumSq3_ ) {
                    rho2 += 0.5 * TMath::Sqrt( 1.0 - mSumSq3_ / s );
                } else {
                    // Continue analytically below higher channel thresholds
                    // This contributes to the real part of the amplitude denominator
                    dMSq += g2Val * resMass * 0.5 * TMath::Sqrt( mSumSq3_ / s - 1.0 );
                }
            } else {
                // Continue analytically below higher channel thresholds
                // This contributes to the real part of the amplitude denominator
                rho2  = 0.0;
                dMSq += g2Val * resMass *
                        ( 0.5 * TMath::Sqrt( mSumSq2_ / s - 1.0 ) +
                          0.5 * TMath::Sqrt( mSumSq3_ / s - 1.0 ) );
            }
        } else {
            // Continue analytically below higher channel thresholds
            // This contributes to the real part of the amplitude denominator
            dMSq += g1Val * resMass * 2.0 * TMath::Sqrt( mSumSq1_ / s - 1.0 ) / 3.0;
        }
    }

    Double_t massFactor = 1.0;
    if ( ! absorbM0_ ) {
        massFactor = resMass;
    }
    if ( useAdlerTerm_ ) {
        massFactor *= ( s - sA_ ) / ( resMassSq - sA_ );
    }
    const Double_t width1    = g1Val * rho1 * massFactor;
    const Double_t width2    = g2Val * rho2 * massFactor;
    const Double_t widthTerm = width1 + width2;

    LauComplex resAmplitude( dMSq, widthTerm );

    const Double_t denomFactor = dMSq * dMSq + widthTerm * widthTerm;

    Double_t invDenomFactor = 0.0;
    if ( denomFactor > 1e-10 ) {
        invDenomFactor = 1.0 / denomFactor;
    }

    resAmplitude.rescale( spinTerm * invDenomFactor );

    return resAmplitude;
}

const std::vector<LauParameter*>& LauFlatteRes::getFloatingParameters()
{
    this->clearFloatingParameters();

    if ( ! this->fixMass() ) {
        this->addFloatingParameter( this->getMassPar() );
    }

    // NB the width is given in terms of g1 and g2 so the normal width
    // parameter should be ignored, hence it is not added to the list

    if ( ! this->fixg1Parameter() ) {
        this->addFloatingParameter( g1_ );
    }

    if ( ! this->fixg2Parameter() ) {
        this->addFloatingParameter( g2_ );
    }

    return this->getParameters();
}

void LauFlatteRes::setResonanceParameter( const TString& name, const Double_t value )
{
    if ( name == "g1" ) {
        this->setg1Parameter( value );
        std::cout << "INFO in LauFlatteRes::setResonanceParameter : Setting g1 parameter to "
                  << this->getg1Parameter() << std::endl;
    } else if ( name == "g2" ) {
        this->setg2Parameter( value );
        std::cout << "INFO in LauFlatteRes::setResonanceParameter : Setting g2 parameter to "
                  << this->getg2Parameter() << std::endl;
    } else {
        std::cerr << "WARNING in LauFlatteRes::setResonanceParameter : Parameter name \"" << name
                  << "\" not recognised." << std::endl;
    }
}

void LauFlatteRes::floatResonanceParameter( const TString& name )
{
    if ( name == "g1" ) {
        if ( g1_->fixed() ) {
            g1_->fixed( kFALSE );
            this->addFloatingParameter( g1_ );
        } else {
            std::cerr << "WARNING in LauFlatteRes::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "g2" ) {
        if ( g2_->fixed() ) {
            g2_->fixed( kFALSE );
            this->addFloatingParameter( g2_ );
        } else {
            std::cerr << "WARNING in LauFlatteRes::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else {
        std::cerr << "WARNING in LauFlatteRes::fixResonanceParameter: Parameter name not reconised.  No parameter changes made."
                  << std::endl;
    }
}

LauParameter* LauFlatteRes::getResonanceParameter( const TString& name )
{
    if ( name == "g1" ) {
        return g1_;
    } else if ( name == "g2" ) {
        return g2_;
    } else {
        std::cerr << "WARNING in LauFlatteRes::getResonanceParameter: Parameter name not reconised."
                  << std::endl;
        return 0;
    }
}

void LauFlatteRes::setg1Parameter( const Double_t g1 )
{
    g1_->value( g1 );
    g1_->genValue( g1 );
    g1_->initValue( g1 );
}

void LauFlatteRes::setg2Parameter( const Double_t g2 )
{
    g2_->value( g2 );
    g2_->genValue( g2 );
    g2_->initValue( g2 );
}
