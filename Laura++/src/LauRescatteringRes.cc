
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

/*! \file LauRescatteringRes.cc
    \brief File containing implementation of LauRescatteringRes class.
*/

#include "LauRescatteringRes.hh"

#include "LauConstants.hh"
#include "LauDaughters.hh"
#include "LauParameter.hh"
#include "LauResonanceInfo.hh"

#include "TMath.h"
#include "TSystem.h"

#include <iostream>

LauRescatteringRes::LauRescatteringRes( LauResonanceInfo* resInfo,
                                        const LauAbsResonance::LauResonanceModel resType,
                                        const Int_t resPairAmpInt,
                                        const LauDaughters* daughters ) :
    LauAbsResonance( resInfo, resPairAmpInt, daughters ),
    lambdaPiPi_( 0 ),
    lambdaKK_( 0 ),
    Mf_( 0 ),
    Ms_( 0 ),
    Mprime_( 0 ),
    Eps1_( 0 ),
    Eps2_( 0 ),
    C0_( 0 ),
    model_( resType )
{
    TString parNameBase = this->getSanitisedName();
    TString lambdaPiPiName( parNameBase );
    lambdaPiPiName += "_lambdaPiPi";
    lambdaPiPi_     = resInfo->getExtraParameter( lambdaPiPiName );
    if ( lambdaPiPi_ == 0 ) {
        lambdaPiPi_ = new LauParameter( lambdaPiPiName, 1.0, 0.0, 10.0, kTRUE );
        lambdaPiPi_->secondStage( kTRUE );
        resInfo->addExtraParameter( lambdaPiPi_ );
    }
    TString lambdaKKName( parNameBase );
    lambdaKKName += "_lambdaKK";
    lambdaKK_     = resInfo->getExtraParameter( lambdaKKName );
    if ( lambdaKK_ == 0 ) {
        lambdaKK_ = new LauParameter( lambdaKKName, 2.8, 0.0, 10.0, kTRUE );
        lambdaKK_->secondStage( kTRUE );
        resInfo->addExtraParameter( lambdaKK_ );
    }
    TString MfName( parNameBase );
    MfName += "_Mf";
    Mf_     = resInfo->getExtraParameter( MfName );
    if ( Mf_ == 0 ) {
        Mf_ = new LauParameter( MfName, 1.32, 0.0, 10.0, kTRUE );
        Mf_->secondStage( kTRUE );
        resInfo->addExtraParameter( Mf_ );
    }
    TString MsName( parNameBase );
    MsName += "_Ms";
    Ms_     = resInfo->getExtraParameter( MsName );
    if ( Ms_ == 0 ) {
        Ms_ = new LauParameter( MsName, 0.92, 0.0, 10.0, kTRUE );
        Ms_->secondStage( kTRUE );
        resInfo->addExtraParameter( Ms_ );
    }
    TString MprimeName( parNameBase );
    MprimeName += "_Mprime";
    Mprime_     = resInfo->getExtraParameter( MprimeName );
    if ( Mprime_ == 0 ) {
        Mprime_ = new LauParameter( MprimeName, 1.5, 0.0, 10.0, kTRUE );
        Mprime_->secondStage( kTRUE );
        resInfo->addExtraParameter( Mprime_ );
    }

    TString Eps1Name( parNameBase );
    Eps1Name += "_Eps1";
    Eps1_     = resInfo->getExtraParameter( Eps1Name );
    if ( Eps1_ == 0 ) {
        Eps1_ = new LauParameter( Eps1Name, 2.4, 0.0, 10.0, kTRUE );
        Eps1_->secondStage( kTRUE );
        resInfo->addExtraParameter( Eps1_ );
    }

    TString Eps2Name( parNameBase );
    Eps2Name += "_Eps2";
    Eps2_     = resInfo->getExtraParameter( Eps2Name );
    if ( Eps2_ == 0 ) {
        Eps2_ = new LauParameter( Eps2Name, -5.5, -10.0, 10.0, kTRUE );
        Eps2_->secondStage( kTRUE );
        resInfo->addExtraParameter( Eps2_ );
    }

    TString C0Name( parNameBase );
    C0Name += "_C0";
    C0_     = resInfo->getExtraParameter( C0Name );
    if ( C0_ == 0 ) {
        C0_ = new LauParameter( C0Name, 1.3, 0.0, 10.0, kTRUE );
        C0_->secondStage( kTRUE );
        resInfo->addExtraParameter( C0_ );
    }
}

LauRescatteringRes::~LauRescatteringRes()
{
}

void LauRescatteringRes::initialise()
{
    const LauDaughters* daughters = this->getDaughters();
    Int_t resPairAmpInt           = this->getPairInt();
    if ( daughters->gotSymmetricalDP() && resPairAmpInt != 3 ) {
        std::cerr << "WARNING in LauRescatteringRes::initialise : Dalitz plot is symmetric - this lineshape is not appropriate."
                  << std::endl;
        std::cerr << "WARNING I think that this warning is not appropiate because Laura Simetrize at LauIsobarModel level."
                  << std::endl;
    }

    if ( ( model_ != LauAbsResonance::Rescattering ) &&
         ( model_ != LauAbsResonance::RescatteringNoInter ) ) {
        std::cerr << "WARNING in LauRescatteringRes::initialise : Unknown model requested, defaulting to exponential."
                  << std::endl;
        model_ = LauAbsResonance::Rescattering;
    }

    if ( ( model_ != LauAbsResonance::RescatteringNoInter ) && ( this->getSpin() != 0 ) ) {
        std::cerr << "WARNING in LauRescatteringRes::initialise : Non-zero spin will be ignored for this model - perhaps you should use LauAbsResonance::BelleSymNRNoInter instead"
                  << std::endl;
    }
}

LauComplex LauRescatteringRes::amplitude( const LauKinematics* kinematics )
{
    // This function returns the complex dynamical amplitude for a Reescatering distribution  o original
    // pelaez paper parameters Eq 2.15a [Pelaez et Yndúrain: arXiv:hep-ph/0411334v2 Mar 2005]
    Double_t Mprime   = this->getMprime();
    Double_t Mf       = this->getMf();
    Double_t Ms       = this->getMs();
    Double_t eps1     = this->getEps1();
    Double_t eps2     = this->getEps2();
    Double_t c0       = this->getC0();
    Double_t lambPiPi = this->getLambdaPiPi();
    Double_t lambKK   = this->getLambdaKK();

    Double_t mk = LauConstants::mK;

    // Calculate Mandelstam variables: s = m_13^2, t = m_23^2, u = m_12^2.
    Double_t s          = 0;
    Double_t t          = 0;
    Int_t resPairAmpInt = getPairInt();

    if ( resPairAmpInt == 1 ) {
        s = kinematics->getm23Sq();
        t = kinematics->getm12Sq();
    } else if ( resPairAmpInt == 2 ) {
        s = kinematics->getm13Sq();
        t = kinematics->getm23Sq();
    } else if ( resPairAmpInt == 3 ) {
        s = kinematics->getm12Sq();
        t = kinematics->getm13Sq();
    } else {
        std::cerr << "ERROR in LauAbsResonance::amplitude : Nonsense setup of resPairAmp array."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Calculate amplitude for s variable.
    Double_t mass_s     = TMath::Sqrt( s );
    Double_t k2Square_s = ( s - 4.0 * mk * mk ) / 4.0;
    Double_t k2Abs_s    = 0;
    if ( k2Square_s > 0 )
        k2Abs_s = TMath::Sqrt( k2Square_s );
    else
        k2Abs_s = TMath::Sqrt( -1.0 * k2Square_s );
    Double_t cotdelta0_s = c0 * ( s - Ms * Ms ) * ( Mf * Mf - s ) * k2Abs_s /
                           ( Mf * Mf * mass_s * k2Square_s );    // Eq 2.15a
    Double_t delta0_s = TMath::ATan( 1.0 / cotdelta0_s );
    Double_t eta0_s   = 1.0 - ( eps1 * k2Abs_s / mass_s + eps2 * k2Square_s / s ) *
                                ( Mprime * Mprime - s ) / s;    // Eq 2.15a

    if ( ( mass_s < 2.0 * mk ) || ( mass_s > Mprime ) )
        eta0_s = 1;
    Double_t mag_s = TMath::Sqrt( 1 - eta0_s * eta0_s );

    Double_t tauRe_s = mag_s * TMath::Cos( 2.0 * delta0_s );
    Double_t tauIm_s = mag_s * TMath::Sin( 2.0 * delta0_s );

    Double_t NR1_s = 1.0 / ( 1.0 + s / ( lambPiPi * lambPiPi ) );
    Double_t NR2_s = 1.0 / ( 1.0 + s / ( lambKK * lambKK ) );

    //LauComplex resAmplitude(-tauIm_s*NR1_s*NR2_s - tauIm_t*NR1_t*NR2_t, tauRe_s*NR1_s*NR2_s + tauRe_t*NR1_t*NR2_t );
    if ( ( model_ == LauAbsResonance::RescatteringNoInter ) && ( t <= s ) ) {
        NR1_s   = 0.0;
        NR2_s   = 0.0;
        tauRe_s = 0.0;
        tauIm_s = 0.0;
    }

    LauComplex resAmplitude( -tauIm_s * NR1_s * NR2_s, tauRe_s * NR1_s * NR2_s );

    return resAmplitude;
}

LauComplex LauRescatteringRes::resAmp( Double_t mass, Double_t spinTerm )
{
    std::cerr << "ERROR in LauRescatteringRes : This method should never be called." << std::endl;
    std::cerr << "                          : Returning zero amplitude for mass = " << mass
              << " and spinTerm = " << spinTerm << "." << std::endl;
    return LauComplex( 0.0, 0.0 );
}

const std::vector<LauParameter*>& LauRescatteringRes::getFloatingParameters()
{
    this->clearFloatingParameters();

    if ( ! this->fixLambdaPiPi() ) {
        this->addFloatingParameter( lambdaPiPi_ );
    }
    if ( ! this->fixLambdaKK() ) {
        this->addFloatingParameter( lambdaKK_ );
    }
    if ( ! this->fixMf() ) {
        this->addFloatingParameter( Mf_ );
    }
    if ( ! this->fixMs() ) {
        this->addFloatingParameter( Ms_ );
    }
    if ( ! this->fixMprime() ) {
        this->addFloatingParameter( Mprime_ );
    }
    if ( ! this->fixEps1() ) {
        this->addFloatingParameter( Eps1_ );
    }
    if ( ! this->fixEps2() ) {
        this->addFloatingParameter( Eps2_ );
    }
    if ( ! this->fixC0() ) {
        this->addFloatingParameter( C0_ );
    }

    return this->getParameters();
}

void LauRescatteringRes::setResonanceParameter( const TString& name, const Double_t value )
{
    // Set various parameters for the lineshape
    if ( name == "lambdaPiPi" ) {
        this->setLambdaPiPi( value );
        std::cout << "INFO in LauRescatteringRes::setResonanceParameter : Setting parameter lambdaPiPi = "
                  << this->getLambdaPiPi() << std::endl;
    } else if ( name == "lambdaKK" ) {
        this->setLambdaKK( value );
        std::cout << "INFO in LauRescatteringRes::setResonanceParameter : Setting parameter lambdaKK = "
                  << this->getLambdaKK() << std::endl;
    } else if ( name == "Mf" ) {
        this->setMf( value );
        std::cout << "INFO in LauRescatteringRes::setResonanceParameter : Setting parameter Mf = "
                  << this->getMf() << std::endl;
    } else if ( name == "Ms" ) {
        this->setMs( value );
        std::cout << "INFO in LauRescatteringRes::setResonanceParameter : Setting parameter Ms = "
                  << this->getMs() << std::endl;
    } else if ( name == "Mprime" ) {
        this->setMprime( value );
        std::cout << "INFO in LauRescatteringRes::setResonanceParameter : Setting parameter Mprime = "
                  << this->getMprime() << std::endl;
    } else if ( name == "Eps1" ) {
        this->setEps1( value );
        std::cout << "INFO in LauRescatteringRes::setResonanceParameter : Setting parameter eps1 = "
                  << this->getEps1() << std::endl;
    } else if ( name == "Eps2" ) {
        this->setEps2( value );
        std::cout << "INFO in LauRescatteringRes::setResonanceParameter : Setting parameter eps2 = "
                  << this->getEps2() << std::endl;
    } else if ( name == "C0" ) {
        this->setC0( value );
        std::cout << "INFO in LauRescatteringRes::setResonanceParameter : Setting parameter eps2 = "
                  << this->getC0() << std::endl;
    } else {
        std::cerr << "WARNING in LauRescatteringRes::setResonanceParameter: Parameter name not reconised.  No parameter changes made."
                  << std::endl;
    }
}

void LauRescatteringRes::floatResonanceParameter( const TString& name )
{
    if ( name == "lambdaPiPi" ) {
        if ( lambdaPiPi_->fixed() ) {
            lambdaPiPi_->fixed( kFALSE );
            this->addFloatingParameter( lambdaPiPi_ );
        } else {
            std::cerr << "WARNING in LauRescatteringRes::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "lambdaKK" ) {
        if ( lambdaKK_->fixed() ) {
            lambdaKK_->fixed( kFALSE );
            this->addFloatingParameter( lambdaKK_ );
        } else {
            std::cerr << "WARNING in LauRescatteringRes::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "Mf" ) {
        if ( Mf_->fixed() ) {
            Mf_->fixed( kFALSE );
            this->addFloatingParameter( Mf_ );
        } else {
            std::cerr << "WARNING in LauRescatteringRes::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "Ms" ) {
        if ( Ms_->fixed() ) {
            Ms_->fixed( kFALSE );
            this->addFloatingParameter( Ms_ );
        } else {
            std::cerr << "WARNING in LauRescatteringRes::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "Mprime" ) {
        if ( Mprime_->fixed() ) {
            Mprime_->fixed( kFALSE );
            this->addFloatingParameter( Mprime_ );
        } else {
            std::cerr << "WARNING in LauRescatteringRes::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "Eps1" ) {
        if ( Eps1_->fixed() ) {
            Eps1_->fixed( kFALSE );
            this->addFloatingParameter( Eps1_ );
        } else {
            std::cerr << "WARNING in LauRescatteringRes::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "Eps2" ) {
        if ( Eps2_->fixed() ) {
            Eps2_->fixed( kFALSE );
            this->addFloatingParameter( Eps2_ );
        } else {
            std::cerr << "WARNING in LauRescatteringRes::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "C0" ) {
        if ( C0_->fixed() ) {
            C0_->fixed( kFALSE );
            this->addFloatingParameter( C0_ );
        } else {
            std::cerr << "WARNING in LauRescatteringRes::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    }

    else {
        std::cerr << "WARNING in LauRescatteringRes::fixResonanceParameter: Parameter name not reconised.  No parameter changes made."
                  << std::endl;
    }
}

LauParameter* LauRescatteringRes::getResonanceParameter( const TString& name )
{
    if ( name == "lambdaPiPi" ) {
        return lambdaPiPi_;
    } else if ( name == "lambdaKK" ) {
        return lambdaKK_;
    }
    if ( name == "Mf" ) {
        return Mf_;
    } else if ( name == "Ms" ) {
        return Ms_;
    }
    if ( name == "Mprime" ) {
        return Mprime_;
    }
    if ( name == "Eps1" ) {
        return Eps1_;
    }
    if ( name == "Eps2" ) {
        return Eps2_;
    }
    if ( name == "C0" ) {
        return C0_;
    } else {
        std::cerr << "WARNING in LauRescatteringRes::getResonanceParameter: Parameter name not reconised."
                  << std::endl;
        return 0;
    }
}

void LauRescatteringRes::setLambdaPiPi( const Double_t lambda )
{
    lambdaPiPi_->value( lambda );
    lambdaPiPi_->genValue( lambda );
    lambdaPiPi_->initValue( lambda );
}

void LauRescatteringRes::setLambdaKK( const Double_t lambda )
{
    lambdaKK_->value( lambda );
    lambdaKK_->genValue( lambda );
    lambdaKK_->initValue( lambda );
}

void LauRescatteringRes::setMf( const Double_t Mf )
{
    Mf_->value( Mf );
    Mf_->genValue( Mf );
    Mf_->initValue( Mf );
}

void LauRescatteringRes::setMs( const Double_t Ms )
{
    Ms_->value( Ms );
    Ms_->genValue( Ms );
    Ms_->initValue( Ms );
}

void LauRescatteringRes::setMprime( const Double_t Mprime )
{
    Mprime_->value( Mprime );
    Mprime_->genValue( Mprime );
    Mprime_->initValue( Mprime );
}

void LauRescatteringRes::setEps1( const Double_t Eps1 )
{
    Eps1_->value( Eps1 );
    Eps1_->genValue( Eps1 );
    Eps1_->initValue( Eps1 );
}

void LauRescatteringRes::setEps2( const Double_t Eps2 )
{
    Eps2_->value( Eps2 );
    Eps2_->genValue( Eps2 );
    Eps2_->initValue( Eps2 );
}

void LauRescatteringRes::setC0( const Double_t C0 )
{
    C0_->value( C0 );
    C0_->genValue( C0 );
    C0_->initValue( C0 );
}
