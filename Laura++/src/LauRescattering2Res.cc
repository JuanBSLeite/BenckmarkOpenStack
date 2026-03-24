
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

/*! \file LauRescattering2Res.cc
    \brief File containing implementation of LauRescattering2Res class.
*/

#include "LauRescattering2Res.hh"

#include "LauConstants.hh"
#include "LauResonanceInfo.hh"

#include "TMath.h"

#include <iostream>

LauRescattering2Res::LauRescattering2Res( LauResonanceInfo* resInfo,
                                          const Int_t resPairAmpInt,
                                          const LauDaughters* daughters ) :
    LauAbsResonance( resInfo, resPairAmpInt, daughters ),
    B1_( 0 ),
    B2_( 0 ),
    B3_( 0 ),
    C1_( 0 ),
    C2_( 0 ),
    C3_( 0 ),
    C4_( 0 ),
    C5_( 0 ),
    D0_( 0 ),
    D1_( 0 ),
    D2_( 0 ),
    D3_( 0 ),
    F1_( 0 ),
    F2_( 0 ),
    F3_( 0 ),
    F4_( 0 )
{
    // Default values for parameters, taken from
    const Double_t B1Val( 23.6 );
    const Double_t B2Val( 29.4 );
    const Double_t B3Val( 0.6 );
    const Double_t C1Val( 34.39 );
    const Double_t C2Val( 4.4 );
    const Double_t C3Val( -32.9 );
    const Double_t C4Val( -16. );
    const Double_t C5Val( 7.4 );
    const Double_t D0Val( 0.59 );
    const Double_t D1Val( -0.38 );
    const Double_t D2Val( 0.12 );
    const Double_t D3Val( -0.09 );
    const Double_t F1Val( -0.043 );
    const Double_t F2Val( -0.008 );
    const Double_t F3Val( -0.28 );
    const Double_t F4Val( 0.026 );
    const TString& parNameBase = this->getSanitisedName();

    TString B1Name( parNameBase );
    B1Name += "_B1";
    B1_     = resInfo->getExtraParameter( B1Name );
    if ( B1_ == 0 ) {
        B1_ = new LauParameter( B1Name, B1Val, -100., 100., kTRUE );
        B1_->secondStage( kTRUE );
        resInfo->addExtraParameter( B1_ );
    }
    TString B2Name( parNameBase );
    B2Name += "_B2";
    B2_     = resInfo->getExtraParameter( B2Name );
    if ( B2_ == 0 ) {
        B2_ = new LauParameter( B2Name, B2Val, -100., 100., kTRUE );
        B2_->secondStage( kTRUE );
        resInfo->addExtraParameter( B2_ );
    }
    TString B3Name( parNameBase );
    B3Name += "_B3";
    B3_     = resInfo->getExtraParameter( B3Name );
    if ( B3_ == 0 ) {
        B3_ = new LauParameter( B3Name, B3Val, -100., 100., kTRUE );
        B3_->secondStage( kTRUE );
        resInfo->addExtraParameter( B3_ );
    }

    TString C1Name( parNameBase );
    C1Name += "_C1";
    C1_     = resInfo->getExtraParameter( C1Name );
    if ( C1_ == 0 ) {
        C1_ = new LauParameter( C1Name, C1Val, -100., 100., kTRUE );
        C1_->secondStage( kTRUE );
        resInfo->addExtraParameter( C1_ );
    }
    TString C2Name( parNameBase );
    C2Name += "_C2";
    C2_     = resInfo->getExtraParameter( C2Name );
    if ( C2_ == 0 ) {
        C2_ = new LauParameter( C2Name, C2Val, -100., 100., kTRUE );
        C2_->secondStage( kTRUE );
        resInfo->addExtraParameter( C2_ );
    }
    TString C3Name( parNameBase );
    C3Name += "_C3";
    C3_     = resInfo->getExtraParameter( C3Name );
    if ( C3_ == 0 ) {
        C3_ = new LauParameter( C3Name, C3Val, -100., 100., kTRUE );
        C3_->secondStage( kTRUE );
        resInfo->addExtraParameter( C3_ );
    }
    TString C4Name( parNameBase );
    C4Name += "_C4";
    C4_     = resInfo->getExtraParameter( C4Name );
    if ( C4_ == 0 ) {
        C4_ = new LauParameter( C4Name, C4Val, -100., 100., kTRUE );
        C4_->secondStage( kTRUE );
        resInfo->addExtraParameter( C4_ );
    }
    TString C5Name( parNameBase );
    C5Name += "_C5";
    C5_     = resInfo->getExtraParameter( C5Name );
    if ( C5_ == 0 ) {
        C5_ = new LauParameter( C5Name, C5Val, -100., 100., kTRUE );
        C5_->secondStage( kTRUE );
        resInfo->addExtraParameter( C5_ );
    }

    TString D0Name( parNameBase );
    D0Name += "_D0";
    D0_     = resInfo->getExtraParameter( D0Name );
    if ( D0_ == 0 ) {
        D0_ = new LauParameter( D0Name, D0Val, -100., 100., kTRUE );
        D0_->secondStage( kTRUE );
        resInfo->addExtraParameter( D0_ );
    }
    TString D1Name( parNameBase );
    D1Name += "_D1";
    D1_     = resInfo->getExtraParameter( D1Name );
    if ( D1_ == 0 ) {
        D1_ = new LauParameter( D1Name, D1Val, -100., 100., kTRUE );
        D1_->secondStage( kTRUE );
        resInfo->addExtraParameter( D1_ );
    }
    TString D2Name( parNameBase );
    D2Name += "_D2";
    D2_     = resInfo->getExtraParameter( D2Name );
    if ( D2_ == 0 ) {
        D2_ = new LauParameter( D2Name, D2Val, -100., 100., kTRUE );
        D2_->secondStage( kTRUE );
        resInfo->addExtraParameter( D2_ );
    }
    TString D3Name( parNameBase );
    D3Name += "_D3";
    D3_     = resInfo->getExtraParameter( D3Name );
    if ( D3_ == 0 ) {
        D3_ = new LauParameter( D3Name, D3Val, -100., 100., kTRUE );
        D3_->secondStage( kTRUE );
        resInfo->addExtraParameter( D3_ );
    }
    TString F1Name( parNameBase );
    F1Name += "_F1";
    F1_     = resInfo->getExtraParameter( F1Name );
    if ( F1_ == 0 ) {
        F1_ = new LauParameter( F1Name, F1Val, -100., 100., kTRUE );
        F1_->secondStage( kTRUE );
        resInfo->addExtraParameter( F1_ );
    }
    TString F2Name( parNameBase );
    F2Name += "_F2";
    F2_     = resInfo->getExtraParameter( F2Name );
    if ( F2_ == 0 ) {
        F2_ = new LauParameter( F2Name, F2Val, -100., 100., kTRUE );
        F2_->secondStage( kTRUE );
        resInfo->addExtraParameter( F2_ );
    }
    TString F3Name( parNameBase );
    F3Name += "_F3";
    F3_     = resInfo->getExtraParameter( F3Name );
    if ( F3_ == 0 ) {
        F3_ = new LauParameter( F3Name, F3Val, -100., 100., kTRUE );
        F3_->secondStage( kTRUE );
        resInfo->addExtraParameter( F3_ );
    }
    TString F4Name( parNameBase );
    F4Name += "_F4";
    F4_     = resInfo->getExtraParameter( F4Name );
    if ( F4_ == 0 ) {
        F4_ = new LauParameter( F4Name, F4Val, -100., 100., kTRUE );
        F4_->secondStage( kTRUE );
        resInfo->addExtraParameter( F4_ );
    }
}

LauRescattering2Res::~LauRescattering2Res()
{
    delete B1_;
    delete B2_;
    delete B3_;
    delete C1_;
    delete C2_;
    delete C3_;
    delete C4_;
    delete C5_;
    delete D0_;
    delete D1_;
    delete D2_;
    delete D3_;
    delete F1_;
    delete F2_;
    delete F3_;
    delete F4_;
}

void LauRescattering2Res::initialise()
{
    sqr_tmin[1] = 2.00 * LauConstants::mK;
    sqr_tmin[2] = 1.47;
    sqr_tmax[1] = 1.47;
    sqr_tmax[2] = 2.00;

    B0_ = 226.5 * TMath::Pi() / 180.0 + B1_->value() - B2_->value() + B3_->value();
    C0_ = phi00( sqr_tmax[1] * sqr_tmax[1], 1 ) + C1_->value() - C2_->value() + C3_->value() -
          C4_->value() + C5_->value();
    F0_ = g00( sqr_tmax[1] * sqr_tmax[1], 1 ) + F1_->value() - F2_->value() + F3_->value() -
          F4_->value();

    std::cout << "##### B = " << this->B1_->value() << ", " << this->B2_->value() << ", "
              << this->B3_->value() << std::endl;
    std::cout << "##### C = " << this->C1_->value() << ", " << this->C2_->value() << ", "
              << this->C3_->value() << ", " << this->C4_->value() << ", " << this->C5_->value()
              << std::endl;
    std::cout << "##### D = " << this->D0_->value() << ", " << this->D1_->value() << ", "
              << this->D2_->value() << ", " << this->D3_->value() << std::endl;
    std::cout << "##### F = " << this->F1_->value() << ", " << this->F2_->value() << ", "
              << this->F3_->value() << ", " << this->F4_->value() << std::endl;
}

LauComplex LauRescattering2Res::resAmp( Double_t mass, Double_t )
{
    Int_t i = 0;
    if ( mass < sqr_tmax[1] )
        i = 1;
    if ( mass > sqr_tmax[1] && mass < sqr_tmax[2] )
        i = 2;
    if ( i == 0 ) {
        std::cout << " ERROR MASS  = " << mass << " out of limmits mI, mII" << std::endl;
        return LauComplex( 0, 0 );
    }

    Double_t sqr_t = mass;
    Double_t mag   = this->g00( sqr_t, i );
    Double_t phase = this->phi00( sqr_t, i );

    return LauComplex( mag * TMath::Cos( phase ), mag * TMath::Sin( phase ) );
}

Double_t LauRescattering2Res::pn( const Double_t x_, const Double_t n ) const
{
    if ( n == 0 )
        return 1.0;
    if ( n == 1 )
        return x_;
    return 2 * x_ * pn( x_, n - 1 ) - pn( x_, n - 2 );
}

Double_t LauRescattering2Res::x( const Double_t sqr_t, const Int_t i ) const
{
    return 2.0 * ( sqr_t - sqr_tmin[i] ) / ( sqr_tmax[i] - sqr_tmin[i] ) - 1.0;
}

Double_t LauRescattering2Res::phi00( const Double_t sqr_t, const Int_t i ) const
{
    Double_t x_t = x( sqr_t, i );

    if ( i == 1 )
        return B0_ * this->pn( x_t, 0 ) + B1_->value() * this->pn( x_t, 1 ) +
               B2_->value() * this->pn( x_t, 2 ) + B3_->value() * this->pn( x_t, 3 );

    if ( i == 2 )
        return C0_ * this->pn( x_t, 0 ) + C1_->value() * this->pn( x_t, 1 ) +
               C2_->value() * this->pn( x_t, 2 ) + C3_->value() * this->pn( x_t, 3 ) +
               C4_->value() * this->pn( x_t, 4 ) + C5_->value() * this->pn( x_t, 5 );

    return 0;
}

Double_t LauRescattering2Res::g00( const Double_t sqr_t, const Int_t i ) const
{
    Double_t x_t = x( sqr_t, i );

    if ( i == 1 )
        return D0_->value() * this->pn( x_t, 0 ) + D1_->value() * this->pn( x_t, 1 ) +
               D2_->value() * this->pn( x_t, 2 ) + D3_->value() * this->pn( x_t, 3 );

    if ( i == 2 )
        return F0_ * this->pn( x_t, 0 ) + F1_->value() * this->pn( x_t, 1 ) +
               F2_->value() * this->pn( x_t, 2 ) + F3_->value() * this->pn( x_t, 3 ) +
               F4_->value() * this->pn( x_t, 4 );

    return 0;
}

// TODO up to here!
const std::vector<LauParameter*>& LauRescattering2Res::getFloatingParameters()
{
    this->clearFloatingParameters();

    if ( ! this->fixB1Parameter() ) {
        this->addFloatingParameter( B1_ );
    }
    if ( ! this->fixB2Parameter() ) {
        this->addFloatingParameter( B2_ );
    }
    if ( ! this->fixB3Parameter() ) {
        this->addFloatingParameter( B3_ );
    }

    if ( ! this->fixC1Parameter() ) {
        this->addFloatingParameter( C1_ );
    }
    if ( ! this->fixC2Parameter() ) {
        this->addFloatingParameter( C2_ );
    }
    if ( ! this->fixC3Parameter() ) {
        this->addFloatingParameter( C3_ );
    }
    if ( ! this->fixC4Parameter() ) {
        this->addFloatingParameter( C4_ );
    }
    if ( ! this->fixC5Parameter() ) {
        this->addFloatingParameter( C5_ );
    }

    if ( ! this->fixD0Parameter() ) {
        this->addFloatingParameter( D0_ );
    }
    if ( ! this->fixD1Parameter() ) {
        this->addFloatingParameter( D1_ );
    }
    if ( ! this->fixD2Parameter() ) {
        this->addFloatingParameter( D2_ );
    }
    if ( ! this->fixD3Parameter() ) {
        this->addFloatingParameter( D3_ );
    }

    if ( ! this->fixF1Parameter() ) {
        this->addFloatingParameter( F1_ );
    }

    if ( ! this->fixF2Parameter() ) {
        this->addFloatingParameter( F2_ );
    }

    if ( ! this->fixF3Parameter() ) {
        this->addFloatingParameter( F3_ );
    }

    if ( ! this->fixF4Parameter() ) {
        this->addFloatingParameter( F4_ );
    }

    return this->getParameters();
}

void LauRescattering2Res::setResonanceParameter( const TString& name, const Double_t value )
{
    // Set various parameters for the NRAmplitude lineshape dynamics
    if ( name == "B1" ) {
        this->setB1Parameter( value );
        std::cout << "INFO in LauRescattering2Res::setResonanceParameter : Setting NRAmplitude B1 = "
                  << this->getB1Parameter() << std::endl;
    } else if ( name == "B2" ) {
        this->setB2Parameter( value );
        std::cout << "INFO in LauRescattering2Res::setResonanceParameter : Setting NRAmplitude B2 = "
                  << this->getB2Parameter() << std::endl;
    } else if ( name == "B3" ) {
        this->setB3Parameter( value );
        std::cout << "INFO in LauRescattering2Res::setResonanceParameter : Setting NRAmplitude B3 = "
                  << this->getB3Parameter() << std::endl;
    } else if ( name == "C1" ) {
        this->setC1Parameter( value );
        std::cout << "INFO in LauRescattering2Res::setResonanceParameter : Setting NRAmplitude C1 = "
                  << this->getC1Parameter() << std::endl;
    } else if ( name == "C2" ) {
        this->setC2Parameter( value );
        std::cout << "INFO in LauRescattering2Res::setResonanceParameter : Setting NRAmplitude C2 = "
                  << this->getC2Parameter() << std::endl;
    } else if ( name == "C3" ) {
        this->setC3Parameter( value );
        std::cout << "INFO in LauRescattering2Res::setResonanceParameter : Setting NRAmplitude C3 = "
                  << this->getC3Parameter() << std::endl;
    } else if ( name == "C4" ) {
        this->setC4Parameter( value );
        std::cout << "INFO in LauRescattering2Res::setResonanceParameter : Setting NRAmplitude C4 = "
                  << this->getC4Parameter() << std::endl;
    } else if ( name == "C5" ) {
        this->setC5Parameter( value );
        std::cout << "INFO in LauRescattering2Res::setResonanceParameter : Setting NRAmplitude C5 = "
                  << this->getC5Parameter() << std::endl;
    } else if ( name == "D0" ) {
        this->setD0Parameter( value );
        std::cout << "INFO in LauRescattering2Res::setResonanceParameter : Setting NRAmplitude D0 = "
                  << this->getD0Parameter() << std::endl;
    } else if ( name == "D1" ) {
        this->setD1Parameter( value );
        std::cout << "INFO in LauRescattering2Res::setResonanceParameter : Setting NRAmplitude D1 = "
                  << this->getD1Parameter() << std::endl;
    } else if ( name == "D2" ) {
        this->setD2Parameter( value );
        std::cout << "INFO in LauRescattering2Res::setResonanceParameter : Setting NRAmplitude D2 = "
                  << this->getD2Parameter() << std::endl;
    } else if ( name == "D3" ) {
        this->setD3Parameter( value );
        std::cout << "INFO in LauRescattering2Res::setResonanceParameter : Setting NRAmplitude D3 = "
                  << this->getD3Parameter() << std::endl;
    } else if ( name == "F1" ) {
        this->setF1Parameter( value );
        std::cout << "INFO in LauRescattering2Res::setResonanceParameter : Setting NRAmplitude F1 = "
                  << this->getF1Parameter() << std::endl;
    } else if ( name == "F2" ) {
        this->setF2Parameter( value );
        std::cout << "INFO in LauRescattering2Res::setResonanceParameter : Setting NRAmplitude F2 = "
                  << this->getF2Parameter() << std::endl;
    } else if ( name == "F3" ) {
        this->setF3Parameter( value );
        std::cout << "INFO in LauRescattering2Res::setResonanceParameter : Setting NRAmplitude F3 = "
                  << this->getF3Parameter() << std::endl;
    } else if ( name == "F4" ) {
        this->setF4Parameter( value );
        std::cout << "INFO in LauRescattering2Res::setResonanceParameter : Setting NRAmplitude F4 = "
                  << this->getF4Parameter() << std::endl;
    } else {
        std::cerr << "WARNING in LauRescattering2Res::setResonanceParameter: Parameter name not reconised.  No parameter changes made."
                  << std::endl;
    }
}

void LauRescattering2Res::floatResonanceParameter( const TString& name )
{
    if ( name == "B1" ) {
        if ( B1_->fixed() ) {
            B1_->fixed( kFALSE );
            this->addFloatingParameter( B1_ );
        } else {
            std::cerr << "WARNING in LauRescattering2Res::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "B2" ) {
        if ( B2_->fixed() ) {
            B2_->fixed( kFALSE );
            this->addFloatingParameter( B2_ );
        } else {
            std::cerr << "WARNING in LauRescattering2Res::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "B3" ) {
        if ( B3_->fixed() ) {
            B3_->fixed( kFALSE );
            this->addFloatingParameter( B3_ );
        } else {
            std::cerr << "WARNING in LauRescattering2Res::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "C1" ) {
        if ( C1_->fixed() ) {
            C1_->fixed( kFALSE );
            this->addFloatingParameter( C1_ );
        } else {
            std::cerr << "WARNING in LauRescattering2Res::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "C2" ) {
        if ( C2_->fixed() ) {
            C2_->fixed( kFALSE );
            this->addFloatingParameter( C2_ );
        } else {
            std::cerr << "WARNING in LauRescattering2Res::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "C3" ) {
        if ( C3_->fixed() ) {
            C3_->fixed( kFALSE );
            this->addFloatingParameter( C3_ );
        } else {
            std::cerr << "WARNING in LauRescattering2Res::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "C4" ) {
        if ( C4_->fixed() ) {
            C4_->fixed( kFALSE );
            this->addFloatingParameter( C4_ );
        } else {
            std::cerr << "WARNING in LauRescattering2Res::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "C5" ) {
        if ( C5_->fixed() ) {
            C5_->fixed( kFALSE );
            this->addFloatingParameter( C5_ );
        } else {
            std::cerr << "WARNING in LauRescattering2Res::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "D0" ) {
        if ( D0_->fixed() ) {
            D0_->fixed( kFALSE );
            this->addFloatingParameter( D0_ );
        } else {
            std::cerr << "WARNING in LauRescattering2Res::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "D1" ) {
        if ( D1_->fixed() ) {
            D1_->fixed( kFALSE );
            this->addFloatingParameter( D1_ );
        } else {
            std::cerr << "WARNING in LauRescattering2Res::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "D2" ) {
        if ( D2_->fixed() ) {
            D2_->fixed( kFALSE );
            this->addFloatingParameter( D2_ );
        } else {
            std::cerr << "WARNING in LauRescattering2Res::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "D3" ) {
        if ( D3_->fixed() ) {
            D3_->fixed( kFALSE );
            this->addFloatingParameter( D3_ );
        } else {
            std::cerr << "WARNING in LauRescattering2Res::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "F1" ) {
        if ( F1_->fixed() ) {
            F1_->fixed( kFALSE );
            this->addFloatingParameter( F1_ );
        } else {
            std::cerr << "WARNING in LauRescattering2Res::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "F2" ) {
        if ( F2_->fixed() ) {
            F2_->fixed( kFALSE );
            this->addFloatingParameter( F2_ );
        } else {
            std::cerr << "WARNING in LauRescattering2Res::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "F3" ) {
        if ( F3_->fixed() ) {
            F3_->fixed( kFALSE );
            this->addFloatingParameter( F3_ );
        } else {
            std::cerr << "WARNING in LauRescattering2Res::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else if ( name == "4F" ) {
        if ( F4_->fixed() ) {
            F4_->fixed( kFALSE );
            this->addFloatingParameter( F4_ );
        } else {
            std::cerr << "WARNING in LauRescattering2Res::floatResonanceParameter: Parameter already floating.  No parameter changes made."
                      << std::endl;
        }
    } else {
        std::cerr << "WARNING in LauRescattering2Res::fixResonanceParameter: Parameter name not reconised.  No parameter changes made."
                  << std::endl;
    }
}

LauParameter* LauRescattering2Res::getResonanceParameter( const TString& name )
{
    if ( name == "B1" ) {
        return B1_;
    } else if ( name == "B2" ) {
        return B2_;
    } else if ( name == "B3" ) {
        return B3_;
    } else if ( name == "C1" ) {
        return C1_;
    } else if ( name == "C2" ) {
        return C2_;
    } else if ( name == "C3" ) {
        return C3_;
    } else if ( name == "C4" ) {
        return C4_;
    } else if ( name == "C5" ) {
        return C5_;
    } else if ( name == "D0" ) {
        return D0_;
    } else if ( name == "D1" ) {
        return D1_;
    } else if ( name == "D2" ) {
        return D2_;
    } else if ( name == "D3" ) {
        return D3_;
    } else if ( name == "F1" ) {
        return F1_;
    } else if ( name == "F2" ) {
        return F2_;
    } else if ( name == "F3" ) {
        return F3_;
    } else if ( name == "F4" ) {
        return F4_;
    } else {
        std::cerr << "WARNING in LauRescattering2Res::getResonanceParameter: Parameter name not reconised."
                  << std::endl;
        return 0;
    }
}

void LauRescattering2Res::setB1Parameter( const Double_t B1 )
{
    B1_->value( B1 );
    B1_->genValue( B1 );
    B1_->initValue( B1 );
}

void LauRescattering2Res::setB2Parameter( const Double_t B2 )
{
    B2_->value( B2 );
    B2_->genValue( B2 );
    B2_->initValue( B2 );
}

void LauRescattering2Res::setB3Parameter( const Double_t B3 )
{
    B3_->value( B3 );
    B3_->genValue( B3 );
    B3_->initValue( B3 );
}

void LauRescattering2Res::setC1Parameter( const Double_t C1 )
{
    C1_->value( C1 );
    C1_->genValue( C1 );
    C1_->initValue( C1 );
}

void LauRescattering2Res::setC2Parameter( const Double_t C2 )
{
    C2_->value( C2 );
    C2_->genValue( C2 );
    C2_->initValue( C2 );
}

void LauRescattering2Res::setC3Parameter( const Double_t C3 )
{
    C3_->value( C3 );
    C3_->genValue( C3 );
    C3_->initValue( C3 );
}

void LauRescattering2Res::setC4Parameter( const Double_t C4 )
{
    C4_->value( C4 );
    C4_->genValue( C4 );
    C4_->initValue( C4 );
}

void LauRescattering2Res::setC5Parameter( const Double_t C5 )
{
    C5_->value( C5 );
    C5_->genValue( C5 );
    C5_->initValue( C5 );
}

void LauRescattering2Res::setD0Parameter( const Double_t D0 )
{
    D0_->value( D0 );
    D0_->genValue( D0 );
    D0_->initValue( D0 );
}

void LauRescattering2Res::setD1Parameter( const Double_t D1 )
{
    D1_->value( D1 );
    D1_->genValue( D1 );
    D1_->initValue( D1 );
}

void LauRescattering2Res::setD2Parameter( const Double_t D2 )
{
    D2_->value( D2 );
    D2_->genValue( D2 );
    D2_->initValue( D2 );
}

void LauRescattering2Res::setD3Parameter( const Double_t D3 )
{
    D3_->value( D3 );
    D3_->genValue( D3 );
    D3_->initValue( D3 );
}

void LauRescattering2Res::setF1Parameter( const Double_t F1 )
{
    F1_->value( F1 );
    F1_->genValue( F1 );
    F1_->initValue( F1 );
}

void LauRescattering2Res::setF2Parameter( const Double_t F2 )
{
    F2_->value( F2 );
    F2_->genValue( F2 );
    F2_->initValue( F2 );
}

void LauRescattering2Res::setF3Parameter( const Double_t F3 )
{
    F3_->value( F3 );
    F3_->genValue( F3 );
    F3_->initValue( F3 );
}

void LauRescattering2Res::setF4Parameter( const Double_t F4 )
{
    F4_->value( F4 );
    F4_->genValue( F4 );
    F4_->initValue( F4 );
}
