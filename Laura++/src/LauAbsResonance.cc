
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

/*! \file LauAbsResonance.cc
    \brief File containing implementation of LauAbsResonance class.
*/

#include "LauAbsResonance.hh"

#include "LauConstants.hh"
#include "LauDaughters.hh"
#include "LauKinematics.hh"
#include "LauParameter.hh"
#include "LauResonanceInfo.hh"

#include "TSystem.h"

#include <iostream>

bool LauAbsResonance::isIncoherentModel( LauResonanceModel model )
{
    switch ( model ) {
        case BW :
        case RelBW :
        case GS :
        case Flatte :
        case Sigma :
        case Kappa :
        case Dabba :
        case LASS :
        case LASS_BW :
        case LASS_NR :
        case EFKLLM :
        case KMatrix :
        case FlatNR :
        case NRModel :
        case BelleNR :
        case PowerLawNR :
        case BelleSymNR :
        case BelleSymNRNoInter :
        case TaylorNR :
        case PolNR :
        case Pole :
        case PolarFFNR :
        case PolarFFSymNR :
        case PolarFFSymNRNoInter :
        case Rescattering :
        case Rescattering2 :
        case RescatteringNoInter :
        case MIPW_MagPhase :
        case MIPW_RealImag :
        case RhoOmegaMix_GS :
        case RhoOmegaMix_RBW :
        case RhoOmegaMix_GS_1 :
        case RhoOmegaMix_RBW_1 :
            break;
        case GaussIncoh :
            return true;
    }
    return false;
}

// Constructor
LauAbsResonance::LauAbsResonance( LauResonanceInfo* resInfo,
                                  const Int_t resPairAmpInt,
                                  const LauDaughters* daughters ) :
    resInfo_( resInfo ),
    daughters_( daughters ),
    resName_( ( resInfo != 0 ) ? resInfo->getName() : "" ),
    sanitisedName_( ( resInfo != 0 ) ? resInfo->getSanitisedName() : "" ),
    resMass_( ( resInfo != 0 ) ? resInfo->getMass() : 0 ),
    resWidth_( ( resInfo != 0 ) ? resInfo->getWidth() : 0 ),
    resSpin_( ( resInfo != 0 ) ? resInfo->getSpin() : 0 ),
    resCharge_( ( resInfo != 0 ) ? resInfo->getCharge() : 0 ),
    resPairAmpInt_( resPairAmpInt )
{
    if ( resInfo == 0 ) {
        std::cerr << "ERROR in LauAbsResonance constructor : null LauResonanceInfo object provided"
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    if ( daughters_ == 0 ) {
        std::cerr << "ERROR in LauAbsResonance constructor : null LauDaughters object provided"
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    nameParent_     = this->getNameParent();
    nameDaug1_      = this->getNameDaug1();
    nameDaug2_      = this->getNameDaug2();
    nameBachelor_   = this->getNameBachelor();
    massParent_     = this->getMassParent();
    massDaug1_      = this->getMassDaug1();
    massDaug2_      = this->getMassDaug2();
    massBachelor_   = this->getMassBachelor();
    chargeParent_   = this->getChargeParent();
    chargeDaug1_    = this->getChargeDaug1();
    chargeDaug2_    = this->getChargeDaug2();
    chargeBachelor_ = this->getChargeBachelor();

    // check that the total charge adds up to that of the resonance
    Int_t totalCharge = chargeDaug1_ + chargeDaug2_;
    if ( ( totalCharge != resCharge_ ) && ( resPairAmpInt_ != 0 ) ) {
        std::cerr << "ERROR in LauAbsResonance : Total charge of daughters = " << totalCharge
                  << ". Resonance charge = " << resCharge_ << "." << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }
}

// Constructor
LauAbsResonance::LauAbsResonance( const TString& resName,
                                  const Int_t resPairAmpInt,
                                  const LauDaughters* daughters,
                                  const Int_t resSpin ) :
    daughters_( daughters ),
    resName_( resName ),
    sanitisedName_( resName ),
    resSpin_( resSpin ),
    resPairAmpInt_( resPairAmpInt )
{
    if ( daughters_ == 0 ) {
        std::cerr << "ERROR in LauAbsResonance constructor : null LauDaughters object provided"
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    nameParent_     = this->getNameParent();
    nameDaug1_      = this->getNameDaug1();
    nameDaug2_      = this->getNameDaug2();
    nameBachelor_   = this->getNameBachelor();
    massParent_     = this->getMassParent();
    massDaug1_      = this->getMassDaug1();
    massDaug2_      = this->getMassDaug2();
    massBachelor_   = this->getMassBachelor();
    chargeParent_   = this->getChargeParent();
    chargeDaug1_    = this->getChargeDaug1();
    chargeDaug2_    = this->getChargeDaug2();
    chargeBachelor_ = this->getChargeBachelor();

    // Since we haven't been provided with a LauResonanceInfo object we can just
    // set the change of the resonance to be the sum of the daughter charges
    resCharge_ = chargeDaug1_ + chargeDaug2_;
}

// Destructor
LauAbsResonance::~LauAbsResonance()
{
}

LauComplex LauAbsResonance::amplitude( const LauKinematics* kinematics )
{
    // Use LauKinematics interface for amplitude

    // For resonance made from tracks i, j, we need the momenta
    // of tracks i and k in the i-j rest frame for spin helicity calculations
    // in the Zemach tensor formalism.
    // Also need the momentum of track k in the parent rest-frame for
    // calculation of the Blatt-Weisskopf factors.
    mass_      = 0.0;
    cosHel_    = 0.0;
    q_         = 0.0;
    p_         = 0.0;
    pstar_     = 0.0;
    erm_       = 1.0;
    covFactor_ = 1.0;

    if ( resPairAmpInt_ == 1 ) {

        mass_   = kinematics->getm23();
        cosHel_ = kinematics->getc23();
        q_      = kinematics->getp2_23();
        p_      = kinematics->getp1_23();
        pstar_  = kinematics->getp1_Parent();
        erm_    = kinematics->getcov23();

    } else if ( resPairAmpInt_ == 2 ) {

        mass_   = kinematics->getm13();
        cosHel_ = kinematics->getc13();
        q_      = kinematics->getp1_13();
        p_      = kinematics->getp2_13();
        pstar_  = kinematics->getp2_Parent();
        erm_    = kinematics->getcov13();

    } else if ( resPairAmpInt_ == 3 ) {

        mass_   = kinematics->getm12();
        cosHel_ = kinematics->getc12();
        q_      = kinematics->getp1_12();
        p_      = kinematics->getp3_12();
        pstar_  = kinematics->getp3_Parent();
        erm_    = kinematics->getcov12();

    } else {
        std::cerr << "ERROR in LauAbsResonance::amplitude : Nonsense setup of resPairAmp array."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    if ( this->flipHelicity() ) {
        cosHel_ *= -1.0;
    }

    if ( this->ignoreMomenta() ) {
        q_     = 1.0;
        p_     = 1.0;
        pstar_ = 1.0;
        erm_   = 1.0;
    }

    // Calculate the spin factors
    Double_t spinTerm( 1.0 );
    Double_t pProd( 1.0 );

    if ( ! this->ignoreSpin() ) {
        switch ( this->getSpinType() ) {

            case Zemach_P :
                pProd    = q_ * p_;
                spinTerm = this->calcZemachSpinFactor( pProd );
                break;

            case Zemach_Pstar :
                pProd    = q_ * pstar_;
                spinTerm = this->calcZemachSpinFactor( pProd );
                break;

            case Covariant :
                pProd    = q_ * pstar_;
                spinTerm = this->calcCovSpinFactor( pProd );
                break;

            case Legendre :
                spinTerm = this->calcLegendrePoly();
                break;
        }
    }

    // Calculate the full amplitude
    LauComplex resAmplitude = this->resAmp( mass_, spinTerm );

    return resAmplitude;
}

void LauAbsResonance::calcCovFactor( const Double_t erm )
{
    if ( resSpin_ == 0 ) {
        covFactor_ = 1.0;
    } else if ( resSpin_ == 1 ) {
        covFactor_ = erm;
    } else if ( resSpin_ == 2 ) {
        covFactor_ = erm * erm + 0.5;
    } else if ( resSpin_ == 3 ) {
        covFactor_ = erm * ( erm * erm + 1.5 );
    } else if ( resSpin_ == 4 ) {
        covFactor_ = ( 8. * erm * erm * erm * erm + 24. * erm * erm + 3. ) / 35.;
    } else if ( resSpin_ > 4 ) {
        std::cerr << "WARNING in LauAbsResonance::calcCovFactor : covariant spin factor cannot (yet) be fully calculated for spin >= 5"
                  << std::endl;
        std::cerr << "                                          : the function of sqrt(1 + (p/mParent)^2) part will be missing"
                  << std::endl;
        covFactor_ = 1.0;
    }
}

Double_t LauAbsResonance::calcCovSpinFactor( const Double_t pProd )
{
    if ( resSpin_ == 0 ) {
        covFactor_ = 1.0;
        return 1.0;
    }

    // Covariant spin factor is (p* q)^L * f_L(erm) * P_L(cosHel)
    Double_t spinFactor( pProd );
    for ( Int_t i( 1 ); i < resSpin_; ++i ) {
        spinFactor *= pProd;
    }

    this->calcCovFactor( erm_ );

    spinFactor *= covFactor_;

    spinFactor *= this->calcLegendrePoly();

    return spinFactor;
}

Double_t LauAbsResonance::calcZemachSpinFactor( const Double_t pProd ) const
{
    // Calculate the spin factors
    //
    // These are calculated as follows
    //
    // -2^j * (q*p)^j * cj * Pj(cosHel)
    //
    // where Pj(coshHel) is the jth order Legendre polynomial and
    //
    // cj = j! / (2j-1)!!

    if ( resSpin_ == 0 ) {
        return 1.0;
    }

    Double_t spinFactor( pProd );
    for ( Int_t i( 1 ); i < resSpin_; ++i ) {
        spinFactor *= pProd;
    }

    spinFactor *= this->calcLegendrePoly();

    return spinFactor;
}

Double_t LauAbsResonance::calcLegendrePoly( const Double_t cosHel )
{
    cosHel_ = cosHel;
    return this->calcLegendrePoly();
}

Double_t LauAbsResonance::calcLegendrePoly() const
{
    Double_t legPol = 1.0;

    if ( resSpin_ == 1 ) {
        // Calculate vector resonance Legendre polynomial
        legPol = -2.0 * cosHel_;
    } else if ( resSpin_ == 2 ) {
        // Calculate tensor resonance Legendre polynomial
        legPol = 4.0 * ( 3.0 * cosHel_ * cosHel_ - 1.0 ) / 3.0;
    } else if ( resSpin_ == 3 ) {
        // Calculate spin 3 resonance Legendre polynomial
        legPol = -8.0 * ( 5.0 * cosHel_ * cosHel_ * cosHel_ - 3.0 * cosHel_ ) / 5.0;
    } else if ( resSpin_ == 4 ) {
        // Calculate spin 4 resonance Legendre polynomial
        legPol = 16.0 *
                 ( 35.0 * cosHel_ * cosHel_ * cosHel_ * cosHel_ - 30.0 * cosHel_ * cosHel_ + 3.0 ) /
                 35.0;
    } else if ( resSpin_ == 5 ) {
        // Calculate spin 5 resonance Legendre polynomial
        legPol = -32.0 *
                 ( 63.0 * cosHel_ * cosHel_ * cosHel_ * cosHel_ * cosHel_ -
                   70.0 * cosHel_ * cosHel_ * cosHel_ + 15.0 * cosHel_ ) /
                 63.0;
    } else if ( resSpin_ > 5 ) {
        std::cerr << "WARNING in LauAbsResonance::calcLegendrePoly : Legendre polynomials not (yet) implemented for spin > 5"
                  << std::endl;
    }

    return legPol;
}

void LauAbsResonance::changeResonance( const Double_t newMass,
                                       const Double_t newWidth,
                                       const Int_t newSpin )
{
    if ( newMass > 0.0 ) {
        resMass_->valueAndRange( newMass, 0.0, 3.0 * newMass );
        resMass_->initValue( newMass );
        resMass_->genValue( newMass );
        std::cout << "INFO in LauAbsResonance::changeResonance : Setting mass to "
                  << resMass_->value() << std::endl;
    }
    if ( newWidth > 0.0 ) {
        resWidth_->valueAndRange( newWidth, 0.0, 3.0 * newWidth );
        resWidth_->initValue( newWidth );
        resWidth_->genValue( newWidth );
        std::cout << "INFO in LauAbsResonance::changeResonance : Setting width to "
                  << resWidth_->value() << std::endl;
    }
    if ( newSpin > -1 ) {
        resSpin_ = newSpin;
        std::cout << "INFO in LauAbsResonance::changeResonance : Setting spin to " << resSpin_
                  << std::endl;
    }
}

void LauAbsResonance::changeBWBarrierRadii( const Double_t resRadius, const Double_t parRadius )
{
    if ( resRadius >= 0.0 && resBWFactor_ != 0 ) {
        LauParameter* resBWRadius = resBWFactor_->getRadiusParameter();
        resBWRadius->value( resRadius );
        resBWRadius->initValue( resRadius );
        resBWRadius->genValue( resRadius );
        std::cout << "INFO in LauAbsResonance::changeBWBarrierRadii : Setting resonance factor radius to "
                  << resBWRadius->value() << std::endl;
    }
    if ( parRadius >= 0.0 && parBWFactor_ != 0 ) {
        LauParameter* parBWRadius = parBWFactor_->getRadiusParameter();
        parBWRadius->value( parRadius );
        parBWRadius->initValue( parRadius );
        parBWRadius->genValue( parRadius );
        std::cout << "INFO in LauAbsResonance::changeBWBarrierRadii : Setting parent factor radius to "
                  << parBWRadius->value() << std::endl;
    }
}

void LauAbsResonance::setResonanceParameter( const TString& name, const Double_t value )
{
    //This function should always be overwritten if needed in classes inheriting from LauAbsResonance.
    std::cerr << "WARNING in LauAbsResonance::setResonanceParameter : Unable to set parameter \""
              << name << "\" to value: " << value << "." << std::endl;
}

void LauAbsResonance::floatResonanceParameter( const TString& name )
{
    //This function should always be overwritten if needed in classes inheriting from LauAbsResonance.
    std::cerr << "WARNING in LauAbsResonance::floatResonanceParameter : Unable to release parameter \""
              << name << "\"." << std::endl;
}

LauParameter* LauAbsResonance::getResonanceParameter( const TString& name )
{
    //This function should always be overwritten if needed in classes inheriting from LauAbsResonance.
    std::cerr << "WARNING in LauAbsResonance::getResonanceParameter : Unable to get parameter \""
              << name << "\"." << std::endl;
    return 0;
}

void LauAbsResonance::addFloatingParameter( LauParameter* param )
{
    if ( param == 0 ) {
        return;
    }

    if ( param->clone() ) {
        resParameters_.push_back( param->parent() );
    } else {
        resParameters_.push_back( param );
    }
}

void LauAbsResonance::fixBarrierRadii( const Bool_t fixResRad, const Bool_t fixParRad )
{
    if ( resBWFactor_ == 0 ) {
        std::cerr << "WARNING in LauAbsResonance::fixBarrierRadii : resonance barrier factor not present, cannot fix/float it"
                  << std::endl;
        return;
    }

    if ( parBWFactor_ == 0 ) {
        std::cerr << "WARNING in LauAbsResonance::fixBarrierRadii : parent barrier factor not present, cannot fix/float it"
                  << std::endl;
        return;
    }

    LauParameter* resBWRadius = resBWFactor_->getRadiusParameter();
    resBWRadius->fixed( fixResRad );

    LauParameter* parBWRadius = parBWFactor_->getRadiusParameter();
    parBWRadius->fixed( fixParRad );
}

Bool_t LauAbsResonance::fixResRadius() const
{
    if ( resBWFactor_ == 0 ) {
        std::cerr << "WARNING in LauAbsResonance::fixResRadius : resonance barrier factor not present"
                  << std::endl;
        return kTRUE;
    }

    LauParameter* bwRadius = resBWFactor_->getRadiusParameter();
    return bwRadius->fixed();
}

Bool_t LauAbsResonance::fixParRadius() const
{
    if ( parBWFactor_ == 0 ) {
        std::cerr << "WARNING in LauAbsResonance::fixParRadius : parent barrier factor not present"
                  << std::endl;
        return kTRUE;
    }

    LauParameter* bwRadius = parBWFactor_->getRadiusParameter();
    return bwRadius->fixed();
}

Double_t LauAbsResonance::getResRadius() const
{
    if ( resBWFactor_ == 0 ) {
        std::cerr << "WARNING in LauAbsResonance::getResRadius : resonance barrier factor not present"
                  << std::endl;
        return -1.0;
    }

    LauParameter* bwRadius = resBWFactor_->getRadiusParameter();
    return bwRadius->unblindValue();
}

Double_t LauAbsResonance::getParRadius() const
{
    if ( parBWFactor_ == 0 ) {
        std::cerr << "WARNING in LauAbsResonance::getParRadius : parent barrier factor not present"
                  << std::endl;
        return -1.0;
    }

    LauParameter* bwRadius = parBWFactor_->getRadiusParameter();
    return bwRadius->unblindValue();
}

Double_t LauAbsResonance::getMassParent() const
{
    // Get the parent mass
    Double_t mass( LauConstants::mB );

    if ( daughters_ ) {
        mass = daughters_->getMassParent();
    }

    return mass;
}

Double_t LauAbsResonance::getMassDaug1() const
{
    // Get the daughter mass
    Double_t mass( LauConstants::mPi );

    if ( daughters_ ) {
        if ( resPairAmpInt_ == 1 ) {
            mass = daughters_->getMassDaug2();
        } else if ( resPairAmpInt_ == 2 ) {
            mass = daughters_->getMassDaug1();
        } else if ( resPairAmpInt_ == 3 ) {
            mass = daughters_->getMassDaug1();
        }
    }

    return mass;
}

Double_t LauAbsResonance::getMassDaug2() const
{
    // Get the daughter mass
    Double_t mass( LauConstants::mPi );

    if ( daughters_ ) {
        if ( resPairAmpInt_ == 1 ) {
            mass = daughters_->getMassDaug3();
        } else if ( resPairAmpInt_ == 2 ) {
            mass = daughters_->getMassDaug3();
        } else if ( resPairAmpInt_ == 3 ) {
            mass = daughters_->getMassDaug2();
        }
    }

    return mass;
}

Double_t LauAbsResonance::getMassBachelor() const
{
    // Get the bachelor mass
    Double_t mass( LauConstants::mPi );

    if ( daughters_ ) {
        if ( resPairAmpInt_ == 1 ) {
            mass = daughters_->getMassDaug1();
        } else if ( resPairAmpInt_ == 2 ) {
            mass = daughters_->getMassDaug2();
        } else if ( resPairAmpInt_ == 3 ) {
            mass = daughters_->getMassDaug3();
        }
    }

    return mass;
}

Int_t LauAbsResonance::getChargeParent() const
{
    // Get the parent charge
    Int_t charge( 0 );

    if ( daughters_ ) {
        charge = daughters_->getChargeParent();
    }

    return charge;
}

Int_t LauAbsResonance::getChargeDaug1() const
{
    // Get the daughter charge
    Int_t charge( 0 );

    if ( daughters_ ) {
        if ( resPairAmpInt_ == 1 ) {
            charge = daughters_->getChargeDaug2();
        } else if ( resPairAmpInt_ == 2 ) {
            charge = daughters_->getChargeDaug1();
        } else if ( resPairAmpInt_ == 3 ) {
            charge = daughters_->getChargeDaug1();
        }
    }

    return charge;
}

Int_t LauAbsResonance::getChargeDaug2() const
{
    // Get the daughter charge
    Int_t charge( 0 );

    if ( daughters_ ) {
        if ( resPairAmpInt_ == 1 ) {
            charge = daughters_->getChargeDaug3();
        } else if ( resPairAmpInt_ == 2 ) {
            charge = daughters_->getChargeDaug3();
        } else if ( resPairAmpInt_ == 3 ) {
            charge = daughters_->getChargeDaug2();
        }
    }

    return charge;
}

Int_t LauAbsResonance::getChargeBachelor() const
{
    // Get the bachelor charge
    Int_t charge( 0 );

    if ( daughters_ ) {
        if ( resPairAmpInt_ == 1 ) {
            charge = daughters_->getChargeDaug1();
        } else if ( resPairAmpInt_ == 2 ) {
            charge = daughters_->getChargeDaug2();
        } else if ( resPairAmpInt_ == 3 ) {
            charge = daughters_->getChargeDaug3();
        }
    }

    return charge;
}

TString LauAbsResonance::getNameParent() const
{
    // Get the parent name
    TString name( "" );

    if ( daughters_ ) {
        name = daughters_->getNameParent();
    }

    return name;
}

TString LauAbsResonance::getNameDaug1() const
{
    // Get the daughter name
    TString name( "" );

    if ( daughters_ ) {
        if ( resPairAmpInt_ == 1 ) {
            name = daughters_->getNameDaug2();
        } else if ( resPairAmpInt_ == 2 ) {
            name = daughters_->getNameDaug1();
        } else if ( resPairAmpInt_ == 3 ) {
            name = daughters_->getNameDaug1();
        }
    }

    return name;
}

TString LauAbsResonance::getNameDaug2() const
{
    // Get the daughter name
    TString name( "" );

    if ( daughters_ ) {
        if ( resPairAmpInt_ == 1 ) {
            name = daughters_->getNameDaug3();
        } else if ( resPairAmpInt_ == 2 ) {
            name = daughters_->getNameDaug3();
        } else if ( resPairAmpInt_ == 3 ) {
            name = daughters_->getNameDaug2();
        }
    }

    return name;
}

TString LauAbsResonance::getNameBachelor() const
{
    // Get the bachelor name
    TString name( "" );

    if ( daughters_ ) {
        if ( resPairAmpInt_ == 1 ) {
            name = daughters_->getNameDaug1();
        } else if ( resPairAmpInt_ == 2 ) {
            name = daughters_->getNameDaug2();
        } else if ( resPairAmpInt_ == 3 ) {
            name = daughters_->getNameDaug3();
        }
    }

    return name;
}
