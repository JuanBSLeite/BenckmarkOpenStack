
/*
Copyright 2016 University of Warwick

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

/*! \file LauRhoOmegaMix.cc
    \brief File containing implementation of LauRhoOmegaMix class.
*/

#include "LauRhoOmegaMix.hh"

#include "LauBlattWeisskopfFactor.hh"
#include "LauConstants.hh"
#include "LauDaughters.hh"
#include "LauParameter.hh"
#include "LauResonanceInfo.hh"
#include "LauResonanceMaker.hh"

#include <iostream>

LauRhoOmegaMix::LauRhoOmegaMix( LauResonanceInfo* resInfo,
                                const LauAbsResonance::LauResonanceModel resType,
                                const Int_t resPairAmpInt,
                                const LauDaughters* daughters ) :
    LauAbsResonance( resInfo, resPairAmpInt, daughters ),
    model_( resType ),
    rhoMass_( 0.0 ),
    rhoResRadius_( 0.0 ),
    rhoParRadius_( 0.0 ),
    mOmega_( 0 ),
    mOmega0_( 0.0 ),
    mOmegaCur_( 0.0 ),
    wOmega_( 0 ),
    wOmega0_( 0.0 ),
    wOmegaCur_( 0.0 ),
    magB_( 0 ),
    phiB_( 0 ),
    delta_( 0 ),
    useDenom_( kTRUE ),
    doneFirstInit_( kFALSE ),
    rhoRes_( 0 ),
    omegaRes_( 0 ),
    whichAmpSq_( 0 )
{

    // Create the rho and omega lineshapes depending on the resonance type enumeration.
    // The narrow omega is always a relativistic Breit-Wigner (RBW), but the broader rho
    // can either be a RBW or a Gounaris-Sakurai lineshape. The" _1" form for the resonance type
    // specifies that we want to assume the amplitude denominator correction term is set to 1.
    // In principle, the rho lineshape can be any resonance defined by the resInfo pointer, so we
    // use that to extract all relevant information about the first resonance.

    // We do not need the barrier nor spin factor terms for the components separately, since
    // they will be added after the lineshapes have been combined

    LauResonanceMaker& resMaker = LauResonanceMaker::get();

    LauAbsResonance::LauResonanceModel rhoType = LauAbsResonance::RelBW;

    if ( resType == LauAbsResonance::RhoOmegaMix_GS || resType == LauAbsResonance::RhoOmegaMix_GS_1 ) {
        rhoType = LauAbsResonance::GS;
    }

    rhoRes_ = resMaker.getResonance( daughters,
                                     resInfo->getName(),
                                     resPairAmpInt,
                                     rhoType,
                                     LauBlattWeisskopfFactor::Light );

    LauResonanceInfo* omegaInfo = resMaker.getResInfo( "omega(782)" );
    omegaRes_                   = resMaker.getResonance( daughters,
                                       omegaInfo->getName(),
                                       resPairAmpInt,
                                       LauAbsResonance::RelBW,
                                       LauBlattWeisskopfFactor::Light );

    // Check to see if we want to set the denominator factor to unity
    if ( resType == LauAbsResonance::RhoOmegaMix_RBW_1 ||
         resType == LauAbsResonance::RhoOmegaMix_GS_1 ) {
        useDenom_ = kFALSE;
    }

    // Initialise various parameters that can be used in the model
    const TString& parNameBase = this->getSanitisedName();

    // Pole mass of the omega (2nd) resonance. Initialise using the resonance maker (PDG) info
    mOmega0_ = omegaInfo->getMass()->unblindValue();
    // Also set the current internal value of the omega mass for initialisation logic
    mOmegaCur_ = mOmega0_;

    TString mOmegaName( parNameBase );
    mOmegaName += "_mOmega";
    mOmega_     = resInfo->getExtraParameter( mOmegaName );
    if ( ! mOmega_ ) {
        mOmega_ = new LauParameter( mOmegaName, mOmega0_, 0.0, 100.0, kTRUE );
        mOmega_->secondStage( kTRUE );
        resInfo->addExtraParameter( mOmega_ );
    }

    // Pole width of the omega (2nd) resonance. Initialise using the resonance maker (PDG) info
    wOmega0_ = omegaInfo->getWidth()->unblindValue();
    // Also set the current internal value of the omega width for initialisation logic
    wOmegaCur_ = wOmega0_;

    TString wOmegaName( parNameBase );
    wOmegaName += "_wOmega";
    wOmega_     = resInfo->getExtraParameter( wOmegaName );
    if ( ! wOmega_ ) {
        wOmega_ = new LauParameter( wOmegaName, wOmega0_, 0.0, 100.0, kTRUE );
        wOmega_->secondStage( kTRUE );
        resInfo->addExtraParameter( wOmega_ );
    }

    // Set the magnitude and phase of the omega amplitude mixing term.
    // These should be fitted
    const Double_t magBVal = 1.0;
    const Double_t phiBVal = 0.0;

    TString magBName( parNameBase );
    magBName += "_magB";
    magB_     = resInfo->getExtraParameter( magBName );
    if ( ! magB_ ) {
        magB_ = new LauParameter( magBName, magBVal, 0.0, 100.0, kTRUE );
        magB_->secondStage( kTRUE );
        resInfo->addExtraParameter(
            magB_,
            kTRUE );    // the kTRUE here allows this value to vary between B and Bbar - TODO: maybe make this configurable?
    }

    TString phiBName( parNameBase );
    phiBName += "_phiB";
    phiB_     = resInfo->getExtraParameter( phiBName );
    if ( ! phiB_ ) {
        phiB_ = new LauParameter( phiBName, phiBVal, -10.0, 10.0, kTRUE );
        phiB_->secondStage( kTRUE );
        resInfo->addExtraParameter(
            phiB_,
            kTRUE );    // the kTRUE here allows this value to vary between B and Bbar - TODO: maybe make this configurable?
    }

    // Set the delta parameter for the omega amplitude mixing term. This
    // is usually fixed but should be varied for systematic error checks.
    // In theory, this parameter can be complex, but we only use its magnitude
    // in the mixing amplitude, since its usually very small: (2.15 +- 0.35) MeV/c^2
    const Double_t deltaVal = 2.15e-3;

    TString deltaName( parNameBase );
    deltaName += "_delta";
    delta_     = resInfo->getExtraParameter( deltaName );
    if ( ! delta_ ) {
        delta_ = new LauParameter( deltaName, deltaVal, 0.0, 100.0, kTRUE );
        delta_->secondStage( kTRUE );
        resInfo->addExtraParameter( delta_ );
    }
}

LauRhoOmegaMix::~LauRhoOmegaMix()
{
}

void LauRhoOmegaMix::initialise()
{

    // Initialise the two resonances. This is done within each amplitude() function
    // call and so floating parameters are checked every time, although secondary
    // initialisation checks will be "skipped" since the parameters will be unchanged
    // for the given set of kinematic variables/parameters
    this->initialiseRho();
    this->initialiseOmega();
}

void LauRhoOmegaMix::initialiseRho()
{
    // Propagate settings
    rhoRes_->setSpinType( this->getSpinType() );
    rhoRes_->flipHelicity( this->flipHelicity() );
    rhoRes_->ignoreMomenta( this->ignoreMomenta() );
    rhoRes_->ignoreSpin( this->ignoreSpin() );
    rhoRes_->ignoreBarrierScaling( this->ignoreBarrierScaling() );

    // Do the initialisation
    rhoRes_->initialise();

    // Keep track of the current pole mass and barrier factor terms so that
    // we can reinitialise the rho resonance if they change
    rhoMass_      = rhoRes_->getMass();
    rhoResRadius_ = rhoRes_->getResRadius();
    rhoParRadius_ = rhoRes_->getParRadius();
}

void LauRhoOmegaMix::initialiseOmega()
{
    // Set the pole mass and width of the omega resonance if this has changed
    // using the parameters mOmega_ and wOmega_

    Double_t newOmegaM( -1.0 ), newOmegaW( -1.0 );
    const Int_t newOmegaSpin( -1 );

    // See if the new pole mass is different from the current value
    Double_t tmpOmegaM = mOmega_->unblindValue();
    if ( fabs( tmpOmegaM - mOmegaCur_ ) > 1e-10 ) {
        newOmegaM = tmpOmegaM;
    }

    // See if the new pole width is different from the current value
    Double_t tmpOmegaW = wOmega_->unblindValue();
    if ( fabs( tmpOmegaW - wOmegaCur_ ) > 1e-10 ) {
        newOmegaW = tmpOmegaW;
    }

    // If any parameter is negative, they are unchanged
    omegaRes_->changeResonance( newOmegaM, newOmegaW, newOmegaSpin );

    Bool_t changedOmegaM( kFALSE );
    if ( newOmegaM > -1.0 ) {
        changedOmegaM = kTRUE;
    }

    if ( doneFirstInit_ == kFALSE ) {

        // Let the omega resonance pointer know if the mass or width are fixed or floating
        omegaRes_->fixMass( this->fixmOmegaValue() );
        omegaRes_->fixWidth( this->fixwOmegaValue() );

        // We do not need to use the spin terms for the omega lineshape, since we
        // use those from the rho for the full amplitude later on
        // (as such we do not need to propagate the flip helicity setting to the omega, just to the rho)
        omegaRes_->ignoreSpin( kTRUE );

        // We want to ignore momentum-dependent width effects: just use the constant pole width
        omegaRes_->ignoreMomenta( kTRUE );

        // And we also need to ignore barrier scaling
        omegaRes_->ignoreBarrierScaling( kTRUE );

        // Initialise the omega resonance pointer
        omegaRes_->initialise();

        doneFirstInit_ = kTRUE;

    } else {

        // Reinitialise the omega resonance pointer only if we have changed
        // its pole mass. It has no barrier factor

        if ( changedOmegaM == kTRUE ) {
            omegaRes_->initialise();
        }
    }

    // Keep track of the current values of the mass and width of the omega (floating/fixed)
    mOmegaCur_ = tmpOmegaM;
    wOmegaCur_ = tmpOmegaW;
}

LauComplex LauRhoOmegaMix::amplitude( const LauKinematics* kinematics )
{

    // This function overrides and returns the complex dynamical amplitude for the
    // rho-omega mass mixing amplitude given the kinematics

    // Check to see if we need to reinitialise the rho resonance pointer
    const Double_t resMass   = rhoRes_->getMass();
    const Double_t resRadius = rhoRes_->getResRadius();
    const Double_t parRadius = rhoRes_->getParRadius();

    if ( ( ( ! this->fixMass() ) && fabs( resMass - rhoMass_ ) > 1e-10 ) ||
         ( ( ! this->fixResRadius() ) && fabs( resRadius - rhoResRadius_ ) > 1e-10 ) ||
         ( ( ! this->fixParRadius() ) && fabs( parRadius - rhoParRadius_ ) > 1e-10 ) ) {

        this->initialiseRho();
    }

    // Always check the initialisaton of the omega resonance in case we have varied
    // its mass/width via the fit parameters
    this->initialiseOmega();

    // First, get the amplitude of the first (rho) resonance.
    // This will include the full barrier and spin terms
    const LauComplex rhoAmp = rhoRes_->amplitude( kinematics );

    // Next, get the amplitude of the second (omega) resonance. This ignores barrier
    // and spin terms, and uses the pole width only (no momentum dependence)
    const LauComplex omegaAmp = omegaRes_->amplitude( kinematics );

    // The Delta parameter, which we assume is purely real. Theoretically, delta can
    // be complex, but in practice we only use its (usually small) magnitude
    const Double_t Delta = ( resMass + mOmegaCur_ ) * this->getdeltaValue();

    // The B amplitude term
    const Double_t magBVal = this->getmagBValue() * Delta;
    const Double_t phiBVal = this->getphiBValue();
    const LauComplex BTerm = LauComplex( magBVal * cos( phiBVal ), magBVal * sin( phiBVal ) );

    // The mass mixing term
    const LauComplex unity( 1.0, 0.0 );
    const LauComplex mixingTerm = BTerm * omegaAmp + unity;

    // Now form the full amplitude
    LauComplex resAmplitude = rhoAmp * mixingTerm;

    // Add the mixing correction denominator term if required
    if ( useDenom_ ) {

        // Here, we need to disable the rho barrier & spin factors, since they are
        // only needed for the numerator term of the full amplitude. Note that we still
        // need to use the momentum-dependent width (with its resonance barrier term)

        // Disable barrier scaling factors for the amplitude (not width)
        rhoRes_->ignoreBarrierScaling( kTRUE );
        // Also ignore spin terms for now
        rhoRes_->ignoreSpin( kTRUE );

        const LauComplex rhoAmp2 = rhoRes_->amplitude( kinematics );

        // Reinstate barrier scaling and spin term flags
        rhoRes_->ignoreBarrierScaling( kFALSE );
        rhoRes_->ignoreSpin( kFALSE );

        // Denominator term
        const LauComplex DeltaSq   = LauComplex( Delta * Delta, 0.0 );
        const LauComplex denomTerm = unity - DeltaSq * rhoAmp2 * omegaAmp;

        // Modify the full amplitude
        resAmplitude = resAmplitude / denomTerm;

        if ( this->whichAmpSq_ == 1 ) {
            // For omega fit fraction
            return LauComplex( sqrt( omegaAmp.abs() * Delta * magBVal ), 0.0 );
        } else if ( this->whichAmpSq_ == 2 ) {
            // For rho fit fraction
            return rhoAmp;
        }
    }

    return resAmplitude;
}

LauComplex LauRhoOmegaMix::resAmp( Double_t mass, Double_t spinTerm )
{
    std::cerr << "ERROR in LauRhoOmegaMix : This method should never be called." << std::endl;
    std::cerr << "                        : Returning zero amplitude for mass = " << mass
              << " and spinTerm = " << spinTerm << "." << std::endl;
    return LauComplex( 0.0, 0.0 );
}

const std::vector<LauParameter*>& LauRhoOmegaMix::getFloatingParameters()
{

    this->clearFloatingParameters();

    if ( ! this->fixmOmegaValue() ) {
        this->addFloatingParameter( mOmega_ );
    }

    if ( ! this->fixwOmegaValue() ) {
        this->addFloatingParameter( wOmega_ );
    }

    if ( ! this->fixmagBValue() ) {
        this->addFloatingParameter( magB_ );
    }

    if ( ! this->fixphiBValue() ) {
        this->addFloatingParameter( phiB_ );
    }

    if ( ! this->fixdeltaValue() ) {
        this->addFloatingParameter( delta_ );
    }

    if ( ! this->fixMass() ) {
        this->addFloatingParameter( this->getMassPar() );
    }

    if ( ! this->fixWidth() ) {
        this->addFloatingParameter( this->getWidthPar() );
    }

    if ( ! this->fixResRadius() ) {
        this->addFloatingParameter( this->getResBWFactor()->getRadiusParameter() );
    }

    if ( ! this->fixParRadius() ) {
        this->addFloatingParameter( this->getParBWFactor()->getRadiusParameter() );
    }

    return this->getParameters();
}

void LauRhoOmegaMix::setResonanceParameter( const TString& name, const Double_t value )
{

    // Set various parameters for the lineshape
    if ( name == "mOmega" ) {
        this->setmOmegaValue( value );
        std::cout << "INFO in LauRhoOmegaMix::setResonanceParameter : Setting parameter mOmega = "
                  << this->getmOmegaValue() << std::endl;
    } else if ( name == "wOmega" ) {
        this->setwOmegaValue( value );
        std::cout << "INFO in LauRhoOmegaMix::setResonanceParameter : Setting parameter wOmega = "
                  << this->getwOmegaValue() << std::endl;
    } else if ( name == "magB" ) {
        this->setmagBValue( value );
        std::cout << "INFO in LauRhoOmegaMix::setResonanceParameter : Setting parameter magB = "
                  << this->getmagBValue() << std::endl;
    } else if ( name == "phiB" ) {
        this->setphiBValue( value );
        std::cout << "INFO in LauRhoOmegaMix::setResonanceParameter : Setting parameter phiB = "
                  << this->getphiBValue() << std::endl;
    } else if ( name == "delta" ) {
        this->setdeltaValue( value );
        std::cout << "INFO in LauRhoOmegaMix::setResonanceParameter : Setting parameter delta = "
                  << this->getdeltaValue() << std::endl;
    } else {
        std::cerr << "WARNING in LauSigmaRes::setResonanceParameter: Parameter name " << name
                  << " not recognised.  No parameter changes made." << std::endl;
    }
}

void LauRhoOmegaMix::floatResonanceParameter( const TString& name )
{

    if ( name == "mOmega" ) {
        if ( mOmega_->fixed() ) {
            mOmega_->fixed( kFALSE );
            this->addFloatingParameter( mOmega_ );
        } else {
            std::cerr << "WARNING in LauRhoOmegaMix::floatResonanceParameter: Parameter " << name
                      << " already floating.  No parameter changes made." << std::endl;
        }
    } else if ( name == "wOmega" ) {
        if ( wOmega_->fixed() ) {
            wOmega_->fixed( kFALSE );
            this->addFloatingParameter( wOmega_ );
        } else {
            std::cerr << "WARNING in LauRhoOmegaMix::floatResonanceParameter: Parameter " << name
                      << " already floating.  No parameter changes made." << std::endl;
        }
    } else if ( name == "magB" ) {
        if ( magB_->fixed() ) {
            magB_->fixed( kFALSE );
            this->addFloatingParameter( magB_ );
        } else {
            std::cerr << "WARNING in LauRhoOmegaMix::floatResonanceParameter: Parameter " << name
                      << " already floating.  No parameter changes made." << std::endl;
        }
    } else if ( name == "phiB" ) {
        if ( phiB_->fixed() ) {
            phiB_->fixed( kFALSE );
            this->addFloatingParameter( phiB_ );
        } else {
            std::cerr << "WARNING in LauRhoOmegaMix::floatResonanceParameter: Parameter " << name
                      << " already floating.  No parameter changes made." << std::endl;
        }
    } else if ( name == "delta" ) {
        if ( delta_->fixed() ) {
            delta_->fixed( kFALSE );
            this->addFloatingParameter( delta_ );
        } else {
            std::cerr << "WARNING in LauRhoOmegaMix::floatResonanceParameter: Parameter " << name
                      << " already floating.  No parameter changes made." << std::endl;
        }
    }
}

LauParameter* LauRhoOmegaMix::getResonanceParameter( const TString& name )
{

    if ( name == "mOmega" ) {
        return mOmega_;
    } else if ( name == "wOmega" ) {
        return wOmega_;
    } else if ( name == "magB" ) {
        return magB_;
    } else if ( name == "phiB" ) {
        return phiB_;
    } else if ( name == "delta" ) {
        return delta_;
    } else {
        std::cerr << "WARNING in LauRhoOmegaMix::getResonanceParameter: Parameter name " << name
                  << " not reconised." << std::endl;
        return 0;
    }
}

void LauRhoOmegaMix::setmOmegaValue( const Double_t mOmega )
{
    mOmega_->value( mOmega );
    mOmega_->genValue( mOmega );
    mOmega_->initValue( mOmega );
}

void LauRhoOmegaMix::setwOmegaValue( const Double_t wOmega )
{
    wOmega_->value( wOmega );
    wOmega_->genValue( wOmega );
    wOmega_->initValue( wOmega );
}

void LauRhoOmegaMix::setmagBValue( const Double_t magB )
{
    magB_->value( magB );
    magB_->genValue( magB );
    magB_->initValue( magB );
}

void LauRhoOmegaMix::setphiBValue( const Double_t phiB )
{
    phiB_->value( phiB );
    phiB_->genValue( phiB );
    phiB_->initValue( phiB );
}

void LauRhoOmegaMix::setdeltaValue( const Double_t delta )
{
    delta_->value( delta );
    delta_->genValue( delta );
    delta_->initValue( delta );
}
