
/*
Copyright 2014 University of Warwick

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

/*! \file LauBlattWeisskopfFactor.hh
    \brief File containing declaration of LauBlattWeisskopfFactor class.
*/

/*! \class LauBlattWeisskopfFactor
    \brief Class that implements the Blatt-Weisskopf barrier factor

    This factor is used in resonance models such as the Relativistic Breit-Wigner.
*/

#ifndef LAU_BLATTWEISSKOPFFACTOR
#define LAU_BLATTWEISSKOPFFACTOR

#include "Rtypes.h"
#include "TString.h"

class LauParameter;
class LauResonanceInfo;

class LauBlattWeisskopfFactor final {

  public:
    //! Define the allowed types of barrier factors
    enum BarrierType {
        BWBarrier, /*!< Blatt-Weisskopf barrier factor (for use when momentum terms not used in angular term) */
        BWPrimeBarrier, /*!< Blatt-Weisskopf barrier factor (for use when momentum terms are used in angular term) - the default */
        ExpBarrier      /*!< expoential barrier factor (mostly used for virtual contributions) */
    };

    //! Define the rest frame in which the momentum should be calculated (only relevant for bachelor)
    enum RestFrame {
        Covariant,     /*!< use expression from Covariant spin factor */
        ParentFrame,   /*!< momentum calculated in parent rest frame */
        ResonanceFrame /*!< momentum calculated in resonance rest frame */
    };

    //! Define resonance categories that will share common barrier factor radii
    enum BlattWeisskopfCategory {
        Default,    //*!< indicates that LauResonanceMaker should use the appropriate category for the given resonance */
        Parent,     //*!< indicates that this is the factor for the decay of the parent particle */
        Indep,      //*!< indicates that this resonance should be independent of all others */
        Light,    //*!< default category for light unflavoured states, e.g. rho(77), f0(980), etc. */
        Kstar,            //*!< default category for K* states */
        Charm,            //*!< default category for D* states */
        StrangeCharm,     //*!< default category for Ds* states */
        Charmonium,       //*!< default category for ccbar states */
        Beauty,           //*!< default category for B* states */
        StrangeBeauty,    //*!< default category for Bs* states */
        CharmBeauty,      //*!< default category for Bc* states */
        Custom1,          //*!< user-defined category */
        Custom2,          //*!< user-defined category */
        Custom3,          //*!< user-defined category */
        Custom4           //*!< user-defined category */
    };

    //! Constructor
    LauBlattWeisskopfFactor( const LauResonanceInfo& resInfo,
                             const BarrierType barrierType,
                             const RestFrame restFrame,
                             const BlattWeisskopfCategory category );

    //! Constructor
    LauBlattWeisskopfFactor( const LauResonanceInfo& resInfo,
                             const Double_t resRadius,
                             const BarrierType barrierType,
                             const RestFrame restFrame,
                             const BlattWeisskopfCategory category );

    //! Constructor
    LauBlattWeisskopfFactor( const Int_t spin,
                             const Double_t resRadius,
                             const BarrierType barrierType,
                             const RestFrame restFrame,
                             const BlattWeisskopfCategory category );

    //! Method to create a new factor with cloned radius parameter
    /*!
        \param newSpin the value of the spin to use for the created instance
        \param newBarrierType the BarrierType to use for the created instance
    */
    LauBlattWeisskopfFactor* createClone( const UInt_t newSpin, const BarrierType newBarrierType );

    //! Retrieve the radius parameter
    const LauParameter* getRadiusParameter() const { return radius_; }

    //! Retrieve the radius parameter
    LauParameter* getRadiusParameter() { return radius_; }

    //! Retrieve the barrier type
    BarrierType getBarrierType() const { return barrierType_; }

    //! Retrieve the rest frame information
    RestFrame getRestFrame() const { return restFrame_; }

    //! Calculate form factor value
    /*!
        \param p the value of the momentum
    */
    Double_t calcFormFactor( const Double_t p ) const;

  private:
    //! Copy constructor
    LauBlattWeisskopfFactor( const LauBlattWeisskopfFactor& other,
                             const UInt_t newSpin,
                             const BarrierType newBarrierType );

    //! Copy assignment operator (not implemented)
    LauBlattWeisskopfFactor& operator=( const LauBlattWeisskopfFactor& other );

    //! Set the name of the radius parameter
    TString setRadiusName( const LauResonanceInfo& resInfo, const BlattWeisskopfCategory category );

    //! Set the name of the radius parameter
    TString setRadiusName( const BlattWeisskopfCategory category );

    //! Set the name of the radius parameter
    TString setRadiusName( const TString& categoryName );

    //! Resonance spin
    const UInt_t spin_;

    //! Radius parameter
    LauParameter* radius_;

    //! Barrier type
    const BarrierType barrierType_;

    //! Rest frame
    const RestFrame restFrame_;

    ClassDef( LauBlattWeisskopfFactor, 0 )
};

#endif
