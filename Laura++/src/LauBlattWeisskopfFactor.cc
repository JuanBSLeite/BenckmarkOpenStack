
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

/*! \file LauBlattWeisskopfFactor.cc
    \brief File containing implementation of LauBlattWeisskopfFactor class.
*/

#include "LauBlattWeisskopfFactor.hh"

#include "LauParameter.hh"
#include "LauResonanceInfo.hh"

#include "TMath.h"

#include <iostream>

// Constructor
LauBlattWeisskopfFactor::LauBlattWeisskopfFactor( const LauResonanceInfo& resInfo,
                                                  const BarrierType barrierType,
                                                  const RestFrame restFrame,
                                                  const BlattWeisskopfCategory category ) :
    spin_( resInfo.getSpin() ),
    radius_( new LauParameter( "NEED_A_GOOD_NAME", resInfo.getBWRadius(), 0.0, 10.0, kTRUE ) ),
    barrierType_( barrierType ),
    restFrame_( restFrame )
{
    TString categoryName = this->setRadiusName( resInfo, category );
    std::cout << "INFO in LauBlattWeisskopfFactor constructor : creating radius parameter for category \""
              << categoryName << "\", with initial value " << resInfo.getBWRadius() << std::endl;
}

LauBlattWeisskopfFactor::LauBlattWeisskopfFactor( const LauResonanceInfo& resInfo,
                                                  const Double_t resRadius,
                                                  const BarrierType barrierType,
                                                  const RestFrame restFrame,
                                                  const BlattWeisskopfCategory category ) :
    spin_( resInfo.getSpin() ),
    radius_( new LauParameter( "NEED_A_GOOD_NAME", resRadius, 0.0, 10.0, kTRUE ) ),
    barrierType_( barrierType ),
    restFrame_( restFrame )
{
    TString categoryName = this->setRadiusName( resInfo, category );
    std::cout << "INFO in LauBlattWeisskopfFactor constructor : creating radius parameter for category \""
              << categoryName << "\", with initial value " << resRadius << std::endl;
}

LauBlattWeisskopfFactor::LauBlattWeisskopfFactor( const Int_t spin,
                                                  const Double_t resRadius,
                                                  const BarrierType barrierType,
                                                  const RestFrame restFrame,
                                                  const BlattWeisskopfCategory category ) :
    spin_( spin ),
    radius_( new LauParameter( "NEED_A_GOOD_NAME", resRadius, 0.0, 10.0, kTRUE ) ),
    barrierType_( barrierType ),
    restFrame_( restFrame )
{
    TString categoryName = this->setRadiusName( category );
    std::cout << "INFO in LauBlattWeisskopfFactor constructor : creating radius parameter for category \""
              << categoryName << "\", with initial value " << resRadius << std::endl;
}

LauBlattWeisskopfFactor::LauBlattWeisskopfFactor( const LauBlattWeisskopfFactor& other,
                                                  const UInt_t newSpin,
                                                  const BarrierType newBarrierType ) :
    spin_( newSpin ),
    radius_( other.radius_->createClone() ),
    barrierType_( newBarrierType ),
    restFrame_( other.restFrame_ )
{
}

TString LauBlattWeisskopfFactor::setRadiusName( const LauResonanceInfo& resInfo,
                                                const BlattWeisskopfCategory category )
{
    switch ( category ) {
        case Indep :
            return this->setRadiusName( resInfo.getSanitisedName() );
        default :
            return this->setRadiusName( category );
    }
}

TString LauBlattWeisskopfFactor::setRadiusName( const BlattWeisskopfCategory category )
{
    switch ( category ) {
        case Default :
            return this->setRadiusName( "Unknown" );
        case Indep :
            // We shouldn't ever end up here
            return this->setRadiusName( "Unknown" );
        case Parent :
            return this->setRadiusName( "Parent" );
        case Light :
            return this->setRadiusName( "Light" );
        case Kstar :
            return this->setRadiusName( "Kstar" );
        case Charm :
            return this->setRadiusName( "Charm" );
        case StrangeCharm :
            return this->setRadiusName( "StrangeCharm" );
        case Charmonium :
            return this->setRadiusName( "Charmonium" );
        case Beauty :
            return this->setRadiusName( "Beauty" );
        case StrangeBeauty :
            return this->setRadiusName( "StrangeBeauty" );
        case CharmBeauty :
            return this->setRadiusName( "CharmBeauty" );
        case Custom1 :
            return this->setRadiusName( "Custom1" );
        case Custom2 :
            return this->setRadiusName( "Custom2" );
        case Custom3 :
            return this->setRadiusName( "Custom3" );
        case Custom4 :
            return this->setRadiusName( "Custom4" );
    }
    // We should never get here but gcc seems to think we can
    return this->setRadiusName( "Unknown" );
}

TString LauBlattWeisskopfFactor::setRadiusName( const TString& categoryName )
{
    TString name = "BarrierRadius_";
    name.Append( categoryName );

    radius_->name( name );

    return categoryName;
}

LauBlattWeisskopfFactor* LauBlattWeisskopfFactor::createClone( const UInt_t newSpin,
                                                               const BarrierType newBarrierType )
{
    LauBlattWeisskopfFactor* clone = new LauBlattWeisskopfFactor( *this, newSpin, newBarrierType );
    return clone;
}

Double_t LauBlattWeisskopfFactor::calcFormFactor( const Double_t p ) const
{
    // Calculate the requested form factor for the resonance, given the momentum value
    Double_t fFactor( 1.0 );

    // For scalars the form factor is always unity
    // TODO: and we currently don't have formulae for spin > 5
    if ( ( spin_ == 0 ) || ( spin_ > 5 ) ) {
        return fFactor;
    }

    const Double_t radius = radius_->unblindValue();
    const Double_t z      = radius * radius * p * p;

    if ( barrierType_ == BWBarrier ) {
        if ( spin_ == 1 ) {
            fFactor = TMath::Sqrt( 2.0 * z / ( z + 1.0 ) );
        } else if ( spin_ == 2 ) {
            fFactor = TMath::Sqrt( 13.0 * z * z / ( z * z + 3.0 * z + 9.0 ) );
        } else if ( spin_ == 3 ) {
            fFactor = TMath::Sqrt( 277.0 * z * z * z /
                                   ( z * z * z + 6.0 * z * z + 45.0 * z + 225.0 ) );
        } else if ( spin_ == 4 ) {
            fFactor = TMath::Sqrt(
                12746.0 * z * z * z * z /
                ( z * z * z * z + 10.0 * z * z * z + 135.0 * z * z + 1575.0 * z + 11025.0 ) );
        } else if ( spin_ == 5 ) {
            fFactor = TMath::Sqrt( 998881.0 * z * z * z * z * z /
                                   ( z * z * z * z * z + 15.0 * z * z * z * z + 315.0 * z * z * z +
                                     6300.0 * z * z + 99225.0 * z + 893025.0 ) );
        }
    } else if ( barrierType_ == BWPrimeBarrier ) {
        if ( spin_ == 1 ) {
            fFactor = TMath::Sqrt( 1.0 / ( z + 1.0 ) );
        } else if ( spin_ == 2 ) {
            fFactor = TMath::Sqrt( 1.0 / ( z * z + 3.0 * z + 9.0 ) );
        } else if ( spin_ == 3 ) {
            fFactor = TMath::Sqrt( 1.0 / ( z * z * z + 6.0 * z * z + 45.0 * z + 225.0 ) );
        } else if ( spin_ == 4 ) {
            fFactor = TMath::Sqrt(
                1.0 / ( z * z * z * z + 10.0 * z * z * z + 135.0 * z * z + 1575.0 * z + 11025.0 ) );
        } else if ( spin_ == 5 ) {
            fFactor = TMath::Sqrt( 1.0 /
                                   ( z * z * z * z * z + 15.0 * z * z * z * z + 315.0 * z * z * z +
                                     6300.0 * z * z + 99225.0 * z + 893025.0 ) );
        }
    } else if ( barrierType_ == ExpBarrier ) {
        if ( spin_ == 1 ) {
            fFactor = TMath::Exp( -TMath::Sqrt( z ) );
        } else if ( spin_ == 2 ) {
            fFactor = TMath::Exp( -z );
        } else if ( spin_ == 3 ) {
            fFactor = TMath::Exp( -TMath::Sqrt( z * z * z ) );
        } else if ( spin_ == 4 ) {
            fFactor = TMath::Exp( -z * z );
        } else if ( spin_ == 5 ) {
            fFactor = TMath::Exp( -TMath::Sqrt( z * z * z * z * z ) );
        }
    }

    return fFactor;
}
