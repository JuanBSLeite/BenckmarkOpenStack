
/*
Copyright 2006 University of Warwick

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

/*! \file LauResonanceInfo.cc
    \brief File containing implementation of LauResonanceInfo class.
*/

#include "LauResonanceInfo.hh"

#include "LauParameter.hh"

#include <iostream>

LauResonanceInfo::LauResonanceInfo( const TString& name,
                                    const Double_t mass,
                                    const Double_t width,
                                    const Int_t spin,
                                    const Int_t charge,
                                    const LauBlattWeisskopfFactor::BlattWeisskopfCategory bwCategory,
                                    const Double_t bwRadius ) :
    name_( name ),
    sanitisedName_( "" ),
    mass_( 0 ),
    width_( 0 ),
    spin_( spin ),
    charge_( charge ),
    bwCategory_( bwCategory ),
    bwRadius_( bwRadius ),
    conjugate_( 0 ),
    extraPars_()
{
    this->sanitiseName();
    mass_ = new LauParameter( sanitisedName_ + "_MASS", mass, 0.0, 3.0 * mass, kTRUE );
    mass_->secondStage( kTRUE );
    width_ = new LauParameter( sanitisedName_ + "_WIDTH", width, 0.0, 3.0 * width, kTRUE );
    width_->secondStage( kTRUE );
}

LauResonanceInfo::~LauResonanceInfo()
{
    delete mass_;
    mass_ = 0;
    delete width_;
    width_ = 0;

    for ( std::set<LauParameter*>::iterator iter = extraPars_.begin(); iter != extraPars_.end();
          ++iter ) {
        delete ( *iter );
    }
    extraPars_.clear();
}

LauResonanceInfo::LauResonanceInfo( const LauResonanceInfo& other,
                                    const TString& newName,
                                    const Int_t newCharge ) :
    name_( newName ),
    sanitisedName_( "" ),
    mass_( 0 ),
    width_( 0 ),
    spin_( other.spin_ ),
    charge_( newCharge ),
    bwCategory_( other.bwCategory_ ),
    bwRadius_( other.bwRadius_ ),
    conjugate_( 0 ),
    extraPars_()
{
    this->sanitiseName();
    mass_  = other.mass_->createClone( sanitisedName_ + "_MASS" );
    width_ = other.width_->createClone( sanitisedName_ + "_WIDTH" );
    for ( std::set<LauParameter*>::iterator iter = other.extraPars_.begin();
          iter != other.extraPars_.end();
          ++iter ) {
        TString parName = ( *iter )->name();
        parName.Remove( 0, parName.Last( '_' ) );
        parName.Prepend( sanitisedName_ );
        LauParameter* par = ( *iter )->createClone( parName );
        extraPars_.insert( par );
    }
}

LauResonanceInfo* LauResonanceInfo::createChargeConjugate()
{
    Int_t newCharge = -charge_;

    TString newName( name_ );
    Ssiz_t index = newName.Index( "+" );
    if ( index != -1 ) {
        newName.Replace( index, 1, "-" );
    } else {
        index = newName.Index( "-" );
        if ( index != -1 ) {
            newName.Replace( index, 1, "+" );
        }
    }

    LauResonanceInfo* conjugate = new LauResonanceInfo( *this, newName, newCharge );
    conjugate->conjugate_       = this;
    this->conjugate_            = conjugate;

    return conjugate;
}

LauResonanceInfo* LauResonanceInfo::createSharedParameterRecord( const TString& name )
{
    LauResonanceInfo* newinfo = new LauResonanceInfo( *this, name, charge_ );

    sharedParRecords_.push_back( newinfo );

    return newinfo;
}

LauParameter* LauResonanceInfo::getExtraParameter( const TString& parName )
{
    LauParameter* par( 0 );
    for ( std::set<LauParameter*>::iterator iter = extraPars_.begin(); iter != extraPars_.end();
          ++iter ) {
        if ( ( *iter )->name() == parName ) {
            par = ( *iter );
        }
    }
    return par;
}

void LauResonanceInfo::addExtraParameter( LauParameter* param, const Bool_t independentPar )
{
    bool ok = extraPars_.insert( param ).second;
    if ( ! ok ) {
        std::cerr << "WARNING in LauResonanceInfo::addExtraParameter : parameter already present, not adding again"
                  << std::endl;
        return;
    }

    if ( conjugate_ != 0 ) {
        conjugate_->addCloneOfExtraParameter( param, independentPar );
    }

    for ( std::vector<LauResonanceInfo*>::iterator iter = sharedParRecords_.begin();
          iter != sharedParRecords_.end();
          ++iter ) {
        ( *iter )->addCloneOfExtraParameter( param, independentPar );
    }
}

void LauResonanceInfo::addCloneOfExtraParameter( LauParameter* param, const Bool_t copyNotClone )
{
    TString parName = param->name();
    parName.Remove( 0, parName.Last( '_' ) );
    parName.Prepend( sanitisedName_ );

    LauParameter* cloneParam( 0 );
    if ( copyNotClone ) {
        cloneParam = new LauParameter( parName,
                                       param->unblindValue(),
                                       param->minValue(),
                                       param->maxValue(),
                                       param->fixed() );
        cloneParam->secondStage( kTRUE );
    } else {
        cloneParam = param->createClone( parName );
    }
    extraPars_.insert( cloneParam );
}

std::ostream& operator<<( std::ostream& stream, const LauResonanceInfo& infoRecord )
{
    stream << infoRecord.getName() << ": ";
    stream << "mass = " << infoRecord.getMass()->value() << ", ";
    stream << "width = " << infoRecord.getWidth()->value() << ", ";
    stream << "spin = " << infoRecord.getSpin() << ", ";
    Int_t charge = infoRecord.getCharge();
    if ( charge < 0 ) {
        stream << "charge = " << infoRecord.getCharge() << ", ";
    } else {
        stream << "charge =  " << infoRecord.getCharge() << ", ";
    }
    stream << "BW radius = " << infoRecord.getBWRadius();

    return stream;
}

void LauResonanceInfo::sanitiseName()
{
    sanitisedName_ = name_;
    sanitisedName_ = sanitisedName_.ReplaceAll( "+", "p" );
    sanitisedName_ = sanitisedName_.ReplaceAll( "-", "m" );
    sanitisedName_ = sanitisedName_.ReplaceAll( "*", "st" );
    sanitisedName_ = sanitisedName_.ReplaceAll( "(", "_" );
    sanitisedName_ = sanitisedName_.ReplaceAll( ")", "_" );
    sanitisedName_ = sanitisedName_.ReplaceAll( "[", "_" );
    sanitisedName_ = sanitisedName_.ReplaceAll( "]", "_" );
    sanitisedName_ = sanitisedName_.ReplaceAll( "<", "_" );
    sanitisedName_ = sanitisedName_.ReplaceAll( ">", "_" );
    sanitisedName_ = sanitisedName_.ReplaceAll( "{", "_" );
    sanitisedName_ = sanitisedName_.ReplaceAll( "}", "_" );
    sanitisedName_ = sanitisedName_.ReplaceAll( " ", "_" );
    sanitisedName_ = sanitisedName_.ReplaceAll( "$", "" );
    sanitisedName_ = sanitisedName_.ReplaceAll( "%", "" );
    sanitisedName_ = sanitisedName_.ReplaceAll( "&", "" );
    sanitisedName_ = sanitisedName_.ReplaceAll( "/", "" );
    sanitisedName_ = sanitisedName_.ReplaceAll( ":", "" );
    sanitisedName_ = sanitisedName_.ReplaceAll( ";", "" );
    sanitisedName_ = sanitisedName_.ReplaceAll( "=", "" );
    sanitisedName_ = sanitisedName_.ReplaceAll( "\\", "" );
    sanitisedName_ = sanitisedName_.ReplaceAll( "^", "" );
    sanitisedName_ = sanitisedName_.ReplaceAll( "|", "" );
    sanitisedName_ = sanitisedName_.ReplaceAll( ",", "" );
    sanitisedName_ = sanitisedName_.ReplaceAll( ".", "" );
    sanitisedName_.Remove( TString::kBoth, '_' );
}
