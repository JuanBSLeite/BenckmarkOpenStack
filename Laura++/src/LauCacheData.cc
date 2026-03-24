
/*
Copyright 2005 University of Warwick

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

/*! \file LauCacheData.cc
    \brief File containing implementation of LauCacheData class.
*/

#include "LauCacheData.hh"

LauCacheData::LauCacheData() :
    m13Sq_( 0.0 ),
    m23Sq_( 0.0 ),
    mPrime_( 0.0 ),
    thPrime_( 0.0 ),
    eff_( 0.0 ),
    scfFraction_( 0.0 ),
    jacobian_( 0.0 )
{
}

LauCacheData::~LauCacheData()
{
}

LauCacheData::LauCacheData( const LauCacheData& rhs )
{
    this->storem13Sq( rhs.retrievem13Sq() );
    this->storem23Sq( rhs.retrievem23Sq() );
    this->storemPrime( rhs.retrievemPrime() );
    this->storethPrime( rhs.retrievethPrime() );
    this->storeEff( rhs.retrieveEff() );
    this->storeScfFraction( rhs.retrieveScfFraction() );
    this->storeJacobian( rhs.retrieveJacobian() );
    this->storeRealAmp( rhs.retrieveRealAmp() );
    this->storeImagAmp( rhs.retrieveImagAmp() );
    this->storeIncohIntensities( rhs.retrieveIncohIntensities() );
}

LauCacheData& LauCacheData::operator=( const LauCacheData& rhs )
{
    if ( &rhs != this ) {
        this->storem13Sq( rhs.retrievem13Sq() );
        this->storem23Sq( rhs.retrievem23Sq() );
        this->storemPrime( rhs.retrievemPrime() );
        this->storethPrime( rhs.retrievethPrime() );
        this->storeEff( rhs.retrieveEff() );
        this->storeScfFraction( rhs.retrieveScfFraction() );
        this->storeJacobian( rhs.retrieveJacobian() );
        this->storeRealAmp( rhs.retrieveRealAmp() );
        this->storeImagAmp( rhs.retrieveImagAmp() );
        this->storeIncohIntensities( rhs.retrieveIncohIntensities() );
    }
    return *this;
}
