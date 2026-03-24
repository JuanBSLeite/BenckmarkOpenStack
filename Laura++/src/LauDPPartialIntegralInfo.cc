
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

/*! \file LauDPPartialIntegralInfo.cc
    \brief File containing implementation of LauDPPartialIntegralInfo class.
*/

#include "LauDPPartialIntegralInfo.hh"

#include "LauIntegrals.hh"
#include "LauKinematics.hh"

#include <iostream>

LauDPPartialIntegralInfo::LauDPPartialIntegralInfo( const Double_t minm13,
                                                    const Double_t maxm13,
                                                    const Double_t minm23,
                                                    const Double_t maxm23,
                                                    const Double_t m13BinWidth,
                                                    const Double_t m23BinWidth,
                                                    const Double_t precision,
                                                    const UInt_t nAmp,
                                                    const UInt_t nIncohAmp,
                                                    const Bool_t squareDP,
                                                    const LauKinematics* kinematics ) :
    minm13_( minm13 ),
    maxm13_( maxm13 ),
    minm23_( minm23 ),
    maxm23_( maxm23 ),
    m13BinWidth_( m13BinWidth ),
    m23BinWidth_( m23BinWidth ),
    nm13Points_( static_cast<UInt_t>( ( maxm13 - minm13 ) / m13BinWidth ) ),
    nm23Points_( static_cast<UInt_t>( ( maxm23 - minm23 ) / m23BinWidth ) ),
    nAmp_( nAmp ),
    nIncohAmp_( nIncohAmp ),
    squareDP_( squareDP )
{
    const Double_t meanm13      = 0.5 * ( minm13 + maxm13 );
    const Double_t rangem13     = maxm13 - minm13;
    const Double_t halfRangem13 = 0.5 * rangem13;

    const Double_t meanm23      = 0.5 * ( minm23 + maxm23 );
    const Double_t rangem23     = maxm23 - minm23;
    const Double_t halfRangem23 = 0.5 * rangem23;

    const Double_t intFactor = halfRangem13 * halfRangem23;

    // Raise error if squareDP is true but the kinematics object is not provided
    if ( squareDP_ && kinematics == 0 ) {
        std::cerr << "ERROR in LauDPPartialIntegralInfo constructor : Integration in the square DP has been specified but no valid kinematics object has been provided!"
                  << std::endl;
        return;
    }

    // Avoid integral if we have no points in either x or y space
    if ( nm13Points_ == 0 || nm23Points_ == 0 ) {
        std::cerr << "ERROR in LauDPPartialIntegralInfo constructor : Range has zero grid points in one or both of the dimensions!"
                  << std::endl;
        return;
    }

    // Print a warning if we have a very large number of points
    if ( ( nm13Points_ * nm23Points_ ) > 8000000 ) {
        std::cerr << "WARNING in LauDPPartialIntegralInfo constructor : The integration binning scheme has a very large number of bins, this could cause high memory consumption!"
                  << std::endl;
        std::cerr << "                                                : In case of problems, consider using LauIsobarDynamics::setNarrowResonanceThreshold and/or LauIsobarDynamics::setIntegralBinningFactor to tune the binning behaviour."
                  << std::endl;
    }

    LauIntegrals dpIntegrals( precision );
    dpIntegrals.calcGaussLegendreWeights( nm13Points_, m13Points_, m13Weights_ );
    dpIntegrals.calcGaussLegendreWeights( nm23Points_, m23Points_, m23Weights_ );

    // Print out total weights for the integration
    Double_t totm13Weight( 0.0 ), totm23Weight( 0.0 );
    for ( UInt_t i = 0; i < nm13Points_; ++i ) {
        totm13Weight += m13Weights_[i];
    }
    for ( UInt_t i = 0; i < nm23Points_; ++i ) {
        totm23Weight += m23Weights_[i];
    }

    if ( squareDP_ ) {
        std::cout << "INFO in LauDPPartialIntegralInfo constructor : nmPrimePoints = " << nm13Points_
                  << ", nthPrimePoints = " << nm23Points_ << std::endl;
        std::cout << "                                             : mPrimeBinWidth = "
                  << m13BinWidth_ << ", thPrimeBinWidth = " << m23BinWidth_ << std::endl;
        std::cout << "                                             : Integrating over mPrime = "
                  << minm13_ << " to " << maxm13_ << ", thPrime = " << minm23_ << " to " << maxm23_
                  << std::endl;
        std::cout << "                                             : totmPrimeWeight = "
                  << totm13Weight << ", totthPrimeWeight = " << totm23Weight << std::endl;
    } else {
        std::cout << "INFO in LauDPPartialIntegralInfo constructor : nm13Points = " << nm13Points_
                  << ", nm23Points = " << nm23Points_ << std::endl;
        std::cout << "                                             : m13BinWidth = " << m13BinWidth_
                  << ", m23BinWidth = " << m23BinWidth_ << std::endl;
        std::cout << "                                             : Integrating over m13 = " << minm13_
                  << " to " << maxm13_ << ", m23 = " << minm23_ << " to " << maxm23_ << std::endl;
        std::cout << "                                             : totm13Weight = " << totm13Weight
                  << ", totm23Weight = " << totm23Weight << std::endl;
    }

    // Calculate the m13 and m23 values at the grid points

    UInt_t midpoint = ( nm13Points_ + 1 ) / 2;
    for ( UInt_t i = 0; i < midpoint; ++i ) {

        UInt_t ii = nm13Points_ - 1 - i;    // symmetric i index

        Double_t dm13   = halfRangem13 * m13Points_[i];
        Double_t m13Val = meanm13 - dm13;
        m13Points_[i]   = m13Val;

        m13Val         = meanm13 + dm13;
        m13Points_[ii] = m13Val;
    }

    midpoint = ( nm23Points_ + 1 ) / 2;
    for ( UInt_t i = 0; i < midpoint; i++ ) {

        UInt_t ii = nm23Points_ - 1 - i;    // symmetric i index

        Double_t dm23   = halfRangem23 * m23Points_[i];
        Double_t m23Val = meanm23 - dm23;
        m23Points_[i]   = m23Val;

        m23Val         = meanm23 + dm23;
        m23Points_[ii] = m23Val;
    }

    // Now compute the combined weights at each grid point
    weights_.resize( nm13Points_ );
    efficiencies_.resize( nm13Points_ );
    amplitudes_.resize( nm13Points_ );
    incohIntensities_.resize( nm13Points_ );
    for ( UInt_t i = 0; i < nm13Points_; ++i ) {

        weights_[i].resize( nm23Points_ );
        efficiencies_[i].resize( nm23Points_ );
        amplitudes_[i].resize( nm23Points_ );
        incohIntensities_[i].resize( nm23Points_ );

        for ( UInt_t j = 0; j < nm23Points_; ++j ) {

            Double_t weight = m13Weights_[i] * m23Weights_[j];

            Double_t jacobian( 0.0 );
            if ( squareDP_ ) {
                jacobian = kinematics->calcSqDPJacobian( m13Points_[i], m23Points_[j] );
            } else {
                jacobian = 4.0 * m13Points_[i] * m23Points_[j];
            }
            weight *= ( jacobian * intFactor );

            weights_[i][j] = weight;

            amplitudes_[i][j].resize( nAmp_ );
            incohIntensities_[i][j].resize( nIncohAmp_ );

        }    // j weights loop
    }        // i weights loop
}

LauDPPartialIntegralInfo::~LauDPPartialIntegralInfo()
{
}

std::ostream& operator<<( std::ostream& stream, const LauDPPartialIntegralInfo& infoRecord )
{
    stream << "minm13 = " << infoRecord.getMinm13() << ", ";
    stream << "maxm13 = " << infoRecord.getMaxm13() << ", ";
    stream << "minm23 = " << infoRecord.getMinm23() << ", ";
    stream << "maxm23 = " << infoRecord.getMaxm23() << ", ";
    stream << "m13BinWidth = " << infoRecord.getM13BinWidth() << ", ";
    stream << "m23BinWidth = " << infoRecord.getM23BinWidth() << std::endl;
    return stream;
}
