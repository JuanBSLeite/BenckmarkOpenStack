
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

/*! \file LauIntegrals.cc
    \brief File containing implementation of LauIntegrals class.
*/

#include <iostream>
using std::cout;
using std::endl;

#include "LauConstants.hh"
#include "LauIntegrals.hh"

#include "TMath.h"

LauIntegrals::LauIntegrals( Double_t weightsPrecision ) :
    weightsPrecision_( weightsPrecision )
{
}

LauIntegrals::~LauIntegrals()
{
}

LauIntegrals::LauIntegrals( const LauIntegrals& rhs ) :
    weightsPrecision_( rhs.weightsPrecision_ )
{
}

LauIntegrals& LauIntegrals::operator=( const LauIntegrals& rhs )
{
    if ( &rhs != this ) {
        weightsPrecision_ = rhs.weightsPrecision_;
    }
    return *this;
}

void LauIntegrals::calcGaussLegendreWeights( const Int_t numPoints,
                                             std::vector<Double_t>& abscissas,
                                             std::vector<Double_t>& weights )
{
    // Calculate the Gauss-Legendre weights that will be used for the
    // simple Gaussian integration method given the number of points

    abscissas.clear();
    weights.clear();
    abscissas.resize( numPoints );
    weights.resize( numPoints );

    Int_t m = ( numPoints + 1 ) / 2;
    Double_t dnumPoints( numPoints );
    Double_t dnumPointsPlusHalf = static_cast<Double_t>( numPoints + 0.5 );
    Int_t i( 0 ), j( 0 );
    Double_t p1( 0.0 ), p2( 0.0 ), p3( 0.0 ), pp( 0.0 ), z( 0.0 ), zSq( 0.0 ), z1( 0.0 ), di( 0.0 );

    for ( i = 1; i <= m; i++ ) {

        di  = static_cast<Double_t>( i );
        z   = TMath::Cos( LauConstants::pi * ( di - 0.25 ) / dnumPointsPlusHalf );
        zSq = z * z;

        // Starting with the above approximation for the ith root, we enter the
        // main loop of refinement by Newton's method
        do {
            p1 = 1.0;
            p2 = 0.0;

            // Calculate the Legendre polynomial at z - recurrence relation
            for ( j = 1; j <= numPoints; j++ ) {
                p3 = p2;
                p2 = p1;
                p1 = ( ( 2.0 * j - 1.0 ) * z * p2 - ( j - 1.0 ) * p3 ) / static_cast<Double_t>( j );
            }
            // p1 = Legendre polynomial. Compute its derivative, pp.

            pp = dnumPoints * ( z * p1 - p2 ) / ( zSq - 1.0 );
            z1 = z;
            z  = z1 - p1 / pp;    // Newton's method
        } while ( TMath::Abs( z - z1 ) > weightsPrecision_ );

        // Scale the root to the desired interval
        // Remember that the vector entries start with 0, hence i-1 (where first i value is 1)
        abscissas[i - 1]         = z;
        abscissas[numPoints - i] = abscissas[i - 1];    // Symmetric abscissa
        weights[i - 1]           = 2.0 / ( ( 1.0 - zSq ) * pp * pp );
        weights[numPoints - i]   = weights[i - 1];    // Symmetric weight
    }
}

void LauIntegrals::calcGaussHermiteWeights( const Int_t numPoints,
                                            std::vector<Double_t>& abscissas,
                                            std::vector<Double_t>& weights )
{
    // Calculate the Gauss-Hermite weights that will be used for the
    // simple Gaussian integration method for a function weighted by an exp(-x**2) term
    // These weights are common to all integration intervals - therefore they
    // should only be calculated once - when the constructor is called, for example.

    abscissas.clear();
    weights.clear();
    abscissas.resize( numPoints );
    weights.resize( numPoints );

    Int_t m = ( numPoints + 1 ) / 2;
    Double_t dnumPoints( numPoints );
    Int_t i( 0 ), j( 0 ), its( 0 );
    Double_t p1( 0.0 ), p2( 0.0 ), p3( 0.0 ), pp( 0.0 ), z( 0.0 ), z1( 0.0 );
    Double_t numPointsTerm( 2 * numPoints + 1 );

    // The roots are symmetric about the origin. Therefore, we only need to find half of them
    // Loop over the desired roots.
    for ( i = 1; i <= m; i++ ) {

        if ( i == 1 ) {
            // Initial guess for largest root
            z = TMath::Sqrt( numPointsTerm ) - 1.85575 * TMath::Power( numPointsTerm, -0.16667 );
        } else if ( i == 2 ) {
            // Initial guess for second largest root
            z -= 1.14 * TMath::Power( dnumPoints, 0.426 ) / z;
        } else if ( i == 3 ) {
            // Initial guess for third largest root
            z = 1.86 * z - 0.86 * abscissas[0];
        } else if ( i == 4 ) {
            // Initial guess for fourth largest root
            z = 1.91 * z - 0.91 * abscissas[1];
        } else {
            // Initial guess for other roots
            z = 2.0 * z - abscissas[i - 3];
        }

        for ( its = 1; its <= 10; its++ ) {
            p1 = LauConstants::pim4;
            p2 = 0.0;

            for ( j = 1; j <= numPoints; j++ ) {
                Double_t dj( j );
                p3 = p2;
                p2 = p1;
                p1 = z * TMath::Sqrt( 2.0 / dj ) * p2 - TMath::Sqrt( ( dj - 1.0 ) / dj ) * p3;
            }
            // p1 is now the desired Hermite polynomial. Compute its derivative, pp
            pp = TMath::Sqrt( 2.0 * dnumPoints ) * p2;
            z1 = z;
            if ( pp > 1e-10 ) {
                z = z1 - p1 / pp;
            }
            if ( TMath::Abs( z - z1 ) < weightsPrecision_ )
                break;
        }

        // Scale the root to the desired interval
        // Remember that the vector entries start with 0, hence i-1 (where first i value is 1)
        abscissas[i - 1]         = z;
        abscissas[numPoints - i] = -z;    // Symmetric abscissa
        weights[i - 1]           = 2.0 / ( pp * pp );
        weights[numPoints - i]   = weights[i - 1];    // Symmetric weight
    }
}
