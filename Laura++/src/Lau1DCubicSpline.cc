
/*
Copyright 2015 University of Warwick

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

/*! \file Lau1DCubicSpline.cc
    \brief File containing implementation of Lau1DCubicSpline class.
*/

#include "Lau1DCubicSpline.hh"

#include "TH2.h"
#include "TMath.h"
#include "TSystem.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

Lau1DCubicSpline::Lau1DCubicSpline( const std::vector<Double_t>& xs,
                                    const std::vector<Double_t>& ys,
                                    LauSplineType type,
                                    LauSplineBoundaryType leftBound,
                                    LauSplineBoundaryType rightBound,
                                    Double_t dydx0,
                                    Double_t dydxn ) :
    nKnots_( xs.size() ),
    x_( xs ),
    y_( ys ),
    type_( type ),
    leftBound_( leftBound ),
    rightBound_( rightBound ),
    dydx0_( dydx0 ),
    dydxn_( dydxn )
{
    init();
}

Lau1DCubicSpline::~Lau1DCubicSpline()
{
}

Double_t Lau1DCubicSpline::evaluate( Double_t x ) const
{
    // do not attempt to extrapolate the spline
    if ( x < x_[0] || x > x_[nKnots_ - 1] ) {
        std::cout << "WARNING in Lau1DCubicSpline::evaluate : function is only defined between "
                  << x_[0] << " and " << x_[nKnots_ - 1] << std::endl;
        std::cout << "                                        value at " << x << " returned as 0"
                  << std::endl;
        return 0.;
    }

    // first determine which 'cell' of the spline x is in
    // cell i runs from knot i to knot i+1
    Int_t cell( 0 );

    while ( x > x_[cell + 1] ) {
        ++cell;
    }

    // obtain x- and y-values of the neighbouring knots
    Double_t xLow  = x_[cell];
    Double_t xHigh = x_[cell + 1];
    Double_t yLow  = y_[cell];
    Double_t yHigh = y_[cell + 1];

    if ( type_ == Lau1DCubicSpline::LinearInterpolation ) {
        return yHigh * ( x - xLow ) / ( xHigh - xLow ) + yLow * ( xHigh - x ) / ( xHigh - xLow );
    }

    // obtain t, the normalised x-coordinate within the cell,
    // and the coefficients a and b, which are defined in cell i as:
    //
    //      a_i =  k_i  *(x_i+1 - x_i) - (y_i+1 - y_i),
    //      b_i = -k_i+1*(x_i+1 - x_i) + (y_i+1 - y_i)
    //
    // where k_i is (by construction) the first derivative at knot i
    Double_t t = ( x - xLow ) / ( xHigh - xLow );
    Double_t a = dydx_[cell] * ( xHigh - xLow ) - ( yHigh - yLow );
    Double_t b = -1. * dydx_[cell + 1] * ( xHigh - xLow ) + ( yHigh - yLow );

    Double_t retVal = ( 1 - t ) * yLow + t * yHigh + t * ( 1 - t ) * ( a * ( 1 - t ) + b * t );

    return retVal;
}

void Lau1DCubicSpline::updateYValues( const std::vector<Double_t>& ys )
{
    y_ = ys;
    this->calcDerivatives();
}

void Lau1DCubicSpline::updateType( LauSplineType type )
{
    if ( type_ != type ) {
        type_ = type;
        this->calcDerivatives();
    }
}

void Lau1DCubicSpline::updateBoundaryConditions( LauSplineBoundaryType leftBound,
                                                 LauSplineBoundaryType rightBound,
                                                 Double_t dydx0,
                                                 Double_t dydxn )
{
    Bool_t updateDerivatives( kFALSE );

    if ( leftBound_ != leftBound || rightBound_ != rightBound ) {
        leftBound_        = leftBound;
        rightBound_       = rightBound;
        updateDerivatives = kTRUE;
    }

    if ( dydx0_ != dydx0 ) {
        dydx0_ = dydx0;
        if ( leftBound_ == Lau1DCubicSpline::Clamped )
            updateDerivatives = kTRUE;
    }

    if ( dydxn_ != dydxn ) {
        dydxn_ = dydxn;
        if ( rightBound_ == Lau1DCubicSpline::Clamped )
            updateDerivatives = kTRUE;
    }

    if ( updateDerivatives ) {
        this->calcDerivatives();
    }
}

void Lau1DCubicSpline::init()
{
    if ( y_.size() != x_.size() ) {
        std::cout << "ERROR in Lau1DCubicSpline::init : The number of y-values given does not match the number of x-values"
                  << std::endl;
        std::cout << "                                  Found " << y_.size() << ", expected "
                  << x_.size() << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    dydx_.insert( dydx_.begin(), nKnots_, 0. );

    a_.insert( a_.begin(), nKnots_, 0. );
    b_.insert( b_.begin(), nKnots_, 0. );
    c_.insert( c_.begin(), nKnots_, 0. );
    d_.insert( d_.begin(), nKnots_, 0. );

    this->calcDerivatives();
}

void Lau1DCubicSpline::calcDerivatives()
{
    switch ( type_ ) {
        case Lau1DCubicSpline::StandardSpline :
            this->calcDerivativesStandard();
            break;
        case Lau1DCubicSpline::AkimaSpline :
            this->calcDerivativesAkima();
            break;
        case Lau1DCubicSpline::LinearInterpolation :
            //derivatives not needed for linear interpolation
            break;
    }
}

void Lau1DCubicSpline::calcDerivativesStandard()
{
    // derivatives are determined such that the second derivative is continuous at internal knots

    // derivatives, k_i, are the solutions to a set of linear equations of the form:
    // a_i * k_i-1  +  b_i * k+i  +  c_i * k_i+1  =  d_i with a_0 = 0, c_n-1 = 0

    // this is solved using the tridiagonal matrix algorithm as on en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm

    // first and last equations give boundary conditions
    // - for natural boundary, require f''(x) = 0 at end knot
    // - for 'not a knot' boundary, require f'''(x) continuous at second knot
    // - for clamped boundary, require predefined value of f'(x) at end knot

    // non-zero values of a_0 and c_n-1 would give cyclic boundary conditions
    a_[0]           = 0.;
    c_[nKnots_ - 1] = 0.;

    // set left boundary condition
    if ( leftBound_ == Lau1DCubicSpline::Natural ) {
        b_[0] = 2. / ( x_[1] - x_[0] );
        c_[0] = 1. / ( x_[1] - x_[0] );
        d_[0] = 3. * ( y_[1] - y_[0] ) / ( ( x_[1] - x_[0] ) * ( x_[1] - x_[0] ) );

    } else if ( leftBound_ == Lau1DCubicSpline::NotAKnot ) {
        // define the width, h, and the 'slope', delta, of the first cell
        Double_t h1( x_[1] - x_[0] ), h2( x_[2] - x_[0] );
        Double_t delta1( ( y_[1] - y_[0] ) / h1 ), delta2( ( y_[2] - y_[1] ) / h2 );

        // these coefficients can be determined by requiring f'''_0(x_1) = f'''_1(x_1)
        // the requirement f''_0(x_1) = f''_1(x_1) has been used to remove the dependence on k_2
        b_[0] = h2;
        c_[0] = h1 + h2;
        d_[0] = delta1 * ( 2. * h2 * h2 + 3. * h1 * h2 ) / ( h1 + h2 ) +
                delta2 * 5. * h1 * h1 / ( h1 + h2 );

    } else {    //Clamped
        b_[0] = 1.;
        c_[0] = 0.;
        d_[0] = dydx0_;
    }

    // set right boundary condition
    if ( rightBound_ == Lau1DCubicSpline::Natural ) {
        a_[nKnots_ - 1] = 1. / ( x_[nKnots_ - 1] - x_[nKnots_ - 2] );
        b_[nKnots_ - 1] = 2. / ( x_[nKnots_ - 1] - x_[nKnots_ - 2] );
        d_[nKnots_ - 1] = 3. * ( y_[nKnots_ - 1] - y_[nKnots_ - 2] ) /
                          ( ( x_[nKnots_ - 1] - x_[nKnots_ - 2] ) *
                            ( x_[nKnots_ - 1] - x_[nKnots_ - 2] ) );

    } else if ( rightBound_ == Lau1DCubicSpline::NotAKnot ) {
        // define the width, h, and the 'slope', delta, of the last cell
        Double_t hnm1( x_[nKnots_ - 1] - x_[nKnots_ - 2] ), hnm2( x_[nKnots_ - 2] - x_[nKnots_ - 3] );
        Double_t deltanm1( ( y_[nKnots_ - 1] - y_[nKnots_ - 2] ) / hnm1 ),
            deltanm2( ( y_[nKnots_ - 2] - y_[nKnots_ - 3] ) / hnm2 );

        // these coefficients can be determined by requiring f'''_n-3(x_n-2) = f'''_n-2(x_n-2)
        // the requirement f''_n-3(x_n-2) = f''_n-2(x_n-2) has been used to remove
        // the dependence on k_n-3
        a_[nKnots_ - 1] = hnm2 + hnm1;
        b_[nKnots_ - 1] = hnm1;
        d_[nKnots_ - 1] = deltanm2 * hnm1 * hnm1 / ( hnm2 + hnm1 ) +
                          deltanm1 * ( 2. * hnm2 * hnm2 + 3. * hnm2 * hnm1 ) / ( hnm2 + hnm1 );

    } else {    //Clamped
        a_[nKnots_ - 1] = 0.;
        b_[nKnots_ - 1] = 1.;
        d_[nKnots_ - 1] = dydxn_;
    }

    // the remaining equations ensure that f_i-1''(x_i) = f''_i(x_i) for all internal knots
    for ( UInt_t i = 1; i < nKnots_ - 1; ++i ) {
        a_[i] = 1. / ( x_[i] - x_[i - 1] );
        b_[i] = 2. / ( x_[i] - x_[i - 1] ) + 2. / ( x_[i + 1] - x_[i] );
        c_[i] = 1. / ( x_[i + 1] - x_[i] );
        d_[i] = 3. * ( y_[i] - y_[i - 1] ) / ( ( x_[i] - x_[i - 1] ) * ( x_[i] - x_[i - 1] ) ) +
                3. * ( y_[i + 1] - y_[i] ) / ( ( x_[i + 1] - x_[i] ) * ( x_[i + 1] - x_[i] ) );
    }

    // forward sweep - replace c_i and d_i with
    //
    // c'_i =  c_i /  b_i                 for i = 0
    //         c_i / (b_i - a_i * c_i-1)  otherwise
    //
    // and
    //
    // d'_i =   d_i                /  b_i                 for i = 0
    //         (d_i - a_i * d_i-1) / (b_i - a_i * c_i-1)  otherwise

    c_[0] /= b_[0];
    d_[0] /= b_[0];

    for ( UInt_t i = 1; i < nKnots_ - 1; ++i ) {
        c_[i] = c_[i] / ( b_[i] - a_[i] * c_[i - 1] );
        d_[i] = ( d_[i] - a_[i] * d_[i - 1] ) / ( b_[i] - a_[i] * c_[i - 1] );
    }

    d_[nKnots_ - 1] = ( d_[nKnots_ - 1] - a_[nKnots_ - 1] * d_[nKnots_ - 2] ) /
                      ( b_[nKnots_ - 1] - a_[nKnots_ - 1] * c_[nKnots_ - 2] );

    // back substitution - calculate k_i = dy/dx at each knot
    //
    // k_i = d'_i                   for i = n-1
    //       d'_i - c'_i * k_i+1    otherwise

    dydx_[nKnots_ - 1] = d_[nKnots_ - 1];

    for ( Int_t i = nKnots_ - 2; i >= 0; --i ) {
        dydx_[i] = d_[i] - c_[i] * dydx_[i + 1];
    }
}

void Lau1DCubicSpline::calcDerivativesAkima()
{
    //derivatives are calculated according to the Akima method
    // J.ACM vol. 17 no. 4 pp 589-602

    Double_t am1( 0. ), an( 0. ), anp1( 0. );

    // a[i] is the slope of the segment from point i-1 to point i
    //
    // n.b. segment 0 is before point 0 and segment n is after point n-1
    //      internal segments are numbered 1 - n-1
    for ( UInt_t i = 1; i < nKnots_; ++i ) {
        a_[i] = ( y_[i] - y_[i - 1] ) / ( x_[i] - x_[i - 1] );
    }

    // calculate slopes between additional points on each end according to the Akima method
    // method assumes that the additional points follow a quadratic defined by the last three points
    // this leads to the relations a[2] - a[1] = a[1] - a[0] = a[0] - a[-1] and a[n-1] - a[n-2] = a[n] - a[n-1] = a[n+1] - a[n]
    a_[0] = 2 * a_[1] - a_[2];
    am1   = 2 * a_[0] - a_[1];

    an   = 2 * a_[nKnots_ - 1] - a_[nKnots_ - 2];
    anp1 = 2 * an - a_[nKnots_ - 1];

    // b[i] is the weight of a[i]   towards dydx[i]
    // c[i] is the weight of a[i+1] towards dydx[i]
    // See Appendix A of J.ACM vol. 17 no. 4 pp 589-602 for a justification of these weights
    b_[0] = TMath::Abs( a_[1] - a_[0] );
    b_[1] = TMath::Abs( a_[2] - a_[1] );

    c_[0] = TMath::Abs( a_[0] - am1 );
    c_[1] = TMath::Abs( a_[1] - a_[0] );

    for ( UInt_t i = 2; i < nKnots_ - 2; ++i ) {
        b_[i] = TMath::Abs( a_[i + 2] - a_[i + 1] );
        c_[i] = TMath::Abs( a_[i] - a_[i - 1] );
    }

    b_[nKnots_ - 2] = TMath::Abs( an - a_[nKnots_ - 1] );
    b_[nKnots_ - 1] = TMath::Abs( anp1 - an );

    c_[nKnots_ - 2] = TMath::Abs( a_[nKnots_ - 2] - a_[nKnots_ - 3] );
    c_[nKnots_ - 1] = TMath::Abs( a_[nKnots_ - 1] - a_[nKnots_ - 2] );

    // dy/dx calculated as the weighted average of a[i] and a[i+1]:
    // dy/dx_i = ( | a_i+2 - a_i+1 | a_i  +  | a_i - a_i-1 | a_i+1 ) / ( | a_i+2 - a_i+1 |  +  | a_i - a_i-1 | )
    // in the special case a_i-1 == a_i != a_i+1 == a_i+2 this function is undefined so dy/dx is then defined as (a_i + a_i+1) / 2
    for ( UInt_t i = 0; i < nKnots_ - 2; ++i ) {
        if ( b_[i] == 0 && c_[i] == 0 ) {
            dydx_[i] = ( a_[i] + a_[i + 1] ) / 2.;
        } else {
            dydx_[i] = ( b_[i] * a_[i] + c_[i] * a_[i + 1] ) / ( b_[i] + c_[i] );
        }
    }

    if ( b_[nKnots_ - 1] == 0 && c_[nKnots_ - 1] == 0 ) {
        dydx_[nKnots_ - 1] = ( a_[nKnots_ - 1] + an ) / 2.;
    } else {
        dydx_[nKnots_ - 1] = ( b_[nKnots_ - 1] * a_[nKnots_ - 1] + c_[nKnots_ - 1] * an ) /
                             ( b_[nKnots_ - 1] + c_[nKnots_ - 1] );
    }
}
