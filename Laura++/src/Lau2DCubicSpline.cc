
/*
Copyright 2013 University of Warwick

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

/*! \file Lau2DCubicSpline.cc
    \brief File containing implementation of Lau2DCubicSpline class.

    Class for defining a 2D cubic spline based on RooBinned2DBicubicBase by
    Manuel Tobias Schiller <manuel.schiller@nikhef.nl> (2012-08-29).
*/

#include "Lau2DCubicSpline.hh"

#include "TH2.h"
#include "TSystem.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

Lau2DCubicSpline::~Lau2DCubicSpline()
{
}

inline Double_t Lau2DCubicSpline::histcont( const TH2& h, Int_t xbin, Int_t ybin ) const
{
    //reflect until we're in range
    while ( xbin < 0 || xbin >= nBinsX - 1 ) {
        if ( xbin < 0 )
            xbin = -xbin - 1;
        if ( xbin >= nBinsX - 1 )
            xbin = 2 * ( nBinsX - 1 ) - xbin - 1;
    }
    while ( ybin < 0 || ybin >= nBinsY - 1 ) {
        if ( ybin < 0 )
            ybin = -ybin - 1;
        if ( ybin >= nBinsY - 1 )
            ybin = 2 * ( nBinsY - 1 ) - ybin - 1;
    }

    return h.GetBinContent( 1 + xbin, 1 + ybin );
}

inline Double_t Lau2DCubicSpline::dhistdx( const TH2& h, Int_t xbin, Int_t ybin ) const
{
    return 0.5 * ( histcont( h, xbin + 1, ybin ) - histcont( h, xbin - 1, ybin ) );
}

inline Double_t Lau2DCubicSpline::dhistdy( const TH2& h, Int_t xbin, Int_t ybin ) const
{
    return 0.5 * ( histcont( h, xbin, ybin + 1 ) - histcont( h, xbin, ybin - 1 ) );
}

inline Double_t Lau2DCubicSpline::d2histdxdy( const TH2& h, Int_t xbin, Int_t ybin ) const
{
    return 0.5 * ( histcont( h, xbin - 1, ybin - 1 ) - histcont( h, xbin + 1, ybin - 1 ) +
                   histcont( h, xbin + 1, ybin + 1 ) - histcont( h, xbin - 1, ybin + 1 ) );
}

Lau2DCubicSpline::Lau2DCubicSpline( const TH2& h ) :
    nBinsX( 1 + h.GetNbinsX() ),
    nBinsY( 1 + h.GetNbinsY() ),
    binSizeX( h.GetXaxis()->GetBinWidth( 1 ) ),
    binSizeY( h.GetYaxis()->GetBinWidth( 1 ) ),
    xmin( h.GetXaxis()->GetBinCenter( 1 ) - binSizeX ),
    xmax( h.GetXaxis()->GetBinCenter( nBinsX - 1 ) + binSizeX ),
    ymin( h.GetYaxis()->GetBinCenter( 1 ) - binSizeY ),
    ymax( h.GetYaxis()->GetBinCenter( nBinsY - 1 ) + binSizeY ),
    coeffs( CoeffRecLen * nBinsX * nBinsY )
{
    const TAxis* xaxis = h.GetXaxis();
    const TAxis* yaxis = h.GetYaxis();
    // verify that all bins have same size
    for ( Int_t i = 1; i < nBinsX; ++i ) {
        if ( std::abs( xaxis->GetBinWidth( i ) / binSizeX - 1. ) > 1e-9 ) {
            std::cerr << "ERROR in Lau2DCubicSpline constructor : the histogram has variable bin sizes."
                      << std::endl;
            gSystem->Exit( EXIT_FAILURE );
        }
    }
    for ( Int_t i = 1; i < nBinsY; ++i ) {
        if ( std::abs( yaxis->GetBinWidth( i ) / binSizeY - 1. ) > 1e-9 ) {
            std::cerr << "ERROR in Lau2DCubicSpline constructor : the histogram has variable bin sizes."
                      << std::endl;
            gSystem->Exit( EXIT_FAILURE );
        }
    }
    // ok, go through histogram to precalculate the interpolation coefficients
    // in rectangles between bin centres
    //
    // for that purpose, we map each of those rectangles to the unit square
    for ( Int_t j = -1; j < nBinsY - 1; ++j ) {
        for ( Int_t i = -1; i < nBinsX - 1; ++i ) {
            const Double_t rhs[NCoeff] = { // function values in bin centres
                                           histcont( h, i, j ),
                                           histcont( h, i + 1, j ),
                                           histcont( h, i, j + 1 ),
                                           histcont( h, i + 1, j + 1 ),
                                           // df/dx in bin centres (finite difference approximation)
                                           dhistdx( h, i, j ),
                                           dhistdx( h, i + 1, j ),
                                           dhistdx( h, i, j + 1 ),
                                           dhistdx( h, i + 1, j + 1 ),
                                           // df/dy in bin centres (finite difference approximation)
                                           dhistdy( h, i, j ),
                                           dhistdy( h, i + 1, j ),
                                           dhistdy( h, i, j + 1 ),
                                           dhistdy( h, i + 1, j + 1 ),
                                           // d^2f/dxdy in bin centres (finite difference approximation)
                                           d2histdxdy( h, i, j ),
                                           d2histdxdy( h, i + 1, j ),
                                           d2histdxdy( h, i, j + 1 ),
                                           d2histdxdy( h, i + 1, j + 1 ) };
            // work out solution - strange array placement is due to the fact
            // that terms with x/y to high powers can be small, so they should
            // be added up first during evaluation to avoid cancellation
            // issues; at the same time you want to access them in order to
            // not confuse the CPU cache, so they're stored back to front
            //
            // a_00 ... a_30
            coeff( 1 + i, 1 + j, 15 ) = rhs[0];
            coeff( 1 + i, 1 + j, 14 ) = rhs[4];
            coeff( 1 + i, 1 + j, 13 ) = 3. * ( -rhs[0] + rhs[1] ) - 2. * rhs[4] - rhs[5];
            coeff( 1 + i, 1 + j, 12 ) = 2. * ( rhs[0] - rhs[1] ) + rhs[4] + rhs[5];
            // a_31 ... a_31
            coeff( 1 + i, 1 + j, 11 ) = rhs[8];
            coeff( 1 + i, 1 + j, 10 ) = rhs[12];
            coeff( 1 + i, 1 + j, 9 )  = 3. * ( -rhs[8] + rhs[9] ) - 2. * rhs[12] - rhs[13];
            coeff( 1 + i, 1 + j, 8 )  = 2. * ( rhs[8] - rhs[9] ) + rhs[12] + rhs[13];
            // a_02 ... a_32
            coeff( 1 + i, 1 + j, 7 ) = 3. * ( -rhs[0] + rhs[2] ) - 2. * rhs[8] - rhs[10];
            coeff( 1 + i, 1 + j, 6 ) = 3. * ( -rhs[4] + rhs[6] ) - 2. * rhs[12] - rhs[14];
            coeff( 1 + i, 1 + j, 5 ) = 9. * ( rhs[0] - rhs[1] - rhs[2] + rhs[3] ) +
                                       6. * ( rhs[4] - rhs[6] + rhs[8] - rhs[9] ) + 4. * rhs[12] +
                                       3. * ( rhs[5] - rhs[7] + rhs[10] - rhs[11] ) +
                                       2. * ( rhs[13] + rhs[14] ) + rhs[15];
            coeff( 1 + i, 1 + j, 4 ) = 6. * ( -rhs[0] + rhs[1] + rhs[2] - rhs[3] ) +
                                       4. * ( -rhs[8] + rhs[9] ) +
                                       3. * ( -rhs[4] - rhs[5] + rhs[6] + rhs[7] ) +
                                       2. * ( -rhs[10] + rhs[11] - rhs[12] - rhs[13] ) - rhs[14] -
                                       rhs[15];
            // a_03 ... a_33
            coeff( 1 + i, 1 + j, 3 ) = 2. * ( rhs[0] - rhs[2] ) + rhs[8] + rhs[10];
            coeff( 1 + i, 1 + j, 2 ) = 2. * ( rhs[4] - rhs[6] ) + rhs[12] + rhs[14];
            coeff( 1 + i, 1 + j, 1 ) = 6. * ( -rhs[0] + rhs[1] + rhs[2] - rhs[3] ) +
                                       4. * ( -rhs[4] + rhs[6] ) +
                                       3. * ( -rhs[8] + rhs[9] - rhs[10] + rhs[11] ) +
                                       2. * ( -rhs[5] + rhs[7] - rhs[12] - rhs[14] ) - rhs[13] -
                                       rhs[15];
            coeff( 1 + i, 1 + j, 0 ) = 4. * ( rhs[0] - rhs[1] - rhs[2] + rhs[3] ) +
                                       2. * ( rhs[4] + rhs[5] - rhs[6] - rhs[7] + rhs[8] - rhs[9] +
                                              rhs[10] - rhs[11] ) +
                                       rhs[12] + rhs[13] + rhs[14] + rhs[15];
            // coeff(1 + i, 1 + j, 17) contains integral of function over the
            // square in "unit square coordinates", i.e. neglecting bin widths
            // this is done to help speed up calculations of 2D integrals
            Double_t sum = 0.;
            for ( Int_t k = 0; k < NCoeff; ++k )
                sum += coeff( 1 + i, 1 + j, k ) / Double_t( ( 4 - ( k % 4 ) ) * ( 4 - ( k / 4 ) ) );
            coeff( 1 + i, 1 + j, NCoeff ) = sum;
        }
    }
}

Double_t Lau2DCubicSpline::evaluate( Double_t x, Double_t y ) const
{
    // protect against NaN and out of range
    if ( x <= xmin || x >= xmax || y <= ymin || y >= ymax || x != x || y != y )
        return 0.;
    // find the bin in question
    const Int_t binx = Int_t( Double_t( nBinsX ) * ( x - xmin ) / ( xmax - xmin ) );
    const Int_t biny = Int_t( Double_t( nBinsY ) * ( y - ymin ) / ( ymax - ymin ) );
    // get low edge of bin
    const Double_t xlo = Double_t( nBinsX - binx ) / Double_t( nBinsX ) * xmin +
                         Double_t( binx ) / Double_t( nBinsX ) * xmax;
    const Double_t ylo = Double_t( nBinsY - biny ) / Double_t( nBinsY ) * ymin +
                         Double_t( biny ) / Double_t( nBinsY ) * ymax;
    // normalise to coordinates in unit sqare
    const Double_t hx = ( x - xlo ) / binSizeX;
    const Double_t hy = ( y - ylo ) / binSizeY;
    // monomials
    const Double_t hxton[4] = { hx * hx * hx, hx * hx, hx, 1. };
    const Double_t hyton[4] = { hy * hy * hy, hy * hy, hy, 1. };
    // sum up
    Double_t retVal = 0.;
    for ( Int_t k = 0; k < NCoeff; ++k )
        retVal += coeff( binx, biny, k ) * hxton[k % 4] * hyton[k / 4];

    return retVal;
}

Double_t Lau2DCubicSpline::analyticalIntegral() const
{
    return evalXY( xmin, xmax, ymin, ymax );
}

Double_t Lau2DCubicSpline::analyticalIntegral( Double_t x1, Double_t x2, Double_t y1, Double_t y2 ) const
{
    if ( y1 == y2 )
        return evalX( x1, x2, y1 );
    if ( x1 == x2 )
        return evalY( x1, y1, y2 );
    return evalXY( x1, x2, y1, y2 );
}

Double_t Lau2DCubicSpline::evalX( Double_t x1, Double_t x2, Double_t y ) const
{
    // protect against NaN
    if ( x1 != x1 || x2 != x2 || y != y )
        return 0.;
    // find the bin in question
    const Int_t biny = Int_t( Double_t( nBinsY ) * ( y - ymin ) / ( ymax - ymin ) );
    // get low edge of bin
    const Double_t ylo = Double_t( nBinsY - biny ) / Double_t( nBinsY ) * ymin +
                         Double_t( biny ) / Double_t( nBinsY ) * ymax;
    // normalise to coordinates in unit sqare
    const Double_t hy = ( y - ylo ) / binSizeY;
    // monomials
    const Double_t hyton[4] = { hy * hy * hy, hy * hy, hy, 1. };
    // integral
    Double_t sum = 0.;
    for ( Int_t binx = 0; binx < nBinsX; ++binx ) {
        // get low/high edge of bin
        const Double_t xhi = Double_t( nBinsX - binx - 1 ) / Double_t( nBinsX ) * xmin +
                             Double_t( binx + 1 ) / Double_t( nBinsX ) * xmax;
        if ( xhi < x1 )
            continue;
        const Double_t xlo = Double_t( nBinsX - binx ) / Double_t( nBinsX ) * xmin +
                             Double_t( binx ) / Double_t( nBinsX ) * xmax;
        if ( xlo > x2 )
            break;
        // work out integration range
        const Double_t a = ( ( xlo > x1 ) ? 0. : ( x1 - xlo ) ) / binSizeX;
        const Double_t b = ( ( xhi < x2 ) ? binSizeX : ( x2 - xlo ) ) / binSizeX;
        // integrated monomials
        const Double_t hxton[4] = { 0.25 * ( b * b * b * b - a * a * a * a ),
                                    ( b * b * b - a * a * a ) / 3.,
                                    0.5 * ( b * b - a * a ),
                                    b - a };
        Double_t lsum           = 0.;
        for ( Int_t k = 0; k < NCoeff; ++k )
            lsum += coeff( binx, biny, k ) * hxton[k % 4] * hyton[k / 4];
        sum += lsum;
    }
    // move from unit square coordinates to user coordinates
    return sum * binSizeX;
}

Double_t Lau2DCubicSpline::evalY( Double_t x, Double_t y1, Double_t y2 ) const
{
    // protect against NaN
    if ( x != x || y1 != y1 || y2 != y2 )
        return 0.;
    // find the bin in question
    const Int_t binx = Int_t( Double_t( nBinsX ) * ( x - xmin ) / ( xmax - xmin ) );
    // get low edge of bin
    const Double_t xlo = Double_t( nBinsX - binx ) / Double_t( nBinsX ) * xmin +
                         Double_t( binx ) / Double_t( nBinsX ) * xmax;
    // normalise to coordinates in unit sqare
    const Double_t hx = ( x - xlo ) / binSizeX;
    // monomials
    const Double_t hxton[4] = { hx * hx * hx, hx * hx, hx, 1. };
    // integral
    Double_t sum = 0.;
    for ( Int_t biny = 0; biny < nBinsY; ++biny ) {
        // get low/high edge of bin
        const Double_t yhi = Double_t( nBinsY - biny - 1 ) / Double_t( nBinsY ) * ymin +
                             Double_t( biny + 1 ) / Double_t( nBinsY ) * ymax;
        if ( yhi < y1 )
            continue;
        const Double_t ylo = Double_t( nBinsY - biny ) / Double_t( nBinsY ) * ymin +
                             Double_t( biny ) / Double_t( nBinsY ) * ymax;
        if ( ylo > y2 )
            break;
        // work out integration range
        const Double_t a = ( ( ylo > y1 ) ? 0. : ( y1 - ylo ) ) / binSizeY;
        const Double_t b = ( ( yhi < y2 ) ? binSizeY : ( y2 - ylo ) ) / binSizeY;
        // integrated monomials
        const Double_t hyton[4] = { 0.25 * ( b * b * b * b - a * a * a * a ),
                                    ( b * b * b - a * a * a ) / 3.,
                                    0.5 * ( b * b - a * a ),
                                    b - a };
        Double_t lsum           = 0.;
        for ( Int_t k = 0; k < NCoeff; ++k )
            lsum += coeff( binx, biny, k ) * hxton[k % 4] * hyton[k / 4];
        sum += lsum;
    }
    // move from unit square coordinates to user coordinates
    return sum * binSizeY;
}

Double_t Lau2DCubicSpline::evalXY( Double_t x1, Double_t x2, Double_t y1, Double_t y2 ) const
{
    // protect against NaN
    if ( x1 != x1 || x2 != x2 || y1 != y1 || y2 != y2 )
        return 0.;
    // integral
    Double_t sum = 0.;
    for ( Int_t biny = 0; biny < nBinsY; ++biny ) {
        // get low/high edge of bin
        const Double_t yhi = Double_t( nBinsY - biny - 1 ) / Double_t( nBinsY ) * ymin +
                             Double_t( biny + 1 ) / Double_t( nBinsY ) * ymax;
        if ( yhi < y1 )
            continue;
        const Double_t ylo = Double_t( nBinsY - biny ) / Double_t( nBinsY ) * ymin +
                             Double_t( biny ) / Double_t( nBinsY ) * ymax;
        if ( ylo > y2 )
            break;
        // work out integration range
        const Double_t ay            = ( ( ylo > y1 ) ? 0. : ( y1 - ylo ) ) / binSizeY;
        const Double_t by            = ( ( yhi < y2 ) ? binSizeY : ( y2 - ylo ) ) / binSizeY;
        const Bool_t yFullyContained = std::abs( by - ay - 1.0 ) < 1e-15;
        // integrated monomials
        const Double_t hyton[4] = { 0.25 * ( by * by * by * by - ay * ay * ay * ay ),
                                    ( by * by * by - ay * ay * ay ) / 3.,
                                    0.5 * ( by * by - ay * ay ),
                                    by - ay };
        for ( Int_t binx = 0; binx < nBinsX; ++binx ) {
            // get low/high edge of bin
            const Double_t xhi = Double_t( nBinsX - binx - 1 ) / Double_t( nBinsX ) * xmin +
                                 Double_t( binx + 1 ) / Double_t( nBinsX ) * xmax;
            if ( xhi < x1 )
                continue;
            const Double_t xlo = Double_t( nBinsX - binx ) / Double_t( nBinsX ) * xmin +
                                 Double_t( binx ) / Double_t( nBinsX ) * xmax;
            if ( xlo > x2 )
                break;
            // work out integration range
            const Double_t ax            = ( ( xlo > x1 ) ? 0. : ( x1 - xlo ) ) / binSizeX;
            const Double_t bx            = ( ( xhi < x2 ) ? binSizeX : ( x2 - xlo ) ) / binSizeX;
            const Bool_t xFullyContained = std::abs( bx - ax - 1.0 ) < 1e-15;
            if ( xFullyContained && yFullyContained ) {
                // for fully contained bins, we have cached the integral
                sum += coeff( binx, biny, NCoeff );
                continue;
            }
            // integrated monomials
            const Double_t hxton[4] = { 0.25 * ( bx * bx * bx * bx - ax * ax * ax * ax ),
                                        ( bx * bx * bx - ax * ax * ax ) / 3.,
                                        0.5 * ( bx * bx - ax * ax ),
                                        bx - ax };
            // integrate over bin

            Double_t lsum = 0.;
            for ( Int_t k = 0; k < NCoeff; ++k )
                lsum += coeff( binx, biny, k ) * hxton[k % 4] * hyton[k / 4];
            sum += lsum;
        }
    }
    // move from unit square coordinates to user coordinates
    return sum * binSizeX * binSizeY;
}
