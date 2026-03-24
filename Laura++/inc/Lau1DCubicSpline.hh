
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

/*! \file Lau1DCubicSpline.hh
    \brief File containing declaration of Lau1DCubicSpline class.
*/

/*! \class Lau1DCubicSpline
    \brief Class for defining a 1D cubic spline based on a set of knots.

    Class for defining a 1D cubic spline based on a set of knots.
    Interpolation between the knots is performed either by one of
    two types of cubic spline (standard or Akima) or by linear interpolation.
    The splines are defined by a piecewise cubic function which, between knots i and i+1, has the form

    f_i(x) = (1 - t)*y_i + t*y_i+1 + t*(1 - t)(a*(1 - t) + b*t)

    where t = (x - x_i)/(x_i+1 - x_i),
    a =  k_i  *(x_i+1 - x_i) - (y_i+1 - y_i),
    b = -k_i+1*(x_i+1 - x_i) + (y_i+1 - y_i)

    and k_i is (by construction) the first derivative at knot i.
    f(x) and f'(x) are continuous at the internal knots by construction.

    For the standard splines, f''(x) is required to be continuous
    at all internal knots placing n-2 constraints on the n parameters, k_i.
    The final two constraints are set by the boundary conditions.
    At each boundary, the function may be:

    (i)   Clamped :     f'(x) = C at the last knot
    (ii)  Natural :     f''(x) = 0 at the last knot
    (iii) Not a knot :  f'''(x) continuous at the second last knot

    The algorithms used in these splines can be found on:
    http://en.wikipedia.org/wiki/Spline_interpolation#Algorithm_to_find_the_interpolating_cubic_spline
    http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm

    For the Akima splines, the values of k_i are determined from the slopes of the four nearest segments (a_i-1, a_i, a_i+1 and a_i+2) as

    k_i = ( | a_i+2 - a_i+1 | a_i  +  | a_i - a_i-1 | a_i+1 ) / ( | a_i+2 - a_i+1 |  +  | a_i - a_i-1 | )

    and as

    k_i = ( a_i + a_i+1 ) / 2

    in the special case a_i-1 == a_i != a_i+1 == a_i+2.

    Boundary conditions are specified by the relations
    a_2 - a_1 = a_1 - a_0 = a_0 - a_-1 and
    a_n-1 - a_n-2 = a_n - a_n-1 = a_n+1 - a_n

    The algorithms used in these splines can be found in:
    J.ACM vol. 17 no. 4 pp 589-602
*/

#ifndef LAU_1DCUBICSPLINE
#define LAU_1DCUBICSPLINE

#include "Rtypes.h"

#include <vector>

class Lau1DCubicSpline {

  public:
    //! Define the allowed interpolation types
    enum LauSplineType {
        StandardSpline, /*!< standard cubic splines with f''(x) continuous at all internal knots */
        AkimaSpline, /*!< Akima cubic splines with f'(x) at each knot defined locally by the positions of only five knots */
        LinearInterpolation /*! Linear interpolation between each pair of knots */
    };

    //! Define the allowed boundary condition types
    /*!
        These are only supported by standard splines
    */
    enum LauSplineBoundaryType {
        Clamped, /*!< clamped boundary - f'(x) = C */
        Natural, /*!< natural boundary - f''(x) = 0 */
        NotAKnot /*!< 'not a knot' boundary - f'''(x) continuous at second last knot */
    };

    //! Constructor
    /*!
        \param [in] xs the x-values of the knots
        \param [in] ys the y-values of the knots
        \param [in] type the type of the spline (e.g. StandardSpline, AkimaSpline, LinearInterpolation)
        \param [in] leftBound the left-hand boundary condition
        \param [in] rightBound the right-hand boundary condition
        \param [in] dydx0 the gradient at the left-hand boundary of a clamped spline
        \param [in] dydxn the gradient at the right-hand boundary of a clamped spline
    */
    Lau1DCubicSpline( const std::vector<Double_t>& xs,
                      const std::vector<Double_t>& ys,
                      LauSplineType type               = Lau1DCubicSpline::StandardSpline,
                      LauSplineBoundaryType leftBound  = Lau1DCubicSpline::NotAKnot,
                      LauSplineBoundaryType rightBound = Lau1DCubicSpline::NotAKnot,
                      Double_t dydx0                   = 0.0,
                      Double_t dydxn                   = 0.0 );

    //! Destructor
    virtual ~Lau1DCubicSpline();

    //! Evaluate the function at given point
    /*!
        \param [in] x the x-coordinate
        \return the value of the spline at x
    */
    Double_t evaluate( Double_t x ) const;

    //! Update the y-values of the knots
    /*!
        \param [in] ys the y-values of the knots
    */
    void updateYValues( const std::vector<Double_t>& ys );

    //! Update the type of interpolation to perform
    /*!
        \param [in] type the type of interpolation
    */
    void updateType( LauSplineType type );

    //! Update the boundary conditions for the spline
    /*!
        \param [in] leftBound the left-hand boundary condition
        \param [in] rightBound the right-hand boundary condition
        \param [in] dydx0 the gradient at the left-hand boundary of a clamped spline
        \param [in] dydxn the gradient at the right-hand boundary of a clamped spline
    */
    void updateBoundaryConditions( LauSplineBoundaryType leftBound,
                                   LauSplineBoundaryType rightBound,
                                   Double_t dydx0 = 0.0,
                                   Double_t dydxn = 0.0 );

  private:
    //! Copy constructor - not implemented
    Lau1DCubicSpline( const Lau1DCubicSpline& rhs );

    //! Copy assignment operator - not implemented
    Lau1DCubicSpline& operator=( const Lau1DCubicSpline& rhs );

    //! Initialise the class
    void init();

    //! Calculate the first derivative at each knot
    void calcDerivatives();

    //! Calculate the first derivatives according to the standard method
    void calcDerivativesStandard();
    //! Calculate the first derivatives according to the Akima method
    void calcDerivativesAkima();

    //! The number of knots in the spline
    const UInt_t nKnots_;

    //! The x-value at each knot
    std::vector<Double_t> x_;
    //! The y-value at each knot
    std::vector<Double_t> y_;
    //! The first derivative at each knot
    std::vector<Double_t> dydx_;

    //! The 'a' coefficients used to determine the derivatives
    std::vector<Double_t> a_;
    //! The 'b' coefficients used to determine the derivatives
    std::vector<Double_t> b_;
    //! The 'c' coefficients used to determine the derivatives
    std::vector<Double_t> c_;
    //! The 'd' coefficients used to determine the derivatives
    std::vector<Double_t> d_;

    //! The type of interpolation to be performed
    LauSplineType type_;

    //! The left-hand boundary condition on the spline
    LauSplineBoundaryType leftBound_;
    //! The right-hand boundary condition on the spline
    LauSplineBoundaryType rightBound_;

    //! The gradient at the left boundary for a clamped spline
    Double_t dydx0_;
    //! The gradient at the right boundary for a clamped spline
    Double_t dydxn_;

    ClassDef( Lau1DCubicSpline, 0 );    // Class for defining a 1D cubic spline
};

#endif
