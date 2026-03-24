
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

/*! \file Lau2DCubicSpline.hh
    \brief File containing declaration of Lau2DCubicSpline class.
*/

/*! \class Lau2DCubicSpline
    \brief Class for defining a 2D cubic spline based on an input histogram

    Class for defining a 2D cubic spline based on RooBinned2DBicubicBase by
    Manuel Tobias Schiller <manuel.schiller@nikhef.nl> (2012-08-29).

    Smoothly interpolate between bin midpoints of a 2D histogram with
    equal-sized bins

    The interpolation function coincides with the saved histogram values in bin
    midpoints; a cubic interpolation polynomial is used. Bins must all have the
    same size. The interpolation polynomial in a "cell" between four bin
    midpoints is:

    \f[ p(x, y) = \sum_{i=0}^3 \sum_{j=0}^3 a_{ij} x^i y^j \f]

    The coefficients \f$a_{ij}\f$ are determined by the requirement that
    function value, first derivatives and the mixed second derivative must agree
    with that of the binned function at the four bin midpoints at the edges of
    the cell. Symmetric finite differences are used to approximate these
    derivatives.

    For each cell, the coefficients are determined at constuction time. The
    object also keeps a cache of 2D integrals over complete cells, such that 2D
    integrations can be done analytically in a reasonably short amount of time.
*/

#ifndef LAU_2DCUBICSPLINE
#define LAU_2DCUBICSPLINE

#include "Rtypes.h"

#include <vector>

class TH2;

class Lau2DCubicSpline {
  private:
    //! Length of coefficient record in array
    enum {
        NCoeff      = 16,
        CoeffRecLen = 17
    };

  public:
    //! Constructor from histogram
    /*!
        \param h the histogram
    */
    Lau2DCubicSpline( const TH2& h );

    //! Destructor
    virtual ~Lau2DCubicSpline();

    //! Evaluate the function at given point
    /*!
        \param [in] x the x co-ordinate
        \param [in] y the y co-ordinate
    */
    virtual Double_t evaluate( Double_t x, Double_t y ) const;

    //! Evaluate analytical integral in x, y, or x and y
    /*!
        \param [in] x1 the lower x limit
        \param [in] x2 the upper x limit
        \param [in] y1 the lower y limit
        \param [in] y2 the upper y limit
    */
    virtual Double_t analyticalIntegral( Double_t x1, Double_t x2, Double_t y1, Double_t y2 ) const;

    //! Evaluate analytical integral across the whole function
    virtual Double_t analyticalIntegral() const;

  private:
    //! Copy constructor - not implemented
    Lau2DCubicSpline( const Lau2DCubicSpline& rhs );

    //! Copy assignment operator - not implemented
    Lau2DCubicSpline& operator=( const Lau2DCubicSpline& rhs );

    //! Number of bins in x axis
    Int_t nBinsX;

    //! Number of bins in y axis
    Int_t nBinsY;

    //! Bin size in x
    Double_t binSizeX;

    //! Bin size in y
    Double_t binSizeY;

    //! Minimum x value
    Double_t xmin;

    //! Maximum x value
    Double_t xmax;

    //! Minimum y value
    Double_t ymin;

    //! Maximum y value
    Double_t ymax;

    //! Coefficients of interpolation polynomials
    std::vector<Double_t> coeffs;

    //! Get contents of a given bin from a histogram
    /*!
        \param [in] h the histogram
        \param [in] xbin the x bin index
        \param [in] ybin the y bin index
        \return the bin contents
    */
    Double_t histcont( const TH2& h, Int_t xbin, Int_t ybin ) const;

    //! Get d/dx finite difference in a given bin from a histogram
    /*!
        \param [in] h the histogram
        \param [in] xbin the x bin index
        \param [in] ybin the y bin index
        \return the d/dx finite difference
    */
    Double_t dhistdx( const TH2& h, Int_t xbin, Int_t ybin ) const;

    //! Get d/dy finite difference in a given bin from a histogram
    /*!
        \param [in] h the histogram
        \param [in] xbin the x bin index
        \param [in] ybin the y bin index
        \return the d/dy finite difference
    */
    Double_t dhistdy( const TH2& h, Int_t xbin, Int_t ybin ) const;

    //! Get d^2/dydx finite difference in a given bin from a histogram
    /*!
        \param [in] h the histogram
        \param [in] xbin the x bin index
        \param [in] ybin the y bin index
        \return the d^2/dydx finite difference
    */
    Double_t d2histdxdy( const TH2& h, Int_t xbin, Int_t ybin ) const;

    //! Const access to coefficients
    /*!
        \param [in] binx the x bin index
        \param [in] biny the y bin index
        \param [in] theCoeff the coefficient index
        \return the coefficient
    */
    inline const Double_t& coeff( Int_t binx, Int_t biny, Int_t theCoeff ) const
    {
        return coeffs[theCoeff + CoeffRecLen * ( binx + nBinsX * biny )];
    }

    //! Access to coefficients
    /*!
        \param [in] binx the x bin index
        \param [in] biny the y bin index
        \param [in] theCoeff the coefficient index
        \return the coefficient
    */
    inline Double_t& coeff( Int_t binx, Int_t biny, Int_t theCoeff )
    {
        return coeffs[theCoeff + CoeffRecLen * ( binx + nBinsX * biny )];
    }

    //! Evaluate integral over x at a given y from (x1, y) to (x2, y)
    /*!
        \param [in] x1 lower x limit
        \param [in] x2 upper x limit
        \param [in] y y value
        \return integral over x
    */
    Double_t evalX( Double_t x1, Double_t x2, Double_t y ) const;

    //! Evaluate integral over y at a given x from (x, y1) to (x, y2)
    /*!
        \param [in] x x value
        \param [in] y1 lower y limit
        \param [in] y2 upper y limit
        \return integral over y
    */
    Double_t evalY( Double_t x, Double_t y1, Double_t y2 ) const;

    //! Evaluate integral over x and y from (x1, y1) to (x2, y2)
    /*!
        \param [in] x1 lower x limit
        \param [in] x2 upper x limit
        \param [in] y1 lower y limit
        \param [in] y2 upper y limit
        \return integral over x and y
    */
    Double_t evalXY( Double_t x1, Double_t x2, Double_t y1, Double_t y2 ) const;

    ClassDef( Lau2DCubicSpline,
              0 );    // Class for defining a 2D cubic spline based on an input histogram
};

#endif
