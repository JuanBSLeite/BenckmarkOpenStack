
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

/*! \file Lau2DAbsDPPdf.hh
    \brief File containing declaration of Lau2DAbsDPPdf class.
*/

/*! \class Lau2DAbsDPPdf
    \brief Pure abstract base class for defining a variation across a 2D DP.

    Pure abstract base class for defining a normalised variation across a 2D DP.
*/

#ifndef LAU_2DABS_DP_PDF
#define LAU_2DABS_DP_PDF

#include "Rtypes.h"

class Lau2DAbsDPPdf {

  public:
    //! Constructor
    Lau2DAbsDPPdf() {}

    //! Destructor
    virtual ~Lau2DAbsDPPdf() {}

    //! Retrieve maximum height
    virtual Double_t getMaxHeight() const = 0;

    //! Retrieve PDF normalisation
    virtual Double_t getHistNorm() const = 0;

    //! Perform the interpolation (unnormalised)
    /*!
        \param [in] x the x-axis value
        \param [in] y the y-aixs value
        \return the unnormalised PDF value
    */
    virtual Double_t interpolateXY( Double_t x, Double_t y ) const = 0;

    //! Perform the interpolation and divide by the normalisation
    /*!
        \param [in] x the x-axis abscissa value
        \param [in] y the y-axis abscissa value
        \return the normalised PDF value
    */
    virtual Double_t interpolateXYNorm( Double_t x, Double_t y ) const = 0;

  private:
    //! Copy constructor - not implemented
    Lau2DAbsDPPdf( const Lau2DAbsDPPdf& rhs );

    //! Copy assignment operator - not implemented
    Lau2DAbsDPPdf& operator=( const Lau2DAbsDPPdf& rhs );

    ClassDef( Lau2DAbsDPPdf, 0 )    // Abstract base class for 2D DP variation
};

#endif
