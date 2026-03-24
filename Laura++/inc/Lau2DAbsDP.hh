
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

/*! \file Lau2DAbsDP.hh
    \brief File containing declaration of Lau2DAbsDP class.
*/

/*! \class Lau2DAbsDP
    \brief Pure abstract base class for defining a variation across a 2D DP.

    Pure abstract base class for defining an unnormalised variation across a 2D DP.
*/

#ifndef LAU_2DABS_DP
#define LAU_2DABS_DP

#include "Rtypes.h"

class Lau2DAbsDP {

  public:
    //!Constructor
    Lau2DAbsDP() {}

    //! Destructor
    virtual ~Lau2DAbsDP() {}

    //! Perform the interpolation
    /*!
        \param [in] x the x-axis value
        \param [in] y the y-aixs value
        \return the interpolated histogram value
    */
    virtual Double_t interpolateXY( Double_t x, Double_t y ) const = 0;

    //! Determine whether the histogram is defined in the square DP
    /*
            \return kTRUE if the square DP is being used, kFALSE otherwise
        */
    virtual Bool_t usingSquareDP() const = 0;

  private:
    //! Copy constructor - not implemented
    Lau2DAbsDP( const Lau2DAbsDP& rhs );

    //! Copy assignment operator - not implemented
    Lau2DAbsDP& operator=( const Lau2DAbsDP& rhs );

    ClassDef( Lau2DAbsDP, 0 )    // Abstract base class for 2D DP variation
};

#endif
