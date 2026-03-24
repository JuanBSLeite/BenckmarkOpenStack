
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

/*! \file Lau2DAbsHistDP.hh
    \brief File containing declaration of Lau2DAbsHistDP class.
*/

/*! \class Lau2DAbsHistDP
    \brief Abstract base class for defining a variation across a 2D DP based on a histogram.

    Abstract base class for defining an unnormalised variation across a 2D DP based on a histogram.
    Contains helper methods to vary bin contents in the input histogram within uncertainties and
    scale the input histogram to match a given average value.
    The histogram can be defined in the conventional DP (m13Sq vs m23Sq) or in the square DP and
    one half can be used to describe symmetric DPs.
*/

#ifndef LAU_2DABSHIST_DP
#define LAU_2DABSHIST_DP

#include "Lau2DAbsDP.hh"

class TH2;
class LauDaughters;
class LauKinematics;

class Lau2DAbsHistDP : public Lau2DAbsDP {

  public:
    //! Constructor
    /*!
        \param [in] daughters the daughter particles
        \param [in] useUpperHalfOnly boolean flag to specify that the supplied histogram contains only the upper half of a symmetric DP (or lower half if using square DP coordinates)
        \param [in] squareDP boolean flag to specify whether the supplied histogram is in square DP coordinates
    */
    Lau2DAbsHistDP( const LauDaughters* daughters,
                    Bool_t useUpperHalfOnly = kFALSE,
                    Bool_t squareDP         = kFALSE );

    //! Copy constructor
    //Lau2DAbsHistDP( const Lau2DAbsHistDP& rhs );

    //! Destructor
    virtual ~Lau2DAbsHistDP();

    //! Perform the interpolation
    /*!
        \param [in] x the x-axis value
        \param [in] y the y-axis value
        \return the interpolated histogram value
    */
    virtual Double_t interpolateXY( Double_t x, Double_t y ) const = 0;

    //! Determine whether the histogram is defined in the square DP
    /*
            \return kTRUE if the square DP is being used, kFALSE otherwise
        */
    virtual Bool_t usingSquareDP() const { return squareDP_; };

  protected:
    //! Fluctuate the contents of each histogram bin independently, in accordance with their errors
    /*!
        \param [in,out] hist the histogram
    */
    void doBinFluctuation( TH2* hist );

    //! Fluctuate the contents of each histogram bin independently, in accordance with their errors
    /*!
        \param [in,out] hist the histogram
        \param [in] errorHi the histogram containing the upper uncertainties on the bins
        \param [in] errorLo the histogram containing the lower uncertainties on the bins
    */
    void doBinFluctuation( TH2* hist, const TH2* errorHi, const TH2* errorLo );

    //! Rescale the histogram bin contents based on the desired average efficiency and its uncertainty
    /*!
        The new average is sampled from a Gaussian distribution G(x;avEff,avEffError).
        The histogram is then rescaled according to newAvg/oldAvg.

        \param [in,out] hist the histogram
        \param [in] avEff the desired average efficiency
        \param [in] avEffError the error on that efficiency
    */
    void raiseOrLowerBins( TH2* hist, const Double_t avEff, const Double_t avEffError );

    //! Compute the average bin content for bins within the kinematic boundary
    /*!
        This method just uses the raw bin contents with no interpolation

        \param [in] hist the histogram
        \return the average value over the DP
    */
    Double_t computeAverageContents( const TH2* hist ) const;

    //! Check whether the given co-ordinates are within the kinematic boundary
    /*!
        \param [in] x the x co-ordinate
        \param [in] y the y co-ordinate
        \return true if the co-ordinates are within the kinematic boundary, otherwise false
    */
    Bool_t withinDPBoundaries( Double_t x, Double_t y ) const;

    //! If only using the upper half of the (symmetric) DP then transform into the correct half
    /*!
        \param [in,out] x the x co-ordinate
        \param [in,out] y the y co-ordinate
    */
    void getUpperHalf( Double_t& x, Double_t& y ) const;

  private:
    //! Copy constructor - not implemented
    Lau2DAbsHistDP( const Lau2DAbsHistDP& rhs );

    //! Copy assignment operator - not implemented
    Lau2DAbsHistDP& operator=( const Lau2DAbsHistDP& rhs );

    //! Kinematics used to check events are within DP boundary
    const LauKinematics* kinematics_;

    //! Boolean for using the upper half of DP
    Bool_t upperHalf_;

    //! Boolean for using square DP variables
    Bool_t squareDP_;

    ClassDef( Lau2DAbsHistDP, 0 )    // Abstract base class for 2D DP variations based on a histogram
};

#endif
