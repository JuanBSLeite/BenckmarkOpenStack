
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

/*! \file Lau2DAbsHistDPPdf.hh
    \brief File containing declaration of Lau2DAbsHistDPPdf class.
*/

/*! \class Lau2DAbsHistDPPdf
    \brief Abstract base class for defining a variation across a 2D DP based on a histogram.

    Abstract base class for defining a normalised variation across a 2D DP based on a histogram.
    The returned values are normalised to the total area.
    The histogram can be defined in the conventional DP (m13Sq vs m23Sq) or in the square DP and
    one half can be used to describe symmetric DPs.
*/

#ifndef LAU_2DABSHIST_DP_PDF
#define LAU_2DABSHIST_DP_PDF

#include "Lau2DAbsDPPdf.hh"

class TH2;
class LauKinematics;
class LauVetoes;

class Lau2DAbsHistDPPdf : public Lau2DAbsDPPdf {

  public:
    //! Constructor
    /*!
        \param [in] kinematics the current DP kinematics
        \param [in] vetoes the vetoes within the DP
        \param [in] useUpperHalfOnly boolean flag to specify that the supplied histogram contains only the upper half of a symmetric DP (or lower half if using square DP coordinates)
        \param [in] squareDP boolean flag to specify whether the supplied histogram is in square DP coordinates
    */
    Lau2DAbsHistDPPdf( LauKinematics* kinematics,
                       const LauVetoes* vetoes,
                       Bool_t useUpperHalfOnly = kFALSE,
                       Bool_t squareDP         = kFALSE );

    //! Destructor
    virtual ~Lau2DAbsHistDPPdf();

    //! Retrieve maximum height of the PDF
    /*!
        \return the maximum height
    */
    Double_t getMaxHeight() const { return maxHeight_; }

    //! Retrieve PDF normalisation
    /*!
        \return the normalisation factor
    */
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

  protected:
    //! Get the kinematics object
    /*!
        \return the kinematics
    */
    const LauKinematics* getKinematics() const { return kinematics_; }

    //! Get the vetoes object
    /*!
        \return the vetoes
    */
    const LauVetoes* getVetoes() const { return vetoes_; }

    //! Calculate maximum height
    /*!
        \param [in,out] hist the histogram
    */
    void calcMaxHeight( TH2* hist );

    //! Fluctuate the histogram bin contents in accordance with their errors
    /*!
        \param [in,out] hist the histogram
    */
    void doBinFluctuation( TH2* hist );

    //! Check whether the given co-ordinates are within the kinematic boundary
    /*!
        \param [in] x the x co-ordinate
        \param [in] y the y co-ordinate
        \return true if the co-ordinates are within the kinematic boundary, otherwise false
    */
    Bool_t withinDPBoundaries( Double_t x, Double_t y ) const;

    //! Update the current co-ordinates in the kinematic space
    /*!
        \param [in] x the x co-ordinate
        \param [in] y the y co-ordinate
    */
    void updateKinematics( Double_t x, Double_t y ) const;

    //! If only using the upper half of the (symmetric) DP then transform into the correct half
    /*!
        \param [in,out] x the x co-ordinate
        \param [in,out] y the y co-ordinate
    */
    void getUpperHalf( Double_t& x, Double_t& y ) const;

  private:
    //! Copy constructor - not implemented
    Lau2DAbsHistDPPdf( const Lau2DAbsHistDPPdf& rhs );

    //! Copy assignment operator - not implemented
    Lau2DAbsHistDPPdf& operator=( const Lau2DAbsHistDPPdf& rhs );

    //! DP kinematics
    LauKinematics* kinematics_;

    //! Vetos within DP
    const LauVetoes* vetoes_;

    //! The maximum height of 2D histogram
    Double_t maxHeight_;

    //! Boolean for using the upper half of DP
    Bool_t upperHalf_;
    //! Boolean for using square DP variables
    Bool_t squareDP_;

    ClassDef( Lau2DAbsHistDPPdf, 0 )    // Abstract base class for 2D DP variations based on a histogram
};

#endif
