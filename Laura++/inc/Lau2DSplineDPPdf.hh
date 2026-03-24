
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

/*! \file Lau2DSplineDPPdf.hh
    \brief File containing declaration of Lau2DSplineDPPdf class.
*/

/*! \class Lau2DSplineDPPdf
    \brief Class for defining a 2D DP spline PDF.

    Class for defining a 2D DP spline PDF.
    Employs a 2D cubic spline to get the histogram value based on an input histogram.
    The returned values are normalised to the total area.
    The histogram can be defined in the conventional DP (m13Sq vs m23Sq) or in the square DP.
*/

#ifndef LAU_2DSPLINE_DP_PDF
#define LAU_2DSPLINE_DP_PDF

#include "Lau2DAbsHistDPPdf.hh"

class TH2;
class Lau2DCubicSpline;
class LauKinematics;
class LauVetoes;

class Lau2DSplineDPPdf : public Lau2DAbsHistDPPdf {

  public:
    //! Constructor
    /*!
        \param [in] hist the 2D DP histogram
        \param [in] kinematics the current DP kinematics
        \param [in] vetoes the vetoes within the DP
        \param [in] fluctuateBins boolean flag to determine whether the bin contents should be fluctuated in accordance with their errors (useful for systematic error evaluation).
        The seed for the random number generator used to fluctuate the bins should first be set using LauRandom::setSeed.

        \param [in] useUpperHalfOnly boolean flag to specify that the supplied histogram contains only the upper half of a symmetric DP (or lower half if using square DP coordinates)
        \param [in] squareDP boolean flag to specify whether the supplied histogram is in square DP coordinates
    */
    Lau2DSplineDPPdf( const TH2* hist,
                      LauKinematics* kinematics,
                      const LauVetoes* vetoes,
                      Bool_t fluctuateBins    = kFALSE,
                      Bool_t useUpperHalfOnly = kFALSE,
                      Bool_t squareDP         = kFALSE );

    //! Destructor
    virtual ~Lau2DSplineDPPdf();

    //! Retrieve PDF normalisation
    /*!
        \return the normalisation factor
    */
    Double_t getHistNorm() const { return norm_; }

    //! Perform the interpolation (unnormalised)
    /*!
        \param [in] x the x-axis value
        \param [in] y the y-aixs value
        \return the unnormalised PDF value
    */
    Double_t interpolateXY( Double_t x, Double_t y ) const;

    //! Perform the interpolation and divide by the normalisation
    /*!
        \param [in] x the x-axis abscissa value
        \param [in] y the y-axis abscissa value
        \return the normalised PDF value
    */
    Double_t interpolateXYNorm( Double_t x, Double_t y ) const;

  protected:
    //! Calculate the PDF normalisation
    void calcHistNorm();

  private:
    //! Copy constructor - not implemented
    Lau2DSplineDPPdf( const Lau2DSplineDPPdf& rhs );

    //! Copy assignment operator - not implemented
    Lau2DSplineDPPdf& operator=( const Lau2DSplineDPPdf& rhs );

    //! The maximum height of 2D histogram
    Double_t maxHeight_;

    //! The histogram normalisation
    Double_t norm_;

    //! A 2D cubic spline generated from the histogram
    Lau2DCubicSpline* spline_;

    ClassDef( Lau2DSplineDPPdf, 0 )    // 2D Spline utility class for DP analyses
};

#endif
