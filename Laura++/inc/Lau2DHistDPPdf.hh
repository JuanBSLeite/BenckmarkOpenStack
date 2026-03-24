
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

/*! \file Lau2DHistDPPdf.hh
    \brief File containing declaration of Lau2DHistDPPdf class.
*/

/*! \class Lau2DHistDPPdf
    \brief Class for defining a 2D DP histogram PDF.

    Class for defining a 2D DP histogram PDF.
    Employs linear interpolation to get the histogram value based on how far away a point in (x,y)
    is to nearby bin centres. The returned values are normalised to the total area.
    The histogram can be defined in the conventional DP (m13Sq vs m23Sq) or in the square DP.
*/

#ifndef LAU_2DHIST_DP_PDF
#define LAU_2DHIST_DP_PDF

#include "Lau2DAbsHistDPPdf.hh"

class TH2;
class LauKinematics;
class LauVetoes;

class Lau2DHistDPPdf : public Lau2DAbsHistDPPdf {

  public:
    //! Constructor
    /*!
        \param [in] hist the 2D DP histogram
        \param [in] kinematics the current DP kinematics
        \param [in] vetoes the vetoes within the DP
        \param [in] useInterpolation boolean flag to determine whether linear interpolation between bins should be used or simply the raw bin values
        \param [in] fluctuateBins boolean flag to determine whether the bin contents should be fluctuated in accordance with their errors (useful for systematic error evaluation).
        The seed for the random number generator used to fluctuate the bins should first be set using LauRandom::setSeed.

        \param [in] useUpperHalfOnly boolean flag to specify that the supplied histogram contains only the upper half of a symmetric DP (or lower half if using square DP coordinates)
        \param [in] squareDP boolean flag to specify whether the supplied histogram is in square DP coordinates
    */
    Lau2DHistDPPdf( const TH2* hist,
                    LauKinematics* kinematics,
                    const LauVetoes* vetoes,
                    Bool_t useInterpolation = kTRUE,
                    Bool_t fluctuateBins    = kFALSE,
                    Bool_t useUpperHalfOnly = kFALSE,
                    Bool_t squareDP         = kFALSE );

    //! Destructor
    virtual ~Lau2DHistDPPdf();

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
    //! Get the bin content from the histogram
    /*!
        \param [in] xBinNo the x-axis bin number
        \param [in] yBinNo the y-axis bin number
        \return the bin conent
    */
    Double_t getBinHistValue( Int_t xBinNo, Int_t yBinNo ) const;

    //! Calculate the PDF normalisation
    void calcHistNorm();

    //! Check the normalisation calculation
    void checkNormalisation();

  private:
    //! Copy constructor - not implemented
    Lau2DHistDPPdf( const Lau2DHistDPPdf& rhs );

    //! Copy assignment operator - not implemented
    Lau2DHistDPPdf& operator=( const Lau2DHistDPPdf& rhs );

    //! The underlying histogram
    TH2* hist_;

    //! The histogram x-axis minimum
    Double_t minX_;
    //! The histogram x-axis maximum
    Double_t maxX_;
    //! The histogram y-axis minimum
    Double_t minY_;
    //! The histogram y-axis maximum
    Double_t maxY_;
    //! The histogram x-axis range
    Double_t rangeX_;
    //! The histogram y-axis range
    Double_t rangeY_;

    //! The histogram x-axis bin width
    Double_t binXWidth_;
    //! The histogram y-axis bin width
    Double_t binYWidth_;
    //! The histogram x-axis inverse bin width
    Double_t invBinXWidth_;
    //! The histogram y-axis inverse bin width
    Double_t invBinYWidth_;

    //! The maximum height of 2D histogram
    Double_t maxHeight_;

    //! The number of bins on the x-axis of the histogram
    Int_t nBinsX_;
    //! The number of bins on the y-axis of the histogram
    Int_t nBinsY_;

    //! The histogram normalisation
    Double_t norm_;

    //! Control boolean for using the linear interpolation
    Bool_t useInterpolation_;

    ClassDef( Lau2DHistDPPdf, 0 )    // 2D Histogram utility class for DP analyses
};

#endif
