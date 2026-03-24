
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

/*! \file Lau2DSplineDPPdf.cc
    \brief File containing implementation of Lau2DSplineDPPdf class.
*/

#include "Lau2DSplineDPPdf.hh"

#include "Lau2DCubicSpline.hh"
#include "LauDaughters.hh"
#include "LauKinematics.hh"
#include "LauRandom.hh"

#include "TAxis.h"
#include "TH2.h"
#include "TRandom.h"
#include "TSystem.h"

#include <iostream>

Lau2DSplineDPPdf::Lau2DSplineDPPdf( const TH2* hist,
                                    LauKinematics* kinematics,
                                    const LauVetoes* vetoes,
                                    Bool_t fluctuateBins,
                                    Bool_t useUpperHalfOnly,
                                    Bool_t squareDP ) :
    Lau2DAbsHistDPPdf( kinematics, vetoes, useUpperHalfOnly, squareDP ),
    spline_( 0 )
{
    //We may need to modify the histogram so clone it
    TH2* tempHist( hist ? dynamic_cast<TH2*>( hist->Clone() ) : 0 );

    if ( ! tempHist ) {
        std::cerr << "ERROR in Lau2DSplineDPPdf constructor : the histogram pointer is null."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    if ( fluctuateBins ) {
        this->doBinFluctuation( tempHist );
    }

    spline_ = new Lau2DCubicSpline( *tempHist );

    // Calculate the PDF normalisation.
    this->calcHistNorm();

    // Also obtain the maximum height
    this->calcMaxHeight( tempHist );

    delete tempHist;
}

Lau2DSplineDPPdf::~Lau2DSplineDPPdf()
{
    delete spline_;
    spline_ = 0;
}

Double_t Lau2DSplineDPPdf::interpolateXYNorm( Double_t x, Double_t y ) const
{
    // Get the normalised interpolated value.
    Double_t value = this->interpolateXY( x, y );
    return value / norm_;
}

Double_t Lau2DSplineDPPdf::interpolateXY( Double_t x, Double_t y ) const
{
    // This function returns the interpolated value of the histogram function
    // for the given values of x and y by finding the adjacent bins and extrapolating
    // using weights based on the inverse distance of the point from the adajcent
    // bin centres.
    // Here, x = m13^2, y = m23^2, or m', theta' for square DP co-ordinates

    // If we're only using one half then flip co-ordinates
    // appropriately for conventional or square DP
    getUpperHalf( x, y );

    // First ask whether the point is inside the kinematic region.
    if ( withinDPBoundaries( x, y ) == kFALSE ) {
        std::cerr << "WARNING in Lau2DSplineDPPdf::interpolateXY : Given position is outside the DP boundary, returning 0.0."
                  << std::endl;
        return 0.0;
    }

    return spline_->evaluate( x, y );
}

void Lau2DSplineDPPdf::calcHistNorm()
{
    norm_ = spline_->analyticalIntegral();
}
