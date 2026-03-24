
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

/*! \file LauModIndPartWaveMagPhase.hh
    \brief File containing declaration of LauModIndPartWaveMagPhase class.
*/

/*! \class LauModIndPartWaveMagPhase
    \brief Class for defining a model independent partial wave component where the amplitudes are parameterised in terms of magnitude and phase

    This model uses splines to produce a partial wave from the magnitude and phase values of the amplitude at a series of points in the phase space.
    The magnitude and phase at each point can be floated in the fit.
*/

#ifndef LAU_MODINDPARTWAVE_MAGPHASE
#define LAU_MODINDPARTWAVE_MAGPHASE

#include "LauAbsModIndPartWave.hh"
#include "LauComplex.hh"

class TSpline3;

class LauModIndPartWaveMagPhase : public LauAbsModIndPartWave {

  public:
    //! Constructor
    /*!
        \param [in] resInfo the object containing information on the resonance name, mass, width, spin, charge, etc.
        \param [in] resPairAmpInt the number of the daughter not produced by the resonance
        \param [in] daughters the daughter particles
    */
    LauModIndPartWaveMagPhase( LauResonanceInfo* resInfo,
                               Int_t resPairAmpInt,
                               const LauDaughters* daughters );

    //! Destructor
    virtual ~LauModIndPartWaveMagPhase();

    //! Set the values of the two real parameters that define the amplitude at a given knot
    /*!
        \param [in] knot the knot to be updated
        \param [in] magVal the value of the magnitude at the knot
        \param [in] phaseVal the value of the phase at the knot
        \param [in] fixMagnitude whether the magnitude should be fixed
        \param [in] fixPhase whether the phase should be fixed
    */
    virtual void setKnotAmp( const UInt_t knot,
                             const Double_t magVal,
                             const Double_t phaseVal,
                             const Bool_t fixMagnitude,
                             const Bool_t fixPhase );

    //! Method to set the type of interpolation used for the splines
    /*!
        \param [in] magType the type of interpolation for the magnitude spline
        \param [in] phaseType the type of interpolation for the phase spline
    */
    void setType( Lau1DCubicSpline::LauSplineType magType, Lau1DCubicSpline::LauSplineType phaseType )
    {
        this->setSplineType( magType, phaseType );
    }

    //! Method to set the boundary conditions of the splines
    /*!
        \param [in] magLeftBound the type of boundary condition for the left edge of the magnitude spline
        \param [in] magRightBound the type of boundary condition for the right edge of the magnitude spline
        \param [in] phaseLeftBound the type of boundary condition for the left edge of the phase spline
        \param [in] phaseRightBound the type of boundary condition for the right edge of the phase spline
        \param [in] magLeftGrad the gradient at the left edge of the magnitude spline if clamped
        \param [in] magRightGrad the gradient at the right edge of the magnitude spline if clamped
        \param [in] phaseLeftGrad the gradient at the left edge of the phase spline if clamped
        \param [in] phaseRightGrad the gradient at the right edge of the phase spline if clamped
    */
    void setBoundaryConditions( Lau1DCubicSpline::LauSplineBoundaryType magLeftBound,
                                Lau1DCubicSpline::LauSplineBoundaryType magRightBound,
                                Lau1DCubicSpline::LauSplineBoundaryType phaseLeftBound,
                                Lau1DCubicSpline::LauSplineBoundaryType phaseRightBound,
                                Double_t magLeftGrad    = 0.0,
                                Double_t magRightGrad   = 0.0,
                                Double_t phaseLeftGrad  = 0.0,
                                Double_t phaseRightGrad = 0.0 )
    {
        this->setSplineBoundaryConditions( magLeftBound,
                                           magRightBound,
                                           phaseLeftBound,
                                           phaseRightBound,
                                           magLeftGrad,
                                           magRightGrad,
                                           phaseLeftGrad,
                                           phaseRightGrad );
    }

    //! Get the resonance model type
    /*!
        \return the resonance model type
    */
    virtual LauAbsResonance::LauResonanceModel getResonanceModel() const
    {
        return LauAbsResonance::MIPW_MagPhase;
    }

  protected:
    //! Evaluate the amplitude at the given point from the splines
    /*!
        \param [in] mass appropriate invariant mass for the resonance
    */
    virtual void evaluateAmplitude( const Double_t mass );

    //! Method to create the parameter objects for the given knot
    /*!
        \param [in] iKnot the index of the knot
    */
    virtual void createAmpParameters( const UInt_t iKnot );

  private:
    ClassDef( LauModIndPartWaveMagPhase, 0 )    // model independent partial wave
};

#endif
