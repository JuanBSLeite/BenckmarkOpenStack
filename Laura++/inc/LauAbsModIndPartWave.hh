
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

/*! \file LauAbsModIndPartWave.hh
    \brief File containing declaration of LauAbsModIndPartWave class.
*/

/*! \class LauAbsModIndPartWave
    \brief Abstract base class for defining a model independent partial wave component

    Abstract base class for defining a model independent partial wave component.
    This model uses splines to produce a partial wave from two sets of real numbers that represent the amplitude at a series of points in the phase space.
    These real numbers at each point can be floated in the fit.
    Classes inheriting from this define whether these real numbers are e.g. the magnitude and phase or e.g. the real and imaginary part of the amplitude.
*/

#ifndef LAU_ABSMODINDPARTWAVE
#define LAU_ABSMODINDPARTWAVE

#include "Lau1DCubicSpline.hh"
#include "LauAbsResonance.hh"
#include "LauComplex.hh"

#include <set>

class LauAbsModIndPartWave : public LauAbsResonance {

  public:
    //! Constructor
    /*!
        \param [in] resInfo the object containing information on the resonance name, mass, width, spin, charge, etc.
        \param [in] resPairAmpInt the number of the daughter not produced by the resonance
        \param [in] daughters the daughter particles
    */
    LauAbsModIndPartWave( LauResonanceInfo* resInfo,
                          Int_t resPairAmpInt,
                          const LauDaughters* daughters );

    //! Destructor
    virtual ~LauAbsModIndPartWave();

    //! Initialise the model
    virtual void initialise();

    //! Define the knot positions
    /*!
        If absent from the set provided, knots are added automatically at the upper and lower kinematic limits

        \param [in] masses the mass values at which the knots are placed
    */
    void defineKnots( const std::set<Double_t>& masses );

    //! Return the number of knots that have been defined (including those at the upper and lower kinematic limits)
    /*!
        \return the number of knots
    */
    UInt_t nKnots() const { return nKnots_; }

    //! Set the values of the two real parameters that define the amplitude at a given knot
    /*!
        \param [in] knot the knot to be updated
        \param [in] ampVal1 the value of first real parameter representing the amplitude at the knot
        \param [in] ampVal2 the value of second real parameter representing the amplitude at the knot
        \param [in] fixAmpVal1 whether the first real parameter should be fixed
        \param [in] fixAmpVal2 whether the second real parameter should be fixed
    */
    virtual void setKnotAmp( const UInt_t knot,
                             const Double_t ampVal1,
                             const Double_t ampVal2,
                             const Bool_t fixAmpVal1,
                             const Bool_t fixAmpVal2 ) = 0;

    //! Set whether the parameters should be floated only in the second-stage of a two-stage fit
    /*!
        By default, the parameters describing the amplitude at each
        knot will float from the outset of a fit.
        If, however, a good estimate of these is already known, it
        can be more efficient to initially fix them and then to
        float them only in a second stage (as is done for other
        resonance lineshape parameters).
        This function allows the toggling of this behaviour.

        \param secondStage whether the parameters should float only in the second stage
    */
    void floatKnotsSecondStage( const Bool_t secondStage );

    //! Retrieve the value of the second stage flag
    Bool_t floatKnotsSecondStage() const { return secondStage_; };

    //! Retrieve the resonance parameters, e.g. so that they can be loaded into a fit
    /*!
        \return floating parameters of the resonance
    */
    virtual const std::vector<LauParameter*>& getFloatingParameters();

  protected:
    //! Complex resonant amplitude
    /*!
        \param [in] mass appropriate invariant mass for the resonance
        \param [in] spinTerm spin term
    */
    virtual LauComplex resAmp( Double_t mass, Double_t spinTerm );

    //! Evaluate the amplitude at the given point from the splines
    /*!
        \param [in] mass appropriate invariant mass for the resonance
    */
    virtual void evaluateAmplitude( const Double_t mass ) = 0;

    //! Method to check that the supplied knot positions are valid
    /*!
        \param [in] masses the mass values at which the knots are placed
    */
    virtual std::set<Double_t> checkKnots( const std::set<Double_t>& masses );

    //! Method to create the parameter objects for the given knot
    /*!
        \param [in] iKnot the index of the knot
    */
    virtual void createAmpParameters( const UInt_t iKnot ) = 0;

    //! Method to set the type of interpolation used for the splines
    /*!
        \param [in] type1 the type of interpolation for the first spline
        \param [in] type2 the type of interpolation for the second spline
    */
    void setSplineType( Lau1DCubicSpline::LauSplineType type1, Lau1DCubicSpline::LauSplineType type2 );

    //! Method to set the boundary conditions of the splines
    /*!
        \param [in] leftBound1 the type of boundary condition for the left edge of the first spline
        \param [in] rightBound1 the type of boundary condition for the right edge of the first spline
        \param [in] leftBound2 the type of boundary condition for the left edge of the second spline
        \param [in] rightBound2 the type of boundary condition for the right edge of the second spline
        \param [in] leftGrad1 the gradient at the left edge of the first spline if clamped
        \param [in] rightGrad1 the gradient at the right edge of the first spline if clamped
        \param [in] leftGrad2 the gradient at the left edge of the second spline if clamped
        \param [in] rightGrad2 the gradient at the right edge of the second spline if clamped
    */
    void setSplineBoundaryConditions( Lau1DCubicSpline::LauSplineBoundaryType leftBound1,
                                      Lau1DCubicSpline::LauSplineBoundaryType rightBound1,
                                      Lau1DCubicSpline::LauSplineBoundaryType leftBound2,
                                      Lau1DCubicSpline::LauSplineBoundaryType rightBound2,
                                      Double_t leftGrad1  = 0.0,
                                      Double_t rightGrad1 = 0.0,
                                      Double_t leftGrad2  = 0.0,
                                      Double_t rightGrad2 = 0.0 );

    //! Helper function to set the current amplitude value
    /*!
        \param [in] realPart the real part of the amplitude
        \param [in] imagPart the imaginary part of the amplitude
    */
    void setAmp( const Double_t realPart, const Double_t imagPart )
    {
        amp_.setRealImagPart( realPart, imagPart );
    }

    //! Helper function to access the masses
    const std::vector<Double_t>& getMasses() { return masses_; }

    //! Helper function to access the 1st parameter set
    std::vector<Double_t>& getAmp1Vals() { return amp1Vals_; }

    //! Helper function to access the 2nd parameter set
    std::vector<Double_t>& getAmp2Vals() { return amp2Vals_; }

    //! Helper function to access the 1st parameter set
    std::vector<LauParameter*>& getAmp1Pars() { return amp1Pars_; }

    //! Helper function to access the 2nd parameter set
    std::vector<LauParameter*>& getAmp2Pars() { return amp2Pars_; }

    //! Helper function to access the 1st spline
    const Lau1DCubicSpline* getSpline1() const { return spline1_; }

    //! Helper function to access the 1st spline
    const Lau1DCubicSpline* getSpline2() const { return spline2_; }

  private:
    //! The number of knots
    UInt_t nKnots_;

    //! The masses at which knots are defined in the magnitude and phase splines
    std::vector<Double_t> masses_;

    //! The values of the first real parameter at each knot
    std::vector<Double_t> amp1Vals_;
    //! The values of the second real parameter at each knot
    std::vector<Double_t> amp2Vals_;

    //! The parameters for the first real value at the knots
    std::vector<LauParameter*> amp1Pars_;
    //! The parameters for the second real value at the knots
    std::vector<LauParameter*> amp2Pars_;

    //! The spline used to interpolate the values of the first real parameter
    Lau1DCubicSpline* spline1_;
    //! The spline used to interpolate the values of the second real parameter
    Lau1DCubicSpline* spline2_;

    //! The type of interpolation used for the first spline
    Lau1DCubicSpline::LauSplineType type1_;
    //! The type of interpolation used for the second spline
    Lau1DCubicSpline::LauSplineType type2_;

    //! The lower boundary condition type for the first spline
    Lau1DCubicSpline::LauSplineBoundaryType leftBound1_;
    //! The upper boundary condition type for the first spline
    Lau1DCubicSpline::LauSplineBoundaryType rightBound1_;
    //! The lower boundary condition type for the second spline
    Lau1DCubicSpline::LauSplineBoundaryType leftBound2_;
    //! The upper boundary condition type for the second spline
    Lau1DCubicSpline::LauSplineBoundaryType rightBound2_;

    //! The gradient at the left boundary for the first spline if clamped
    Double_t leftGrad1_;
    //! The gradient at the right boundary for the first spline if clamped
    Double_t rightGrad1_;
    //! The gradient at the left boundary for the second spline if clamped
    Double_t leftGrad2_;
    //! The gradient at the right boundary for the second spline if clamped
    Double_t rightGrad2_;

    //! Flag to determine if the parameters should be floated only in the second stage of the fit
    Bool_t secondStage_;

    //! The current value of the amplitude
    LauComplex amp_;

    ClassDef( LauAbsModIndPartWave, 0 )    // model independent partial wave
};

#endif
