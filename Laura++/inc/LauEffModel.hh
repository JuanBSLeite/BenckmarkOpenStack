
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

/*! \file LauEffModel.hh
    \brief File containing declaration of LauEffModel class.
*/

/*! \class LauEffModel
    \brief Class that implements the efficiency description across the signal Dalitz plot.

    Class that defines the efficiency model variation across the signal Dalitz plot.
    The method is based in a predetermined two-dimensional histogram to characterize the phase space acceptance.
    The efficiency variation is defined in terms of x = m_13^2, y = m_23^2 for the Dalitz plot (default)
    or x = m', y = theta' for the square Dalitz plot
*/

#ifndef LAUEFFMODEL
#define LAUEFFMODEL

#include "LauAbsEffModel.hh"

class TH2;

class LauDaughters;
class LauKinematics;
class LauVetoes;
class Lau2DAbsDP;

class LauEffModel : public LauAbsEffModel {

  public:
    //! Constructor
    /*!
        \param [in] daughters the daughters particles of the Dalitz plot model
        \param [in] vetoes the object describing the vetoes applied in the phase space
    */
    LauEffModel( const LauDaughters* daughters, const LauVetoes* vetoes );

    //! Destructor
    virtual ~LauEffModel();

    //! Set the efficiency variation across the phase space using a predetermined 2D histogram.
    /*!
        The efficiency is defined in terms of x = m_13^2, y = m_23^2 or x = m', y = theta' for the square Dalitz plot

        \param [in] effHisto the 2-dimensional histogram that describes the efficiency variation
        \param [in] useInterpolation boolean flag decision to switch on/off linear interpolation between bins should be used or simply the raw bin values.
        \param [in] fluctuateBins boolean flag to determine whether the bin contents should be fluctuated in accordance with their errors.
        The seed for the random number generator used to fluctuate the bins should first be set using LauRandom::setSeed.

        \param [in] avEff the desired average efficiency - see Lau2DHistDP::raiseOrLowerBins, values less than zero switch off this behaviour.
        The seed for the random number generator used to raise or lower the bins should first be set using LauRandom::setSeed.

        \param [in] absError the error on that efficiency - see Lau2DHistDP::raiseOrLowerBins, values less than zero switch off this behaviour
        \param [in] useUpperHalfOnly boolean flag to specify that the supplied histogram contains only the upper half of a symmetric DP (or lower half if using square DP coordinates)
        \param [in] squareDP boolean flag to determine whether the supplied histogram is given in square DP coordinates
    */
    void setEffHisto( const TH2* effHisto,
                      Bool_t useInterpolation = kTRUE,
                      Bool_t fluctuateBins    = kFALSE,
                      Double_t avEff          = -1.0,
                      Double_t absError       = -1.0,
                      Bool_t useUpperHalfOnly = kFALSE,
                      Bool_t squareDP         = kFALSE );

    //! Set the efficiency variation across the phase space using a predetermined 2D histogram.
    /*!
        The efficiency is defined in terms of x = m_13^2, y = m_23^2 or x = m', y = theta' for the square Dalitz plot

        \param [in] effHisto the 2-dimensional histogram that describes the efficiency variation
        \param [in] errorHi the 2-dimensional histogram that describes the upper uncertainty on the efficiency variation
        \param [in] errorLo the 2-dimensional histogram that describes the lower uncertainty on the efficiency variation
        \param [in] useInterpolation boolean flag decision to switch on/off linear interpolation between bins should be used or simply the raw bin values.
        \param [in] fluctuateBins boolean flag to determine whether the bin contents should be fluctuated in accordance with their errors.
        The seed for the random number generator used to fluctuate the bins should first be set using LauRandom::setSeed.

        \param [in] avEff the desired average efficiency - see Lau2DHistDP::raiseOrLowerBins, values less than zero switch off this behaviour.
        The seed for the random number generator used to raise or lower the bins should first be set using LauRandom::setSeed.

        \param [in] absError the error on that efficiency - see Lau2DHistDP::raiseOrLowerBins, values less than zero switch off this behaviour
        \param [in] useUpperHalfOnly boolean flag to specify that the supplied histogram contains only the upper half of a symmetric DP (or lower half if using square DP coordinates)
        \param [in] squareDP boolean flag to determine whether the supplied histogram is given in square DP coordinates
    */
    void setEffHisto( const TH2* effHisto,
                      const TH2* errorHi,
                      const TH2* errorLo,
                      Bool_t useInterpolation = kTRUE,
                      Bool_t fluctuateBins    = kFALSE,
                      Double_t avEff          = -1.0,
                      Double_t absError       = -1.0,
                      Bool_t useUpperHalfOnly = kFALSE,
                      Bool_t squareDP         = kFALSE );

    //! Set the efficiency variation across the phase space using a spline based on a predetermined 2D histogram.
    /*!
        The efficiency is defined in terms of x = m_13^2, y = m_23^2 or x = m', y = theta' for the square Dalitz plot

        \param [in] effHisto the 2-dimensional histogram that describes the efficiency variation
        \param [in] fluctuateBins boolean flag to determine whether the bin contents should be fluctuated in accordance with their errors.
        The seed for the random number generator used to fluctuate the bins should first be set using LauRandom::setSeed.

        \param [in] avEff the desired average efficiency - see Lau2DHistDP::raiseOrLowerBins, values less than zero switch off this behaviour.
        The seed for the random number generator used to raise or lower the bins should first be set using LauRandom::setSeed.

        \param [in] absError the error on that efficiency - see Lau2DHistDP::raiseOrLowerBins, values less than zero switch off this behaviour
        \param [in] useUpperHalfOnly boolean flag to specify that the supplied histogram contains only the upper half of a symmetric DP (or lower half if using square DP coordinates)
        \param [in] squareDP boolean flag to determine whether the supplied histogram is given in square DP coordinates
    */
    void setEffSpline( const TH2* effHisto,
                       Bool_t fluctuateBins    = kFALSE,
                       Double_t avEff          = -1.0,
                       Double_t absError       = -1.0,
                       Bool_t useUpperHalfOnly = kFALSE,
                       Bool_t squareDP         = kFALSE );

    //! Set the efficiency variation across the phase space using a spline based on a predetermined 2D histogram.
    /*!
        The efficiency is defined in terms of x = m_13^2, y = m_23^2 or x = m', y = theta' for the square Dalitz plot

        \param [in] effHisto the 2-dimensional histogram that describes the efficiency variation
        \param [in] errorHi the 2-dimensional histogram that describes the upper uncertainty on the efficiency variation
        \param [in] errorLo the 2-dimensional histogram that describes the lower uncertainty on the efficiency variation
        \param [in] fluctuateBins boolean flag to determine whether the bin contents should be fluctuated in accordance with their errors.
        The seed for the random number generator used to fluctuate the bins should first be set using LauRandom::setSeed.

        \param [in] avEff the desired average efficiency - see Lau2DHistDP::raiseOrLowerBins, values less than zero switch off this behaviour.
        The seed for the random number generator used to raise or lower the bins should first be set using LauRandom::setSeed.

        \param [in] absError the error on that efficiency - see Lau2DHistDP::raiseOrLowerBins, values less than zero switch off this behaviour
        \param [in] useUpperHalfOnly boolean flag to specify that the supplied histogram contains only the upper half of a symmetric DP (or lower half if using square DP coordinates)
        \param [in] squareDP boolean flag to determine whether the supplied histogram is given in square DP coordinates
    */
    void setEffSpline( const TH2* effHisto,
                       const TH2* errorHi,
                       const TH2* errorLo,
                       Bool_t fluctuateBins    = kFALSE,
                       Double_t avEff          = -1.0,
                       Double_t absError       = -1.0,
                       Bool_t useUpperHalfOnly = kFALSE,
                       Bool_t squareDP         = kFALSE );

    //! Add a multiplicative efficiency variation across the phase space using a predetermined 2D histogram.
    /*!
        The efficiency is defined in terms of x = m_13^2, y = m_23^2 or x = m', y = theta' for the square Dalitz plot

        \param [in] effHisto the 2-dimensional histogram that describes the efficiency variation
        \param [in] useInterpolation boolean flag decision to switch on/off linear interpolation between bins should be used or simply the raw bin values.
        \param [in] avEff the desired average efficiency - see Lau2DHistDP::raiseOrLowerBins, values less than zero switch off this behaviour.
        The seed for the random number generator used to raise or lower the bins should first be set using LauRandom::setSeed.

        \param [in] absError the error on that efficiency - see Lau2DHistDP::raiseOrLowerBins, values less than zero switch off this behaviour
        \param [in] useUpperHalfOnly boolean flag to specify that the supplied histogram contains only the upper half of a symmetric DP (or lower half if using square DP coordinates)
        \param [in] squareDP boolean flag to determine whether the supplied histogram is given in square DP coordinates
    */
    void addEffHisto( const TH2* effHisto,
                      Bool_t useInterpolation = kTRUE,
                      Double_t avEff          = -1.0,
                      Double_t absError       = -1.0,
                      Bool_t useUpperHalfOnly = kFALSE,
                      Bool_t squareDP         = kFALSE );

    //! Add a multiplicative efficiency variation across the phase space using a predetermined 2D histogram.
    /*!
        The efficiency is defined in terms of x = m_13^2, y = m_23^2 or x = m', y = theta' for the square Dalitz plot

        \param [in] effHisto the 2-dimensional histogram that describes the efficiency variation
        \param [in] errorHi the 2-dimensional histogram that describes the upper uncertainty on the efficiency variation
        \param [in] errorLo the 2-dimensional histogram that describes the lower uncertainty on the efficiency variation
        \param [in] useInterpolation boolean flag decision to switch on/off linear interpolation between bins should be used or simply the raw bin values.
        \param [in] avEff the desired average efficiency - see Lau2DHistDP::raiseOrLowerBins, values less than zero switch off this behaviour.
        The seed for the random number generator used to raise or lower the bins should first be set using LauRandom::setSeed.

        \param [in] absError the error on that efficiency - see Lau2DHistDP::raiseOrLowerBins, values less than zero switch off this behaviour
        \param [in] useUpperHalfOnly boolean flag to specify that the supplied histogram contains only the upper half of a symmetric DP (or lower half if using square DP coordinates)
        \param [in] squareDP boolean flag to determine whether the supplied histogram is given in square DP coordinates
    */
    void addEffHisto( const TH2* effHisto,
                      const TH2* errorHi,
                      const TH2* errorLo,
                      Bool_t useInterpolation = kTRUE,
                      Double_t avEff          = -1.0,
                      Double_t absError       = -1.0,
                      Bool_t useUpperHalfOnly = kFALSE,
                      Bool_t squareDP         = kFALSE );

    //! Add a multiplicative efficiency variation across the phase space using a spline based on a predetermined 2D histogram.
    /*!
        The efficiency is defined in terms of x = m_13^2, y = m_23^2 or x = m', y = theta' for the square Dalitz plot

        \param [in] effHisto the 2-dimensional histogram that describes the efficiency variation
        \param [in] avEff the desired average efficiency - see Lau2DHistDP::raiseOrLowerBins, values less than zero switch off this behaviour.
        The seed for the random number generator used to raise or lower the bins should first be set using LauRandom::setSeed.

        \param [in] absError the error on that efficiency - see Lau2DHistDP::raiseOrLowerBins, values less than zero switch off this behaviour
        \param [in] useUpperHalfOnly boolean flag to specify that the supplied histogram contains only the upper half of a symmetric DP (or lower half if using square DP coordinates)
        \param [in] squareDP boolean flag to determine whether the supplied histogram is given in square DP coordinates
    */
    void addEffSpline( const TH2* effHisto,
                       Double_t avEff          = -1.0,
                       Double_t absError       = -1.0,
                       Bool_t useUpperHalfOnly = kFALSE,
                       Bool_t squareDP         = kFALSE );

    //! Add a multiplicative efficiency variation across the phase space using a spline based on a predetermined 2D histogram.
    /*!
        The efficiency is defined in terms of x = m_13^2, y = m_23^2 or x = m', y = theta' for the square Dalitz plot

        \param [in] effHisto the 2-dimensional histogram that describes the efficiency variation
        \param [in] errorHi the 2-dimensional histogram that describes the upper uncertainty on the efficiency variation
        \param [in] errorLo the 2-dimensional histogram that describes the lower uncertainty on the efficiency variation
        \param [in] avEff the desired average efficiency - see Lau2DHistDP::raiseOrLowerBins, values less than zero switch off this behaviour.
        The seed for the random number generator used to raise or lower the bins should first be set using LauRandom::setSeed.

        \param [in] absError the error on that efficiency - see Lau2DHistDP::raiseOrLowerBins, values less than zero switch off this behaviour
        \param [in] useUpperHalfOnly boolean flag to specify that the supplied histogram contains only the upper half of a symmetric DP (or lower half if using square DP coordinates)
        \param [in] squareDP boolean flag to determine whether the supplied histogram is given in square DP coordinates
    */
    void addEffSpline( const TH2* effHisto,
                       const TH2* errorHi,
                       const TH2* errorLo,
                       Double_t avEff          = -1.0,
                       Double_t absError       = -1.0,
                       Bool_t useUpperHalfOnly = kFALSE,
                       Bool_t squareDP         = kFALSE );

    //! Determine the efficiency for a given point in the Dalitz plot.
    /*!
        The method uses the 2D histogram set by the setEffHisto() function and the vetoes information.

        \param [in] kinematics the object that defines the DP position
        \return the efficiency value at the given point in the DP
    */
    Double_t calcEfficiency( const LauKinematics* kinematics ) const;

    //! Determine whether the given DP position is outside the vetoes
    /*!
        \param [in] kinematics the object that defines the DP position
        \return kTRUE if the DP position is outside all veto regions, kFALSE otherwise
    */
    Bool_t passVeto( const LauKinematics* kinematics ) const;

    //! Determine whether the efficiency histogram has had its bins fluctuated within their errors
    Bool_t fluctuateEffHisto() const { return fluctuateEffHisto_; }

    //! Return the daughters object
    /*
            \return the LauDaughters object associated with the DP
        */
    const LauDaughters* getDaughters() const { return daughters_; }

  private:
    //! Copy constructor - not implemented
    LauEffModel( const LauEffModel& rhs );

    //! Copy assignment operator - not implemented
    LauEffModel& operator=( const LauEffModel& rhs );

    //! Get the efficiency from a two-dimensional histogram
    /*!
        \param [in] kinematics the object that defines the DP position
    */
    Double_t getEffHistValue( const LauKinematics* kinematics ) const;

    //! The daughters object
    const LauDaughters* daughters_;

    //! The vetoes object
    const LauVetoes* vetoes_;

    //! The efficiency histogram objects
    std::vector<Lau2DAbsDP*> effHisto_;

    //! Fluctuate histogram within the error
    Bool_t fluctuateEffHisto_;

    //! Flag to track whether a warning has been issued for bin values less than zero
    mutable Bool_t lowBinWarningIssued_;

    //! Flag to track whether a warning has been issued for bin values greater than one
    mutable Bool_t highBinWarningIssued_;

    ClassDef( LauEffModel, 0 )    // Implement the signal efficiency across the DP
};

#endif
