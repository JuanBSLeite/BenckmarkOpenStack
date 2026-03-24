
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

/*! \file LauBkgndDPModel.hh
    \brief File containing declaration of LauBkgndDPModel class.
*/

/*! \class LauBkgndDPModel
    \brief Class for defining a histogram-based background Dalitz plot model

    Class for defining a histogram-based background Dalitz plot model
*/

#ifndef LAU_BKGND_DP_MODEL
#define LAU_BKGND_DP_MODEL

#include "LauAbsBkgndDPModel.hh"

#include <vector>

class TH2;
class Lau2DAbsDPPdf;
class LauDaughters;
class LauFitDataTree;
class LauVetoes;

class LauBkgndDPModel : public LauAbsBkgndDPModel {

  public:
    //! Constructor
    /*!
        \param [in] daughters the daughters in the decay
        \param [in] vetoes the vetoes in the Datliz plot
    */
    LauBkgndDPModel( LauDaughters* daughters, LauVetoes* vetoes );

    //! Destructor
    virtual ~LauBkgndDPModel();

    //! Set background histogram
    /*!
        \param [in] histo the 2D histogram describing the DP distribution
        \param [in] useInterpolation boolean flag to determine whether linear interpolation between bins should be used or simply the raw bin values
        \param [in] fluctuateBins boolean flag to determine whether the bin contents should be fluctuated in accordance with their errors.
        The seed for the random number generator used to fluctuate the bins should first be set using LauRandom::setSeed.

        \param [in] useUpperHalfOnly boolean flag to specify that the supplied histogram contains only the upper half of a symmetric DP (or lower half if using square DP coordinates)
        \param [in] squareDP boolean flag to determine whether the supplied histogram is given in square DP coordinates
    */
    void setBkgndHisto( const TH2* histo,
                        Bool_t useInterpolation,
                        Bool_t fluctuateBins,
                        Bool_t useUpperHalfOnly,
                        Bool_t squareDP = kFALSE );

    //! Set the background histogram and generate a spline
    /*!
        \param [in] histo the 2D histogram describing the DP distribution
        \param [in] fluctuateBins boolean flag to determine whether the bin contents should be fluctuated in accordance with their errors.
        The seed for the random number generator used to fluctuate the bins should first be set using LauRandom::setSeed.

        \param [in] useUpperHalfOnly boolean flag to specify that the supplied histogram contains only the upper half of a symmetric DP (or lower half if using square DP coordinates)
        \param [in] squareDP boolean flag to determine whether the supplied histogram is given in square DP coordinates
    */
    void setBkgndSpline( const TH2* histo,
                         Bool_t fluctuateBins,
                         Bool_t useUpperHalfOnly,
                         Bool_t squareDP );

    //! Initialise the model
    virtual void initialise();

    //! Generate a toy MC event from the model
    /*!
        \return success/failure flag
    */
    virtual Bool_t generate();

    //! Get likelihood for a given DP position
    /*!
        \param [in] m13Sq the invariant mass squared of children 1 and 3
        \param [in] m23Sq the invariant mass squared of children 2 and 3
        \return the likelihood value
    */
    virtual Double_t getLikelihood( const Double_t m13Sq, const Double_t m23Sq );

    //! Get unnormalised likelihood for a given DP position
    /*!
        \param [in] m13Sq the invariant mass squared of children 1 and 3
        \param [in] m23Sq the invariant mass squared of children 2 and 3
        \return the unnormalised likelihood value
    */
    virtual Double_t getUnNormValue( const Double_t m13Sq, const Double_t m23Sq );

    //! Get likelihood for a given event
    /*!
        \param [in] iEvt the event number
        \return the likelihood value
    */
    virtual Double_t getLikelihood( UInt_t iEvt );

    //! Get unnormalised likelihood for a given event
    /*!
        \param [in] iEvt the event number
        \return the unnormalised likelihood value
    */
    virtual Double_t getUnNormValue( UInt_t iEvt );

    //! Get PDF normalisation constant
    /*!
        \return the PDF normalisation constant
    */
    virtual Double_t getPdfNorm() const { return pdfNorm_; }

    //! Cache the input data and (if appropriate) the per-event likelihood values
    /*!
        \param [in] fitDataTree the input data
    */
    virtual void fillDataTree( const LauFitDataTree& fitDataTree );

  protected:
    //! Calulate histogram value at a given point
    /*!
        \param [in] xVal the x-value
        \param [in] yVal the y-value
        \return the histogram value
    */
    Double_t calcHistValue( Double_t xVal, Double_t yVal );

    //! Set data event number
    /*!
        \param [in] iEvt the event number
    */
    virtual void setDataEventNo( UInt_t iEvt );

  private:
    //! Copy constructor (not implemented)
    LauBkgndDPModel( const LauBkgndDPModel& rhs );

    //! Copy assignment operator (not implemented)
    LauBkgndDPModel& operator=( const LauBkgndDPModel& rhs );

    //! Flags whether the DP is symmetrical or not
    Bool_t symmetricalDP_;

    //! Flags whether or not to work in square DP coordinates
    Bool_t squareDP_;

    //! PDF of Dalitz plot background, from a 2D histogram
    Lau2DAbsDPPdf* bgHistDPPdf_;

    //! Cached histogram values for each event
    std::vector<Double_t> bgData_;

    //! Histogram value for the current event
    Double_t curEvtHistVal_;

    //! Maximum height of PDF
    Double_t maxPdfHeight_;

    //! Normalisation of PDF
    Double_t pdfNorm_;

    //! Boolean to indicate if the warning that there is no histogram has already been issued
    Bool_t doneGenWarning_;

    //! Flag to track whether a warning has been issued for bin values less than zero
    mutable Bool_t lowBinWarningIssued_;

    ClassDef( LauBkgndDPModel, 0 )    // DP background model
};

#endif
