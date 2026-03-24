
/*
Copyright 2008 University of Warwick

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

/*! \file Lau2DHistPdf.hh
    \brief File containing declaration of Lau2DHistPdf class.
*/

/*! \class Lau2DHistPdf
    \brief Class for defining a 2D histogram PDF.

    Class for defining a 2D histogram PDF.
    Employs linear interpolation to get the histogram value based on how far away a
    point is to nearby bin centres.  The returned values are normalised to
    the total area.
*/

#ifndef LAU_2DHIST_PDF
#define LAU_2DHIST_PDF

#include "LauAbsPdf.hh"

class TH2;
class TH1;
class Lau1DHistPdf;

class Lau2DHistPdf : public LauAbsPdf {

  public:
    //! Constructor
    /*!
        \param [in] theVarNames the names of the abscissa variables
        \param [in] hist the 2D histogram from which the PDF should be constructed
        \param [in] minVals the minimum values of the abscissas
        \param [in] maxVals the maximum values of the abscissas
        \param [in] useInterpolation boolean flag to determine whether linear interpolation between bins should be used or simply the raw bin values
        \param [in] fluctuateBins boolean flag to determine whether the bin contents should be fluctuated in accordance with their errors.
        The seed for the random number generator used to fluctuate the bins should first be set using LauRandom::setSeed.
    */
    Lau2DHistPdf( const std::vector<TString>& theVarNames,
                  const TH2* hist,
                  const LauFitData& minVals,
                  const LauFitData& maxVals,
                  Bool_t useInterpolation = kTRUE,
                  Bool_t fluctuateBins    = kFALSE );

    //! Destructor
    virtual ~Lau2DHistPdf();

    //! Get the minimum value of x-axis abscissa
    Double_t getMinX() const { return minX_; }

    //! Get the maximum value of x-axis abscissa
    Double_t getMaxX() const { return maxX_; }

    //! Get the range of x-axis abscissa values
    Double_t getRangeX() const { return rangeX_; }

    //! Get the minimum value of y-axis abscissa
    Double_t getMinY() const { return minY_; }

    //! Get the maximum value of y-axis abscissa
    Double_t getMaxY() const { return maxY_; }

    //! Get the range of y-axis abscissa values
    Double_t getRangeY() const { return rangeY_; }

    //! Cache information from data
    /*!
        \param [in] inputData the input data
    */
    virtual void cacheInfo( const LauFitDataTree& inputData );

    //! Calculate the likelihood (and intermediate info) for a given value of the abscissas
    /*!
        \param [in] abscissas the abscissa values
    */
    virtual void calcLikelihoodInfo( const LauAbscissas& abscissas );

    //! Calculate the likelihood (and intermediate info) for a given event number
    /*!
        \param [in] iEvt event number
    */
    virtual void calcLikelihoodInfo( UInt_t iEvt );

    //! Get likelihood only for the given variable
    /*!
        \param [in] theVarName the variable name (must be one of the two dimensions of this histogram!)
        \return likelihood value
    */
    virtual Double_t getLikelihood( const TString& theVarName ) const;
    using LauAbsPdf::getLikelihood;

    //! Calculate the normalisation
    virtual void calcNorm();

    //! Calculate the PDF height
    /*!
        \param [in] kinematics the current DP kinematics (not strictly required in this case since PDF has no DP dependence)
    */
    virtual void calcPDFHeight( const LauKinematics* kinematics );

    //! Generate an event from the PDF
    /*!
        \param [in] kinematics the current DP kinematics (used to determine the DP position by some PDFs that have dependence on it)
        \return generated data
    */
    virtual LauFitData generate( const LauKinematics* kinematics );

  protected:
    //! Fluctuate the histogram bin contents in accordance with their errors
    void doBinFluctuation();

    //! Check the normalisation calculation
    void checkNormalisation();

    //! Get the bin content from the histogram
    /*!
        \param [in] i the x-axis bin number
        \param [in] j the y-axis bin number
        \return the bin conent
    */
    Double_t getBinHistValue( Int_t i, Int_t j ) const;

    //! Perform the interpolation (unnormalised)
    /*!
        \param [in] x the x-axis abscissa value
        \param [in] y the y-aixs abscissa value
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

  private:
    //! Copy constructor (not implemented)
    Lau2DHistPdf( const Lau2DHistPdf& other );

    //! Copy assignment operator (not implemented)
    Lau2DHistPdf& operator=( const Lau2DHistPdf& other );

    //! The underlying histogram
    TH2* hist_;
    //! Projection of histogram x-axis
    TH1* xProj_;
    //! Projection of histogram y-axis
    TH1* yProj_;

    //! 1D PDF for x variable
    Lau1DHistPdf* xVarPdf_;
    //! 1D PDF for y variable
    Lau1DHistPdf* yVarPdf_;

    //! x variable name
    TString xName_;
    //! y variable name
    TString yName_;

    //! The number of bins on the x-axis of the histogram
    Int_t nBinsX_;
    //! The number of bins on the y-axis of the histogram
    Int_t nBinsY_;
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

    //! Control boolean for using the linear interpolation
    Bool_t useInterpolation_;
    //! Control boolean for performing the fluctuation of the histogram bin contents
    Bool_t fluctuateBins_;

    ClassDef( Lau2DHistPdf, 0 )    // 2D histogram pdf class
};

#endif
