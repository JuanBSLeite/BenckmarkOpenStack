
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

/*! \file LauAbsPdf.hh
    \brief File containing declaration of LauAbsPdf class.
*/

/*! \class LauAbsPdf
    \brief Class for defining the abstract interface for PDF classes.

    Class for defining the abstract interface for PDF classes.
    Much common code is implemented but the following three methods are not:
    evaluatePDF, calcNorm, calcPDFHeight
*/

#ifndef LAU_ABS_PDF
#define LAU_ABS_PDF

#include "LauAbsRValue.hh"
#include "LauFitDataTree.hh"
#include "LauParamFixed.hh"
#include "LauParameter.hh"

#include "Rtypes.h"
#include "TRandom.h"
#include "TString.h"

#include <map>
#include <vector>

class LauKinematics;

class LauAbsPdf {

  public:
    //! The type used for containing multiple abscissa values.
    typedef std::vector<Double_t> LauAbscissas;

    //! The possible numerical intergration methods.
    enum IntMethod {
        GaussLegendre, /*!< Gauss-Legendre integration */
        Trapezoid      /*!< Simple trapezoid integration */
    };

    //! Constructor for a 1D PDF
    /*!
        \param [in] theVarName the name of the abscissa variable
        \param [in] params the parameters of the PDF
        \param [in] minAbscissa the minimum value of the abscissa
        \param [in] maxAbscissa the maximum value of the abscissa
    */
    LauAbsPdf( const TString& theVarName,
               const std::vector<LauAbsRValue*>& params,
               Double_t minAbscissa,
               Double_t maxAbscissa );

    //! Constructor for a multidimensional PDF
    /*!
        \param [in] theVarNames the names of the abscissa variables
        \param [in] params the parameters of the PDF
        \param [in] minAbscissas the minimum values of the abscissas
        \param [in] maxAbscissas the maximum values of the abscissas
    */
    LauAbsPdf( const std::vector<TString>& theVarNames,
               const std::vector<LauAbsRValue*>& params,
               const LauFitData& minAbscissas,
               const LauFitData& maxAbscissas );

    //! Destructor
    virtual ~LauAbsPdf() {}

    //! Retrieve the name of the abscissa
    /*!
        \return the abscissa's name
    */
    virtual const TString& varName() const { return varNames_.find( 0 )->second; }

    //! Retrieve the names of the abscissas
    /*!
        \return the abscissas' names
    */
    virtual std::vector<TString> varNames() const;

    //! Retrieve the number of PDF parameters
    /*!
        \return the number of PDF parameters
    */
    virtual UInt_t nParameters() const { return param_.size(); }

    //! Retrieve the number of fixed PDF parameters
    /*!
        \return the number of fixed PDF parameters
    */
    virtual UInt_t nFixedParameters() const;

    //! Retrieve the number of abscissas
    /*!
        \return the number of abscissas
    */
    virtual UInt_t nInputVars() const { return varNames_.size(); }

    //! Specifies whether or not the PDF is DP dependent.
    /*!
        DP-dependent PDFs should override this method.

        \return true if the PDF is DP-dependent, false otherwise (the default)
    */
    virtual Bool_t isDPDependent() const { return kFALSE; }

    //! Retrieve the minimum value of the (primary) abscissa
    /*!
        \return the minimum value of the (primary) abscissa
    */
    virtual Double_t getMinAbscissa() const { return minAbscissas_.front(); }

    //! Retrieve the maximum value of the (primary) abscissa
    /*!
        \return the maximum value of the (primary) abscissa
    */
    virtual Double_t getMaxAbscissa() const { return maxAbscissas_.front(); }

    //! Retrieve the range of the (primary) abscissa
    /*!
        \return the range of the (primary) abscissa
    */
    virtual Double_t getRange() const { return this->getMaxAbscissa() - this->getMinAbscissa(); }

    //! Retrieve the minimum value of the specified abscissa
    /*!
        \param [in] theVarName the name of the abscissa
        \return the minimum value of the specified abscissa
    */
    virtual Double_t getMinAbscissa( const TString& theVarName ) const;

    //! Retrieve the maximum value of the specified abscissa
    /*!
        \param [in] theVarName the name of the abscissa
        \return the maximum value of the specified abscissa
    */
    virtual Double_t getMaxAbscissa( const TString& theVarName ) const;

    //! Retrieve the range of the specified abscissa
    /*!
        \param [in] theVarName the name of the abscissa
        \return the range of the specified abscissa
    */
    virtual Double_t getRange( const TString& theVarName ) const
    {
        return this->getMaxAbscissa( theVarName ) - this->getMinAbscissa( theVarName );
    }

    //! Retrieve the minimum values of all the abscissas
    /*!
        \return the minimum values of the abscissas
    */
    virtual LauFitData getMinAbscissas() const;

    //! Retrieve the maximum values of all the abscissas
    /*!
        \return the maximum values of the abscissas
    */
    virtual LauFitData getMaxAbscissas() const;

    //! Retrieve the ranges of all the abscissas
    /*!
        \return the ranges of the abscissas
    */
    virtual LauFitData getRanges() const;

    //! Update the pulls for all parameters
    virtual void updatePulls();

    //! Cache information from data
    /*!
        Will cache, for every event, the abscissa values and, if all parameters are fixed, the PDF value.

        \param [in] inputData the data to be used to calculate everything
    */
    virtual void cacheInfo( const LauFitDataTree& inputData );

    //! Calculate the likelihood (and all associated information) given value(s) of the abscissa(s)
    /*!
        \param [in] abscissas the value(s) of the abscissa(s)
    */
    virtual void calcLikelihoodInfo( const LauAbscissas& abscissas ) = 0;

    //! Retrieve the likelihood (and all associated information) given the event number
    /*!
        \param [in] iEvt the event number
    */
    virtual void calcLikelihoodInfo( UInt_t iEvt );

    //! Retrieve the unnormalised likelihood value
    /*!
        \return the unnormalised likelihood value
    */
    virtual Double_t getUnNormLikelihood() const { return unNormPDFVal_; }

    //! Retrieve the normalisation factor
    /*!
        \return the normalisation factor
    */
    virtual Double_t getNorm() const { return norm_; }

    //! Retrieve the normalised likelihood value
    /*!
        \return the normalised likelihood value
    */
    virtual Double_t getLikelihood() const;

    //! For multidimentional PDFs, retrieve the normalised likelihood value of a named variable
    /*!
        Defaults to simply returning the complete likelihood.  Should be overridden by multidimensional PDFs.

        \return the normalised likelihood value of the named variable
    */
    virtual Double_t getLikelihood( const TString& theVarName ) const;

    //! Retrieve the maximum height
    /*!
        \return the maximum height
    */
    virtual Double_t getMaxHeight() const { return maxHeight_; }

    //! Generate an event from the PDF
    /*!
        \param [in] kinematics used by some PDFs to determine the DP position, on which they have dependence
    */
    virtual LauFitData generate( const LauKinematics* kinematics );

    //! Set the random function used for toy MC generation
    /*!
        \param [in] randomFun the random function to be used
    */
    virtual void setRandomFun( TRandom* randomFun ) { randomFun_ = randomFun; }

    //! Retrieve the parameters of the PDF, e.g. so that they can be loaded into a fit
    /*!
        \return the parameters of the PDF
    */
    virtual const std::vector<LauAbsRValue*>& getParameters() const { return param_; }

    //! Retrieve the parameters of the PDF, e.g. so that they can be loaded into a fit
    /*!
        \return the parameters of the PDF
    */
    virtual std::vector<LauAbsRValue*>& getParameters() { return param_; }

    //! Calculate the normalisation factor of the PDF
    /*!
        By default uses numerical integration.
        PDFs with analytical integrals should override this method.
    */
    virtual void calcNorm();

    //! Calculate the maximum height of the PDF
    /*!
        \param [in] kinematics used by some PDFs to determine the DP position, on which they have dependence
    */
    virtual void calcPDFHeight( const LauKinematics* kinematics ) = 0;

    //! Check if the maximum height of the PDF is up to date
    /*!
        \return true if the maximum height is up to date
    */
    virtual Bool_t heightUpToDate() const { return heightUpToDate_; }

    //! Set whether the height is up to date
    /*!
        \param [in] hutd whether the height is up to date
    */
    virtual void heightUpToDate( Bool_t hutd ) { heightUpToDate_ = hutd; }

    //! Check if the PDF is to be cached
    /*!
        \return true if the PDF is to be cached
    */
    virtual Bool_t cachePDF() const { return cachePDF_; }

    //! Retrieve the number of points to integrate over when normalising
    /*!
        \return the number of points to integrate over when normalising
    */
    virtual Int_t nNormPoints() const { return nNormPoints_; }

    //! Set the number of points to integrate over when normalising
    /*!
        \param [in] nPoints the number of points
    */
    virtual void nNormPoints( Int_t nPoints ) { nNormPoints_ = nPoints; }

    //! Retrieve the integration method used to normalise the PDF
    /*!
        \return the integration method used to normalise the PDF
    */
    virtual IntMethod integMethod() const { return integMethod_; }

    //! Set the integration method used to normalise the PDF
    /*!
        \param [in] method the integration method to be used
    */
    virtual void integMethod( IntMethod method ) { integMethod_ = method; }

  protected:
    //! Set whether the PDF is to be cached
    /*!
        \param [in] doCachePDF whether the PDF is to be cached
    */
    virtual void cachePDF( Bool_t doCachePDF ) { cachePDF_ = doCachePDF; }

    //! Integrate the PDF using the Gauss-Legendre method
    /*!
        \return the integral of the PDF
    */
    virtual Double_t integrGaussLegendre();

    //! Integrate the PDF using the simple trapezoid method
    /*!
        \return the integral of the PDF
    */
    virtual Double_t integTrapezoid();

    //! Set the normalisation factor
    /*!
        \param [in] norm the normalisation factor
    */
    virtual void setNorm( Double_t norm ) { norm_ = norm; }

    //! Set the maximum height
    /*!
        \param [in] maxHeight the new maximum height
    */
    virtual void setMaxHeight( Double_t maxHeight )
    {
        maxHeight_ = maxHeight;
        this->heightUpToDate( kTRUE );
    }

    //! Set the minimum value of the specified abscissa
    /*!
        \param [in] theVarName the name of the abscissa
        \param [in] minAbscissa value to set as minimum
    */
    virtual void setMinAbscissa( const TString& theVarName, Double_t minAbscissa );

    //! Set the maximum value of the specified abscissa
    /*!
        \param [in] theVarName the name of the abscissa
        \param [in] maxAbscissa value to set as maximum
    */
    virtual void setMaxAbscissa( const TString& theVarName, Double_t maxAbscissa );

    //! Set the range of the specified abscissa
    /*!
        \param [in] theVarName the name of the abscissa
        \param [in] minAbscissa value to set as minimum
        \param [in] maxAbscissa value to set as maximum
    */
    virtual void setRange( const TString& theVarName, Double_t minAbscissa, Double_t maxAbscissa )
    {
        this->setMinAbscissa( theVarName, minAbscissa );
        this->setMaxAbscissa( theVarName, maxAbscissa );
    }

    //! Check that all abscissas are within their allowed ranges
    /*!
        \param [in] abscissas the value(s) of the abscissa(s)
        \return true if all abscissas are within their allowed ranges
    */
    virtual Bool_t checkRange( const LauAbscissas& abscissas ) const;

    //! Set the unnormalised likelihood
    virtual void setUnNormPDFVal( Double_t unNormPDFVal ) { unNormPDFVal_ = unNormPDFVal; }

    //! Retrieve the specified parameter
    /*!
        \param [in] parName the parameter to retrieve
    */
    virtual LauAbsRValue* findParameter( const TString& parName );

    //! Retrieve the specified parameter
    /*!
        \param [in] parName the parameter to retrieve
    */
    virtual const LauAbsRValue* findParameter( const TString& parName ) const;

    //! Retrieve the random function used for MC generation
    /*!
        \return the random function used for MC generation
    */
    virtual TRandom* getRandomFun() const { return randomFun_; }

    //! Retrieve the abscissa(s)
    /*!
        \return the abscissa(s)
    */
    virtual std::vector<LauAbscissas>& getAbscissas() { return abscissas_; }

    //! Retrieve the abscissa(s)
    /*!
        \return the abscissa(s)
    */
    virtual const std::vector<LauAbscissas>& getAbscissas() const { return abscissas_; }

    //! Retrieve the cached unnormalised likelihood values
    /*!
        \return the cached unnormalised likelihood values
    */
    virtual std::vector<Double_t>& getUnNormPDFValues() { return unNormPDFValues_; }

    //! Retrieve the cached unnormalised likelihood values
    /*!
        \return the cached unnormalised likelihood values
    */
    virtual const std::vector<Double_t>& getUnNormPDFValues() const { return unNormPDFValues_; }

    //! Add parameters to the PDF
    /*!
        \param [in] params the parameters to add
    */
    virtual void addParameters( std::vector<LauAbsRValue*>& params );

    //! Check whether the calcNorm method is running
    /*!
        \return true if the calcNorm method is running
    */
    virtual Bool_t withinNormCalc() const { return withinNormCalc_; }

    //! Set flag to track whether the calcNorm method is running
    /*!
        \param [in] yorn whether the calcNorm method is running
    */
    virtual void withinNormCalc( Bool_t yorn ) { withinNormCalc_ = yorn; }

    //! Check whether the generate method is running
    /*!
        \return true if the generate method is running
    */
    virtual Bool_t withinGeneration() const { return withinGeneration_; }

    //! Set flag to track whether the generate method is running
    /*!
        \param [in] yorn whether the generate method is running
    */
    virtual void withinGeneration( Bool_t yorn ) { withinGeneration_ = yorn; }

    //! Check whether the normalisation weights have been calculated
    /*!
        \return true if the normalisation weights have been calculated
    */
    virtual Bool_t normWeightsDone() const { return normWeightsDone_; }

    //! Set whether the normalisation weights have been calculated
    /*!
        \param [in] yorn whether the normalisation weights have been calculated
    */
    virtual void normWeightsDone( Bool_t yorn ) { normWeightsDone_ = yorn; }

    //! Calculate the weights and abscissas used for normalisation
    virtual void getNormWeights();

    //! Retrieve the abscissa points used for normalisation
    /*!
        \return the abscissa points used for normalisation
    */
    virtual const std::vector<LauAbscissas>& normAbscissas() const { return normAbscissas_; }

    //! Retrieve the weights used for normalisation
    /*!
        \return the weights used for normalisation
    */
    virtual const std::vector<Double_t>& normWeights() const { return normWeights_; }

  private:
    //! Copy constructor (not implemented)
    LauAbsPdf( const LauAbsPdf& rhs );

    //! Copy assignment operator (not implemented)
    LauAbsPdf& operator=( const LauAbsPdf& rhs );

    //! The names of the PDF variables
    std::map<UInt_t, TString> varNames_;

    //! The parameters of the PDF (if any)
    std::vector<LauAbsRValue*> param_;

    //! Normalisation factor of the PDF
    Double_t norm_;

    //! Maximum height of the PDF
    Double_t maxHeight_;

    //! Track whether the height is up to date
    Bool_t heightUpToDate_;

    //! The minimum value(s) of the abscissa(s)
    LauAbscissas minAbscissas_;

    //! The maximum value(s) of the abscissa(s)
    LauAbscissas maxAbscissas_;

    //! The random function used for MC generation
    TRandom* randomFun_;

    //! Cached values of the abscissas
    std::vector<LauAbscissas> abscissas_;

    //! Cached unnormalised likelihood values
    std::vector<Double_t> unNormPDFValues_;

    //! Whether the unnormalised PDF values are cached
    /*! They are only cached if all parameters are fixed. */
    Bool_t cachePDF_;

    //! The unnormalised liklihood value
    Double_t unNormPDFVal_;

    //! number of points to integrate over when normalising
    Int_t nNormPoints_;

    //! The integration method used for normalising the PDF
    IntMethod integMethod_;

    //! Flag to keep track of whether the calcNorm method is running
    Bool_t withinNormCalc_;

    //! Flag to keep track of whether the generate method is running
    Bool_t withinGeneration_;

    //! Whether the normalisation weights have been calculated
    Bool_t normWeightsDone_;

    //! The normalisation abscissas
    std::vector<LauAbscissas> normAbscissas_;

    //! The normalisation weights
    std::vector<Double_t> normWeights_;

    ClassDef( LauAbsPdf, 0 )    // Define the abstract PDF class
};

#endif
