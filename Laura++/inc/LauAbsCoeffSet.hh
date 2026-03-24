
/*
Copyright 2006 University of Warwick

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

/*! \file LauAbsCoeffSet.hh
    \brief File containing declaration of LauAbsCoeffSet class.
*/

/*! \class LauAbsCoeffSet
    \brief Class for defining the abstract interface for complex coefficient classes.

    Class for defining the abstract interface for complex coefficient classes.
    Some common code is implemented but most methods are not.
*/

#ifndef LAU_ABS_COEFF_SET
#define LAU_ABS_COEFF_SET

#include "TString.h"

#include <iosfwd>
#include <vector>

class TRandom;
class LauComplex;
class LauParameter;

class LauAbsCoeffSet {

  public:
    //! Options for cloning operation
    enum CloneOption {
        All,          /*!< no special operation, all parameters cloned */
        TiePhase,     /*!< phase cloned, magnitude free to vary */
        TieMagnitude, /*!< magnitude cloned, phase free to vary */
        TieRealPart,  /*!< real part cloned, imaginary part free to vary */
        TieImagPart,  /*!< imaginary part cloned, real part free to vary */
        TieCPPars     /*!< CP-violating parameters cloned, CP-conserving ones free to vary */
    };

    //! Destructor
    virtual ~LauAbsCoeffSet() {}

    //! Retrieve the parameters of the coefficient so that they can be loaded into a fit
    /*!
        \return the parameters of the coefficient
    */
    virtual std::vector<LauParameter*> getParameters() = 0;

    //! Print the current values of the parameters
    virtual void printParValues() const = 0;

    //! Print the column headings for a results table
    /*!
        \param [out] stream the stream to print to
    */
    virtual void printTableHeading( std::ostream& stream ) const = 0;

    //! Print the parameters of the complex coefficient as a row in the results table
    /*!
        \param [out] stream the stream to print to
    */
    virtual void printTableRow( std::ostream& stream ) const = 0;

    //! Randomise the starting values of the parameters for a fit
    virtual void randomiseInitValues() = 0;

    //! Make sure values are in "standard" ranges, e.g. phases should be between -pi and pi
    virtual void finaliseValues() = 0;

    //! Retrieve the complex coefficient for a particle
    /*!
        \return the complex coefficient for a particle
    */
    virtual const LauComplex& particleCoeff() = 0;

    //! Retrieve the complex coefficient for an antiparticle
    /*!
        \return the complex coefficient for an antiparticle
    */
    virtual const LauComplex& antiparticleCoeff() = 0;

    //! Set the parameters based on the complex coefficients for particles and antiparticles
    /*!
        \param [in] coeff the complex coefficient for a particle
        \param [in] coeffBar the complex coefficient for an antiparticle
        \param [in] init whether or not the initial and generated values should also be adjusted
    */
    virtual void setCoeffValues( const LauComplex& coeff, const LauComplex& coeffBar, Bool_t init ) = 0;

    //! Calculate the CP asymmetry
    /*!
        \return the CP asymmetry
    */
    virtual LauParameter acp() = 0;

    //! Create a clone of the coefficient set
    /*!
        \param [in] newName the clone's name
        \param [in] cloneOption special option for the cloning operation
        \param [in] constFactor a constant factor by which to multiply the cloned parameters
        \return a clone of the coefficient set
    */
    virtual LauAbsCoeffSet* createClone( const TString& newName,
                                         CloneOption cloneOption = All,
                                         Double_t constFactor    = 1.0 ) = 0;

    //! Retrieve the name of the coefficient set
    /*!
        The name should correspond to the name of the resonance in the model.

        \return the name of the coefficient set
    */
    virtual TString name() const { return name_; }

    //! Set the name of the coefficient set
    /*!
        The name should correspond to the name of the resonance in the model.

        \param [in] theName the name to set
    */
    virtual void name( const TString& theName ) { name_ = theName; }

    //! Retrieve the base name of the coefficient set
    /*!
        The base name is generally of the form "Ai", where i is an integer.
        This is used in the fit results ntuple.

        \return the base name of the coefficient set
    */
    virtual const TString& baseName() const { return basename_; }

    //! Set the base name of the coefficient set
    /*!
        The base name is generally of the form "Ai", where i is an integer.
        This is used in the fit results ntuple.

        \param [in] theBasename the base name to set
    */
    virtual void baseName( const TString& theBasename ) { basename_ = theBasename; }

    //! Retrieve the index number of the coefficient set
    /*!
        \return the index number of the coefficient set
    */
    virtual UInt_t index() const { return index_; }

    //! Set the index number of the coefficient set
    /*!
        \param [in] newIndex the new index
    */
    virtual void index( UInt_t newIndex );

    //! Set the value of the named parameter
    /*!
        \param [in] parName the name of the parameter to adjust
        \param [in] value the new value for the parameter to take
        \param [in] init whether or not the initial and generated values should also be adjusted
    */
    virtual void setParameterValue( const TString& parName, Double_t value, Bool_t init );

    //! Set the error of the named parameter
    /*!
        This is particularly useful for tuning the step size used by MINUIT

        \param [in] parName the name of the parameter to adjust
        \param [in] error the new error value for the parameter to take
    */
    virtual void setParameterError( const TString& parName, Double_t error );

    //! Set the named parameter to be fixed in the fit
    /*!
        \param [in] parName the name of the parameter to adjust
    */
    virtual void fixParameter( const TString& parName );

    //! Set the named parameter to float in the fit
    /*!
        \param [in] parName the name of the parameter to adjust
    */
    virtual void floatParameter( const TString& parName );

    //! Blind the named parameter
    /*!
        See LauBlind documentation for details of blinding procedure

        \param [in] parName the name of the parameter to adjust
        \param [in] blindingString the unique blinding string used to seed the random number generator
        \param [in] width the width of the Gaussian from which the offset should be sampled
    */
    virtual void blindParameter( const TString& parName,
                                 const TString& blindingString,
                                 const Double_t width );

    //! Add Gaussian constraint to the named parameter
    /*!
        \param [in] parName the name of the parameter to adjust
        \param [in] mean the mean of the Gaussian constraint
        \param [in] width the width of the Gaussian constraint
    */
    virtual void addGaussianConstraint( const TString& parName,
                                        const Double_t mean,
                                        const Double_t width );

    //! Add suffix to the name of the given parameter
    /*!
        \param [in] parName the name of the parameter to adjust
        \param [in] suffix the suffix to add to the parameter name
    */
    virtual void addSuffixToParameterName( const TString& parName, const TString& suffix );

    //! Set the allowed range for magnitude parameters
    /*!
        \param [in] minMag the lower edge of the range
        \param [in] maxMag the upper edge of the range
    */
    static void setMagnitudeRange( Double_t minMag, Double_t maxMag )
    {
        minMagnitude_ = minMag;
        maxMagnitude_ = maxMag;
    }

    //! Set the allowed range for phase parameters
    /*!
        \param [in] minPhase the lower edge of the range
        \param [in] maxPhase the upper edge of the range
    */
    static void setPhaseRange( Double_t minPhase, Double_t maxPhase )
    {
        minPhase_ = minPhase;
        maxPhase_ = maxPhase;
    }

    //! Set the allowed range for real/imaginary part parameters
    /*!
        \param [in] minPar the lower edge of the range
        \param [in] maxPar the upper edge of the range
    */
    static void setRealImagRange( Double_t minPar, Double_t maxPar )
    {
        minRealImagPart_ = minPar;
        maxRealImagPart_ = maxPar;
    }

    //! Set the allowed range for CP-violating parameters
    /*!
        \param [in] minPar the lower edge of the range
        \param [in] maxPar the upper edge of the range
    */
    static void setCPParRange( Double_t minPar, Double_t maxPar )
    {
        minDelta_ = minPar;
        maxDelta_ = maxPar;
    }

    //! Set the randomiser
    /*!
        Set the random number generator to use for randomising parameter starting values.
        Will default to LauRandom::zeroSeedRandom if not explicitly supplied via this function.

        \param [in] randomiser the random number generator to use for randomising parameter starting values
    */
    static void setRandomiser( TRandom* randomiser ) { randomiser_ = randomiser; }

    //! Access the randomiser
    /*!
        \return the random number generator to use for randomising parameter starting values
    */
    static TRandom* getRandomiser();

  protected:
    //! Constructor
    /*!
        \param [in] theName the name of the coefficient set
        \param [in] theBaseName the single character base for the parameter names
    */
    LauAbsCoeffSet( const TString& theName, const TString& theBaseName = "A" );

    //! Find the parameter with the given name
    /*!
        \param [in] parName the name of the parameter to be found
        return the retrieved parameter
    */
    LauParameter* findParameter( const TString& parName );

    //! Prepend the base name and index to the name of a parameter
    /*!
        \param [in,out] par pointer to the parameter
        \param [in] oldBaseName the old base name, which might need to be removed before adding the new one
    */
    virtual void adjustName( LauParameter* par, const TString& oldBaseName );

    //! Minimum allowed value of magnitude parameters
    static Double_t minMagnitude_;
    //! Maximum allowed value of magnitude parameters
    static Double_t maxMagnitude_;
    //! Minimum allowed value of phase parameters
    static Double_t minPhase_;
    //! Maximum allowed value of phase parameters
    static Double_t maxPhase_;
    //! Minimum allowed value of real/imaginary part parameters
    static Double_t minRealImagPart_;
    //! Maximum allowed value of real/imaginary part parameters
    static Double_t maxRealImagPart_;
    //! Minimum allowed value of CP-violating real/imaginary part parameters
    static Double_t minDelta_;
    //! Maximum allowed value of CP-violating real/imaginary part parameters
    static Double_t maxDelta_;

    //! Copy constructor
    /*!
        \param [in] rhs the coefficient to clone
    */
    LauAbsCoeffSet( const LauAbsCoeffSet& rhs );

    //! Copy assignment operator
    /*!
        \param [in] rhs the coefficient to clone
    */
    LauAbsCoeffSet& operator=( const LauAbsCoeffSet& rhs );

  private:
    //! Random number generator to use for randomising parameter starting values
    static TRandom* randomiser_;

    //! The name of the coefficient set
    TString name_;

    //! The base name of the coefficient set
    TString basename_;

    //! The index number of the coefficient set
    UInt_t index_;

    ClassDef( LauAbsCoeffSet, 0 )
};

#endif
