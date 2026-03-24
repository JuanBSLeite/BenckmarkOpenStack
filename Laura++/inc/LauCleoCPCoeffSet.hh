
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

/*! \file LauCleoCPCoeffSet.hh
    \brief File containing declaration of LauCleoCPCoeffSet class.
*/

/*! \class LauCleoCPCoeffSet
    \brief Class for defining a complex coefficient using the Cleo CP convention.

    Holds a set of real values that define the complex coefficient of an amplitude component.
    The amplitude has the form (a +/- b) * exp( i*(delta +/- phi) ) where
    a is the average magnitude,
    b is the asymmetric magnitude,
    delta is the strong phase and
    phi is the weak phase.
    [Phys.Rev. D70 (2004) 091101]
*/

#ifndef LAU_CLEOCP_COEFF_SET
#define LAU_CLEOCP_COEFF_SET

#include "LauAbsCoeffSet.hh"
#include "LauComplex.hh"
#include "LauParameter.hh"

#include "Rtypes.h"

#include <iosfwd>
#include <vector>

class LauCleoCPCoeffSet : public LauAbsCoeffSet {

  public:
    //! Constructor
    /*!
        \param [in] compName the name of the coefficient set
        \param [in] a the magnitude a
        \param [in] delta the strong phase
        \param [in] b the magnitude b
        \param [in] phi the weak phase
        \param [in] aFixed whether a is fixed
        \param [in] deltaFixed whether delta is fixed
        \param [in] bFixed whether b is fixed
        \param [in] phiFixed whether phi is fixed
        \param [in] bSecondStage whether b should be floated only in the second stage of the fit
        \param [in] phiSecondStage whether phi should be floated only in the second stage of the fit
    */
    LauCleoCPCoeffSet( const TString& compName,
                       Double_t a,
                       Double_t delta,
                       Double_t b,
                       Double_t phi,
                       Bool_t aFixed,
                       Bool_t deltaFixed,
                       Bool_t bFixed,
                       Bool_t phiFixed,
                       Bool_t bSecondStage   = kFALSE,
                       Bool_t phiSecondStage = kFALSE );

    //! Destructor
    virtual ~LauCleoCPCoeffSet() {}

    //! Retrieve the parameters of the coefficient, e.g. so that they can be loaded into a fit
    /*!
        \return the parameters of the coefficient
    */
    virtual std::vector<LauParameter*> getParameters();

    //! Print the current values of the parameters
    virtual void printParValues() const;

    //! Print the column headings for a results table
    /*!
        \param [out] stream the stream to print to
    */
    virtual void printTableHeading( std::ostream& stream ) const;

    //! Print the parameters of the complex coefficient as a row in the results table
    /*!
        \param [out] stream the stream to print to
    */
    virtual void printTableRow( std::ostream& stream ) const;

    //! Randomise the starting values of the parameters for a fit
    virtual void randomiseInitValues();

    //! Make sure values are in "standard" ranges, e.g. phases should be between -pi and pi
    virtual void finaliseValues();

    //! Retrieve the complex coefficient for a particle
    /*!
        \return the complex coefficient for a particle
    */
    virtual const LauComplex& particleCoeff();

    //! Retrieve the complex coefficient for an antiparticle
    /*!
        \return the complex coefficient for an antiparticle
    */
    virtual const LauComplex& antiparticleCoeff();

    //! Set the parameters based on the complex coefficients for particles and antiparticles
    /*!
        \param [in] coeff the complex coefficient for a particle
        \param [in] coeffBar the complex coefficient for an antiparticle
        \param [in] init whether or not the initial and generated values should also be adjusted
    */
    virtual void setCoeffValues( const LauComplex& coeff, const LauComplex& coeffBar, Bool_t init );

    //! Calculate the CP asymmetry
    /*!
        \return the CP asymmetry
    */
    virtual LauParameter acp();

    //! Create a clone of the coefficient set
    /*!
        \param [in] newName the clone's name
        \param [in] cloneOption special option for the cloning operation
        \param [in] constFactor a constant factor to multiply the clone's parameters by
        \return a clone of the coefficient set
    */
    virtual LauAbsCoeffSet* createClone( const TString& newName,
                                         CloneOption cloneOption = All,
                                         Double_t constFactor    = 1.0 );

  private:
    //! Copy constructor
    /*!
        This creates cloned parameters, not copies.

        \param [in] rhs the coefficient to clone
        \param [in] cloneOption special option for the cloning operation
        \param [in] constFactor a constant factor to multiply the clone's parameters by
    */
    LauCleoCPCoeffSet( const LauCleoCPCoeffSet& rhs,
                       CloneOption cloneOption = All,
                       Double_t constFactor    = 1.0 );

    //! Copy assignment operator (not implemented)
    /*!
        \param [in] rhs the coefficient to clone
    */
    LauCleoCPCoeffSet& operator=( const LauCleoCPCoeffSet& rhs );

    // the actual fit parameters
    // (need to be pointers so they can be cloned)
    //! The magnitude a
    LauParameter* a_;
    //! The magnitude b
    LauParameter* b_;
    //! The strong phase
    LauParameter* delta_;
    //! The weak phase
    LauParameter* phi_;

    //! The particle complex coefficient
    LauComplex particleCoeff_;
    //! The antiparticle complex coefficient
    LauComplex antiparticleCoeff_;

    //! The CP asymmetry
    LauParameter acp_;

    ClassDef( LauCleoCPCoeffSet, 0 )
};

#endif
