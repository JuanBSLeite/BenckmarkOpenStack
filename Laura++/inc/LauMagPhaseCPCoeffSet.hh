
/*
Copyright 2011 University of Warwick

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

/*! \file LauMagPhaseCPCoeffSet.hh
    \brief File containing declaration of LauMagPhaseCPCoeffSet class.
*/

/*! \class LauMagPhaseCPCoeffSet
    \brief Class for defining a complex coefficient using seperate magnitudes and phases for particles and antiparticles.

    Holds a set of real values that define the complex coefficient of an amplitude component.
    The amplitudes have the form:
    c    = mag    * exp(i*phase)
    cBar = magBar * exp(i*phaseBar)
    where mag and magBar are the magnitudes for particle and antiparticle and
    phase and phaseBar are the phases for particle and antiparticle.
*/

#ifndef LAU_MAGPHASECP_COEFF_SET
#define LAU_MAGPHASECP_COEFF_SET

#include "LauAbsCoeffSet.hh"
#include "LauComplex.hh"
#include "LauParameter.hh"

#include "Rtypes.h"

#include <iosfwd>
#include <vector>

class LauMagPhaseCPCoeffSet : public LauAbsCoeffSet {

  public:
    //! Constructor
    /*!
        \param [in] compName the name of the coefficient set
        \param [in] mag the magnitude for particles
        \param [in] phase the phase for particles
        \param [in] magBar the magnitude for antiparticles
        \param [in] phaseBar the phase for antiparticles
        \param [in] magFixed whether mag is fixed
        \param [in] phaseFixed whether phase is fixed
        \param [in] magBarFixed whether magBar is fixed
        \param [in] phaseBarFixed whether phaseBar is fixed
    */
    LauMagPhaseCPCoeffSet( const TString& compName,
                           Double_t mag,
                           Double_t phase,
                           Double_t magBar,
                           Double_t phaseBar,
                           Bool_t magFixed,
                           Bool_t phaseFixed,
                           Bool_t magBarFixed,
                           Bool_t phaseBarFixed );

    //! Destructor
    virtual ~LauMagPhaseCPCoeffSet() {}

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
    LauMagPhaseCPCoeffSet( const LauMagPhaseCPCoeffSet& rhs,
                           CloneOption cloneOption = All,
                           Double_t constFactor    = 1.0 );

    //! Copy assignment operator (not implemented)
    /*!
        \param [in] rhs the coefficient to clone
    */
    LauMagPhaseCPCoeffSet& operator=( const LauMagPhaseCPCoeffSet& rhs );

    // the actual fit parameters
    // (need to be pointers so they can be cloned)
    //! The magnitude for particles
    LauParameter* mag_;
    //! The phase for particles
    LauParameter* phase_;
    //! The magnitude for antiparticles
    LauParameter* magBar_;
    //! The phase for antiparticles
    LauParameter* phaseBar_;

    //! The particle complex coefficient
    LauComplex particleCoeff_;
    //! The antiparticle complex coefficient
    LauComplex antiparticleCoeff_;

    //! The CP asymmetry
    LauParameter acp_;

    ClassDef( LauMagPhaseCPCoeffSet, 0 )
};

#endif
