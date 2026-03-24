
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

/*! \file LauCartesianCPCoeffSet.hh
    \brief File containing declaration of LauCartesianCPCoeffSet class.
*/

/*! \class LauCartesianCPCoeffSet
    \brief Class for defining a complex coefficient using the Cartesian CP convention.

    Holds a set of real values that define the complex coefficient of an amplitude component.
    The amplitude has the form x +/- delta_x + i * ( y +/- delta_y ).
    [Phys.Rev. D78 (2008) 012004]
*/

#ifndef LAU_CARTESIANCP_COEFF_SET
#define LAU_CARTESIANCP_COEFF_SET

#include "LauAbsCoeffSet.hh"
#include "LauComplex.hh"
#include "LauParameter.hh"

#include "Rtypes.h"

#include <iosfwd>
#include <vector>

class LauCartesianCPCoeffSet : public LauAbsCoeffSet {

  public:
    //! Constructor
    /*!
        \param [in] compName the name of the coefficient set
        \param [in] x the average real part
        \param [in] y the average imaginary part
        \param [in] deltaX the asymmetric real part
        \param [in] deltaY the asymmetric imaginary part
        \param [in] xFixed whether x is fixed
        \param [in] yFixed whether y is fixed
        \param [in] deltaXFixed whether deltaX is fixed
        \param [in] deltaYFixed whether deltaY is fixed
        \param [in] deltaXSecondStage whether deltaX should be floated only in the second stage of the fit
        \param [in] deltaYSecondStage whether deltaY should be floated only in the second stage of the fit
    */
    LauCartesianCPCoeffSet( const TString& compName,
                            Double_t x,
                            Double_t y,
                            Double_t deltaX,
                            Double_t deltaY,
                            Bool_t xFixed,
                            Bool_t yFixed,
                            Bool_t deltaXFixed,
                            Bool_t deltaYFixed,
                            Bool_t deltaXSecondStage = kFALSE,
                            Bool_t deltaYSecondStage = kFALSE );

    //! Destructor
    virtual ~LauCartesianCPCoeffSet() {}

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
    LauCartesianCPCoeffSet( const LauCartesianCPCoeffSet& rhs,
                            CloneOption cloneOption = All,
                            Double_t constFactor    = 1.0 );

    //! Copy assignment operator (not implemented)
    /*!
        \param [in] rhs the coefficient to clone
    */
    LauCartesianCPCoeffSet& operator=( const LauCartesianCPCoeffSet& rhs );

    // the actual fit parameters
    // (need to be pointers so they can be cloned)
    //! The average real part
    LauParameter* x_;
    //! The average imaginary part
    LauParameter* y_;
    //! The asymmetric real part
    LauParameter* deltaX_;
    //! The asymmetric imaginary part
    LauParameter* deltaY_;

    //! The particle complex coefficient
    LauComplex particleCoeff_;
    //! The antiparticle complex coefficient
    LauComplex antiparticleCoeff_;

    //! The CP asymmetry
    LauParameter acp_;

    ClassDef( LauCartesianCPCoeffSet, 0 )
};

#endif
