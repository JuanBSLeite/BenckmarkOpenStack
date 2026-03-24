
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

/*! \file LauCartesianGammaCPCoeffSet.hh
    \brief File containing declaration of LauCartesianGammaCPCoeffSet class.
*/

/*! \class LauCartesianGammaCPCoeffSet
    \brief Class for defining a complex coefficient using the Cartesian gamma CP convention.

    Holds a set of real values that define the complex coefficient of an amplitude component.
    The amplitude has the form ( x + i * y ) * ( 1 + xCP +/- delta_xCP + i * ( yCP +/- delta_yCP ) ).
    [Phys. Rev. D79, 051301 (2009)]
*/

#ifndef LAU_CARTESIANGAMMACP_COEFF_SET
#define LAU_CARTESIANGAMMACP_COEFF_SET

#include "LauAbsCoeffSet.hh"
#include "LauComplex.hh"
#include "LauParameter.hh"

#include "Rtypes.h"

#include <iosfwd>
#include <vector>

class LauCartesianGammaCPCoeffSet : public LauAbsCoeffSet {

  public:
    //! Constructor
    /*!
        \param [in] compName the name of the coefficient set
        \param [in] x the real nonCP part
        \param [in] y the imaginary nonCP part
        \param [in] xCP the average real CP part
        \param [in] yCP the average imaginary CP part
        \param [in] deltaXCP the asymmetric real CP part
        \param [in] deltaYCP the asymmetric imaginary CP part
        \param [in] xFixed whether x is fixed
        \param [in] yFixed whether y is fixed
        \param [in] xCPFixed whether xCP is fixed
        \param [in] yCPFixed whether yCP is fixed
        \param [in] deltaXCPFixed whether deltaXCP is fixed
        \param [in] deltaYCPFixed whether deltaYCP is fixed
        \param [in] deltaXCPSecondStage whether deltaXCP should be floated only in the second stage of the fit
        \param [in] deltaYCPSecondStage whether deltaYCP should be floated only in the second stage of the fit
    */
    LauCartesianGammaCPCoeffSet( const TString& compName,
                                 const Double_t x,
                                 const Double_t y,
                                 const Double_t xCP,
                                 const Double_t yCP,
                                 const Double_t deltaXCP,
                                 const Double_t deltaYCP,
                                 const Bool_t xFixed,
                                 const Bool_t yFixed,
                                 const Bool_t xCPFixed,
                                 const Bool_t yCPFixed,
                                 const Bool_t deltaXCPFixed,
                                 const Bool_t deltaYCPFixed,
                                 const Bool_t deltaXCPSecondStage = kFALSE,
                                 const Bool_t deltaYCPSecondStage = kFALSE );

    //! Destructor
    virtual ~LauCartesianGammaCPCoeffSet() {}

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
        This method is not supported by this class because there are more than four parameters so there is not a unique solution.

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
    LauCartesianGammaCPCoeffSet( const LauCartesianGammaCPCoeffSet& rhs,
                                 CloneOption cloneOption = All,
                                 Double_t constFactor    = 1.0 );

    //! Copy assignment operator (not implemented)
    /*!
        \param [in] rhs the coefficient to clone
    */
    LauCartesianGammaCPCoeffSet& operator=( const LauCartesianGammaCPCoeffSet& rhs );

    // the actual fit parameters
    // (need to be pointers so they can be cloned)
    //! The nonCP real part
    LauParameter* x_;
    //! The nonCP imaginary part
    LauParameter* y_;
    //! The average CP real part
    LauParameter* xCP_;
    //! The average CP imaginary part
    LauParameter* yCP_;
    //! The asymmetric CP real part
    LauParameter* deltaXCP_;
    //! The asymmetric CP imaginary part
    LauParameter* deltaYCP_;

    //! The nonCP part of the complex coefficient
    LauComplex nonCPPart_;
    //! The CP part of the complex coefficient for the particle
    LauComplex cpPart_;
    //! The CP part of the complex coefficient for the antiparticle
    LauComplex cpAntiPart_;

    //! The particle complex coefficient
    LauComplex particleCoeff_;
    //! The antiparticle complex coefficient
    LauComplex antiparticleCoeff_;

    //! The CP asymmetry
    LauParameter acp_;

    ClassDef( LauCartesianGammaCPCoeffSet, 0 )
};

#endif
