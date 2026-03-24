
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

/*! \file LauPolarGammaCPCoeffSet.hh
    \brief File containing declaration of LauPolarGammaCPCoeffSet class.
*/

/*! \class LauPolarGammaCPCoeffSet
    \brief Class for defining a complex coefficient useful for extracting the CKM angle gamma from B -> D h h Dalitz plots.

    Holds a set of real values that define the complex coefficient of an amplitude component.
    Depending on the type of the D decay, the amplitude has one of the following forms:
    CP-odd eigenstate:    ( x + i * y ) * ( 1 - rB * exp( i * ( deltaB +/- gamma ) ) )
    CP-even eigenstate:   ( x + i * y ) * ( 1 + rB * exp( i * ( deltaB +/- gamma ) ) )
    ADS favoured state:   ( x + i * y ) * ( 1 + rB * rD * exp( i * ( deltaB - deltaD +/- gamma ) ) )
    ADS suppressed state: ( x + i * y ) * ( rD * exp( - i * deltaD ) + rB * exp( i * ( deltaB +/- gamma ) ) )
    [Phys. Rev. D79, 051301 (2009)]
*/

#ifndef LAU_POLARGAMMACP_COEFF_SET
#define LAU_POLARGAMMACP_COEFF_SET

#include "LauAbsCoeffSet.hh"
#include "LauComplex.hh"
#include "LauParameter.hh"

#include "Rtypes.h"

#include <iosfwd>
#include <vector>

class LauPolarGammaCPCoeffSet : public LauAbsCoeffSet {

  public:
    //! The possible D decay modes
    enum DecayType {
        GLW_CPOdd,      /*!< GLW CP-odd, e.g. D0 -> K0 pi0 */
        GLW_CPEven,     /*!< GLW CP-even, e.g. D0 -> K+ K-  */
        ADS_Favoured,   /*!< ADS Favoured, e.g. D0 -> K- pi+ */
        ADS_Suppressed, /*!< ADS Suppressed, e.g. D0 -> K+ pi- */
        GLW_CPOdd_btouOnly, /*!< GLW CP-odd, e.g. D0 -> K0 pi0, where B decay only proceeds via b -> u transition */
        GLW_CPEven_btouOnly, /*!< GLW CP-even, e.g. D0 -> K+ K-, where B decay only proceeds via b -> u transition  */
        ADS_Favoured_btouOnly, /*!< ADS Favoured, e.g. D0 -> K- pi+, where B decay only proceeds via b -> u transition */
        ADS_Suppressed_btouOnly, /*!< ADS Suppressed, e.g. D0 -> K+ pi-, where B decay only proceeds via b -> u transition */
    };

    //! Constructor
    /*!
        \param [in] compName the name of the coefficient set
        \param [in] decayType the type of the D decay
        \param [in] x the real part of the b -> c amplitude
        \param [in] y the imaginary part of the b -> c amplitude
        \param [in] rB the magnitude of the ratio of the b -> u and b -> c amplitudes
        \param [in] deltaB the relative CP-conserving (strong) phase of the b -> u and b -> c amplitudes
        \param [in] gamma the relative CP-violating (weak) phase of the b -> u and b -> c amplitudes
        \param [in] rD the magnitude of the ratio of the favoured and suppressed D-decay amplitudes
        \param [in] deltaD the relative strong phase of the favoured and suppressed D-decay amplitudes
        \param [in] xFixed whether x is fixed
        \param [in] yFixed whether y is fixed
        \param [in] rBFixed whether rB is fixed
        \param [in] deltaBFixed whether deltaB is fixed
        \param [in] gammaFixed whether gamma is fixed
        \param [in] rDFixed whether rD is fixed
        \param [in] deltaDFixed whether deltaD is fixed
        \param [in] rBSecondStage whether rB should be floated only in the second stage of the fit
        \param [in] deltaBSecondStage whether deltaB should be floated only in the second stage of the fit
        \param [in] gammaSecondStage whether gamma should be floated only in the second stage of the fit
        \param [in] rDSecondStage whether rD should be floated only in the second stage of the fit
        \param [in] deltaDSecondStage whether deltaD should be floated only in the second stage of the fit
        \param [in] useGlobalGamma whether gamma should be shared with other resonances
        \param [in] useGlobalADSPars whether rD and deltaD should be shared with other resonances
    */
    LauPolarGammaCPCoeffSet( const TString& compName,
                             const DecayType decayType,
                             const Double_t x,
                             const Double_t y,
                             const Double_t rB,
                             const Double_t deltaB,
                             const Double_t gamma,
                             const Double_t rD,
                             const Double_t deltaD,
                             const Bool_t xFixed,
                             const Bool_t yFixed,
                             const Bool_t rBFixed,
                             const Bool_t deltaBFixed,
                             const Bool_t gammaFixed,
                             const Bool_t rDFixed,
                             const Bool_t deltaDFixed,
                             const Bool_t rBSecondStage     = kFALSE,
                             const Bool_t deltaBSecondStage = kFALSE,
                             const Bool_t gammaSecondStage  = kFALSE,
                             const Bool_t rDSecondStage     = kFALSE,
                             const Bool_t deltaDSecondStage = kFALSE,
                             const Bool_t useGlobalGamma    = kFALSE,
                             const Bool_t useGlobalADSPars  = kFALSE );

    //! Destructor
    virtual ~LauPolarGammaCPCoeffSet() {}

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
        This method is not supported by this class because there are more than four parameters, hence there is not a unique solution.

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
    LauPolarGammaCPCoeffSet( const LauPolarGammaCPCoeffSet& rhs,
                             CloneOption cloneOption = All,
                             Double_t constFactor    = 1.0 );

    //! Copy assignment operator (not implemented)
    /*!
        \param [in] rhs the coefficient to clone
    */
    LauPolarGammaCPCoeffSet& operator=( const LauPolarGammaCPCoeffSet& rhs );

    //! Prepend the base name and index to the name of a parameter
    /*!
        \param [in,out] par pointer to the parameter
        \param [in] oldBaseName the old base name, which might need to be removed before adding the new one
    */
    virtual void adjustName( LauParameter* par, const TString& oldBaseName );

    //! Update the amplitudes based on the new values of the parameters
    void updateAmplitudes();

    //! The type of the D decay
    const DecayType decayType_;

    // the actual fit parameters
    // (need to be pointers so they can be cloned)

    //! The real part of the b -> c amplitude
    LauParameter* x_;

    //! The imaginary part of the b -> c amplitude
    LauParameter* y_;

    //! the magnitude of the ratio of the b -> u and b -> c amplitudes
    LauParameter* rB_;

    //! the relative CP-conserving (strong) phase of the b -> u and b -> c amplitudes
    LauParameter* deltaB_;

    //! the relative CP-violating (weak) phase of the b -> u and b -> c amplitudes
    LauParameter* gamma_;

    //! the magnitude of the ratio of the favoured and suppressed D-decay amplitudes
    LauParameter* rD_;

    //! the relative strong phase of the favoured and suppressed D-decay amplitudes
    LauParameter* deltaD_;

    //! The CP-violating phase (shared by multiple resonances)
    static LauParameter* gammaGlobal_;

    //! the magnitude of the ratio of the favoured and suppressed D-decay amplitudes (shared by multiple resonances)
    static LauParameter* rDGlobal_;

    //! the relative strong phase of the favoured and suppressed D-decay amplitudes (shared by multiple resonances)
    static LauParameter* deltaDGlobal_;

    //! Whether the global gamma is used for this resonance
    const Bool_t useGlobalGamma_;

    //! Whether the global rD and deltaD are used for this resonance
    const Bool_t useGlobalADSPars_;

    //! The b -> c part of the complex coefficient
    LauComplex nonCPPart_;
    //! The b -> u part of the complex coefficient for the particle
    LauComplex cpPart_;
    //! The b -> u part of the complex coefficient for the antiparticle
    LauComplex cpAntiPart_;

    //! The particle complex coefficient
    LauComplex particleCoeff_;
    //! The antiparticle complex coefficient
    LauComplex antiparticleCoeff_;

    //! The CP asymmetry
    LauParameter acp_;

    ClassDef( LauPolarGammaCPCoeffSet, 0 )
};

#endif
