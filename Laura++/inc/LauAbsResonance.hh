
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

/*! \file LauAbsResonance.hh
    \brief File containing declaration of LauAbsResonance class.
*/

/*! \class LauAbsResonance
    \brief Abstract class for defining type for resonance amplitude models (Breit-Wigner, Flatte etc.)

    Abstract Class for defining the type for all classes used to model
    resonances in the Dalitz plot, such as Breit-Wigner functions.
    In addition, some common functionality is implemented, including data such
    as the mass and width of the desired state.
*/

#ifndef LAU_ABS_RESONANCE
#define LAU_ABS_RESONANCE

#include "LauBlattWeisskopfFactor.hh"
#include "LauComplex.hh"
#include "LauParameter.hh"

#include "TString.h"

class LauDaughters;
class LauKinematics;
class LauResonanceInfo;

class LauAbsResonance {

  public:
    //! Define the allowed resonance types
    enum LauResonanceModel {
        BW,         /*!< simple Breit-Wigner */
        RelBW,      /*!< relativistic Breit-Wigner */
        GS,         /*!< a modified Breit-Wigner from Gounaris-Sakurai */
        Flatte,     /*!< Flatte or coupled-channel Breit-Wigner */
        Sigma,      /*!< special shape for the sigma or f_0(600) */
        Kappa,      /*!< special shape for the kappa, a low-mass Kpi scalar */
        Dabba,      /*!< special shape for the dabba, a low-mass Dpi scalar */
        LASS,       /*!< the LASS amplitude to describe the Kpi S-wave */
        LASS_BW,    /*!< the resonant part of the LASS amplitude */
        LASS_NR,    /*!< the nonresonant part of the LASS amplitude */
        EFKLLM,     /*!< a form-factor-based description of the Kpi S-wave */
        KMatrix,    /*!< description using K-matrix and P-vector */
        FlatNR,     /*!< a uniform nonresonant amplitude */
        NRModel,    /*!< a theoretical model nonresonant amplitude */
        BelleNR,    /*!< an empirical exponential nonresonant amplitude */
        PowerLawNR, /*!< an empirical power law nonresonant amplitude */
        BelleSymNR, /*!< an empirical exponential nonresonant amplitude for symmetrised DPs */
        BelleSymNRNoInter, /*!< an empirical exponential nonresonant amplitude for symmetrised DPs without interference */
        TaylorNR,  /*!< an empirical Taylor expansion nonresonant amplitude for symmetrised DPs */
        PolNR,     /*!< an empirical polynomial nonresonant amplitude */
        Pole,      /*!< scalar Pole lineshape */
        PolarFFNR, /*!< Polar Form Factor nonresonant amplitude */
        PolarFFSymNR, /*!< Polar Form Factor nonresonant amplitude for symmetrised DPs */
        PolarFFSymNRNoInter, /*!< Polar Form Factor nonresonant amplitude for symmetrised DPs without interference */
        Rescattering,        /*!< KK-PiPi inelastic scattering amplitude */
        Rescattering2,       /*!< KK-PiPi inelastic scattering amplitude */
        RescatteringNoInter, /*!< KK-PiPi inelastic scattering amplitude */
        MIPW_MagPhase, /*!< a model independent partial wave - magnitude and phase representation */
        MIPW_RealImag, /*!< a model independent partial wave - real and imaginary part representation */
        GaussIncoh,      /*!< an incoherent Gaussian shape */
        RhoOmegaMix_GS,  /*!< mass mixing model using GS for res 1 and RBW for res 2 */
        RhoOmegaMix_RBW, /*!< mass mixing model using two RBWs  */
        RhoOmegaMix_GS_1, /*!< mass mixing model using GS for res 1 and RBW for res 2, with denominator factor = 1 */
        RhoOmegaMix_RBW_1 /*!< mass mixing model using two RBWs, with denominator factor = 1 */
    };

    //! Define the allowed spin formalisms
    enum LauSpinType {
        Zemach_P,     /*!< Zemach tensor formalism, bachelor momentum in resonance rest frame */
        Zemach_Pstar, /*!< Zemach tensor formalism, bachelor momentum in parent rest frame */
        Covariant,    /*!< Covariant tensor formalism */
        Legendre      /*!< Legendre polynomials only */
    };

    //! Is the resonance model incoherent?
    /*!
        \param [in] model the resonance model
        \return true if the model is incoherent
    */
    static bool isIncoherentModel( LauResonanceModel model );

    //! Constructor (for use by standard resonances)
    /*!
        \param [in] resInfo the object containing information on the resonance name, mass, width, spin, charge, etc.
        \param [in] resPairAmpInt the number of the daughter not produced by the resonance
        \param [in] daughters the daughter particles
    */
    LauAbsResonance( LauResonanceInfo* resInfo,
                     const Int_t resPairAmpInt,
                     const LauDaughters* daughters );

    //! Constructor (for use by K-matrix components)
    /*!
        \param [in] resName the name of the component
        \param [in] resPairAmpInt the number of the daughter not produced by the resonance
        \param [in] daughters the daughter particles
        \param [in] resSpin the spin of the final channel into which the K-matrix scatters
    */
    LauAbsResonance( const TString& resName,
                     const Int_t resPairAmpInt,
                     const LauDaughters* daughters,
                     const Int_t resSpin );

    //! Destructor
    virtual ~LauAbsResonance();

    //! Initialise the model
    virtual void initialise() = 0;

    //! Calculate the complex amplitude
    /*!
        \param [in] kinematics the kinematic variables of the current event
        \return the complex amplitude
    */
    virtual LauComplex amplitude( const LauKinematics* kinematics );

    //! Get the resonance model type
    /*!
        \return the resonance model type
    */
    virtual LauResonanceModel getResonanceModel() const = 0;

    //! Get the spin type
    /*!
        \return the spin formalism
    */
    LauSpinType getSpinType() const { return spinType_; }

    //! Get the name of the resonance
    /*!
        \return the resonance name
    */
    const TString& getResonanceName() const { return resName_; }

    //! Get the name of the resonance
    /*!
        \return the resonance name
    */
    const TString& getSanitisedName() const { return sanitisedName_; }

    //! Get the integer to identify which DP axis the resonance belongs to
    /*!
        \return the DP axis identification number, the ID of the bachelor
    */
    Int_t getPairInt() const { return resPairAmpInt_; }

    //! Get the spin of the resonance
    /*!
        \return the resonance spin
    */
    Int_t getSpin() const { return resSpin_; }

    //! Get the charge of the resonance
    /*!
        \return the resonance charge
    */
    Int_t getCharge() const { return resCharge_; }

    //! Get the mass of the resonance
    /*!
        \return the resonance mass
    */
    Double_t getMass() const { return ( resMass_ != 0 ) ? resMass_->unblindValue() : -1.0; }

    //! Get the width of the resonance
    /*!
        \return the resonance width
    */
    Double_t getWidth() const { return ( resWidth_ != 0 ) ? resWidth_->unblindValue() : -1.0; }

    //! Get the mass parameter of the resonance
    /*!
        \return the resonance mass parameter
    */
    LauParameter* getMassPar() { return resMass_; }

    //! Get the width parameter of the resonance
    /*!
        \return the resonance width parameter
    */
    LauParameter* getWidthPar() { return resWidth_; }

    //! Retrieve the resonance parameters, e.g. so that they can be loaded into a fit
    /*!
        \return floating parameters of the resonance
    */
    virtual const std::vector<LauParameter*>& getFloatingParameters()
    {
        return this->getParameters();
    };

    //! Is the amplitude pre-symmetrised?
    /*!
        The default value is kFALSE, so pre-symmetrised lineshapes should override this.

        \return whether the amplitude is already symmetrised
    */
    virtual Bool_t preSymmetrised() const { return kFALSE; }

    //! Get the helicity flip flag
    /*!
        \return the flip helicity flag
    */
    Bool_t flipHelicity() const { return flipHelicity_; }

    //! Set the helicity flip flag
    /*!
        \param [in] boolean the helicity flip status
    */
    void flipHelicity( const Bool_t boolean ) { flipHelicity_ = boolean; }

    //! Get the ignore momenta flag
    /*!
        Whether to ignore the momentum factors in both the spin factor and the mass-dependent width

        \return the ignore momenta flag
    */
    Bool_t ignoreMomenta() const { return ignoreMomenta_; }

    //! Set the ignore momenta flag
    /*!
        Whether to ignore the momentum factors in both the spin factor and the mass-dependent width

        \param [in] boolean the ignore momenta status
    */
    void ignoreMomenta( const Bool_t boolean ) { ignoreMomenta_ = boolean; }

    //! Get the ignore spin flag
    /*!
        Whether to set the spinTerm to unity always

        \return the ignore spin flag
    */
    Bool_t ignoreSpin() const { return ignoreSpin_; }

    //! Set the ignore spin flag
    /*!
        Whether to set the spinTerm to unity always

        \param [in] boolean the ignore spin status
    */
    void ignoreSpin( const Bool_t boolean ) { ignoreSpin_ = boolean; }

    //! Get the ignore barrier factor scaling flag
    /*!
        Whether to ignore barrier factor scaling in the amplitude numerator, they are still used for the mass-dependent width

        \return the ignore barrier amplitude scaling flag
    */
    Bool_t ignoreBarrierScaling() const { return ignoreBarrierScaling_; }

    //! Set the ignore barrier factor scaling flag
    /*!
        Whether to ignore barrier factor scaling in the amplitude numerator, they are still used for the mass-dependent width

        \param [in] boolean the ignore barrier factor scaling status
    */
    void ignoreBarrierScaling( const Bool_t boolean ) { ignoreBarrierScaling_ = boolean; }

    //! Allow the mass, width and spin of the resonance to be changed
    /*!
        Negative values wil be ignored, so if, for example, you
        want to only change the spin you can provide negative
        values for the mass and width

        \param [in] newMass new value of the resonance mass
        \param [in] newWidth new value of the resonance width
        \param [in] newSpin new value of the resonance spin
    */
    void changeResonance( const Double_t newMass, const Double_t newWidth, const Int_t newSpin );

    //! Allow the Blatt-Weisskopf radius for the resonance and parent factors to be changed
    /*!
        Negative values wil be ignored, so if, for example, you
        want to only change the parent radius you can provide a
        negative value for the resonance radius

        \param [in] resRadius new value of the resonance radius
        \param [in] parRadius new value of the parent radius
    */
    void changeBWBarrierRadii( const Double_t resRadius, const Double_t parRadius );

    //! Set value of the various parameters
    /*!
        \param [in] name the name of the parameter to be changed
        \param [in] value the new parameter value
    */
    virtual void setResonanceParameter( const TString& name, const Double_t value );

    //! Allow the various parameters to float in the fit
    /*!
        \param [in] name the name of the parameter to be floated
    */
    virtual void floatResonanceParameter( const TString& name );

    //! Access the given resonance parameter
    /*!
        \param [in] name the name of the parameter
        \return the corresponding parameter
    */
    virtual LauParameter* getResonanceParameter( const TString& name );

    //! Fix or release the resonance mass
    /*!
        \param [in] parFixed new status of mass
    */
    void fixMass( const Bool_t parFixed )
    {
        if ( resMass_ != 0 ) {
            resMass_->fixed( parFixed );
        }
    }

    //! Fix or release the resonance width
    /*!
        \param [in] parFixed new status of width
    */
    void fixWidth( const Bool_t parFixed )
    {
        if ( resWidth_ != 0 ) {
            resWidth_->fixed( parFixed );
        }
    }

    //! Get the status of resonance mass (fixed or released)
    /*!
        \return the status of resonance mass (fixed or released)
    */
    Bool_t fixMass() const { return ( resMass_ != 0 ) ? resMass_->fixed() : kTRUE; }

    //! Get the status of resonance width (fixed or released)
    /*!
        \return the status of resonance width (fixed or released)
    */
    Bool_t fixWidth() const { return ( resWidth_ != 0 ) ? resWidth_->fixed() : kTRUE; }

    //! Set the spin formalism to be used
    /*!
        \param [in] spinType the spin formalism
    */
    void setSpinType( const LauSpinType spinType ) { spinType_ = spinType; }

    //! Set the form factor model and parameters
    /*!
        \param [in] resFactor the barrier factor for the resonance decay
        \param [in] parFactor the barrier factor for the parent decay
    */
    void setBarrierRadii( LauBlattWeisskopfFactor* resFactor, LauBlattWeisskopfFactor* parFactor )
    {
        resBWFactor_ = resFactor;
        parBWFactor_ = parFactor;
    }

    //! Fix or release the Blatt-Weisskopf barrier radii
    void fixBarrierRadii( const Bool_t fixResRadius, const Bool_t fixParRadius );

    //! Get the status of resonance barrier radius (fixed or released)
    Bool_t fixResRadius() const;

    //! Get the status of parent barrier radius (fixed or released)
    Bool_t fixParRadius() const;

    //! Get the radius of the resonance barrier factor
    Double_t getResRadius() const;

    //! Get the radius of the parent barrier factor
    Double_t getParRadius() const;

  protected:
    //! Get the name of the parent particle
    TString getNameParent() const;
    //! Get the name of the first daughter of the resonance
    TString getNameDaug1() const;
    //! Get the name of the second daughter of the resonance
    TString getNameDaug2() const;
    //! Get the name of the daughter that does not originate form the resonance
    TString getNameBachelor() const;
    //! Get the parent particle mass
    Double_t getMassParent() const;
    //! Get the mass of daughter 1
    Double_t getMassDaug1() const;
    //! Get the mass of daughter 2
    Double_t getMassDaug2() const;
    //! Get the mass of the bachelor daughter
    Double_t getMassBachelor() const;
    //! Get the Charge of the parent particle
    Int_t getChargeParent() const;
    //! Get the charge of daughter 1
    Int_t getChargeDaug1() const;
    //! Get the charge of daughter 2
    Int_t getChargeDaug2() const;
    //! Get the charge of the bachelor daughter
    Int_t getChargeBachelor() const;

    //! Get the current value of the daughter momentum in the resonance rest frame
    Double_t getQ() const { return q_; }

    //! Get the current value of the bachelor momentum in the resonance rest frame
    Double_t getP() const { return p_; }

    //! Get the current value of the bachelor momentum in the parent rest frame
    Double_t getPstar() const { return pstar_; }

    //! Get the current value of the full spin-dependent covariant factor
    Double_t getCovFactor() const { return covFactor_; }

    //! Get the centrifugal barrier for the parent decay
    LauBlattWeisskopfFactor* getParBWFactor() { return parBWFactor_; }

    //! Get the centrifugal barrier for the parent decay
    const LauBlattWeisskopfFactor* getParBWFactor() const { return parBWFactor_; }

    //! Get the centrifugal barrier for the resonance decay
    LauBlattWeisskopfFactor* getResBWFactor() { return resBWFactor_; }

    //! Get the centrifugal barrier for the resonance decay
    const LauBlattWeisskopfFactor* getResBWFactor() const { return resBWFactor_; }

    //! Access the resonance info object
    LauResonanceInfo* getResInfo() const { return resInfo_; }

    //! Access the daughters object
    const LauDaughters* getDaughters() const { return daughters_; }

    //! Calculate the amplitude spin term using the Zemach tensor formalism
    /*!
        \param [in] pProd the momentum factor (either q * p or q * pstar)
    */
    Double_t calcZemachSpinFactor( const Double_t pProd ) const;

    //! Calculate the amplitude spin term using the covariant tensor formalism
    /*!
        \param [in] pProd the momentum factor (q * pstar)
    */
    Double_t calcCovSpinFactor( const Double_t pProd );

    //! Calculate the spin-dependent covariant factor
    /*!
        \param [in] erm E_ij in the parent rest-frame divided by m_ij (equivalent to sqrt(1 + p^2/mParent^2))
    */
    void calcCovFactor( const Double_t erm );

    //! Calculate the Legendre polynomial for the spin factor
    /*!
        Uses the current-event value of cosHel_
    */
    Double_t calcLegendrePoly() const;

    //! Calculate the Legendre polynomial for the spin factor (specifying the cosHel value)
    /*!
        \param [in] cosHel the cosine of the helicity angle
    */
    Double_t calcLegendrePoly( const Double_t cosHel );

    //! Complex resonant amplitude
    /*!
        \param [in] mass appropriate invariant mass for the resonance
        \param [in] spinTerm spin term
    */
    virtual LauComplex resAmp( Double_t mass, Double_t spinTerm ) = 0;

    //! Clear list of floating parameters
    void clearFloatingParameters() { resParameters_.clear(); }

    //! Add parameter to the list of floating parameters
    /*!
        \param [in] param the parameter to be added to the list
    */
    void addFloatingParameter( LauParameter* param );

    //! Access the list of floating parameters
    std::vector<LauParameter*>& getParameters() { return resParameters_; }

  private:
    //! Copy constructor (not implemented)
    LauAbsResonance( const LauAbsResonance& rhs );

    //! Copy assignment operator (not implemented)
    LauAbsResonance& operator=( const LauAbsResonance& rhs );

    //! Information on the resonance
    LauResonanceInfo* resInfo_ { 0 };

    //! Information on the particles
    const LauDaughters* daughters_;

    //! Parent name
    TString nameParent_ { "" };
    //! Daughter 1 name
    TString nameDaug1_ { "" };
    //! Daughter 2 name
    TString nameDaug2_ { "" };
    //! Bachelor name
    TString nameBachelor_ { "" };

    //! Parent charge
    Int_t chargeParent_ { 0 };
    //! Daughter 1 charge
    Int_t chargeDaug1_ { 0 };
    //! Daughter 2 charge
    Int_t chargeDaug2_ { 0 };
    //! Bachelor charge
    Int_t chargeBachelor_ { 0 };

    //! Parent mass
    Double_t massParent_ { 0.0 };
    //! Daughter 1 mass
    Double_t massDaug1_ { 0.0 };
    //! Daughter 2 mass
    Double_t massDaug2_ { 0.0 };
    //! Bachelor mass
    Double_t massBachelor_ { 0.0 };

    //! Resonance name
    TString resName_;

    //! Resonance name with illegal characters removed
    TString sanitisedName_;

    //! Resonance mass
    LauParameter* resMass_ { 0 };
    //! Resonance width
    LauParameter* resWidth_ { 0 };

    //! All parameters of the resonance
    std::vector<LauParameter*> resParameters_;

    //! Resonance spin
    Int_t resSpin_;
    //! Resonance charge
    Int_t resCharge_ { 0 };
    //! DP axis identifier
    Int_t resPairAmpInt_;
    //! Blatt Weisskopf barrier for parent decay
    LauBlattWeisskopfFactor* parBWFactor_ { 0 };
    //! Blatt Weisskopf barrier for resonance decay
    LauBlattWeisskopfFactor* resBWFactor_ { 0 };

    //! Spin formalism
    LauSpinType spinType_ { Zemach_P };

    //! Boolean to flip helicity
    Bool_t flipHelicity_ { kFALSE };
    //! Boolean to ignore the momentum factors in both the spin factor and the mass-dependent width
    Bool_t ignoreMomenta_ { kFALSE };
    //! Boolean to set the spinTerm to unity always
    Bool_t ignoreSpin_ { kFALSE };
    //! Boolean to ignore barrier factor scaling in the amplitude numerator, they are still used for the mass-dependent width
    Bool_t ignoreBarrierScaling_ { kFALSE };

    // Event kinematics information

    //! Invariant mass
    Double_t mass_ { 0.0 };
    //! Helicity angle cosine
    Double_t cosHel_ { 0.0 };

    //! Daughter momentum in resonance rest frame
    Double_t q_ { 0.0 };
    //! Bachelor momentum in resonance rest frame
    Double_t p_ { 0.0 };
    //! Bachelor momentum in parent rest frame
    Double_t pstar_ { 0.0 };

    //! Covariant factor
    /*!
        sqrt(1 + z*z), where z = p / mParent

        Can also be expressed as E_ij in the parent rest-frame divided by m_ij - indeed this is how LauKinematics calculates it.

        \see LauKinematics::getcov12
        \see LauKinematics::getcov13
        \see LauKinematics::getcov23
    */
    Double_t erm_ { 1.0 };

    //! Covariant factor (full spin-dependent expression)
    Double_t covFactor_ { 1.0 };

    ClassDef( LauAbsResonance, 0 )    // Abstract resonance class
};

#endif
