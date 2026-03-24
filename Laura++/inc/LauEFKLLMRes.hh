
/*
Copyright 2015 University of Warwick

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

/*! \file LauEFKLLMRes.hh
    \brief File containing declaration of LauEFKLLMRes class.
*/

/*! \class LauEFKLLMRes
    \brief Class for defining the EFKLLM K-pi S-wave model

    Class for defining the EFKLLM form-factor model for the K-pi S-wave.
    The model consists of a tabulated form-factor, which is interpolated using cubic splines (one for magnitude values and one for phase values), multiplied by a mass-dependence (e.g. constant, 1/m^2, etc.).
    The massFactor resonance parameter is the power of the mass dependence - defaults to zero, i.e. constant.
    For more details see B. El-Bennich et al. Phys. Rev. D 79, 094005 (2009), arXiv:0902.3645 [hep-ph].
    (The acronym EFKLLM is constructed from the surnames of the authors of the above paper.)
*/

#ifndef LAU_EFKLLM_RES
#define LAU_EFKLLM_RES

#include "LauAbsResonance.hh"
#include "LauComplex.hh"

#include "TString.h"

class Lau1DCubicSpline;

class LauEFKLLMRes : public LauAbsResonance {

  public:
    //! Constructor
    /*!
        \param [in] resInfo the object containing information on the resonance name, mass, width, spin, charge, etc.
        \param [in] resPairAmpInt the number of the daughter not produced by the resonance
        \param [in] daughters the daughter particles
    */
    LauEFKLLMRes( LauResonanceInfo* resInfo, const Int_t resPairAmpInt, const LauDaughters* daughters );

    //! Destructor
    virtual ~LauEFKLLMRes();

    //! Read the form factor information from text file
    /*!
        Creates the splines from the tabulated form factor data.
        These are shared between all instances of this class.

        \param [in] inputFile the name of the file to be read
    */
    static void setupFormFactor( const TString& inputFile );

    //! Initialise the model
    virtual void initialise();

    //! Get the resonance model type
    /*!
        \return the resonance model type
    */
    virtual LauAbsResonance::LauResonanceModel getResonanceModel() const
    {
        return LauAbsResonance::EFKLLM;
    }

    //! Set value of a resonance parameter
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

    //! Retrieve the resonance parameters, e.g. so that they can be loaded into a fit
    /*!
        \return floating parameters of the resonance
    */
    virtual const std::vector<LauParameter*>& getFloatingParameters();

  protected:
    //! Set the power of the mass dependence
    /*!
        \param [in] massFactor the new power of the mass dependence
    */
    void setMassFactor( const Double_t massFactor );

    //! Get the power of the mass dependence
    /*!
        \return the power of the mass dependence
    */
    Double_t getMassFactor() const
    {
        return ( massFactor_ != 0 ) ? massFactor_->unblindValue() : 0.0;
    }

    //! See if the mass factor parameter is fixed or floating
    /*!
        \return kTRUE if the mass factor parameter is fixed, kFALSE otherwise
    */
    Bool_t fixMassFactor() const { return ( massFactor_ != 0 ) ? massFactor_->fixed() : kTRUE; }

    //! Complex resonant amplitude
    /*!
        \param [in] mass appropriate invariant mass for the resonance
        \param [in] spinTerm spin term
    */
    virtual LauComplex resAmp( Double_t mass, Double_t spinTerm );

  private:
    //! Spline describing the magnitude variation of the form-factor
    static Lau1DCubicSpline* magSpline_;
    //! Spline describing the phase variation of the form-factor
    static Lau1DCubicSpline* phaseSpline_;

    //! The power of the mass dependence
    LauParameter* massFactor_;

    ClassDef( LauEFKLLMRes, 0 )    // EFKLLM resonance model
};

#endif
