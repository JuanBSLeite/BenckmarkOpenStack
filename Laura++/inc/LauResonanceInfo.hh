
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

/*! \file LauResonanceInfo.hh
    \brief File containing declaration of LauResonanceInfo class.
*/

/*! \class LauResonanceInfo
    \brief Class for defining the properties of a resonant particle.
*/

#ifndef LAU_RESONANCE_INFO
#define LAU_RESONANCE_INFO

#include "LauBlattWeisskopfFactor.hh"

#include "TString.h"

#include <iosfwd>
#include <set>

class LauParameter;

class LauResonanceInfo {

  public:
    //! Constructor
    /*!
        \param [in] name the name of the resonant particle
        \param [in] mass the mass of the resonant particle
        \param [in] width the width of the resonant particle
        \param [in] spin the spin of the resonant particle
        \param [in] charge the charge of the resonant particle
        \param [in] bwCategory the Blatt-Weisskopf barrier factor category
        \param [in] bwRadius the Blatt-Weisskopf radius of the resonant particle
    */
    LauResonanceInfo( const TString& name,
                      const Double_t mass,
                      const Double_t width,
                      const Int_t spin,
                      const Int_t charge,
                      const LauBlattWeisskopfFactor::BlattWeisskopfCategory bwCategory,
                      const Double_t bwRadius = 4.0 );

    //! Destructor
    virtual ~LauResonanceInfo();

    //! Retrieve the name of the resonant particle
    /*!
        \return the name of the resonant particle
    */
    TString getName() const { return name_; }

    //! Retrieve the sanitised name of the resonant particle
    /*!
        Removes/replaces characters from the name that cause
        problems when used in TBranch names

        \return the sanitised name of the resonant particle
    */
    TString getSanitisedName() const { return sanitisedName_; }

    //! Retrieve the mass of the resonant particle
    /*!
        \return the mass of the resonant particle
    */
    LauParameter* getMass() const { return mass_; }

    //! Retrieve the width of the resonant particle
    /*!
        \return the width of the resonant particle
    */
    LauParameter* getWidth() const { return width_; }

    //! Retrieve the spin of the resonant particle
    /*!
        \return the spin of the resonant particle
    */
    UInt_t getSpin() const { return spin_; }

    //! Retrieve the charge of the resonant particle
    /*!
        \return the charge of the resonant particle
    */
    Int_t getCharge() const { return charge_; }

    //! Retrieve the BW category of the resonant particle
    /*!
        \return the BW category of the resonant particle
    */
    LauBlattWeisskopfFactor::BlattWeisskopfCategory getBWCategory() const { return bwCategory_; }

    //! Retrieve the BW radius of the resonant particle
    /*!
        \return the BW radius of the resonant particle
    */
    Double_t getBWRadius() const { return bwRadius_; }

    //! Create the charge conjugate particle info record
    /*!
        The mass and width parameters are cloned
    */
    LauResonanceInfo* createChargeConjugate();

    //! Create another record that will share parameters with this one
    /*!
        The name needs to be specified.
        The spin and charge are assumed to be the same.
        The mass, width and other parameters will be cloned.

        \param [in] name the name of the resonant particle
    */
    LauResonanceInfo* createSharedParameterRecord( const TString& name );

    //! Retrieve an extra parameter of the resonance
    /*!
        \return the extra parameter (or null pointer if not found)
    */
    LauParameter* getExtraParameter( const TString& parName );

    //! Add an extra parameter of the resonance
    /*!
        \param [in] param the extra parameter to be added
        \param [in] independentPar governs whether any info record that usually shares parameters with this one should also share this one (the default) or make its own independent version
    */
    void addExtraParameter( LauParameter* param, const Bool_t independentPar = kFALSE );

  protected:
    //! Add a clone of an extra parameter of the resonance
    /*!
        \param [in] param the extra parameter to be added
        \param [in] copyNotClone should we create an unlinked copy instead of cloning - default is to clone
    */
    void addCloneOfExtraParameter( LauParameter* param, const Bool_t copyNotClone = kFALSE );

  private:
    //! Copy constructor (not implemented)
    LauResonanceInfo( const LauResonanceInfo& other );

    //! Copy constructor (with new name and charge)
    LauResonanceInfo( const LauResonanceInfo& other, const TString& newName, const Int_t newCharge );

    //! Copy assignment operator (not implemented)
    LauResonanceInfo& operator=( const LauResonanceInfo& other );

    //! Create the sanitised name by removing characters that are illegal in TBranch names
    void sanitiseName();

    //! The name of the resonant particle
    TString name_;

    //! The name of the resonant particle with illegal characters removed
    TString sanitisedName_;

    //! The mass of the resonant particle
    LauParameter* mass_;

    //! The width of the resonant particle
    LauParameter* width_;

    //! The spin of the resonant particle
    UInt_t spin_;

    //! The charge of the resonant particle
    Int_t charge_;

    //! The Blatt-Weisskopf barrier factor category
    LauBlattWeisskopfFactor::BlattWeisskopfCategory bwCategory_;

    //! The Blatt-Weisskopf radius of the resonant particle
    Double_t bwRadius_;

    //! The conjugate info object
    LauResonanceInfo* conjugate_;

    //! Other info objects that share parameters with this one
    std::vector<LauResonanceInfo*> sharedParRecords_;

    //! Extra parameters
    std::set<LauParameter*> extraPars_;

    ClassDef( LauResonanceInfo, 0 )    // Specify each allowed resonance
};

//! output operator formatting of an info record
std::ostream& operator<<( std::ostream& stream, const LauResonanceInfo& infoRecord );

#endif
