
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

/*! \file LauResonanceMaker.hh
    \brief File containing declaration of LauResonanceMaker class.
*/

/*! \class LauResonanceMaker
    \brief Singleton factory class for creating resonances.

    Singleton factory class for creating resonances. Information records for all known resonances are stored within this class.
*/

#ifndef LAU_RESONANCE_MAKER
#define LAU_RESONANCE_MAKER

#include "LauAbsResonance.hh"

#include "TString.h"

#include <iosfwd>
#include <vector>

class LauDaughters;
class LauResonanceInfo;

class LauResonanceMaker {

  public:
    //! Get the factory instance
    static LauResonanceMaker& get();

    //! Set the type of BW factor (for all categories)
    /*!
        This must be used before creating any resonances

        \param [in] bwType the Blatt-Weisskopf barrier type
    */
    void setBWType( const LauBlattWeisskopfFactor::BarrierType bwType );

    //! Set the rest frame in which the bachelor momentum should be calculated (for all BW categories)
    /*!
        This must be used before creating any resonances

        \param [in] restFrame the rest frame in which the bachelor momentum should be calculated for the Blatt-Weisskopf factors
    */
    void setBWBachelorRestFrame( const LauBlattWeisskopfFactor::RestFrame restFrame );

    //! Set the spin formalism to be used for all resonances
    /*!
        This must be used before creating any resonances

        \param [in] spinType the spin formalism to be used
    */
    void setSpinFormalism( const LauAbsResonance::LauSpinType spinType );

    //! Set the BW radius for the given category
    /*!
        \param [in] bwCategory the Blatt-Weisskopf barrier factor category
        \param [in] bwRadius the radius value to be used for the given category
    */
    void setDefaultBWRadius( const LauBlattWeisskopfFactor::BlattWeisskopfCategory bwCategory,
                             const Double_t bwRadius );

    //! Fix or release the Blatt-Weisskopf barrier radius for the given category
    /*!
        \param [in] bwCategory the Blatt-Weisskopf barrier factor category
        \param [in] fixRadius new status of the radius (kTRUE = fixed, kFALSE = floating)
    */
    void fixBWRadius( const LauBlattWeisskopfFactor::BlattWeisskopfCategory bwCategory,
                      const Bool_t fixRadius );

    //! Create a resonance
    /*!
        \param [in] daughters defines the Dalitz plot in which the resonance should be created
        \param [in] resName the name of the resonant particle
        \param [in] resPairAmpInt the index of the daughter not produced by the resonance
        \param [in] resType the type of the resonance
        \param [in] bwCategory the Blatt-Weisskopf barrier factor category
        \return the resonance
    */
    LauAbsResonance* getResonance( const LauDaughters* daughters,
                                   const TString& resName,
                                   const Int_t resPairAmpInt,
                                   const LauAbsResonance::LauResonanceModel resType,
                                   const LauBlattWeisskopfFactor::BlattWeisskopfCategory bwCategory =
                                       LauBlattWeisskopfFactor::Default );

    //! Retrieve the integer index for the specified resonance
    /*!
        \param [in] name the name of the resonant particle
        \return the index
    */
    Int_t resTypeInt( const TString& name ) const;

    //! Retrieve the number of defined resonances
    /*!
        \return the number of defined resonances
    */
    UInt_t getNResDefMax() const { return nResDefMax_; }

    //! Print the information records, one per line, to the requested stream
    /*!
        \param [in,out] stream the stream to which to print the info
    */
    void printAll( std::ostream& stream ) const;

    //! Get the information for the given resonance name
    /*!
        \param [in] resName the name of the resonant particle
        \return the LauResonanceInfo pointer if we can find the resonance name
    */
    LauResonanceInfo* getResInfo( const TString& resName ) const;

    //! Retrieve parent Blatt-Weisskopf factor (for use by K-matrix pole/SVP which doesn't have a `resInfo')
    LauBlattWeisskopfFactor* getParentBWFactor( Int_t newSpin,
                                                LauBlattWeisskopfFactor::BarrierType barrierType );

  protected:
    //! Create the list of known resonances
    void createResonanceVector();

    //! Retrieve Blatt-Weisskopf factor for the given category
    LauBlattWeisskopfFactor* getBWFactor( const LauBlattWeisskopfFactor::BlattWeisskopfCategory bwCategory,
                                          const LauResonanceInfo* resInfo );

  private:
    /*!
        \struct BlattWeisskopfCategoryInfo
        \brief Data structure to store information on a given Blatt-Weisskopf category
    */
    struct BlattWeisskopfCategoryInfo {
        //! The BW factor object
        LauBlattWeisskopfFactor* bwFactor_;
        //! The default value for the radius in this category
        Double_t defaultRadius_;
        //! Whether or not the radius value for this category should be fixed in the fit
        Bool_t radiusFixed_;
    };

    //! Define a type to hold information on each BW category
    typedef std::map<LauBlattWeisskopfFactor::BlattWeisskopfCategory, BlattWeisskopfCategoryInfo>
        BWFactorCategoryMap;

    //! Constructor
    LauResonanceMaker();

    //! Destructor
    virtual ~LauResonanceMaker();

    //! Copy constructor (not inplemented)
    LauResonanceMaker( const LauResonanceMaker& other );

    //! Copy assignment (not implemented)
    LauResonanceMaker& operator=( const LauResonanceMaker& other );

    //! The singleton instance
    static LauResonanceMaker* resonanceMaker_;

    //! The number of known resonances
    UInt_t nResDefMax_;

    //! The known resonances
    std::vector<LauResonanceInfo*> resInfo_;

    //! The type of the Blatt-Weisskopf barrier to use for all resonances
    LauBlattWeisskopfFactor::BarrierType bwBarrierType_;

    //! The rest frame in which the bachelor momentum used in the Blatt-Weisskopf factors should be calculated
    LauBlattWeisskopfFactor::RestFrame bwRestFrame_;

    //! The spin formalism that should be used for all resonances
    LauAbsResonance::LauSpinType spinFormalism_;

    //! The Blatt-Weisskopf factor objects (and related information) for each category
    BWFactorCategoryMap bwFactors_;

    //! The Blatt-Weisskopf factor objects for resonances in the independent category
    std::vector<LauBlattWeisskopfFactor*> bwIndepFactors_;

    //! Boolean flag to control printing a summary of the formalism to be used when the first resonance is created
    Bool_t summaryPrinted_;

    ClassDef( LauResonanceMaker, 0 )    // Kinematic routines
};

#endif
