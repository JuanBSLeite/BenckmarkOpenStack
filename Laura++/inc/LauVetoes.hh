
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

/*! \file LauVetoes.hh
    \brief File containing declaration of LauVetoes class.
*/

/*! \class LauVetoes
    \brief Class for defining vetoes within the Dalitz plot.

    Each veto is defined by an index corresponding to one of the three invariant mass-squared variables, and maximum and minimum values for the veto's range in the specified variable.

    Each mass-squared variable is numbered according to the index of the daughter not involved in the vetoed mass-squared variable. For example a veto in m12 squared would receive the index 3.

    Since v3r2, in the case of a symmetric DP, vetoes in m13 or m23 will automatically be symmetrised to the other variable (so only one needs to be specified).
    However, if one wants to veto a region from a symmetric DP in the already-symmetrised co-ordinates, the special indices 4 and 5 can be used that will apply the veto to mMin or mMax, respectively.
*/

#ifndef LAU_VETOES
#define LAU_VETOES

#include "Rtypes.h"

#include <vector>

class LauKinematics;

class LauVetoes {

  public:
    //! Constructor
    LauVetoes();

    //! Destructor
    virtual ~LauVetoes();

    //! Copy constructor
    /*!
        \param [in] other the object to be copied
    */
    LauVetoes( const LauVetoes& other );

    //! Copy assignment operator
    /*!
        \param [in] other the object to be copied
        \return the assigned vetoes object
    */
    LauVetoes& operator=( const LauVetoes& other );

    //! Add a veto to the Dalitz plot
    /*!
        \param [in] resPairAmpInt the index of the mass-squared variable to be vetoed (in the case of symmetric DPs the special indices 4 and 5 can be used to indicate that a veto should be applied to mMin or mMax)
        \param [in] minMass the minimum mass of the veto
        \param [in] maxMass the maximum mass of the veto
    */
    void addMassVeto( const Int_t resPairAmpInt, const Double_t minMass, const Double_t maxMass );

    //! Add a veto to the Dalitz plot
    /*!
        \param [in] resPairAmpInt the index of the mass-squared variable to be vetoed (in the case of symmetric DPs the special indices 4 and 5 can be used to indicate that a veto should be applied to mMinSq or mMaxSq)
        \param [in] minMassSq the minimum mass-squared of the veto
        \param [in] maxMassSq the maximum mass-squared of the veto
    */
    void addMassSqVeto( const Int_t resPairAmpInt, const Double_t minMassSq, const Double_t maxMassSq );

    //! Check whether the specified Dalitz plot point passes the vetoes
    /*!
        \param [in] kinematics a point in the Dalitz plot
        \return true if the specified Dalitz plot point is outside all veto regions, false otherwise
    */
    Bool_t passVeto( const LauKinematics* kinematics ) const;

  protected:
    //! Check whether the specified Dalitz plot point passes the vetoes
    /*!
        \param [in] m12Sq the mass-squared of the first and second daughters
        \param [in] m23Sq the mass-squared of the second and third daughters
        \param [in] m13Sq the mass-squared of the first and third daughters
        \param [in] symmetricDP is the DP symmetric
        \param [in] fullySymmetricDP is the DP fully symmetric
        \return true if the specified Dalitz plot point is outside all veto regions, false otherwise
    */
    Bool_t passVeto( const Double_t m12Sq,
                     const Double_t m23Sq,
                     const Double_t m13Sq,
                     const Bool_t symmetricDP,
                     const Bool_t fullySymmetricDP ) const;

    //! Retrieve the number of vetoes
    /*!
        \return the number of vetoes
    */
    UInt_t getNVetoes() const { return nVetoes_; }

    //! Retrieve the index of the vetoed mass-squared variable for each veto
    /*!
        \return the index of the vetoed mass-squared variable for each veto
    */
    const std::vector<Int_t>& getVetoPairs() const { return vetoPair_; }

    //! Retrieve the minimum mass-squared for each veto
    /*!
        \return the minimum mass-squared for each veto
    */
    const std::vector<Double_t>& getVetoMinMass() const { return vetoMinMass_; }

    //! Retrieve the maximum mass-squared for each veto
    /*!
        \return the maximum mass-squared for each veto
    */
    const std::vector<Double_t>& getVetoMaxMass() const { return vetoMaxMass_; }

  private:
    //! The number of vetoes
    UInt_t nVetoes_;

    //! The index of the vetoed mass-squared variable for each veto
    std::vector<Int_t> vetoPair_;

    //! The minimum mass-squared for each veto
    std::vector<Double_t> vetoMinMass_;

    //! The maximum mass-squared for each veto
    std::vector<Double_t> vetoMaxMass_;

    ClassDef( LauVetoes, 0 )    // Vetoes in the Dalitz plot
};

#endif
