
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

/*! \file LauWeightedSumEffModel.hh
    \brief File containing declaration of LauWeightedSumEffModel class.
*/

/*! \class LauWeightedSumEffModel
    \brief Class that implements the efficiency description across the signal Dalitz plot.

    Class that defines the efficiency model variation across the signal Dalitz plot.
    The phase space acceptance is determined from a weighted sum of LauAbsEffModel objects.
    This is useful for samples that contain two different categories of signal with different efficiencies.
    The efficiency variation is defined in terms of x = m_13^2, y = m_23^2 for the Dalitz plot (default)
    or x = m', y = theta' for the square Dalitz plot
*/

#ifndef LAUWEIGHTEDSUMEFFMODEL
#define LAUWEIGHTEDSUMEFFMODEL

#include "LauAbsEffModel.hh"

class LauDaughters;
class LauKinematics;

class LauWeightedSumEffModel : public LauAbsEffModel {

  public:
    //! Constructor
    /*!
        \param [in] daughters the daughters particles of the Dalitz plot model
    */
    LauWeightedSumEffModel( const LauDaughters* daughters );

    //! Destructor
    virtual ~LauWeightedSumEffModel() {}

    //! Add an efficiency variation across the phase space using a predetermined LauAbsEffModel object.
    /*!
        The efficiency is defined in terms of x = m_13^2, y = m_23^2 or x = m', y = theta' for the square Dalitz plot

        \param [in] effModel the LauAbsEffModel object that describes the efficiency variation
        \param [in] coeff the coefficient to multiply this efficiency by
    */
    void addEffModel( const LauAbsEffModel* effModel, Double_t coeff );

    //! Determine the efficiency for a given point in the Dalitz plot.
    /*!
        The method uses the models set by the addEffModel() function and the vetoes information.

        \param [in] kinematics the object that defines the DP position
        \return the efficiency value at the given point in the DP
    */
    Double_t calcEfficiency( const LauKinematics* kinematics ) const;

    //! Determine whether the given DP position is outside the vetoes
    /*!
        \param [in] kinematics the object that defines the DP position
        \return kTRUE if the DP position is outside all veto regions, kFALSE otherwise
    */
    Bool_t passVeto( const LauKinematics* kinematics ) const;

    //! Determine whether the efficiency histogram has had its bins fluctuated within their errors
    Bool_t fluctuateEffHisto() const { return fluctuateEffHisto_; }

    //! Return the daughters object
    /*
            \return the LauDaughters object associated with the DP
        */
    const LauDaughters* getDaughters() const { return daughters_; }

  private:
    //! Copy constructor - not implemented
    LauWeightedSumEffModel( const LauWeightedSumEffModel& rhs );

    //! Copy assignment operator - not implemented
    LauWeightedSumEffModel& operator=( const LauWeightedSumEffModel& rhs );

    //! The daughters object
    const LauDaughters* daughters_;

    //! The efficiency model objects
    std::vector<const LauAbsEffModel*> effModel_;

    //! The efficiency model objects
    std::vector<Double_t> coeff_;

    //! Fluctuate histogram within the error
    Bool_t fluctuateEffHisto_;

    //! Flag to track whether a warning has been issued for bin values less than zero
    mutable Bool_t lowBinWarningIssued_;

    //! Flag to track whether a warning has been issued for bin values greater than one
    mutable Bool_t highBinWarningIssued_;

    ClassDef( LauWeightedSumEffModel, 0 )    // Implement the signal efficiency across the DP
};

#endif
