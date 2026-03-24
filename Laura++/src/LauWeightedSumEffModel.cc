
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

/*! \file LauWeightedSumEffModel.cc
    \brief File containing implementation of LauWeightedSumEffModel class.
*/

#include "LauWeightedSumEffModel.hh"

#include "Lau2DHistDP.hh"
#include "Lau2DSplineDP.hh"
#include "LauDaughters.hh"
#include "LauKinematics.hh"
#include "LauVetoes.hh"

#include "TSystem.h"

#include <cstdlib>
#include <iostream>

LauWeightedSumEffModel::LauWeightedSumEffModel( const LauDaughters* daughters ) :
    daughters_( daughters ),
    effModel_( 0 ),
    fluctuateEffHisto_( kFALSE ),
    lowBinWarningIssued_( kFALSE ),
    highBinWarningIssued_( kFALSE )
{
    if ( daughters_ == 0 ) {
        std::cerr << "ERROR in LauWeightedSumEffModel Constructor : invalid pointer to daughters object supplied."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }
}

void LauWeightedSumEffModel::addEffModel( const LauAbsEffModel* effModel, Double_t coeff )
{
    const LauDaughters* otherDaughters = effModel->getDaughters();
    if ( otherDaughters->getTypeDaug1() != daughters_->getTypeDaug1() ||
         otherDaughters->getTypeDaug2() != daughters_->getTypeDaug2() ||
         otherDaughters->getTypeDaug3() != daughters_->getTypeDaug3() ||
         otherDaughters->getTypeParent() != daughters_->getTypeParent() ) {
        std::cerr << "ERROR in LauWeightedSumEffModel::addEffModel : daughters of provided efficiency model do not match those expected."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    effModel_.push_back( effModel );
    coeff_.push_back( coeff );

    if ( effModel->fluctuateEffHisto() ) {
        fluctuateEffHisto_ = kTRUE;
    }
}

Double_t LauWeightedSumEffModel::calcEfficiency( const LauKinematics* kinematics ) const
{
    // Routine to calculate the efficiency for the given event/point in
    // the Dalitz plot. This routine uses the models set by the
    // addEffModel() function.
    Double_t eff( 0.0 );

    std::vector<const LauAbsEffModel*>::const_iterator it  = effModel_.begin();
    std::vector<const LauAbsEffModel*>::const_iterator end = effModel_.end();

    std::vector<Double_t>::const_iterator coeffIt = coeff_.begin();

    for ( ; it != end; ++it ) {
        eff += ( *coeffIt ) * ( *it )->calcEfficiency( kinematics );
        ++coeffIt;
    }

    // Check that the efficiency is in the allowed range (0-1)
    // Out of range efficiencies could be caused by incorrect coefficients.
    if ( eff < 0.0 ) {
        if ( ! lowBinWarningIssued_ ) {
            std::cerr << "WARNING in LauWeightedSumEffModel::calcEfficiency : Efficiency " << eff
                      << " is less than 0 - setting to 0.  You may want to check your coefficients!"
                      << std::endl
                      << "                                       : Further warnings will be suppressed."
                      << std::endl;
            lowBinWarningIssued_ = kTRUE;
        }
        eff = 0.0;
    } else if ( eff > 1.0 ) {
        if ( ! highBinWarningIssued_ ) {
            std::cerr << "WARNING in LauWeightedSumEffModel::calcEfficiency : Efficiency " << eff
                      << " is greater than 1 - setting to 1.  You may want to check your coefficients!"
                      << std::endl
                      << "                                       : Further warnings will be suppressed."
                      << std::endl;
            highBinWarningIssued_ = kTRUE;
        }
        eff = 1.0;
    }

    return eff;
}

Bool_t LauWeightedSumEffModel::passVeto( const LauKinematics* kinematics ) const
{
    std::vector<const LauAbsEffModel*>::const_iterator it  = effModel_.begin();
    std::vector<const LauAbsEffModel*>::const_iterator end = effModel_.end();

    for ( ; it != end; ++it ) {
        if ( ! ( *it )->passVeto( kinematics ) )
            return kFALSE;
    }
    return kTRUE;
}
