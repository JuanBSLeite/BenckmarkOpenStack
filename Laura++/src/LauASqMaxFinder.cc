
/*
Copyright 2022 University of Warwick

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

/*! \file LauASqMaxFinder.cc
    \brief File containing implementation of LauASqMaxFinder class.
*/

#include "LauASqMaxFinder.hh"

#include "LauAbsResonance.hh"
#include "LauFitter.hh"
#include "LauIsobarDynamics.hh"
#include "LauKinematics.hh"
#include "LauParameter.hh"

#include "TSystem.h"

#include <array>
#include <iostream>
#include <memory>

LauASqMaxFinder::LauASqMaxFinder( LauIsobarDynamics& iso ) :
    iso_ { iso }
{
    const LauKinematics* kinematics { iso_.getKinematics() };
    const Double_t m13SqMin { kinematics->getm13SqMin() };
    const Double_t m13SqMax { kinematics->getm13SqMax() };
    const Double_t m23SqMin { kinematics->getm23SqMin() };
    const Double_t m23SqMax { kinematics->getm23SqMax() };

    ownedParams_.resize( 2 );
    ownedParams_[0] = std::make_unique<LauParameter>( "m13Sq",
                                                      0.5 * ( m13SqMin + m13SqMax ),
                                                      m13SqMin,
                                                      m13SqMax,
                                                      kFALSE );
    ownedParams_[1] = std::make_unique<LauParameter>( "m23Sq",
                                                      0.5 * ( m23SqMin + m23SqMax ),
                                                      m23SqMin,
                                                      m23SqMax,
                                                      kFALSE );

    params_.resize( 2 );
    params_[0] = ownedParams_[0].get();
    params_[1] = ownedParams_[1].get();
}

void LauASqMaxFinder::setParsFromMinuit( Double_t* par, Int_t npar )
{
    if ( static_cast<std::size_t>( npar ) != params_.size() ) {
        std::cerr << "ERROR in LauASqMaxFinder::setParsFromMinuit : Unexpected number of free parameters: "
                  << npar << ".\n";
        std::cerr << "                                              Expected: " << params_.size()
                  << ".\n"
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    for ( Int_t i { 0 }; i < npar; ++i ) {
        params_[i]->value( par[i] );
    }

    return;
}

Double_t LauASqMaxFinder::getTotNegLogLikelihood()
{
    const Double_t m13Sq { params_[0]->unblindValue() };
    const Double_t m23Sq { params_[1]->unblindValue() };

    if ( ! iso_.getKinematics()->withinDPLimits( m13Sq, m23Sq ) ) {
        if ( printMinimisationWarnings_ ) {
            std::cerr << "WARNING in LauASqMaxFinder::getTotNegLogLikelihood : Fit has wandered outside the phase space boundary\n";
            std::cerr << "                                                   : Returning 0.0 value (further warnings of this type will be suppressed)"
                      << std::endl;
            printMinimisationWarnings_ = kFALSE;
        }
        return 0.0;
    }

    iso_.calcLikelihoodInfo( m13Sq, m23Sq );
    return -iso_.getEvtIntensity();
}

Double_t LauASqMaxFinder::find()
{
    Double_t ASqMax { -1.0 };

    const LauKinematics* kinematics { iso_.getKinematics() };
    LauKinematics localKin { kinematics->getm1(),
                             kinematics->getm2(),
                             kinematics->getm3(),
                             kinematics->getmParent() };

    const Double_t m13SqMin { localKin.getm13Min() };
    const Double_t m13SqMax { localKin.getm13Max() };
    const Double_t m23SqMin { localKin.getm23Min() };
    const Double_t m23SqMax { localKin.getm23Max() };

    LauFitter::setFitterMaxPars( 2 );
    LauFitter::setFitterVerbosity( LauOutputLevel::None );

    for ( UInt_t iRes { 0 }; iRes < iso_.getnTotAmp(); ++iRes ) {
        // Find the center of the resonance and set that as the start of the search
        const LauAbsResonance* res {
            const_cast<const LauIsobarDynamics&>( iso_ ).getResonance( iRes ) };
        const Int_t bachelor { res->getPairInt() };
        const Double_t mass { res->getMass() };
        const Double_t width { res->getWidth() };

        const std::array<Double_t, 3> massVals { mass - width, mass, mass + width };
        const std::array<Double_t, 3> cosHelVals { -1.0, 0.0, 1.0 };
        for ( Double_t massVal : massVals ) {
            for ( Double_t cosHel : cosHelVals ) {

                switch ( bachelor ) {
                    case 0 :
                        localKin.updateKinematics( 0.5 * ( m13SqMax + m13SqMin ),
                                                   0.5 * ( m23SqMax + m23SqMin ) );
                        break;
                    case 1 :
                        localKin.updateKinematicsFrom23( massVal, cosHel );
                        break;
                    case 2 :
                        localKin.updateKinematicsFrom13( massVal, cosHel );
                        break;
                    case 3 :
                        localKin.updateKinematicsFrom12( massVal, cosHel );
                        break;

                    default :
                        break;
                }

                params_[0]->initValue( localKin.getm13Sq() );
                params_[1]->initValue( localKin.getm23Sq() );

                params_[0]->error( 0.1 );
                params_[1]->error( 0.1 );

                LauFitter::fitter().initialise( this, params_ );

                const LauAbsFitter::FitStatus fitResult { LauFitter::fitter().minimise() };

                if ( fitResult.status != 0 and -fitResult.NLL > ASqMax ) {
                    ASqMax = -fitResult.NLL;
                }
            }
        }
    }

    LauFitter::destroyFitter();

    if ( ASqMax < 0.0 ) {
        std::cerr << "WARNING in LauASqMaxFinder::find : Could not find ASqMax\n";
        std::cerr << "                                 : Returning -1.0" << std::endl;
    }

    return ASqMax;
}
