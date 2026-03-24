
/*
Copyright 2005 University of Warwick

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

/*! \file LauIsobarDynamics.cc
    \brief File containing implementation of LauIsobarDynamics class.
*/

#include "LauIsobarDynamics.hh"

#include "LauASqMaxFinder.hh"
#include "LauAbsEffModel.hh"
#include "LauAbsIncohRes.hh"
#include "LauAbsResonance.hh"
#include "LauBelleNR.hh"
#include "LauBelleSymNR.hh"
#include "LauCacheData.hh"
#include "LauConstants.hh"
#include "LauDPPartialIntegralInfo.hh"
#include "LauDaughters.hh"
#include "LauFitDataTree.hh"
#include "LauKMatrixProdPole.hh"
#include "LauKMatrixProdSVP.hh"
#include "LauKMatrixPropFactory.hh"
#include "LauKMatrixPropagator.hh"
#include "LauKinematics.hh"
#include "LauNRAmplitude.hh"
#include "LauPrint.hh"
#include "LauRandom.hh"
#include "LauResonanceInfo.hh"
#include "LauResonanceMaker.hh"
#include "LauRhoOmegaMix.hh"

#include "TFile.h"
#include "TMath.h"
#include "TRandom.h"
#include "TSystem.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <vector>

// for Kpipi: only one scfFraction 2D histogram is needed
LauIsobarDynamics::LauIsobarDynamics( LauDaughters* daughters,
                                      LauAbsEffModel* effModel,
                                      LauAbsEffModel* scfFractionModel ) :
    daughters_( daughters ),
    kinematics_( daughters_ ? daughters_->getKinematics() : 0 ),
    effModel_( effModel ),
    nAmp_( 0 ),
    nIncohAmp_( 0 ),
    DPNorm_( 0.0 ),
    DPRate_( "DPRate", 0.0, 0.0, 1000.0 ),
    meanDPEff_( "meanDPEff", 0.0, 0.0, 1.0 ),
    currentEvent_( 0 ),
    symmetricalDP_( kFALSE ),
    fullySymmetricDP_( kFALSE ),
    flavConjDP_( kFALSE ),
    integralsDone_( kFALSE ),
    normalizationSchemeDone_( kFALSE ),
    forceSymmetriseIntegration_( kFALSE ),
    intFileName_( "integ.dat" ),
    m13BinWidth_( 0.005 ),
    m23BinWidth_( 0.005 ),
    mPrimeBinWidth_( 0.001 ),
    thPrimeBinWidth_( 0.001 ),
    narrowWidth_( 0.020 ),
    binningFactor_( 100.0 ),
    m13Sq_( 0.0 ),
    m23Sq_( 0.0 ),
    mPrime_( 0.0 ),
    thPrime_( 0.0 ),
    tagCat_( -1 ),
    eff_( 1.0 ),
    scfFraction_( 0.0 ),
    jacobian_( 0.0 ),
    ASq_( 0.0 ),
    evtLike_( 0.0 ),
    iterationsMax_( 100000 ),
    nSigGenLoop_( 0 ),
    aSqMaxSet_( 1.25 ),
    aSqMaxVar_( 0.0 ),
    flipHelicity_( kTRUE ),
    recalcNormalisation_( kFALSE ),
    calculateRhoOmegaFitFractions_( kFALSE )
{
    if ( daughters != 0 ) {
        symmetricalDP_    = daughters->gotSymmetricalDP();
        fullySymmetricDP_ = daughters->gotFullySymmetricDP();
        flavConjDP_       = daughters->gotFlavourConjugateDP();
        typDaug_.push_back( daughters->getTypeDaug1() );
        typDaug_.push_back( daughters->getTypeDaug2() );
        typDaug_.push_back( daughters->getTypeDaug3() );
    }

    if ( scfFractionModel != 0 ) {
        scfFractionModel_[0] = scfFractionModel;
    }

    sigResonances_.clear();
    sigIncohResonances_.clear();
    kMatrixPropagators_.clear();
    kMatrixPropSet_.clear();
    extraParameters_.clear();
}

// for Kspipi, we need a scfFraction 2D histogram for each tagging category. They are provided by the map.
// Also, we need to know the place that the tagging category of the current event occupies in the data structure inputFitTree
LauIsobarDynamics::LauIsobarDynamics( LauDaughters* daughters,
                                      LauAbsEffModel* effModel,
                                      LauTagCatScfFractionModelMap scfFractionModel ) :
    daughters_( daughters ),
    kinematics_( daughters_ ? daughters_->getKinematics() : 0 ),
    effModel_( effModel ),
    scfFractionModel_( scfFractionModel ),
    nAmp_( 0 ),
    nIncohAmp_( 0 ),
    DPNorm_( 0.0 ),
    DPRate_( "DPRate", 0.0, 0.0, 1000.0 ),
    meanDPEff_( "meanDPEff", 0.0, 0.0, 1.0 ),
    currentEvent_( 0 ),
    symmetricalDP_( kFALSE ),
    fullySymmetricDP_( kFALSE ),
    flavConjDP_( kFALSE ),
    integralsDone_( kFALSE ),
    normalizationSchemeDone_( kFALSE ),
    forceSymmetriseIntegration_( kFALSE ),
    intFileName_( "integ.dat" ),
    m13BinWidth_( 0.005 ),
    m23BinWidth_( 0.005 ),
    mPrimeBinWidth_( 0.001 ),
    thPrimeBinWidth_( 0.001 ),
    narrowWidth_( 0.020 ),
    binningFactor_( 100.0 ),
    m13Sq_( 0.0 ),
    m23Sq_( 0.0 ),
    mPrime_( 0.0 ),
    thPrime_( 0.0 ),
    tagCat_( -1 ),
    eff_( 1.0 ),
    scfFraction_( 0.0 ),
    jacobian_( 0.0 ),
    ASq_( 0.0 ),
    evtLike_( 0.0 ),
    iterationsMax_( 100000 ),
    nSigGenLoop_( 0 ),
    aSqMaxSet_( 1.25 ),
    aSqMaxVar_( 0.0 ),
    flipHelicity_( kTRUE ),
    recalcNormalisation_( kFALSE ),
    calculateRhoOmegaFitFractions_( kFALSE )
{
    // Constructor for the isobar signal model
    if ( daughters != 0 ) {
        symmetricalDP_    = daughters->gotSymmetricalDP();
        fullySymmetricDP_ = daughters->gotFullySymmetricDP();
        flavConjDP_       = daughters->gotFlavourConjugateDP();
        typDaug_.push_back( daughters->getTypeDaug1() );
        typDaug_.push_back( daughters->getTypeDaug2() );
        typDaug_.push_back( daughters->getTypeDaug3() );
    }

    sigResonances_.clear();
    sigIncohResonances_.clear();
    kMatrixPropagators_.clear();
    kMatrixPropSet_.clear();
    extraParameters_.clear();
}

LauIsobarDynamics::~LauIsobarDynamics()
{
    extraParameters_.clear();

    for ( std::vector<LauCacheData*>::iterator iter = data_.begin(); iter != data_.end(); ++iter ) {
        delete ( *iter );
    }
    data_.clear();

    for ( std::vector<LauDPPartialIntegralInfo*>::iterator it = dpPartialIntegralInfo_.begin();
          it != dpPartialIntegralInfo_.end();
          ++it ) {
        delete ( *it );
    }
    dpPartialIntegralInfo_.clear();
}

void LauIsobarDynamics::resetNormVectors()
{
    for ( UInt_t i = 0; i < nAmp_; i++ ) {

        fSqSum_[i]    = 0.0;
        fSqEffSum_[i] = 0.0;
        fNorm_[i]     = 0.0;
        ff_[i].zero();

        for ( UInt_t j = 0; j < nAmp_; j++ ) {

            fifjEffSum_[i][j].zero();
            fifjSum_[i][j].zero();
        }
    }

    for ( UInt_t i = 0; i < nIncohAmp_; i++ ) {
        fSqSum_[i + nAmp_]    = 0.0;
        fSqEffSum_[i + nAmp_] = 0.0;
        fNorm_[i + nAmp_]     = 0.0;
        incohInten_[i]        = 0.0;
    }
}

void LauIsobarDynamics::recalculateNormalisation()
{
    if ( recalcNormalisation_ == kFALSE ) {
        return;
    }

    // We need to calculate the normalisation constants for the
    // Dalitz plot generation/fitting.

    integralsDone_ = kFALSE;

    this->resetNormVectors();
    this->findIntegralsToBeRecalculated();
    this->calcDPNormalisation();

    integralsDone_ = kTRUE;
}

void LauIsobarDynamics::findIntegralsToBeRecalculated()
{
    // Loop through the resonance parameters and see which ones have changed
    // For those that have changed mark the corresponding resonance(s) as needing to be re-evaluated

    integralsToBeCalculated_.clear();

    const UInt_t nResPars = resonancePars_.size();
    for ( UInt_t iPar( 0 ); iPar < nResPars; ++iPar ) {
        const Double_t newValue = resonancePars_[iPar]->value();
        if ( newValue != resonanceParValues_[iPar] ) {
            resonanceParValues_[iPar] = newValue;

            const std::vector<UInt_t>& indices                 = resonanceParResIndex_[iPar];
            std::vector<UInt_t>::const_iterator indexIter      = indices.begin();
            const std::vector<UInt_t>::const_iterator indexEnd = indices.end();
            for ( ; indexIter != indexEnd; ++indexIter ) {
                integralsToBeCalculated_.insert( *indexIter );
            }
        }
    }
}

void LauIsobarDynamics::collateResonanceParameters()
{
    // Initialise all resonance models
    resonancePars_.clear();
    resonanceParValues_.clear();
    resonanceParResIndex_.clear();

    std::set<LauParameter*> uniqueResPars;
    UInt_t resIndex( 0 );
    for ( std::vector<LauAbsResonance*>::iterator resIter = sigResonances_.begin();
          resIter != sigResonances_.end();
          ++resIter ) {
        ( *resIter )->initialise();

        // Check if this resonance has floating parameters
        // Append all unique parameters to our list
        const std::vector<LauParameter*>& resPars = ( *resIter )->getFloatingParameters();

        for ( std::vector<LauParameter*>::const_iterator parIter = resPars.begin();
              parIter != resPars.end();
              ++parIter ) {
            if ( uniqueResPars.insert( *parIter ).second ) {
                // This parameter has not already been added to
                // the list of unique ones.  Add it, its value
                // and its associated resonance ID to the
                // appropriate lists.
                resonancePars_.push_back( *parIter );
                resonanceParValues_.push_back( ( *parIter )->value() );
                std::vector<UInt_t> resIndices( 1, resIndex );
                resonanceParResIndex_.push_back( resIndices );
            } else {
                // This parameter has already been added to the
                // list of unique ones.  However, we still need
                // to indicate that this resonance should be
                // associated with it.
                std::vector<LauParameter*>::const_iterator uniqueParIter = resonancePars_.begin();
                std::vector<std::vector<UInt_t>>::iterator indicesIter = resonanceParResIndex_.begin();
                while ( ( *uniqueParIter ) != ( *parIter ) ) {
                    ++uniqueParIter;
                    ++indicesIter;
                }
                ( *indicesIter ).push_back( resIndex );
            }
        }

        ++resIndex;
    }

    for ( std::vector<LauAbsIncohRes*>::iterator resIter = sigIncohResonances_.begin();
          resIter != sigIncohResonances_.end();
          ++resIter ) {
        ( *resIter )->initialise();

        // Check if this resonance has floating parameters
        // Append all unique parameters to our list
        const std::vector<LauParameter*>& resPars = ( *resIter )->getFloatingParameters();

        for ( std::vector<LauParameter*>::const_iterator parIter = resPars.begin();
              parIter != resPars.end();
              ++parIter ) {
            if ( uniqueResPars.insert( *parIter ).second ) {
                // This parameter has not already been added to
                // the list of unique ones.  Add it, its value
                // and its associated resonance ID to the
                // appropriate lists.
                resonancePars_.push_back( *parIter );
                resonanceParValues_.push_back( ( *parIter )->value() );
                std::vector<UInt_t> resIndices( 1, resIndex );
                resonanceParResIndex_.push_back( resIndices );
            } else {
                // This parameter has already been added to the
                // list of unique ones.  However, we still need
                // to indicate that this resonance should be
                // associated with it.
                std::vector<LauParameter*>::const_iterator uniqueParIter = resonancePars_.begin();
                std::vector<std::vector<UInt_t>>::iterator indicesIter = resonanceParResIndex_.begin();
                while ( ( *uniqueParIter ) != ( *parIter ) ) {
                    ++uniqueParIter;
                    ++indicesIter;
                }
                ( *indicesIter ).push_back( resIndex );
            }
        }

        ++resIndex;
    }
}

void LauIsobarDynamics::initialise( const std::vector<LauComplex>& coeffs )
{
    // Check whether we have a valid set of integration constants for
    // the normalisation of the signal likelihood function.
    this->initialiseVectors();

    // Mark the DP integrals as undetermined
    integralsDone_ = kFALSE;

    this->collateResonanceParameters();

    if ( resonancePars_.empty() ) {
        recalcNormalisation_ = kFALSE;
    } else {
        recalcNormalisation_ = kTRUE;
    }

    // Print summary of what we have so far to screen
    this->initSummary();

    if ( nAmp_ + nIncohAmp_ == 0 ) {
        std::cout << "INFO in LauIsobarDynamics::initialise : No contributions to DP model, not performing normalisation integrals."
                  << std::endl;
    } else {

        // We need to calculate the normalisation constants for the Dalitz plot generation/fitting.

        std::cout << "INFO in LauIsobarDynamics::initialise : Starting special run to generate the integrals for normalising the PDF..."
                  << std::endl;

        // Since this is the initialisation, we need to calculate everything for every resonance
        integralsToBeCalculated_.clear();
        for ( UInt_t i( 0 ); i < nAmp_ + nIncohAmp_; ++i ) {
            integralsToBeCalculated_.insert( i );
        }

        // Calculate and cache the normalisations of each resonance _dynamic_ amplitude
        // (e.g. Breit-Wigner contribution, not from the complex coefficients).
        // These are stored in fNorm_[i].
        // fSqSum[i] is the total of the dynamical amplitude squared for a given resonance, i.
        // We require that:
        // |fNorm_[i]|^2 * |fSqSum[i]|^2 = 1,
        // i.e. fNorm_[i] normalises each resonance contribution to give the same number of
        // events in the DP, accounting for the total DP area and the dynamics of the resonance.
        this->calcDPNormalisation();

        // Write the integrals to a file (mainly for debugging purposes)
        this->writeIntegralsFile();
    }

    integralsDone_ = kTRUE;

    std::cout << std::setprecision( 10 );

    std::cout << "INFO in LauIsobarDynamics::initialise : Summary of the integrals:" << std::endl;

    for ( UInt_t i = 0; i < nAmp_ + nIncohAmp_; i++ ) {
        std::cout << " fNorm[" << i << "] = " << fNorm_[i] << std::endl;
        std::cout << " fSqSum[" << i << "] = " << fSqSum_[i] << std::endl;
        std::cout << " fSqEffSum[" << i << "] = " << fSqEffSum_[i] << std::endl;
    }

    for ( UInt_t i = 0; i < nAmp_; i++ ) {
        for ( UInt_t j = 0; j < nAmp_; j++ ) {
            std::cout << " fifjEffSum[" << i << "][" << j << "] = " << fifjEffSum_[i][j];
        }
        std::cout << std::endl;
    }

    for ( UInt_t i = 0; i < nAmp_; i++ ) {
        for ( UInt_t j = 0; j < nAmp_; j++ ) {
            std::cout << " fifjSum[" << i << "][" << j << "] = " << fifjSum_[i][j];
        }
        std::cout << std::endl;
    }

    // Calculate the initial fit fractions (for later comparison with Toy MC, if required)
    this->updateCoeffs( coeffs );
    this->calcExtraInfo( kTRUE );

    for ( UInt_t i = 0; i < nAmp_; i++ ) {
        for ( UInt_t j = i; j < nAmp_; j++ ) {
            std::cout << "INFO in LauIsobarDynamics::initialise : Initial fit fraction for amplitude ("
                      << i << "," << j << ") = " << fitFrac_[i][j].genValue() << std::endl;
        }
    }

    for ( UInt_t i = 0; i < nIncohAmp_; i++ ) {
        std::cout << "INFO in LauIsobarDynamics::initialise : Initial fit fraction for incoherent amplitude ("
                  << i << "," << i << ") = " << fitFrac_[i + nAmp_][i + nAmp_].genValue()
                  << std::endl;
    }

    std::cout << "INFO in LauIsobarDynamics::initialise : Initial efficiency = "
              << meanDPEff_.initValue() << std::endl;

    std::cout << "INFO in LauIsobarDynamics::initialise : Initial DPRate = " << DPRate_.initValue()
              << std::endl;
}

void LauIsobarDynamics::initSummary()
{

    UInt_t i( 0 );
    TString nameP = daughters_->getNameParent();
    TString name1 = daughters_->getNameDaug1();
    TString name2 = daughters_->getNameDaug2();
    TString name3 = daughters_->getNameDaug3();

    std::cout << "INFO in LauIsobarDynamics::initSummary : We are going to do a DP with " << nameP
              << " going to " << name1 << " " << name2 << " " << name3 << std::endl;

    std::cout << "                                       : For the following resonance combinations:"
              << std::endl;

    std::cout << "                                       : In tracks 2 and 3:" << std::endl;
    for ( i = 0; i < nAmp_; i++ ) {
        if ( resPairAmp_[i] == 1 ) {
            std::cout << "                                       :   A" << i << ": "
                      << resTypAmp_[i] << " to " << name2 << ", " << name3 << std::endl;
        }
    }

    for ( i = 0; i < nIncohAmp_; i++ ) {
        if ( incohResPairAmp_[i] == 1 ) {
            std::cout << "                                       :   A" << nAmp_ + i << ": "
                      << incohResTypAmp_[i] << " (incoherent) to " << name2 << ", " << name3
                      << std::endl;
        }
    }

    std::cout << "                                       : In tracks 1 and 3:" << std::endl;
    for ( i = 0; i < nAmp_; i++ ) {
        if ( resPairAmp_[i] == 2 ) {
            std::cout << "                                       :   A" << i << ": "
                      << resTypAmp_[i] << " to " << name1 << ", " << name3 << std::endl;
        }
    }

    for ( i = 0; i < nIncohAmp_; i++ ) {
        if ( incohResPairAmp_[i] == 2 ) {
            std::cout << "                                       :   A" << nAmp_ + i << ": "
                      << incohResTypAmp_[i] << " (incoherent) to " << name1 << ", " << name3
                      << std::endl;
        }
    }

    std::cout << "                                       : In tracks 1 and 2:" << std::endl;
    for ( i = 0; i < nAmp_; i++ ) {
        if ( resPairAmp_[i] == 3 ) {
            std::cout << "                                       :   A" << i << ": "
                      << resTypAmp_[i] << " to " << name1 << ", " << name2 << std::endl;
        }
    }

    for ( i = 0; i < nIncohAmp_; i++ ) {
        if ( incohResPairAmp_[i] == 3 ) {
            std::cout << "                                       :   A" << nAmp_ + i << ": "
                      << incohResTypAmp_[i] << " (incoherent) to " << name1 << ", " << name2
                      << std::endl;
        }
    }

    for ( i = 0; i < nAmp_; i++ ) {
        if ( resPairAmp_[i] == 0 ) {
            std::cout << "                                       : and a non-resonant amplitude:"
                      << std::endl;
            std::cout << "                                       :   A" << i << std::endl;
        }
    }

    std::cout << std::endl;
}

void LauIsobarDynamics::initialiseVectors()
{
    fSqSum_.clear();
    fSqSum_.resize( nAmp_ + nIncohAmp_ );
    fSqEffSum_.clear();
    fSqEffSum_.resize( nAmp_ + nIncohAmp_ );
    fifjEffSum_.clear();
    fifjEffSum_.resize( nAmp_ );
    fifjSum_.clear();
    fifjSum_.resize( nAmp_ );
    ff_.clear();
    ff_.resize( nAmp_ );
    incohInten_.clear();
    incohInten_.resize( nIncohAmp_ );
    fNorm_.clear();
    fNorm_.resize( nAmp_ + nIncohAmp_ );
    fitFrac_.clear();
    fitFrac_.resize( nAmp_ + nIncohAmp_ );
    fitFracEffUnCorr_.clear();
    fitFracEffUnCorr_.resize( nAmp_ + nIncohAmp_ );

    LauComplex null( 0.0, 0.0 );

    for ( UInt_t i = 0; i < nAmp_; i++ ) {

        fSqSum_[i] = 0.0;
        fNorm_[i]  = 0.0;
        ff_[i]     = null;
        fifjEffSum_[i].resize( nAmp_ );
        fifjSum_[i].resize( nAmp_ );
        fitFrac_[i].resize( nAmp_ + nIncohAmp_ );
        fitFracEffUnCorr_[i].resize( nAmp_ + nIncohAmp_ );

        for ( UInt_t j = 0; j < nAmp_; j++ ) {
            fifjEffSum_[i][j] = null;
            fifjSum_[i][j]    = null;
            fitFrac_[i][j].valueAndRange( 0.0, -100.0, 100.0 );
            fitFracEffUnCorr_[i][j].valueAndRange( 0.0, -100.0, 100.0 );
        }
        for ( UInt_t j = nAmp_; j < nAmp_ + nIncohAmp_; j++ ) {
            fitFrac_[i][j].valueAndRange( 0.0, -100.0, 100.0 );
            fitFracEffUnCorr_[i][j].valueAndRange( 0.0, -100.0, 100.0 );
        }
    }
    for ( UInt_t i = nAmp_; i < nAmp_ + nIncohAmp_; i++ ) {
        fSqSum_[i]             = 0.0;
        fNorm_[i]              = 0.0;
        incohInten_[i - nAmp_] = 0.0;
        fitFrac_[i].resize( nAmp_ + nIncohAmp_ );
        fitFracEffUnCorr_[i].resize( nAmp_ + nIncohAmp_ );

        for ( UInt_t j = 0; j < nAmp_ + nIncohAmp_; j++ ) {
            fitFrac_[i][j].valueAndRange( 0.0, -100.0, 100.0 );
            fitFracEffUnCorr_[i][j].valueAndRange( 0.0, -100.0, 100.0 );
        }
    }

    UInt_t nKMatrixProps = kMatrixPropagators_.size();
    extraParameters_.clear();
    extraParameters_.resize( nKMatrixProps );

    for ( UInt_t i = 0; i < nKMatrixProps; ++i ) {
        extraParameters_[i].valueAndRange( 0.0, -100.0, 100.0 );
    }
}

void LauIsobarDynamics::writeIntegralsFile()
{
    // Routine to end integration calculation for loglike normalisation.
    // This writes out the normalisation integral output into the file named
    // outputFileName.
    std::cout << "INFO in LauIsobarDynamics::writeIntegralsFile : Writing integral output to integrals file "
              << intFileName_.Data() << std::endl;

    UInt_t i( 0 ), j( 0 );
    std::ofstream getChar( intFileName_.Data() );

    getChar << std::setprecision( 10 );

    // Write out daughter types (pi, pi0, K, K0s?)
    for ( i = 0; i < 3; i++ ) {
        getChar << typDaug_[i] << " ";
    }

    // Write out number of resonances in the Dalitz plot model
    getChar << nAmp_ << std::endl;

    // Write out the resonances
    for ( i = 0; i < nAmp_; i++ ) {
        getChar << resTypAmp_[i] << " ";
    }

    getChar << std::endl;

    // Write out the resonance model types (BW, RelBW etc...)
    for ( i = 0; i < nAmp_; i++ ) {
        LauAbsResonance* theResonance = sigResonances_[i];
        Int_t resModelInt             = theResonance->getResonanceModel();
        getChar << resModelInt << " ";
    }

    getChar << std::endl;

    // Write out the track pairings for each resonance. This is specified
    // by the resPairAmpInt integer in the addResonance function.
    for ( i = 0; i < nAmp_; i++ ) {
        getChar << resPairAmp_[i] << " ";
    }

    getChar << std::endl;

    // Write out the fSqSum = |ff|^2, where ff = resAmp()
    for ( i = 0; i < nAmp_; i++ ) {
        getChar << fSqSum_[i] << " ";
    }

    getChar << std::endl;

    // Similar to fSqSum, but with the efficiency term included.
    for ( i = 0; i < nAmp_; i++ ) {
        getChar << fSqEffSum_[i] << " ";
    }

    getChar << std::endl;

    // Write out the f_i*f_j_conj*eff values = resAmp_i*resAmp_j_conj*eff.
    // Note that only the top half of the i*j "matrix" is required, as it
    // is symmetric w.r.t i, j.
    for ( i = 0; i < nAmp_; i++ ) {
        for ( j = i; j < nAmp_; j++ ) {
            getChar << fifjEffSum_[i][j] << " ";
        }
    }

    getChar << std::endl;

    // Similar to fifjEffSum, but without the efficiency term included.
    for ( i = 0; i < nAmp_; i++ ) {
        for ( j = i; j < nAmp_; j++ ) {
            getChar << fifjSum_[i][j] << " ";
        }
    }

    getChar << std::endl;

    // Write out number of incoherent resonances in the Dalitz plot model
    getChar << nIncohAmp_ << std::endl;

    // Write out the incoherent resonances
    for ( i = 0; i < nIncohAmp_; i++ ) {
        getChar << incohResTypAmp_[i] << " ";
    }

    getChar << std::endl;

    // Write out the incoherent resonance model types (BW, RelBW etc...)
    for ( i = 0; i < nIncohAmp_; i++ ) {
        LauAbsResonance* theResonance = sigIncohResonances_[i];
        Int_t resModelInt             = theResonance->getResonanceModel();
        getChar << resModelInt << " ";
    }

    getChar << std::endl;

    // Write out the track pairings for each incoherent resonance. This is specified
    // by the resPairAmpInt integer in the addIncohResonance function.
    for ( i = 0; i < nIncohAmp_; i++ ) {
        getChar << incohResPairAmp_[i] << " ";
    }

    getChar << std::endl;

    // Write out the fSqSum = |ff|^2, where |ff|^2 = incohResAmp()
    for ( i = nAmp_; i < nAmp_ + nIncohAmp_; i++ ) {
        getChar << fSqSum_[i] << " ";
    }

    getChar << std::endl;

    // Similar to fSqSum, but with the efficiency term included.
    for ( i = nAmp_; i < nAmp_ + nIncohAmp_; i++ ) {
        getChar << fSqEffSum_[i] << " ";
    }

    getChar << std::endl;
}

LauAbsResonance* LauIsobarDynamics::addResonance(
    const TString& resName,
    const Int_t resPairAmpInt,
    const LauAbsResonance::LauResonanceModel resType,
    const LauBlattWeisskopfFactor::BlattWeisskopfCategory bwCategory )
{
    // Function to add a resonance in a Dalitz plot.
    // No check is made w.r.t flavour and charge conservation rules, and so
    // the user is responsible for checking the internal consistency of
    // their function statements with these laws. For example, the program
    // will not prevent the user from asking for a rho resonance in a K-pi
    // pair or a K* resonance in a pi-pi pair.
    // However, to assist the user, a summary of the resonant structure requested
    // by the user is printed before the program runs. It is important to check this
    // information when you first define your Dalitz plot model before doing
    // any fitting/generating.
    // Arguments are: resonance name, integer to specify the resonance track pairing
    // (1 => m_23, 2 => m_13, 3 => m_12), i.e. the bachelor track number.
    // The third argument resType specifies whether the resonance is a Breit-Wigner (BW)
    // Relativistic Breit-Wigner (RelBW) or Flatte distribution (Flatte), for example.

    if ( LauAbsResonance::isIncoherentModel( resType ) == true ) {
        std::cerr << "ERROR in LauIsobarDynamics::addResonance : Resonance type \"" << resType
                  << "\" not allowed for a coherent resonance" << std::endl;
        return 0;
    }

    LauResonanceMaker& resonanceMaker = LauResonanceMaker::get();
    LauAbsResonance* theResonance =
        resonanceMaker.getResonance( daughters_, resName, resPairAmpInt, resType, bwCategory );

    if ( theResonance == 0 ) {
        std::cerr << "ERROR in LauIsobarDynamics::addResonance : Couldn't create the resonance \""
                  << resName << "\"" << std::endl;
        return 0;
    }

    // implement the helicity flip here
    if ( flipHelicity_ && daughters_->getCharge( resPairAmpInt ) == 0 &&
         daughters_->getChargeParent() == 0 && daughters_->getTypeParent() > 0 ) {
        if ( ( resPairAmpInt == 1 && TMath::Abs( daughters_->getTypeDaug2() ) ==
                                         TMath::Abs( daughters_->getTypeDaug3() ) ) ||
             ( resPairAmpInt == 2 && TMath::Abs( daughters_->getTypeDaug1() ) ==
                                         TMath::Abs( daughters_->getTypeDaug3() ) ) ||
             ( resPairAmpInt == 3 && TMath::Abs( daughters_->getTypeDaug1() ) ==
                                         TMath::Abs( daughters_->getTypeDaug2() ) ) ) {
            theResonance->flipHelicity( kTRUE );
        }
    }

    // Set the resonance name and what track is the bachelor
    TString resonanceName = theResonance->getResonanceName();
    resTypAmp_.push_back( resonanceName );

    // Always force the non-resonant amplitude pair to have resPairAmp = 0
    // in case the user chooses the wrong number.
    if ( resType == LauAbsResonance::FlatNR || resType == LauAbsResonance::NRModel ) {
        std::cout << "INFO in LauIsobarDynamics::addResonance : Setting resPairAmp to 0 for "
                  << resonanceName << " contribution." << std::endl;
        resPairAmp_.push_back( 0 );
    } else {
        resPairAmp_.push_back( resPairAmpInt );
    }

    // Increment the number of resonance amplitudes we have so far
    ++nAmp_;

    // Finally, add the resonance object to the internal array
    sigResonances_.push_back( theResonance );

    std::cout << "INFO in LauIsobarDynamics::addResonance : Successfully added resonance. Total number of resonances so far = "
              << nAmp_ << std::endl;

    return theResonance;
}

LauAbsResonance* LauIsobarDynamics::addIncoherentResonance(
    const TString& resName,
    const Int_t resPairAmpInt,
    const LauAbsResonance::LauResonanceModel resType )
{
    // Function to add an incoherent resonance in a Dalitz plot.
    // No check is made w.r.t flavour and charge conservation rules, and so
    // the user is responsible for checking the internal consistency of
    // their function statements with these laws. For example, the program
    // will not prevent the user from asking for a rho resonance in a K-pi
    // pair or a K* resonance in a pi-pi pair.
    // However, to assist the user, a summary of the resonant structure requested
    // by the user is printed before the program runs. It is important to check this
    // information when you first define your Dalitz plot model before doing
    // any fitting/generating.
    // Arguments are: resonance name, integer to specify the resonance track pairing
    // (1 => m_23, 2 => m_13, 3 => m_12), i.e. the bachelor track number.
    // The third argument resType specifies the shape of the resonance
    // Gaussian (GaussIncoh), for example.

    if ( LauAbsResonance::isIncoherentModel( resType ) == false ) {
        std::cerr << "ERROR in LauIsobarDynamics::addIncohResonance : Resonance type \"" << resType
                  << "\" not allowed for an incoherent resonance" << std::endl;
        return 0;
    }

    LauResonanceMaker& resonanceMaker = LauResonanceMaker::get();
    LauAbsIncohRes* theResonance      = dynamic_cast<LauAbsIncohRes*>(
        resonanceMaker.getResonance( daughters_, resName, resPairAmpInt, resType ) );

    if ( theResonance == 0 ) {
        std::cerr << "ERROR in LauIsobarDynamics::addIncohResonance : Couldn't create the resonance \""
                  << resName << "\"" << std::endl;
        return 0;
    }

    // Set the resonance name and what track is the bachelor
    TString resonanceName = theResonance->getResonanceName();
    incohResTypAmp_.push_back( resonanceName );
    incohResPairAmp_.push_back( resPairAmpInt );

    // Increment the number of resonance amplitudes we have so far
    ++nIncohAmp_;

    // Finally, add the resonance object to the internal array
    sigIncohResonances_.push_back( theResonance );

    std::cout << "INFO in LauIsobarDynamics::addIncohResonance : Successfully added incoherent resonance. Total number of incoherent resonances so far = "
              << nIncohAmp_ << std::endl;

    return theResonance;
}

LauKMatrixPropagator* LauIsobarDynamics::defineKMatrixPropagator( const TString& propName,
                                                                  const TString& paramFileName,
                                                                  Int_t resPairAmpInt,
                                                                  Int_t nChannels,
                                                                  Int_t nPoles,
                                                                  Int_t rowIndex )
{
    // Define the K-matrix propagator. The resPairAmpInt integer specifies which mass combination should be used
    // for the invariant mass-squared variable "s". The pole masses and coupling constants are defined in the
    // paramFileName parameter file.
    // The number of channels and poles are defined by the nChannels and nPoles integers, respectively.
    // The integer rowIndex specifies which row of the propagator should be used when
    // summing over all amplitude channels: S-wave will be the first row, so rowIndex = 1.

    // Check that the rowIndex is valid
    if ( rowIndex < 1 || rowIndex > nChannels ) {
        std::cerr << "ERROR in LauIsobarDynamics::defineKMatrixPropagator. The rowIndex, which is set to "
                  << rowIndex << ", must be between 1 and the number of channels " << nChannels
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    TString propagatorName( propName ), parameterFile( paramFileName );

    LauKMatrixPropagator* thePropagator =
        LauKMatrixPropFactory::getInstance()->getPropagator( propagatorName,
                                                             parameterFile,
                                                             resPairAmpInt,
                                                             nChannels,
                                                             nPoles,
                                                             rowIndex );
    kMatrixPropagators_[propagatorName] = thePropagator;

    return thePropagator;
}

void LauIsobarDynamics::addKMatrixProdPole( const TString& poleName,
                                            const TString& propName,
                                            Int_t poleIndex,
                                            Bool_t useProdAdler )
{

    // Add a K-matrix production pole term, using the K-matrix propagator given by the propName.
    // Here, poleIndex is the integer specifying the pole number.

    // First, find the K-matrix propagator.
    KMPropMap::iterator mapIter = kMatrixPropagators_.find( propName );
    if ( mapIter != kMatrixPropagators_.end() ) {

        LauKMatrixPropagator* thePropagator = mapIter->second;

        // Make sure the pole index is valid
        Int_t nPoles = thePropagator->getNPoles();
        if ( poleIndex < 1 || poleIndex > nPoles ) {
            std::cerr << "ERROR in LauIsobarDynamics::addKMatrixProdPole : The pole index "
                      << poleIndex << " is not between 1 and " << nPoles
                      << ". Not adding production pole " << poleName << " for K-matrix propagator "
                      << propName << std::endl;
            return;
        }

        // Now add the K-matrix production pole amplitude to the vector of LauAbsResonance pointers.
        Int_t resPairAmpInt       = thePropagator->getResPairAmpInt();
        LauAbsResonance* prodPole = new LauKMatrixProdPole( poleName,
                                                            poleIndex,
                                                            resPairAmpInt,
                                                            thePropagator,
                                                            daughters_,
                                                            useProdAdler );
        prodPole->setSpinType( LauAbsResonance::Legendre );

        resTypAmp_.push_back( poleName );
        resPairAmp_.push_back( resPairAmpInt );

        ++nAmp_;
        sigResonances_.push_back( prodPole );

        // Also store the propName-poleName pair for calculating total fit fractions later on
        // (avoiding the need to use dynamic casts to check which resonances are of the K-matrix type)
        kMatrixPropSet_[poleName] = propName;

        std::cout << "INFO in LauIsobarDynamics::addKMatrixProdPole : Successfully added K-matrix production pole term. Total number of resonances so far = "
                  << nAmp_ << std::endl;

    } else {

        std::cerr << "ERROR in LauIsobarDynamics::addKMatrixProdPole : The propagator of the name "
                  << propName << " could not be found for the production pole " << poleName
                  << std::endl;
    }
}

void LauIsobarDynamics::addKMatrixProdSVP( const TString& SVPName,
                                           const TString& propName,
                                           Int_t channelIndex,
                                           Bool_t useProdAdler )
{

    // Add a K-matrix production "slowly-varying part" (SVP) term, using the K-matrix propagator
    // given by the propName. Here, channelIndex is the integer specifying the channel number.

    // First, find the K-matrix propagator.
    KMPropMap::iterator mapIter = kMatrixPropagators_.find( propName );
    if ( mapIter != kMatrixPropagators_.end() ) {

        LauKMatrixPropagator* thePropagator = mapIter->second;

        // Make sure the channel index is valid
        Int_t nChannels = thePropagator->getNChannels();
        if ( channelIndex < 1 || channelIndex > nChannels ) {
            std::cerr << "ERROR in LauIsobarDynamics::addKMatrixProdSVP : The channel index "
                      << channelIndex << " is not between 1 and " << nChannels
                      << ". Not adding production slowly-varying part " << SVPName
                      << " for K-matrix propagator " << propName << std::endl;
            return;
        }

        // Now add the K-matrix production SVP amplitude to the vector of LauAbsResonance pointers.
        Int_t resPairAmpInt      = thePropagator->getResPairAmpInt();
        LauAbsResonance* prodSVP = new LauKMatrixProdSVP( SVPName,
                                                          channelIndex,
                                                          resPairAmpInt,
                                                          thePropagator,
                                                          daughters_,
                                                          useProdAdler );
        prodSVP->setSpinType( LauAbsResonance::Legendre );

        resTypAmp_.push_back( SVPName );
        resPairAmp_.push_back( resPairAmpInt );

        ++nAmp_;
        sigResonances_.push_back( prodSVP );

        // Also store the SVPName-propName pair for calculating total fit fractions later on
        // (avoiding the need to use dynamic casts to check which resonances are of the K-matrix type)
        kMatrixPropSet_[SVPName] = propName;

        std::cout << "INFO in LauIsobarDynamics::addKMatrixProdSVP : Successfully added K-matrix production slowly-varying (SVP) term. Total number of resonances so far = "
                  << nAmp_ << std::endl;

    } else {

        std::cerr << "ERROR in LauIsobarDynamics::addKMatrixProdSVP : The propagator of the name "
                  << propName << " could not be found for the production slowly-varying part "
                  << SVPName << std::endl;
    }
}

Int_t LauIsobarDynamics::resonanceIndex( const TString& resName ) const
{
    Int_t index( 0 );
    const LauAbsResonance* theResonance( 0 );

    for ( std::vector<LauAbsResonance*>::const_iterator iter = sigResonances_.begin();
          iter != sigResonances_.end();
          ++iter ) {
        theResonance = ( *iter );
        if ( theResonance != 0 ) {
            const TString& resString = theResonance->getResonanceName();
            if ( resString == resName ) {
                return index;
            }
        }
        ++index;
    }

    for ( std::vector<LauAbsIncohRes*>::const_iterator iter = sigIncohResonances_.begin();
          iter != sigIncohResonances_.end();
          ++iter ) {
        theResonance = ( *iter );
        if ( theResonance != 0 ) {
            const TString& resString = theResonance->getResonanceName();
            if ( resString == resName ) {
                return index;
            }
        }
        ++index;
    }

    return -1;
}

Bool_t LauIsobarDynamics::hasResonance( const TString& resName ) const
{
    const Int_t index = this->resonanceIndex( resName );
    if ( index < 0 ) {
        return kFALSE;
    } else {
        return kTRUE;
    }
}

const LauAbsResonance* LauIsobarDynamics::getResonance( const UInt_t resIndex ) const
{
    if ( resIndex < this->getnCohAmp() ) {
        return sigResonances_[resIndex];
    } else if ( resIndex < this->getnTotAmp() ) {
        return sigIncohResonances_[resIndex - nAmp_];
    } else {
        std::cerr << "ERROR in LauIsobarDynamics::getResonance : Couldn't find resonance with index \""
                  << resIndex << "\" in the model." << std::endl;
        return 0;
    }
}

LauAbsResonance* LauIsobarDynamics::getResonance( const UInt_t resIndex )
{
    if ( resIndex < this->getnCohAmp() ) {
        return sigResonances_[resIndex];
    } else if ( resIndex < this->getnTotAmp() ) {
        return sigIncohResonances_[resIndex - nAmp_];
    } else {
        std::cerr << "ERROR in LauIsobarDynamics::getResonance : Couldn't find resonance with index \""
                  << resIndex << "\" in the model." << std::endl;
        return 0;
    }
}

LauAbsResonance* LauIsobarDynamics::findResonance( const TString& resName )
{
    const Int_t index = this->resonanceIndex( resName );
    if ( index < 0 ) {
        std::cerr << "ERROR in LauIsobarDynamics::findResonance : Couldn't find resonance with name \""
                  << resName << "\" in the model." << std::endl;
        return 0;
    } else {
        return this->getResonance( index );
    }
}

const LauAbsResonance* LauIsobarDynamics::findResonance( const TString& resName ) const
{
    const Int_t index = this->resonanceIndex( resName );
    if ( index < 0 ) {
        std::cerr << "ERROR in LauIsobarDynamics::findResonance : Couldn't find resonance with name \""
                  << resName << "\" in the model." << std::endl;
        return 0;
    } else {
        return this->getResonance( index );
    }
}

void LauIsobarDynamics::removeCharge( TString& string ) const
{
    Ssiz_t index = string.Index( "+" );
    if ( index != -1 ) {
        string.Remove( index, 1 );
    }
    index = string.Index( "-" );
    if ( index != -1 ) {
        string.Remove( index, 1 );
    }
}

void LauIsobarDynamics::calcDPNormalisation()
{
    if ( ! normalizationSchemeDone_ ) {
        this->calcDPNormalisationScheme();
    }

    for ( std::vector<LauDPPartialIntegralInfo*>::iterator it = dpPartialIntegralInfo_.begin();
          it != dpPartialIntegralInfo_.end();
          ++it ) {
        this->calcDPPartialIntegral( *it );
    }

    for ( UInt_t i = 0; i < nAmp_ + nIncohAmp_; ++i ) {
        fNorm_[i] = 0.0;
        if ( fSqSum_[i] > 0.0 ) {
            fNorm_[i] = TMath::Sqrt( 1.0 / ( fSqSum_[i] ) );
        }
    }
}

std::vector<std::pair<Double_t, Double_t>> LauIsobarDynamics::formGapsFromRegions(
    const std::vector<std::pair<Double_t, Double_t>>& regions,
    const Double_t min,
    const Double_t max ) const
{
    std::vector<std::pair<Double_t, Double_t>> gaps( regions.size() + 1, std::make_pair( 0., 0. ) );

    // Given some narrow resonance regions, find the regions that correspond to the gaps between them

    gaps[0].first = min;

    for ( UInt_t i = 0; i < regions.size(); ++i ) {
        gaps[i].second    = regions[i].first;
        gaps[i + 1].first = regions[i].second;
    }

    gaps[gaps.size() - 1].second = max;

    return gaps;
}

void LauIsobarDynamics::cullNullRegions( std::vector<LauDPPartialIntegralInfo*>& regions ) const
{
    LauDPPartialIntegralInfo* tmp( 0 );
    regions.erase( std::remove( regions.begin(), regions.end(), tmp ), regions.end() );
}

void LauIsobarDynamics::correctDPOverlap( std::vector<std::pair<Double_t, Double_t>>& regions,
                                          const std::vector<Double_t>& binnings ) const
{
    if ( regions.empty() ) {
        return;
    }

    // If the regions overlap, ensure that the one with the finest binning takes precedence (i.e., extends its full width)

    for ( UInt_t i = 0; i < regions.size() - 1; ++i ) {
        if ( regions[i + 1].first <= regions[i].second ) {
            if ( ( binnings[i] < binnings[i + 1] ) ) {
                regions[i + 1] = std::make_pair( regions[i].second, regions[i + 1].second );
            } else {
                regions[i] = std::make_pair( regions[i].first, regions[i + 1].first );
            }
        }
    }
}

std::vector<LauDPPartialIntegralInfo*> LauIsobarDynamics::m13IntegrationRegions(
    const std::vector<std::pair<Double_t, Double_t>>& m13Regions,
    const std::vector<std::pair<Double_t, Double_t>>& m23Regions,
    const std::vector<Double_t>& m13Binnings,
    const Double_t precision,
    const Double_t defaultBinning ) const
{
    // Create integration regions for all narrow resonances in m13 except for the overlaps with narrow resonances in m23
    std::vector<LauDPPartialIntegralInfo*> integrationRegions;

    const Double_t m23Min = kinematics_->getm23Min();
    const Double_t m23Max = kinematics_->getm23Max();

    // Loop over narrow resonances in m13
    for ( UInt_t m13i = 0; m13i < m13Regions.size(); ++m13i ) {

        const Double_t m13Binning = m13Binnings[m13i];

        const Double_t resMin13 = m13Regions[m13i].first;
        const Double_t resMax13 = m13Regions[m13i].second;

        // Initialise to the full height of the DP in case there are no narrow resonances in m23
        Double_t lastResMax23 = m23Min;

        // Loop over narrow resonances in m23
        for ( UInt_t m23i = 0; m23i < m23Regions.size(); m23i++ ) {

            const Double_t resMin23 = m23Regions[m23i].first;
            const Double_t resMax23 = m23Regions[m23i].second;

            // For the first entry, add the area between m23 threshold and this first entry
            if ( m23i == 0 ) {
                integrationRegions.push_back( this->newDPIntegrationRegion( resMin13,
                                                                            resMax13,
                                                                            m23Min,
                                                                            resMin23,
                                                                            m13Binning,
                                                                            defaultBinning,
                                                                            precision,
                                                                            nAmp_,
                                                                            nIncohAmp_ ) );
            }

            // For all entries except the last one, add the area between this and the next entry
            if ( m23i != ( m23Regions.size() - 1 ) ) {
                const Double_t nextResMin23 = m23Regions[m23i + 1].first;
                integrationRegions.push_back( this->newDPIntegrationRegion( resMin13,
                                                                            resMax13,
                                                                            resMax23,
                                                                            nextResMin23,
                                                                            m13Binning,
                                                                            defaultBinning,
                                                                            precision,
                                                                            nAmp_,
                                                                            nIncohAmp_ ) );
            } else {
                lastResMax23 = resMax23;
            }
        }

        // Add the area between the last entry and the maximum m23 (which could be the whole strip if there are no entries in m23Regions)
        integrationRegions.push_back( this->newDPIntegrationRegion( resMin13,
                                                                    resMax13,
                                                                    lastResMax23,
                                                                    m23Max,
                                                                    m13Binning,
                                                                    defaultBinning,
                                                                    precision,
                                                                    nAmp_,
                                                                    nIncohAmp_ ) );
    }

    return integrationRegions;
}

std::vector<LauDPPartialIntegralInfo*> LauIsobarDynamics::m23IntegrationRegions(
    const std::vector<std::pair<Double_t, Double_t>>& m13Regions,
    const std::vector<std::pair<Double_t, Double_t>>& m23Regions,
    const std::vector<Double_t>& m13Binnings,
    const std::vector<Double_t>& m23Binnings,
    const Double_t precision,
    const Double_t defaultBinning ) const
{
    // Create integration regions for all narrow resonances in m23 (including the overlap regions with m13 narrow resonances)
    std::vector<LauDPPartialIntegralInfo*> integrationRegions;

    const Double_t m13Min = kinematics_->getm13Min();
    const Double_t m13Max = kinematics_->getm13Max();

    // Loop over narrow resonances in m23
    for ( UInt_t m23i = 0; m23i < m23Regions.size(); m23i++ ) {

        const Double_t m23Binning = m23Binnings[m23i];

        const Double_t resMin23 = m23Regions[m23i].first;
        const Double_t resMax23 = m23Regions[m23i].second;

        // Initialise to the full width of the DP in case there are no narrow resonances in m13
        Double_t lastResMax13 = m13Min;

        // Loop over narrow resonances in m13
        for ( UInt_t m13i = 0; m13i < m13Regions.size(); m13i++ ) {

            const Double_t m13Binning = m13Binnings[m23i];

            const Double_t resMin13 = m13Regions[m13i].first;
            const Double_t resMax13 = m13Regions[m13i].second;

            // Overlap region (only needed in m23)
            integrationRegions.push_back( this->newDPIntegrationRegion( resMin13,
                                                                        resMax13,
                                                                        resMin23,
                                                                        resMax23,
                                                                        m13Binning,
                                                                        m23Binning,
                                                                        precision,
                                                                        nAmp_,
                                                                        nIncohAmp_ ) );

            // For the first entry, add the area between m13 threshold and this first entry
            if ( m13i == 0 ) {
                integrationRegions.push_back( this->newDPIntegrationRegion( m13Min,
                                                                            resMin13,
                                                                            resMin23,
                                                                            resMax23,
                                                                            defaultBinning,
                                                                            m23Binning,
                                                                            precision,
                                                                            nAmp_,
                                                                            nIncohAmp_ ) );
            }

            // For all entries except the last one, add the area between this and the next entry
            if ( m13i != m13Regions.size() - 1 ) {
                const Double_t nextResMin13 = m23Regions[m13i + 1].first;
                integrationRegions.push_back( this->newDPIntegrationRegion( resMax13,
                                                                            nextResMin13,
                                                                            resMin23,
                                                                            resMax23,
                                                                            defaultBinning,
                                                                            m23Binning,
                                                                            precision,
                                                                            nAmp_,
                                                                            nIncohAmp_ ) );
            } else {
                lastResMax13 = resMax13;
            }
        }

        // Add the area between the last entry and the maximum m13 (which could be the whole strip if there are no entries in m13Regions)
        integrationRegions.push_back( this->newDPIntegrationRegion( lastResMax13,
                                                                    m13Max,
                                                                    resMin23,
                                                                    resMax23,
                                                                    defaultBinning,
                                                                    m23Binning,
                                                                    precision,
                                                                    nAmp_,
                                                                    nIncohAmp_ ) );
    }

    return integrationRegions;
}

LauDPPartialIntegralInfo* LauIsobarDynamics::newDPIntegrationRegion( const Double_t minm13,
                                                                     const Double_t maxm13,
                                                                     const Double_t minm23,
                                                                     const Double_t maxm23,
                                                                     const Double_t m13BinWidth,
                                                                     const Double_t m23BinWidth,
                                                                     const Double_t precision,
                                                                     const UInt_t nAmp,
                                                                     const UInt_t nIncohAmp ) const
{
    const UInt_t nm13Points = static_cast<UInt_t>( ( maxm13 - minm13 ) / m13BinWidth );
    const UInt_t nm23Points = static_cast<UInt_t>( ( maxm23 - minm23 ) / m23BinWidth );

    // If we would create a region with no interior points, just return a null pointer
    if ( nm13Points == 0 || nm23Points == 0 ) {
        return 0;
    }

    return new LauDPPartialIntegralInfo( minm13,
                                         maxm13,
                                         minm23,
                                         maxm23,
                                         m13BinWidth,
                                         m23BinWidth,
                                         precision,
                                         nAmp,
                                         nIncohAmp );
}

void LauIsobarDynamics::calcDPNormalisationScheme()
{
    if ( ! dpPartialIntegralInfo_.empty() ) {
        std::cerr << "ERROR in LauIsobarDynamics::calcDPNormalisationScheme : Scheme already stored!"
                  << std::endl;
        return;
    }

    // The precision for the Gauss-Legendre weights
    const Double_t precision( 1e-6 );

    // Get the rectangle that encloses the DP
    const Double_t minm13 = kinematics_->getm13Min();
    const Double_t maxm13 = kinematics_->getm13Max();
    const Double_t minm23 = kinematics_->getm23Min();
    const Double_t maxm23 = kinematics_->getm23Max();
    const Double_t minm12 = kinematics_->getm12Min();
    const Double_t maxm12 = kinematics_->getm12Max();

    // Find out whether we have narrow resonances in the DP (defined here as width < 20 MeV).
    std::vector<std::pair<Double_t, Double_t>> m13NarrowRes;
    std::vector<std::pair<Double_t, Double_t>> m23NarrowRes;
    std::vector<std::pair<Double_t, Double_t>> m12NarrowRes;

    // Rho-omega mixing models implicitly contains omega(782) model, but width is of rho(770) - handle as a special case
    LauResonanceMaker& resonanceMaker = LauResonanceMaker::get();
    LauResonanceInfo* omega_info      = resonanceMaker.getResInfo( "omega(782)" );
    const Double_t omegaMass = ( omega_info != 0 ) ? omega_info->getMass()->unblindValue() : 0.78265;
    const Double_t omegaWidth = ( omega_info != 0 ) ? omega_info->getWidth()->unblindValue()
                                                    : 0.00849;

    for ( std::vector<LauAbsResonance*>::const_iterator iter = sigResonances_.begin();
          iter != sigResonances_.end();
          ++iter ) {

        LauAbsResonance::LauResonanceModel model = ( *iter )->getResonanceModel();
        const TString& name                      = ( *iter )->getResonanceName();
        Int_t pair                               = ( *iter )->getPairInt();
        Double_t mass                            = ( *iter )->getMass();
        Double_t width                           = ( *iter )->getWidth();

        if ( model == LauAbsResonance::RhoOmegaMix_GS || model == LauAbsResonance::RhoOmegaMix_GS_1 ||
             model == LauAbsResonance::RhoOmegaMix_RBW ||
             model == LauAbsResonance::RhoOmegaMix_RBW_1 ) {
            mass  = omegaMass;
            width = omegaWidth;
        }

        if ( width > narrowWidth_ || width <= 0.0 ) {
            continue;
        }

        std::cout << "INFO in LauIsobarDynamics::calcDPNormalisationScheme : Found narrow resonance: "
                  << name << ", mass = " << mass << ", width = " << width << ", pair int = " << pair
                  << std::endl;
        if ( pair == 1 ) {
            if ( mass < minm23 || mass > maxm23 ) {
                std::cout << std::string( 53, ' ' )
                          << ": But its pole is outside the kinematically allowed range, so will not consider it narrow for the purposes of integration"
                          << std::endl;
            } else {
                m23NarrowRes.push_back( std::make_pair( mass, width ) );
                if ( fullySymmetricDP_ ) {
                    m13NarrowRes.push_back( std::make_pair( mass, width ) );
                    m12NarrowRes.push_back( std::make_pair( mass, width ) );
                } else if ( symmetricalDP_ || ( flavConjDP_ && forceSymmetriseIntegration_ ) ) {
                    m13NarrowRes.push_back( std::make_pair( mass, width ) );
                }
            }
        } else if ( pair == 2 ) {
            if ( mass < minm13 || mass > maxm13 ) {
                std::cout << std::string( 53, ' ' )
                          << ": But its pole is outside the kinematically allowed range, so will not consider it narrow for the purposes of integration"
                          << std::endl;
            } else {
                m13NarrowRes.push_back( std::make_pair( mass, width ) );
                if ( fullySymmetricDP_ ) {
                    m23NarrowRes.push_back( std::make_pair( mass, width ) );
                    m12NarrowRes.push_back( std::make_pair( mass, width ) );
                } else if ( symmetricalDP_ || ( flavConjDP_ && forceSymmetriseIntegration_ ) ) {
                    m23NarrowRes.push_back( std::make_pair( mass, width ) );
                }
            }
        } else if ( pair == 3 ) {
            if ( mass < minm12 || mass > maxm12 ) {
                std::cout << std::string( 53, ' ' )
                          << ": But its pole is outside the kinematically allowed range, so will not consider it narrow for the purposes of integration"
                          << std::endl;
            } else {
                m12NarrowRes.push_back( std::make_pair( mass, width ) );
                if ( fullySymmetricDP_ ) {
                    m13NarrowRes.push_back( std::make_pair( mass, width ) );
                    m12NarrowRes.push_back( std::make_pair( mass, width ) );
                }
            }
        } else {
            std::cerr << "WARNING in LauIsobarDynamics::calcDPNormalisationScheme : strange pair integer, "
                      << pair << ", for resonance \"" << ( *iter )->getResonanceName() << std::endl;
        }
    }

    for ( std::vector<LauAbsIncohRes*>::const_iterator iter = sigIncohResonances_.begin();
          iter != sigIncohResonances_.end();
          ++iter ) {

        const TString& name = ( *iter )->getResonanceName();
        Int_t pair          = ( *iter )->getPairInt();
        Double_t mass       = ( *iter )->getMass();
        Double_t width      = ( *iter )->getWidth();

        if ( width > narrowWidth_ || width == 0.0 ) {
            continue;
        }

        std::cout << "INFO in LauIsobarDynamics::calcDPNormalisationScheme : Found narrow resonance: "
                  << name << ", mass = " << mass << ", width = " << width << ", pair int = " << pair
                  << std::endl;

        if ( pair == 1 ) {
            if ( mass < minm23 || mass > maxm23 ) {
                std::cout << std::string( 53, ' ' )
                          << ": But its pole is outside the kinematically allowed range, so will not consider it narrow for the purposes of integration"
                          << std::endl;
            } else {
                m23NarrowRes.push_back( std::make_pair( mass, width ) );
                if ( fullySymmetricDP_ ) {
                    m13NarrowRes.push_back( std::make_pair( mass, width ) );
                    m12NarrowRes.push_back( std::make_pair( mass, width ) );
                } else if ( symmetricalDP_ || ( flavConjDP_ && forceSymmetriseIntegration_ ) ) {
                    m13NarrowRes.push_back( std::make_pair( mass, width ) );
                }
            }
        } else if ( pair == 2 ) {
            if ( mass < minm13 || mass > maxm13 ) {
                std::cout << std::string( 53, ' ' )
                          << ": But its pole is outside the kinematically allowed range, so will not consider it narrow for the purposes of integration"
                          << std::endl;
            } else {
                m13NarrowRes.push_back( std::make_pair( mass, width ) );
                if ( fullySymmetricDP_ ) {
                    m23NarrowRes.push_back( std::make_pair( mass, width ) );
                    m12NarrowRes.push_back( std::make_pair( mass, width ) );
                } else if ( symmetricalDP_ || ( flavConjDP_ && forceSymmetriseIntegration_ ) ) {
                    m23NarrowRes.push_back( std::make_pair( mass, width ) );
                }
            }
        } else if ( pair == 3 ) {
            if ( mass < minm12 || mass > maxm12 ) {
                std::cout << std::string( 53, ' ' )
                          << ": But its pole is outside the kinematically allowed range, so will not consider it narrow for the purposes of integration"
                          << std::endl;
            } else {
                m12NarrowRes.push_back( std::make_pair( mass, width ) );
                if ( fullySymmetricDP_ ) {
                    m13NarrowRes.push_back( std::make_pair( mass, width ) );
                    m12NarrowRes.push_back( std::make_pair( mass, width ) );
                }
            }
        } else {
            std::cerr << "WARNING in LauIsobarDynamics::calcDPNormalisationScheme : strange pair integer, "
                      << pair << ", for resonance \"" << ( *iter )->getResonanceName() << std::endl;
        }
    }

    // Depending on how many narrow resonances we have and where they are
    // we adopt different approaches
    if ( ! m12NarrowRes.empty() ) {

        // We have at least one narrow resonance in m12
        // Switch to using the square DP for the integration
        // TODO - for the time being just use a single, reasonably fine by default and tunable, grid
        //      - can later consider whether there's a need to split up the mPrime axis into regions around particularly narrow resonances in m12
        //      - but it seems that this isn't really needed since even the default tune gives a good resolution for most narrow resonances such as phi / chi_c0
        std::cout << "INFO in LauIsobarDynamics::calcDPNormalisationScheme : One or more narrow resonances found in m12, integrating over whole square Dalitz plot with bin widths of "
                  << mPrimeBinWidth_ << " in mPrime and " << thPrimeBinWidth_ << " in thetaPrime..."
                  << std::endl;

        // Make sure that the kinematics will calculate the square DP co-ordinates
        if ( ! kinematics_->squareDP() ) {
            std::cerr << "WARNING in LauIsobarDynamics::calcDPNormalisationScheme  : forcing kinematics to calculate the required square DP co-ordinates"
                      << std::endl;
            kinematics_->squareDP( kTRUE );
        }

        dpPartialIntegralInfo_.push_back( new LauDPPartialIntegralInfo( 0.0,
                                                                        1.0,
                                                                        0.0,
                                                                        1.0,
                                                                        mPrimeBinWidth_,
                                                                        thPrimeBinWidth_,
                                                                        precision,
                                                                        nAmp_,
                                                                        nIncohAmp_,
                                                                        kTRUE,
                                                                        kinematics_ ) );

    } else if ( m13NarrowRes.empty() && m23NarrowRes.empty() ) {

        // There are no narrow resonances, so we just do a single grid over the whole DP
        std::cout << "INFO in LauIsobarDynamics::calcDPNormalisationScheme : No narrow resonances found, integrating over whole Dalitz plot..."
                  << std::endl;
        dpPartialIntegralInfo_.push_back( new LauDPPartialIntegralInfo( minm13,
                                                                        maxm13,
                                                                        minm23,
                                                                        maxm23,
                                                                        m13BinWidth_,
                                                                        m23BinWidth_,
                                                                        precision,
                                                                        nAmp_,
                                                                        nIncohAmp_ ) );

    } else {

        // Get regions in that correspond to narrow resonances in m13 and m23, and correct for overlaps in each dimension (to use the finest binning)

        // Sort resonances by ascending mass to calculate regions properly
        std::sort( m13NarrowRes.begin(), m13NarrowRes.end() );
        std::sort( m23NarrowRes.begin(), m23NarrowRes.end() );

        // For each narrow resonance in m13, determine the corresponding window and its binning
        std::vector<std::pair<Double_t, Double_t>> m13Regions;
        std::vector<Double_t> m13Binnings;

        for ( std::vector<std::pair<Double_t, Double_t>>::const_iterator iter = m13NarrowRes.begin();
              iter != m13NarrowRes.end();
              ++iter ) {
            Double_t mass  = iter->first;
            Double_t width = iter->second;

            Double_t regionBegin = mass - 5.0 * width;
            Double_t regionEnd   = mass + 5.0 * width;
            Double_t binning     = width / binningFactor_;

            // check if we ought to extend the region to the edge of the phase space (in either direction)
            if ( regionBegin < ( minm13 + 50.0 * m13BinWidth_ ) ) {
                std::cout << "INFO in LauIsobarDynamics::calcDPNormalisationScheme : Resonance at m13 = "
                          << mass << " is close to threshold, extending integration region"
                          << std::endl;
                regionBegin = minm13;
            }
            if ( regionEnd > ( maxm13 - 50.0 * m13BinWidth_ ) ) {
                std::cout << "INFO in LauIsobarDynamics::calcDPNormalisationScheme : Resonance at m13 = "
                          << mass
                          << " is close to upper edge of phase space, extending integration region"
                          << std::endl;
                regionEnd = maxm13;
            }

            m13Regions.push_back( std::make_pair( regionBegin, regionEnd ) );
            m13Binnings.push_back( binning );
        }

        // For each narrow resonance in m23, determine the corresponding window and its binning
        std::vector<std::pair<Double_t, Double_t>> m23Regions;
        std::vector<Double_t> m23Binnings;

        for ( std::vector<std::pair<Double_t, Double_t>>::const_iterator iter = m23NarrowRes.begin();
              iter != m23NarrowRes.end();
              ++iter ) {
            Double_t mass  = iter->first;
            Double_t width = iter->second;

            Double_t regionBegin = mass - 5.0 * width;
            Double_t regionEnd   = mass + 5.0 * width;
            Double_t binning     = width / binningFactor_;

            // check if we ought to extend the region to the edge of the phase space (in either direction)
            if ( regionBegin < ( minm23 + 50.0 * m23BinWidth_ ) ) {
                std::cout << "INFO in LauIsobarDynamics::calcDPNormalisationScheme : Resonance at m23 = "
                          << mass << " is close to threshold, extending integration region"
                          << std::endl;
                regionBegin = minm23;
            }
            if ( regionEnd > ( maxm23 - 50.0 * m23BinWidth_ ) ) {
                std::cout << "INFO in LauIsobarDynamics::calcDPNormalisationScheme : Resonance at m23 = "
                          << mass
                          << " is close to upper edge of phase space, extending integration region"
                          << std::endl;
                regionEnd = maxm23;
            }

            m23Regions.push_back( std::make_pair( regionBegin, regionEnd ) );
            m23Binnings.push_back( binning );
        }

        // Sort out overlaps between regions in the same mass pairing
        this->correctDPOverlap( m13Regions, m13Binnings );
        this->correctDPOverlap( m23Regions, m23Binnings );

        // Get the narrow resonance regions plus any overlap region
        std::vector<LauDPPartialIntegralInfo*> fineScheme13 =
            this->m13IntegrationRegions( m13Regions, m23Regions, m13Binnings, precision, m13BinWidth_ );
        std::vector<LauDPPartialIntegralInfo*> fineScheme23 =
            this->m23IntegrationRegions( m13Regions,
                                         m23Regions,
                                         m13Binnings,
                                         m23Binnings,
                                         precision,
                                         m23BinWidth_ );

        // Get coarse regions by calculating the gaps between the
        // narrow resonances and using the same functions to create
        // the integration grid object for each
        std::vector<std::pair<Double_t, Double_t>> coarseRegions =
            this->formGapsFromRegions( m13Regions, minm13, maxm13 );
        std::vector<Double_t> coarseBinning( fineScheme13.size() + 1, m13BinWidth_ );
        std::vector<LauDPPartialIntegralInfo*> coarseScheme =
            this->m13IntegrationRegions( coarseRegions,
                                         m23Regions,
                                         coarseBinning,
                                         precision,
                                         m13BinWidth_ );

        dpPartialIntegralInfo_.insert( dpPartialIntegralInfo_.end(),
                                       fineScheme13.begin(),
                                       fineScheme13.end() );
        dpPartialIntegralInfo_.insert( dpPartialIntegralInfo_.end(),
                                       fineScheme23.begin(),
                                       fineScheme23.end() );
        dpPartialIntegralInfo_.insert( dpPartialIntegralInfo_.end(),
                                       coarseScheme.begin(),
                                       coarseScheme.end() );

        // Remove any null pointer entries in the integral list
        // (that are produced when an integration region with no
        // interior points is defined)
        this->cullNullRegions( dpPartialIntegralInfo_ );
    }

    normalizationSchemeDone_ = kTRUE;
}

void LauIsobarDynamics::setIntegralBinWidths( const Double_t m13BinWidth,
                                              const Double_t m23BinWidth,
                                              const Double_t mPrimeBinWidth,
                                              const Double_t thPrimeBinWidth )
{
    // Set the bin widths for the m13 vs m23 integration grid
    m13BinWidth_ = m13BinWidth;
    m23BinWidth_ = m23BinWidth;

    // Set the bin widths for the m' vs theta' integration grid
    mPrimeBinWidth_  = mPrimeBinWidth;
    thPrimeBinWidth_ = thPrimeBinWidth;
}

void LauIsobarDynamics::calcDPPartialIntegral( LauDPPartialIntegralInfo* intInfo )
{
    // Calculate the integrals for all parts of the amplitude in the given region of the DP

    const Bool_t squareDP   = intInfo->getSquareDP();
    const UInt_t nm13Points = intInfo->getnm13Points();
    const UInt_t nm23Points = intInfo->getnm23Points();

    //Double_t dpArea(0.0);
    for ( UInt_t i = 0; i < nm13Points; ++i ) {

        const Double_t m13   = intInfo->getM13Value( i );
        const Double_t m13Sq = m13 * m13;

        for ( UInt_t j = 0; j < nm23Points; ++j ) {

            const Double_t m23    = intInfo->getM23Value( j );
            const Double_t m23Sq  = m23 * m23;
            const Double_t weight = intInfo->getWeight( i, j );

            // Calculate the integral contributions for each resonance.
            // Only points within the DP area contribute.
            // This also calculates the total DP area as a check.
            // NB if squareDP is true, m13 and m23 are actually mPrime and thetaPrime
            Bool_t withinDP = squareDP ? kinematics_->withinSqDPLimits( m13, m23 )
                                       : kinematics_->withinDPLimits( m13Sq, m23Sq );
            if ( withinDP == kTRUE ) {

                if ( squareDP ) {
                    // NB m13 and m23 are actually mPrime and thetaPrime
                    kinematics_->updateSqDPKinematics( m13, m23 );
                } else {
                    kinematics_->updateKinematics( m13Sq, m23Sq );
                }

                this->calculateAmplitudes( intInfo, i, j );

                this->addGridPointToIntegrals( weight );

                // Increment total DP area
                //dpArea += weight;
            }

        }    // j weights loop
    }        // i weights loop

    // Print out DP area to check whether we have a sensible output
    //std::cout<<"                                                 : dpArea = "<<dpArea<<std::endl;
}

void LauIsobarDynamics::calculateAmplitudes( LauDPPartialIntegralInfo* intInfo,
                                             const UInt_t m13Point,
                                             const UInt_t m23Point )
{
    const std::set<UInt_t>::const_iterator intEnd = integralsToBeCalculated_.end();

    for ( UInt_t iAmp = 0; iAmp < nAmp_; ++iAmp ) {

        if ( integralsToBeCalculated_.find( iAmp ) != intEnd ) {
            // Calculate the dynamics for this resonance
            ff_[iAmp] = this->resAmp( iAmp );

            // Store the new value in the integration info object
            intInfo->storeAmplitude( m13Point, m23Point, iAmp, ff_[iAmp] );
        } else {
            // Retrieve the cached value of the amplitude
            ff_[iAmp] = intInfo->getAmplitude( m13Point, m23Point, iAmp );
        }
    }

    for ( UInt_t iAmp = 0; iAmp < nIncohAmp_; ++iAmp ) {

        if ( integralsToBeCalculated_.find( iAmp + nAmp_ ) != intEnd ) {
            // Calculate the dynamics for this resonance
            incohInten_[iAmp] = this->incohResAmp( iAmp );

            // Store the new value in the integration info object
            intInfo->storeIntensity( m13Point, m23Point, iAmp, incohInten_[iAmp] );
        } else {
            // Retrieve the cached value of the amplitude
            incohInten_[iAmp] = intInfo->getIntensity( m13Point, m23Point, iAmp );
        }
    }

    // If symmetric, do as above with flipped kinematics and add to amplitude
    // (No need to retrive the cache if this was done in the first case)

    if ( symmetricalDP_ == kTRUE ) {
        kinematics_->flipAndUpdateKinematics();

        for ( UInt_t iAmp = 0; iAmp < nAmp_; ++iAmp ) {

            if ( ( integralsToBeCalculated_.find( iAmp ) != intEnd ) &&
                 ! sigResonances_[iAmp]->preSymmetrised() ) {
                // Calculate the dynamics for this resonance
                ff_[iAmp] += this->resAmp( iAmp );

                // Store the new value in the integration info object
                intInfo->storeAmplitude( m13Point, m23Point, iAmp, ff_[iAmp] );
            }
        }

        for ( UInt_t iAmp = 0; iAmp < nIncohAmp_; ++iAmp ) {

            if ( ( integralsToBeCalculated_.find( iAmp + nAmp_ ) != intEnd ) &&
                 ! sigResonances_[iAmp]->preSymmetrised() ) {
                // Calculate the dynamics for this resonance
                incohInten_[iAmp] += this->incohResAmp( iAmp );

                // Store the new value in the integration info object
                intInfo->storeIntensity( m13Point, m23Point, iAmp, incohInten_[iAmp] );
            }
        }

        kinematics_->flipAndUpdateKinematics();
    }

    if ( fullySymmetricDP_ == kTRUE ) {
        // rotate and evaluate
        kinematics_->rotateAndUpdateKinematics();
        for ( UInt_t iAmp = 0; iAmp < nAmp_; ++iAmp ) {
            if ( ( integralsToBeCalculated_.find( iAmp ) != intEnd ) &&
                 ! sigResonances_[iAmp]->preSymmetrised() ) {
                // Calculate the dynamics for this resonance
                ff_[iAmp] += this->resAmp( iAmp );
                // Store the new value in the integration info object
                intInfo->storeAmplitude( m13Point, m23Point, iAmp, ff_[iAmp] );
            }
        }
        for ( UInt_t iAmp = 0; iAmp < nIncohAmp_; ++iAmp ) {
            if ( ( integralsToBeCalculated_.find( iAmp + nAmp_ ) != intEnd ) &&
                 ! sigResonances_[iAmp]->preSymmetrised() ) {
                // Calculate the dynamics for this resonance
                incohInten_[iAmp] += this->incohResAmp( iAmp );
                // Store the new value in the integration info object
                intInfo->storeIntensity( m13Point, m23Point, iAmp, incohInten_[iAmp] );
            }
        }

        // rotate and evaluate
        kinematics_->rotateAndUpdateKinematics();
        for ( UInt_t iAmp = 0; iAmp < nAmp_; ++iAmp ) {
            if ( ( integralsToBeCalculated_.find( iAmp ) != intEnd ) &&
                 ! sigResonances_[iAmp]->preSymmetrised() ) {
                // Calculate the dynamics for this resonance
                ff_[iAmp] += this->resAmp( iAmp );
                // Store the new value in the integration info object
                intInfo->storeAmplitude( m13Point, m23Point, iAmp, ff_[iAmp] );
            }
        }
        for ( UInt_t iAmp = 0; iAmp < nIncohAmp_; ++iAmp ) {
            if ( ( integralsToBeCalculated_.find( iAmp + nAmp_ ) != intEnd ) &&
                 ! sigResonances_[iAmp]->preSymmetrised() ) {
                // Calculate the dynamics for this resonance
                incohInten_[iAmp] += this->incohResAmp( iAmp );
                // Store the new value in the integration info object
                intInfo->storeIntensity( m13Point, m23Point, iAmp, incohInten_[iAmp] );
            }
        }

        // rotate, flip and evaluate
        kinematics_->rotateAndUpdateKinematics();
        kinematics_->flipAndUpdateKinematics();
        for ( UInt_t iAmp = 0; iAmp < nAmp_; ++iAmp ) {
            if ( ( integralsToBeCalculated_.find( iAmp ) != intEnd ) &&
                 ! sigResonances_[iAmp]->preSymmetrised() ) {
                // Calculate the dynamics for this resonance
                ff_[iAmp] += this->resAmp( iAmp );
                // Store the new value in the integration info object
                intInfo->storeAmplitude( m13Point, m23Point, iAmp, ff_[iAmp] );
            }
        }
        for ( UInt_t iAmp = 0; iAmp < nIncohAmp_; ++iAmp ) {
            if ( ( integralsToBeCalculated_.find( iAmp + nAmp_ ) != intEnd ) &&
                 ! sigResonances_[iAmp]->preSymmetrised() ) {
                // Calculate the dynamics for this resonance
                incohInten_[iAmp] += this->incohResAmp( iAmp );
                // Store the new value in the integration info object
                intInfo->storeIntensity( m13Point, m23Point, iAmp, incohInten_[iAmp] );
            }
        }

        // rotate and evaluate
        kinematics_->rotateAndUpdateKinematics();
        for ( UInt_t iAmp = 0; iAmp < nAmp_; ++iAmp ) {
            if ( ( integralsToBeCalculated_.find( iAmp ) != intEnd ) &&
                 ! sigResonances_[iAmp]->preSymmetrised() ) {
                // Calculate the dynamics for this resonance
                ff_[iAmp] += this->resAmp( iAmp );
                // Store the new value in the integration info object
                intInfo->storeAmplitude( m13Point, m23Point, iAmp, ff_[iAmp] );
            }
        }
        for ( UInt_t iAmp = 0; iAmp < nIncohAmp_; ++iAmp ) {
            if ( ( integralsToBeCalculated_.find( iAmp + nAmp_ ) != intEnd ) &&
                 ! sigResonances_[iAmp]->preSymmetrised() ) {
                // Calculate the dynamics for this resonance
                incohInten_[iAmp] += this->incohResAmp( iAmp );
                // Store the new value in the integration info object
                intInfo->storeIntensity( m13Point, m23Point, iAmp, incohInten_[iAmp] );
            }
        }

        // rotate and evaluate
        kinematics_->rotateAndUpdateKinematics();
        for ( UInt_t iAmp = 0; iAmp < nAmp_; ++iAmp ) {
            if ( ( integralsToBeCalculated_.find( iAmp ) != intEnd ) &&
                 ! sigResonances_[iAmp]->preSymmetrised() ) {
                // Calculate the dynamics for this resonance
                ff_[iAmp] += this->resAmp( iAmp );
                // Store the new value in the integration info object
                intInfo->storeAmplitude( m13Point, m23Point, iAmp, ff_[iAmp] );
            }
        }
        for ( UInt_t iAmp = 0; iAmp < nIncohAmp_; ++iAmp ) {
            if ( ( integralsToBeCalculated_.find( iAmp + nAmp_ ) != intEnd ) &&
                 ! sigResonances_[iAmp]->preSymmetrised() ) {
                // Calculate the dynamics for this resonance
                incohInten_[iAmp] += this->incohResAmp( iAmp );
                // Store the new value in the integration info object
                intInfo->storeIntensity( m13Point, m23Point, iAmp, incohInten_[iAmp] );
            }
        }

        // rotate and flip to get us back to where we started
        kinematics_->rotateAndUpdateKinematics();
        kinematics_->flipAndUpdateKinematics();
    }

    // If we haven't cached the data, then we need to find out the efficiency.
    eff_ = this->retrieveEfficiency();

    intInfo->storeEfficiency( m13Point, m23Point, eff_ );
}

void LauIsobarDynamics::calculateAmplitudes()
{
    std::set<UInt_t>::const_iterator iter         = integralsToBeCalculated_.begin();
    const std::set<UInt_t>::const_iterator intEnd = integralsToBeCalculated_.end();

    for ( iter = integralsToBeCalculated_.begin(); iter != intEnd; ++iter ) {

        // Calculate the dynamics for this resonance
        if ( *iter < nAmp_ ) {
            ff_[*iter] = this->resAmp( *iter );
        } else {
            incohInten_[*iter - nAmp_] = this->incohResAmp( *iter - nAmp_ );
        }
    }

    if ( symmetricalDP_ == kTRUE ) {
        kinematics_->flipAndUpdateKinematics();

        for ( iter = integralsToBeCalculated_.begin(); iter != intEnd; ++iter ) {

            // Calculate the dynamics for this resonance
            if ( *iter < nAmp_ && ! sigResonances_[*iter]->preSymmetrised() ) {
                ff_[*iter] += this->resAmp( *iter );
            } else if ( *iter >= nAmp_ && ! sigResonances_[*iter - nAmp_]->preSymmetrised() ) {
                incohInten_[*iter - nAmp_] += this->incohResAmp( *iter - nAmp_ );
            }
        }

        kinematics_->flipAndUpdateKinematics();
    }

    if ( fullySymmetricDP_ == kTRUE ) {
        // rotate and evaluate
        kinematics_->rotateAndUpdateKinematics();
        for ( iter = integralsToBeCalculated_.begin(); iter != intEnd; ++iter ) {
            if ( *iter < nAmp_ && ! sigResonances_[*iter]->preSymmetrised() ) {
                ff_[*iter] += this->resAmp( *iter );
            } else if ( *iter >= nAmp_ && ! sigResonances_[*iter - nAmp_]->preSymmetrised() ) {
                incohInten_[*iter - nAmp_] += this->incohResAmp( *iter - nAmp_ );
            }
        }

        // rotate and evaluate
        kinematics_->rotateAndUpdateKinematics();
        for ( iter = integralsToBeCalculated_.begin(); iter != intEnd; ++iter ) {
            if ( *iter < nAmp_ && ! sigResonances_[*iter]->preSymmetrised() ) {
                ff_[*iter] += this->resAmp( *iter );
            } else if ( *iter >= nAmp_ && ! sigResonances_[*iter - nAmp_]->preSymmetrised() ) {
                incohInten_[*iter - nAmp_] += this->incohResAmp( *iter - nAmp_ );
            }
        }

        // rotate, flip and evaluate
        kinematics_->rotateAndUpdateKinematics();
        kinematics_->flipAndUpdateKinematics();
        for ( iter = integralsToBeCalculated_.begin(); iter != intEnd; ++iter ) {
            if ( *iter < nAmp_ && ! sigResonances_[*iter]->preSymmetrised() ) {
                ff_[*iter] += this->resAmp( *iter );
            } else if ( *iter >= nAmp_ && ! sigResonances_[*iter - nAmp_]->preSymmetrised() ) {
                incohInten_[*iter - nAmp_] += this->incohResAmp( *iter - nAmp_ );
            }
        }

        // rotate and evaluate
        kinematics_->rotateAndUpdateKinematics();
        for ( iter = integralsToBeCalculated_.begin(); iter != intEnd; ++iter ) {
            if ( *iter < nAmp_ && ! sigResonances_[*iter]->preSymmetrised() ) {
                ff_[*iter] += this->resAmp( *iter );
            } else if ( *iter >= nAmp_ && ! sigResonances_[*iter - nAmp_]->preSymmetrised() ) {
                incohInten_[*iter - nAmp_] += this->incohResAmp( *iter - nAmp_ );
            }
        }

        // rotate and evaluate
        kinematics_->rotateAndUpdateKinematics();
        for ( iter = integralsToBeCalculated_.begin(); iter != intEnd; ++iter ) {
            if ( *iter < nAmp_ && ! sigResonances_[*iter]->preSymmetrised() ) {
                ff_[*iter] += this->resAmp( *iter );
            } else if ( *iter >= nAmp_ && ! sigResonances_[*iter - nAmp_]->preSymmetrised() ) {
                incohInten_[*iter - nAmp_] += this->incohResAmp( *iter - nAmp_ );
            }
        }

        // rotate and flip to get us back to where we started
        kinematics_->rotateAndUpdateKinematics();
        kinematics_->flipAndUpdateKinematics();
    }

    // If we haven't cached the data, then we need to find out the efficiency.
    eff_ = this->retrieveEfficiency();
}

void LauIsobarDynamics::calcTotalAmp( const Bool_t useEff )
{
    // Reset the total amplitude to zero
    totAmp_.zero();

    // Loop over all signal amplitudes
    LauComplex ATerm;
    for ( UInt_t i = 0; i < nAmp_; ++i ) {

        // Get the partial complex amplitude - (mag, phase)*(resonance dynamics)
        ATerm = Amp_[i] * ff_[i];

        // Scale this contribution by its relative normalisation w.r.t. the whole dynamics
        ATerm.rescale( fNorm_[i] );

        // Add this partial amplitude to the sum
        totAmp_ += ATerm;

    }    // Loop over amplitudes

    // |Sum of partial amplitudes|^2
    ASq_ = totAmp_.abs2();

    for ( UInt_t i = 0; i < nIncohAmp_; ++i ) {

        // Get the partial complex amplitude - (mag, phase)
        ATerm = Amp_[i + nAmp_];

        // Scale this contribution by its relative normalisation w.r.t. the whole dynamics
        ATerm.rescale( fNorm_[i + nAmp_] );

        // Add this partial amplitude to the sum
        ASq_ += ATerm.abs2() * incohInten_[i];
    }

    // Apply the efficiency correction for this event.
    // Multiply the amplitude squared sum by the DP efficiency
    if ( useEff ) {
        ASq_ *= eff_;
    }
}

void LauIsobarDynamics::addGridPointToIntegrals( const Double_t weight )
{
    // Combine the Gauss-Legendre weight with the efficiency
    const Double_t effWeight = eff_ * weight;

    LauComplex fifjEffSumTerm;
    LauComplex fifjSumTerm;

    // Calculates the half-matrix of amplitude-squared and interference
    // terms (dynamical part only)
    // Add the values at this point on the integration grid to the sums
    // (one weighted only by the integration weights, one also weighted by
    // the efficiency)
    for ( UInt_t i = 0; i < nAmp_; ++i ) {

        // Add the dynamical amplitude squared for this resonance.
        Double_t fSqVal  = ff_[i].abs2();
        fSqSum_[i]      += fSqVal * weight;
        fSqEffSum_[i]   += fSqVal * effWeight;

        for ( UInt_t j = i; j < nAmp_; ++j ) {

            fifjEffSumTerm = fifjSumTerm = ff_[i] * ff_[j].conj();

            fifjEffSumTerm.rescale( effWeight );
            fifjEffSum_[i][j] += fifjEffSumTerm;

            fifjSumTerm.rescale( weight );
            fifjSum_[i][j] += fifjSumTerm;
        }
    }
    for ( UInt_t i = 0; i < nIncohAmp_; ++i ) {

        // Add the dynamical amplitude squared for this resonance.
        Double_t fSqVal        = incohInten_[i];
        fSqSum_[i + nAmp_]    += fSqVal * weight;
        fSqEffSum_[i + nAmp_] += fSqVal * effWeight;
    }
}

LauComplex LauIsobarDynamics::resAmp( const UInt_t index )
{
    // Routine to calculate the resonance dynamics (amplitude)
    // using the appropriate Breit-Wigner/Form Factors.

    LauComplex amp = LauComplex( 0.0, 0.0 );

    if ( index >= nAmp_ ) {
        std::cerr << "ERROR in LauIsobarDynamics::resAmp : index = " << index
                  << " is not within the range 0 to " << nAmp_ - 1 << std::endl;
        return amp;
    }

    // Get the signal resonance from the stored vector
    LauAbsResonance* sigResonance = sigResonances_[index];

    if ( sigResonance == 0 ) {
        std::cerr << "ERROR in LauIsobarDynamics::resAmp : Couldn't retrieve valid resonance with index = "
                  << index << std::endl;
        return amp;
    }

    amp = sigResonance->amplitude( kinematics_ );

    return amp;
}

Double_t LauIsobarDynamics::incohResAmp( const UInt_t index )
{
    // Routine to calculate the resonance dynamics (amplitude)
    // using the appropriate Breit-Wigner/Form Factors.

    Double_t intensity = 0.;

    if ( index >= nIncohAmp_ ) {
        std::cerr << "ERROR in LauIsobarDynamics::incohResAmp : index = " << index
                  << " is not within the range 0 to " << nIncohAmp_ - 1 << std::endl;
        return intensity;
    }

    // Get the signal resonance from the stored vector
    LauAbsIncohRes* sigResonance = sigIncohResonances_[index];

    if ( sigResonance == 0 ) {
        std::cerr << "ERROR in LauIsobarDynamics::incohResAmp : Couldn't retrieve valid incoherent resonance with index = "
                  << index << std::endl;
        return intensity;
    }

    LauComplex ff = sigResonance->amplitude( kinematics_ );
    intensity     = sigResonance->intensityFactor( kinematics_ ) * ff.abs2();

    return intensity;
}

void LauIsobarDynamics::setFFTerm( const UInt_t index, const Double_t realPart, const Double_t imagPart )
{
    // Function to set the internal ff term (normally calculated using resAmp(index).
    if ( index >= nAmp_ ) {
        std::cerr << "ERROR in LauIsobarDynamics::setFFTerm : index = " << index
                  << " is not within the range 0 to " << nAmp_ - 1 << std::endl;
        return;
    }

    ff_[index].setRealImagPart( realPart, imagPart );
}

void LauIsobarDynamics::setIncohIntenTerm( const UInt_t index, const Double_t value )
{
    // Function to set the internal incoherent intensity term (normally calculated using incohResAmp(index).
    if ( index >= nIncohAmp_ ) {
        std::cerr << "ERROR in LauIsobarDynamics::setFFTerm : index = " << index
                  << " is not within the range 0 to " << nIncohAmp_ - 1 << std::endl;
        return;
    }

    incohInten_[index] = value;
}

void LauIsobarDynamics::calcExtraInfo( const Bool_t init )
{
    // This method calculates the fit fractions, mean efficiency and total DP rate

    Double_t fifjEffTot( 0.0 ), fifjTot( 0.0 );
    UInt_t i, j;
    for ( i = 0; i < nAmp_; i++ ) {

        // Calculate the diagonal terms
        TString name  = "A";
        name         += i;
        name         += "Sq_FitFrac";
        fitFrac_[i][i].name( name );

        name += "EffUnCorr";
        fitFracEffUnCorr_[i][i].name( name );

        Double_t fifjSumReal  = fifjSum_[i][i].re();
        Double_t sumTerm      = Amp_[i].abs2() * fifjSumReal * fNorm_[i] * fNorm_[i];
        fifjTot              += sumTerm;

        Double_t fifjEffSumReal  = fifjEffSum_[i][i].re();
        Double_t sumEffTerm      = Amp_[i].abs2() * fifjEffSumReal * fNorm_[i] * fNorm_[i];
        fifjEffTot              += sumEffTerm;

        fitFrac_[i][i].value( sumTerm );
        fitFracEffUnCorr_[i][i].value( sumEffTerm );
    }

    for ( i = 0; i < nAmp_; i++ ) {
        for ( j = i + 1; j < nAmp_ + nIncohAmp_; j++ ) {
            // Calculate the cross-terms
            TString name  = "A";
            name         += i;
            name         += "A";
            name         += j;
            name         += "_FitFrac";
            fitFrac_[i][j].name( name );

            name += "EffUnCorr";
            fitFracEffUnCorr_[i][j].name( name );

            if ( j >= nAmp_ ) {
                //Set off-diagonal incoherent terms to zero
                fitFrac_[i][j].value( 0. );
                fitFracEffUnCorr_[i][j].value( 0. );
                continue;
            }

            LauComplex AmpjConj = Amp_[j].conj();
            LauComplex AmpTerm  = Amp_[i] * AmpjConj;

            Double_t crossTerm  = 2.0 * ( AmpTerm * fifjSum_[i][j] ).re() * fNorm_[i] * fNorm_[j];
            fifjTot            += crossTerm;

            Double_t crossEffTerm = 2.0 * ( AmpTerm * fifjEffSum_[i][j] ).re() * fNorm_[i] *
                                    fNorm_[j];
            fifjEffTot += crossEffTerm;

            fitFrac_[i][j].value( crossTerm );
            fitFracEffUnCorr_[i][j].value( crossEffTerm );
        }
    }

    for ( i = nAmp_; i < nAmp_ + nIncohAmp_; i++ ) {

        // Calculate the incoherent terms
        TString name  = "A";
        name         += i;
        name         += "Sq_FitFrac";
        fitFrac_[i][i].name( name );

        name += "EffUnCorr";
        fitFracEffUnCorr_[i][i].name( name );

        Double_t sumTerm  = Amp_[i].abs2() * fSqSum_[i] * fNorm_[i] * fNorm_[i];
        fifjTot          += sumTerm;

        Double_t sumEffTerm  = Amp_[i].abs2() * fSqEffSum_[i] * fNorm_[i] * fNorm_[i];
        fifjEffTot          += sumEffTerm;

        fitFrac_[i][i].value( sumTerm );
        fitFracEffUnCorr_[i][i].value( sumEffTerm );
    }

    for ( i = nAmp_; i < nAmp_ + nIncohAmp_; i++ ) {
        for ( j = i + 1; j < nAmp_ + nIncohAmp_; j++ ) {
            //Set off-diagonal incoherent terms to zero
            TString name  = "A";
            name         += i;
            name         += "A";
            name         += j;
            name         += "_FitFrac";
            fitFrac_[i][j].name( name );

            name += "EffUnCorr";
            fitFracEffUnCorr_[i][j].name( name );

            fitFrac_[i][j].value( 0. );
            fitFracEffUnCorr_[i][j].value( 0. );
        }
    }

    if ( TMath::Abs( fifjTot ) > 1e-10 ) {
        meanDPEff_.value( fifjEffTot / fifjTot );
        if ( init ) {
            meanDPEff_.genValue( meanDPEff_.value() );
            meanDPEff_.initValue( meanDPEff_.value() );
        }
    }
    DPRate_.value( fifjTot );
    if ( init ) {
        DPRate_.genValue( DPRate_.value() );
        DPRate_.initValue( DPRate_.value() );
    }

    // Now divide the fitFraction sums by the overall integral
    for ( i = 0; i < nAmp_ + nIncohAmp_; i++ ) {
        for ( j = i; j < nAmp_ + nIncohAmp_; j++ ) {
            // Get the actual fractions by dividing by the total DP rate
            Double_t fitFracVal  = fitFrac_[i][j].value();
            fitFracVal          /= fifjTot;
            fitFrac_[i][j].value( fitFracVal );

            Double_t fitFracEffUnCorrVal  = fitFracEffUnCorr_[i][j].value();
            fitFracEffUnCorrVal          /= fifjEffTot;
            fitFracEffUnCorr_[i][j].value( fitFracEffUnCorrVal );

            if ( init ) {
                fitFrac_[i][j].genValue( fitFrac_[i][j].value() );
                fitFrac_[i][j].initValue( fitFrac_[i][j].value() );
                fitFracEffUnCorr_[i][j].genValue( fitFracEffUnCorr_[i][j].value() );
                fitFracEffUnCorr_[i][j].initValue( fitFracEffUnCorr_[i][j].value() );
            }
        }
    }

    // Work out total fit fraction over all K-matrix components (for each propagator)
    KMPropMap::iterator mapIter;
    Int_t propInt( 0 );

    for ( mapIter = kMatrixPropagators_.begin(); mapIter != kMatrixPropagators_.end(); ++mapIter ) {

        LauKMatrixPropagator* thePropagator = mapIter->second;

        TString propName = thePropagator->getName();

        // Now loop over all resonances and find those which are K-matrix components for this propagator
        Double_t kMatrixTotFitFrac( 0.0 );

        for ( i = 0; i < nAmp_; i++ ) {

            Bool_t gotKMRes1 = this->gotKMatrixMatch( i, propName );
            if ( gotKMRes1 == kFALSE ) {
                continue;
            }

            Double_t fifjSumReal = fifjSum_[i][i].re();
            Double_t sumTerm     = Amp_[i].abs2() * fifjSumReal * fNorm_[i] * fNorm_[i];

            //Double_t fifjEffSumReal = fifjEffSum_[i][i].re();
            //Double_t sumEffTerm = Amp_[i].abs2()*fifjEffSumReal*fNorm_[i]*fNorm_[i];

            kMatrixTotFitFrac += sumTerm;

            for ( j = i + 1; j < nAmp_; j++ ) {

                Bool_t gotKMRes2 = this->gotKMatrixMatch( j, propName );
                if ( gotKMRes2 == kFALSE ) {
                    continue;
                }

                LauComplex AmpjConj = Amp_[j].conj();
                LauComplex AmpTerm  = Amp_[i] * AmpjConj;

                Double_t crossTerm = 2.0 * ( AmpTerm * fifjSum_[i][j] ).re() * fNorm_[i] * fNorm_[j];
                //Double_t crossEffTerm = 2.0*(AmpTerm*fifjEffSum_[i][j]).re()*fNorm_[i]*fNorm_[j];

                kMatrixTotFitFrac += crossTerm;
            }
        }

        kMatrixTotFitFrac /= fifjTot;

        TString parName( "KMatrixTotFF_" );
        parName += propInt;
        extraParameters_[propInt].name( parName );
        extraParameters_[propInt].value( kMatrixTotFitFrac );
        if ( init ) {
            extraParameters_[propInt].genValue( kMatrixTotFitFrac );
            extraParameters_[propInt].initValue( kMatrixTotFitFrac );
        }

        std::cout << "INFO in LauIsobarDynamics::calcExtraInfo : Total K-matrix fit fraction for propagator "
                  << propName << " is " << kMatrixTotFitFrac << std::endl;

        ++propInt;
    }

    // Calculate rho and omega fit fractions from LauRhoOmegaMix

    // In principle the presence of LauRhoOmegaMix can be detected, however
    // working out the index of the correct component is hard.

    // To prevent unwanted effects, this has to be manually turned on, with the
    // assumption that the first two components are rho and rho_COPY in a LauCPFitModel

    if ( this->calculateRhoOmegaFitFractions_ && ! init ) {

        int omegaID = 0;
        int storeID = 1;

        // Check which B flavour (and therefore which rho_COPY we are) by whether the FF is non-zero
        // Only for CP fit though - for a 'simple' fit this is more complicated
        if ( fitFrac_[omegaID][omegaID].value() < 1E-4 ) {
            omegaID = 1;
            storeID = 0;
        }

        // Check this is really the correct model

        LauRhoOmegaMix* rhomega = dynamic_cast<LauRhoOmegaMix*>( getResonance( omegaID ) );

        if ( rhomega != NULL ) {    // Bail out

            std::cout << "INFO in LauIsobarDynamics::calcExtraInfo : Calculating omega fit fraction from resonance "
                      << omegaID << std::endl;
            std::cout << "INFO in LauIsobarDynamics::calcExtraInfo : Storing omega fit fraction in resonance "
                      << storeID << std::endl;

            // Tell the RhoOmegaMix model only to give us the omega amplitude-squared

            rhomega->setWhichAmpSq( 1 );

            // Recalculate the integrals for the omega fit-fraction

            integralsDone_ = kFALSE;

            this->resetNormVectors();
            for ( UInt_t k( 0 ); k < nAmp_ + nIncohAmp_; ++k ) {
                integralsToBeCalculated_.insert( k );
            }
            this->calcDPNormalisation();

            integralsDone_ = kTRUE;

            Double_t fifjSumRealOmega = fifjSum_[omegaID][omegaID].re();

            // Recalculate the integrals for the rho fit-fraction

            rhomega->setWhichAmpSq( 2 );

            integralsDone_ = kFALSE;

            this->resetNormVectors();
            for ( UInt_t k( 0 ); k < nAmp_ + nIncohAmp_; ++k ) {
                integralsToBeCalculated_.insert( k );
            }
            this->calcDPNormalisation();

            integralsDone_ = kTRUE;

            Double_t fitFracPartRho = Amp_[omegaID].abs2() * fifjSum_[omegaID][omegaID].re();

            // Reset the RhoOmegaMix model and the integrals

            rhomega->setWhichAmpSq( 0 );

            integralsDone_ = kFALSE;

            this->resetNormVectors();
            for ( UInt_t k( 0 ); k < nAmp_ + nIncohAmp_; ++k ) {
                integralsToBeCalculated_.insert( k );
            }
            this->calcDPNormalisation();

            integralsDone_ = kTRUE;

            // Store the omega fit-fraction in the rho_COPY location (which is otherwise empty)
            // Store the rho fit-fraction in the rho location (overwriting the combined FF)

            Double_t omegaFF = fifjSumRealOmega * fitFrac_[omegaID][omegaID].value();

            fitFrac_[storeID][storeID].value( omegaFF );
            fitFrac_[omegaID][omegaID].value( fitFracPartRho * fNorm_[omegaID] * fNorm_[omegaID] /
                                              DPRate_.value() );

        } else {
            std::cout << "INFO in LauIsobarDynamics::calcExtraInfo : calculateRhoOmegaFitFractions is set, but the RhoOmegaMix model isn't in the right place. Ignoring this option."
                      << std::endl;
        }
    }
}

Bool_t LauIsobarDynamics::gotKMatrixMatch( UInt_t resAmpInt, const TString& propName ) const
{

    Bool_t gotMatch( kFALSE );

    if ( resAmpInt >= nAmp_ ) {
        return kFALSE;
    }

    const LauAbsResonance* theResonance = sigResonances_[resAmpInt];

    if ( theResonance == 0 ) {
        return kFALSE;
    }

    Int_t resModelInt = theResonance->getResonanceModel();

    if ( resModelInt == LauAbsResonance::KMatrix ) {

        TString resName = theResonance->getResonanceName();

        KMStringMap::const_iterator kMPropSetIter = kMatrixPropSet_.find( resName );

        if ( kMPropSetIter != kMatrixPropSet_.end() ) {
            TString kmPropString = kMPropSetIter->second;
            if ( kmPropString == propName ) {
                gotMatch = kTRUE;
            }
        }
    }

    return gotMatch;
}

Double_t LauIsobarDynamics::calcSigDPNorm()
{
    // Calculate the normalisation for the log-likelihood function.
    DPNorm_ = 0.0;

    for ( UInt_t i = 0; i < nAmp_; i++ ) {
        // fifjEffSum is the contribution from the term involving the resonance
        // dynamics (f_i for resonance i) and the efficiency term.
        Double_t fifjEffSumReal = fifjEffSum_[i][i].re();
        // We need to normalise this contribution w.r.t. the complete dynamics in the DP.
        // Hence we scale by the fNorm_i factor (squared), which is calculated by the
        // initialise() function, when the normalisation integrals are calculated and cached.
        // We also include the complex amplitude squared to get the total normalisation
        // contribution from this resonance.
        DPNorm_ += Amp_[i].abs2() * fifjEffSumReal * fNorm_[i] * fNorm_[i];
    }

    // We now come to the cross-terms (between resonances i and j) in the normalisation.
    for ( UInt_t i = 0; i < nAmp_; i++ ) {
        for ( UInt_t j = i + 1; j < nAmp_; j++ ) {
            LauComplex AmpjConj = Amp_[j].conj();
            LauComplex AmpTerm  = Amp_[i] * AmpjConj;
            // Again, fifjEffSum is the contribution from the term involving the resonance
            // dynamics (f_i*f_j_conjugate) and the efficiency cross term.
            // Also include the relative normalisation between these two resonances w.r.t. the
            // total DP dynamical structure (fNorm_i and fNorm_j) and the complex
            // amplitude squared (mag,phase) part.
            DPNorm_ += 2.0 * ( AmpTerm * fifjEffSum_[i][j] ).re() * fNorm_[i] * fNorm_[j];
        }
    }

    for ( UInt_t i = 0; i < nIncohAmp_; i++ ) {
        DPNorm_ += Amp_[i + nAmp_].abs2() * fSqEffSum_[i + nAmp_] * fNorm_[i + nAmp_] *
                   fNorm_[i + nAmp_];
    }

    return DPNorm_;
}

Bool_t LauIsobarDynamics::generate()
{
    // Routine to generate a signal event according to the Dalitz plot
    // model we have defined.

    // We need to make sure to calculate everything for every resonance
    integralsToBeCalculated_.clear();
    for ( UInt_t i { 0 }; i < nAmp_ + nIncohAmp_; ++i ) {
        integralsToBeCalculated_.insert( i );
    }

    if ( aSqMaxAuto_ ) {
        std::cout << "INFO in LauIsobarDynamics::generate : Starting auto-location of |A|^2 maximum"
                  << std::endl;

        LauASqMaxFinder finder { *this };
        aSqMaxSet_  = finder.find();
        aSqMaxAuto_ = kFALSE;

        std::cout << "INFO in LauIsobarDynamics::generate : auto-location of |A|^2 maximum finds: "
                  << aSqMaxSet_ << std::endl;
        aSqMaxSet_ *= 1.10;
        std::cout << "                                    : applying safety factor of 10\% gives:  "
                  << aSqMaxSet_ << std::endl;
    }

    nSigGenLoop_ = 0;
    Bool_t generatedSig( kFALSE );

    while ( generatedSig == kFALSE && nSigGenLoop_ < iterationsMax_ ) {

        // Generates uniform DP phase-space distribution
        Double_t m13Sq( 0.0 ), m23Sq( 0.0 );
        kinematics_->genFlatPhaseSpace( m13Sq, m23Sq );

        // If we're in a symmetrical DP then we should only generate events in one half
        // TODO - what do we do for fully symmetric?
        if ( symmetricalDP_ && ! fullySymmetricDP_ && m13Sq > m23Sq ) {
            Double_t tmpSq = m13Sq;
            m13Sq          = m23Sq;
            m23Sq          = tmpSq;
        }

        // Calculate the amplitudes and total amplitude for the given DP point
        this->calcLikelihoodInfo( m13Sq, m23Sq );

        // Throw the random number and check it against the ratio of ASq and the accept/reject ceiling
        const Double_t randNo = LauRandom::randomFun()->Rndm();
        if ( randNo > ASq_ / aSqMaxSet_ ) {
            ++nSigGenLoop_;
        } else {
            generatedSig = kTRUE;
            nSigGenLoop_ = 0;

            // Keep a note of the maximum ASq that we've found
            if ( ASq_ > aSqMaxVar_ ) {
                aSqMaxVar_ = ASq_;
            }
        }

    }    // while loop

    // Check that all is well with the generation
    Bool_t sigGenOK( kTRUE );
    if ( GenOK != this->checkToyMC( kTRUE, kFALSE ) ) {
        sigGenOK = kFALSE;
    }

    return sigGenOK;
}

LauIsobarDynamics::ToyMCStatus LauIsobarDynamics::checkToyMC( Bool_t printErrorMessages,
                                                              Bool_t printInfoMessages )
{
    // Check whether we have generated the toy MC OK.
    ToyMCStatus ok( GenOK );

    if ( nSigGenLoop_ >= iterationsMax_ ) {
        // Exceeded maximum allowed iterations - the generation is too inefficient
        if ( printErrorMessages ) {
            std::cerr << "WARNING in LauIsobarDynamics::checkToyMC : More than " << iterationsMax_
                      << " iterations performed and no event accepted." << std::endl;
        }

        if ( aSqMaxSet_ > 1.01 * aSqMaxVar_ ) {
            if ( printErrorMessages ) {
                std::cerr << "                                         : |A|^2 maximum was set to "
                          << aSqMaxSet_ << " but this appears to be too high." << std::endl;
                std::cerr << "                                         : Maximum value of |A|^2 found so far = "
                          << aSqMaxVar_ << std::endl;
                std::cerr << "                                         : The value of |A|^2 maximum will be decreased and the generation restarted."
                          << std::endl;
            }
            aSqMaxSet_ = 1.01 * aSqMaxVar_;
            std::cout << "INFO in LauIsobarDynamics::checkToyMC : |A|^2 max reset to " << aSqMaxSet_
                      << std::endl;
            ok = MaxIterError;
        } else {
            if ( printErrorMessages ) {
                std::cerr << "                                         : |A|^2 maximum was set to "
                          << aSqMaxSet_ << ", which seems to be correct for the given model."
                          << std::endl;
                std::cerr << "                                         : However, the generation is very inefficient - please check your model."
                          << std::endl;
                std::cerr << "                                         : The maximum number of iterations will be increased and the generation restarted."
                          << std::endl;
            }
            iterationsMax_ *= 2;
            std::cout << "INFO in LauIsobarDynamics::checkToyMC : max number of iterations reset to "
                      << iterationsMax_ << std::endl;
            ok = MaxIterError;
        }
    } else if ( aSqMaxVar_ > aSqMaxSet_ ) {
        // Found a value of ASq higher than the accept/reject ceiling - the generation is biased
        if ( printErrorMessages ) {
            std::cerr << "WARNING in LauIsobarDynamics::checkToyMC : |A|^2 maximum was set to "
                      << aSqMaxSet_ << " but a value exceeding this was found: " << aSqMaxVar_
                      << std::endl;
            std::cerr << "                                         : This value was found at m13Sq = "
                      << kinematics_->getm13Sq() << " and m23Sq = " << kinematics_->getm23Sq()
                      << std::endl;
            std::cerr << "                                         : Run was invalid, as any generated MC will be biased, according to the accept/reject method!"
                      << std::endl;
            std::cerr << "                                         : The value of |A|^2 maximum be reset to be > "
                      << aSqMaxVar_ << " and the generation restarted." << std::endl;
        }
        aSqMaxSet_ = 1.01 * aSqMaxVar_;
        std::cout << "INFO in LauIsobarDynamics::checkToyMC : |A|^2 max reset to " << aSqMaxSet_
                  << std::endl;
        ok = ASqMaxError;
    } else if ( printInfoMessages ) {
        std::cout << "INFO in LauIsobarDynamics::checkToyMC : aSqMaxSet = " << aSqMaxSet_
                  << " and aSqMaxVar = " << aSqMaxVar_ << std::endl;
    }

    return ok;
}

void LauIsobarDynamics::setDataEventNo( UInt_t iEvt )
{
    // Retrieve the data for event iEvt
    if ( data_.size() > iEvt ) {
        currentEvent_ = data_[iEvt];
    } else {
        std::cerr << "ERROR in LauIsobarDynamics::setDataEventNo : Event index too large: " << iEvt
                  << " >= " << data_.size() << "." << std::endl;
    }

    m13Sq_   = currentEvent_->retrievem13Sq();
    m23Sq_   = currentEvent_->retrievem23Sq();
    mPrime_  = currentEvent_->retrievemPrime();
    thPrime_ = currentEvent_->retrievethPrime();
    tagCat_  = currentEvent_->retrieveTagCat();
    eff_     = currentEvent_->retrieveEff();
    scfFraction_ =
        currentEvent_->retrieveScfFraction();    // These two are necessary, even though the dynamics don't actually use scfFraction_ or jacobian_,
    jacobian_ = currentEvent_->retrieveJacobian();    // since this is at the heart of the caching mechanism.
}

void LauIsobarDynamics::calcLikelihoodInfo( const UInt_t iEvt )
{
    // Calculate the likelihood and associated info
    // for the given event using cached information
    evtLike_ = 0.0;

    // retrieve the cached dynamics from the tree:
    // realAmp, imagAmp for each resonance plus efficiency, scf fraction and jacobian
    this->setDataEventNo( iEvt );

    // use realAmp and imagAmp to create the resonance amplitudes
    const std::vector<Double_t>& realAmp    = currentEvent_->retrieveRealAmp();
    const std::vector<Double_t>& imagAmp    = currentEvent_->retrieveImagAmp();
    const std::vector<Double_t>& incohInten = currentEvent_->retrieveIncohIntensities();
    for ( UInt_t i = 0; i < nAmp_; i++ ) {
        this->setFFTerm( i, realAmp[i], imagAmp[i] );
    }
    for ( UInt_t i = 0; i < nIncohAmp_; i++ ) {
        this->setIncohIntenTerm( i, incohInten[i] );
    }

    // Update the dynamics - calculates totAmp_ and then ASq_ = totAmp_.abs2() * eff_
    // All calculated using cached information on the individual amplitudes and efficiency.
    this->calcTotalAmp( kTRUE );

    // Calculate the normalised matrix element squared value
    if ( DPNorm_ > 1e-10 ) {
        evtLike_ = ASq_ / DPNorm_;
    }
}

void LauIsobarDynamics::calcLikelihoodInfo( const Double_t m13Sq, const Double_t m23Sq )
{
    this->calcLikelihoodInfo( m13Sq, m23Sq, -1 );
}

void LauIsobarDynamics::calcLikelihoodInfo( const Double_t m13Sq,
                                            const Double_t m23Sq,
                                            const Int_t tagCat )
{
    // Calculate the likelihood and associated info
    // for the given point in the Dalitz plot
    // Also retrieves the SCF fraction in the bin where the event lies (done
    // here to cache it along with the the rest of the DP quantities, like eff)
    // The jacobian for the square DP is calculated here for the same reason.
    evtLike_ = 0.0;

    // update the kinematics for the specified DP point
    kinematics_->updateKinematics( m13Sq, m23Sq );

    // calculate the jacobian and the scfFraction to cache them later
    scfFraction_ = this->retrieveScfFraction( tagCat );
    if ( kinematics_->squareDP() == kTRUE ) {
        jacobian_ = kinematics_->calcSqDPJacobian();
    }

    // calculate the ff_ terms and retrieves eff_ from the efficiency model
    this->calculateAmplitudes();
    // then calculate totAmp_ and finally ASq_ = totAmp_.abs2() * eff_
    this->calcTotalAmp( kTRUE );

    // Calculate the normalised matrix element squared value
    if ( DPNorm_ > 1e-10 ) {
        evtLike_ = ASq_ / DPNorm_;
    }
}

void LauIsobarDynamics::modifyDataTree()
{
    if ( recalcNormalisation_ == kFALSE ) {
        return;
    }

    const UInt_t nEvents = data_.size();

    std::set<UInt_t>::const_iterator iter         = integralsToBeCalculated_.begin();
    const std::set<UInt_t>::const_iterator intEnd = integralsToBeCalculated_.end();

    for ( UInt_t iEvt = 0; iEvt < nEvents; ++iEvt ) {

        currentEvent_ = data_[iEvt];

        std::vector<Double_t>& realAmp    = currentEvent_->retrieveRealAmp();
        std::vector<Double_t>& imagAmp    = currentEvent_->retrieveImagAmp();
        std::vector<Double_t>& incohInten = currentEvent_->retrieveIncohIntensities();

        const Double_t m13Sq = currentEvent_->retrievem13Sq();
        const Double_t m23Sq = currentEvent_->retrievem23Sq();
        const Int_t tagCat   = currentEvent_->retrieveTagCat();

        this->calcLikelihoodInfo( m13Sq, m23Sq, tagCat );

        for ( iter = integralsToBeCalculated_.begin(); iter != intEnd; ++iter ) {
            const UInt_t i = *iter;
            if ( *iter < nAmp_ ) {
                realAmp[i] = ff_[i].re();
                imagAmp[i] = ff_[i].im();
            } else {
                incohInten[i - nAmp_] = incohInten_[i - nAmp_];
            }
        }
    }
}

void LauIsobarDynamics::fillDataTree( const LauFitDataTree& inputFitTree )
{
    // In LauFitDataTree, the first two variables should always be m13^2 and m23^2.
    // Other variables follow thus: charge/flavour tag prob, etc.

    // Since this is the first caching, we need to make sure to calculate everything for every resonance
    integralsToBeCalculated_.clear();
    for ( UInt_t i( 0 ); i < nAmp_ + nIncohAmp_; ++i ) {
        integralsToBeCalculated_.insert( i );
    }

    UInt_t nBranches = inputFitTree.nBranches();

    if ( nBranches < 2 ) {
        std::cerr << "ERROR in LauIsobarDynamics::fillDataTree : Expecting at least 2 variables "
                  << "in input data tree, but there are " << nBranches << "!\n";
        std::cerr << "                                         : Make sure you have the right number of variables in your input data file!"
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Data structure that will cache the variables required to
    // calculate the signal likelihood for this experiment
    for ( std::vector<LauCacheData*>::iterator iter = data_.begin(); iter != data_.end(); ++iter ) {
        delete ( *iter );
    }
    data_.clear();

    Double_t m13Sq( 0.0 ), m23Sq( 0.0 );
    Double_t mPrime( 0.0 ), thPrime( 0.0 );
    Int_t tagCat( -1 );
    std::vector<Double_t> realAmp( nAmp_ ), imagAmp( nAmp_ );
    Double_t eff( 0.0 ), scfFraction( 0.0 ), jacobian( 0.0 );

    UInt_t nEvents = inputFitTree.nEvents() + inputFitTree.nFakeEvents();

    data_.reserve( nEvents );

    for ( UInt_t iEvt = 0; iEvt < nEvents; ++iEvt ) {

        const LauFitData& dataValues    = inputFitTree.getData( iEvt );
        LauFitData::const_iterator iter = dataValues.find( "m13Sq" );
        m13Sq                           = iter->second;
        iter                            = dataValues.find( "m23Sq" );
        m23Sq                           = iter->second;

        // is there more than one tagging category?
        // if so then we need to know the category from the data
        if ( scfFractionModel_.size() > 1 ) {
            iter   = dataValues.find( "tagCat" );
            tagCat = static_cast<Int_t>( iter->second );
        }

        // calculates the amplitudes and total amplitude for the given DP point
        // tagging category not needed by dynamics, but to find out the scfFraction
        this->calcLikelihoodInfo( m13Sq, m23Sq, tagCat );

        // extract the real and imaginary parts of the ff_ terms for storage
        for ( UInt_t i = 0; i < nAmp_; i++ ) {
            realAmp[i] = ff_[i].re();
            imagAmp[i] = ff_[i].im();
        }

        if ( kinematics_->squareDP() ) {
            mPrime  = kinematics_->getmPrime();
            thPrime = kinematics_->getThetaPrime();
        }

        eff         = this->getEvtEff();
        scfFraction = this->getEvtScfFraction();
        jacobian    = this->getEvtJacobian();

        // store the data for each event in the list
        data_.push_back( new LauCacheData() );
        data_[iEvt]->storem13Sq( m13Sq );
        data_[iEvt]->storem23Sq( m23Sq );
        data_[iEvt]->storemPrime( mPrime );
        data_[iEvt]->storethPrime( thPrime );
        data_[iEvt]->storeTagCat( tagCat );
        data_[iEvt]->storeEff( eff );
        data_[iEvt]->storeScfFraction( scfFraction );
        data_[iEvt]->storeJacobian( jacobian );
        data_[iEvt]->storeRealAmp( realAmp );
        data_[iEvt]->storeImagAmp( imagAmp );
        data_[iEvt]->storeIncohIntensities( incohInten_ );
    }
}

Bool_t LauIsobarDynamics::gotReweightedEvent()
{
    // Select the event (kinematics_) using an accept/reject method based on the
    // ratio of the current value of ASq to the maximal value.
    Bool_t accepted( kFALSE );

    // calculate the ff_ terms and retrieves eff_ from the efficiency model
    this->calculateAmplitudes();
    // then calculate totAmp_ and finally ASq_ = totAmp_.abs2() (without the efficiency correction!)
    this->calcTotalAmp( kFALSE );

    // Compare the ASq value with the maximal value (set by the user)
    if ( LauRandom::randomFun()->Rndm() < ASq_ / aSqMaxSet_ ) {
        accepted = kTRUE;
    }

    if ( ASq_ > aSqMaxVar_ ) {
        aSqMaxVar_ = ASq_;
    }

    return accepted;
}

Double_t LauIsobarDynamics::getEventWeight()
{
    // calculate the ff_ terms and retrieves eff_ from the efficiency model
    this->calculateAmplitudes();

    // then calculate totAmp_ and finally ASq_ = totAmp_.abs2() (without the efficiency correction!)
    this->calcTotalAmp( kFALSE );

    // return the event weight = the value of the squared amplitude
    return ASq_;
}

void LauIsobarDynamics::updateCoeffs( const std::vector<LauComplex>& coeffs )
{
    // Check that the number of coeffs is correct
    if ( coeffs.size() != this->getnTotAmp() ) {
        std::cerr << "ERROR in LauIsobarDynamics::updateCoeffs : Expected " << this->getnTotAmp()
                  << " but got " << coeffs.size() << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // Now check if the coeffs have changed
    Bool_t changed = ( Amp_ != coeffs );
    if ( changed ) {
        // Copy the coeffs
        Amp_ = coeffs;
    }

    // TODO should perhaps keep track of whether the resonance parameters have changed here and if none of those and none of the coeffs have changed then we don't need to update the norm

    // Update the total normalisation for the signal likelihood
    this->calcSigDPNorm();
}

TString LauIsobarDynamics::getConjResName( const TString& resName ) const
{
    // Get the name of the charge conjugate resonance
    TString conjName( resName );

    Ssiz_t index1 = resName.Index( "+" );
    Ssiz_t index2 = resName.Index( "-" );
    if ( index1 != -1 ) {
        conjName.Replace( index1, 1, "-" );
    } else if ( index2 != -1 ) {
        conjName.Replace( index2, 1, "+" );
    }

    return conjName;
}

Double_t LauIsobarDynamics::retrieveEfficiency()
{
    Double_t eff( 1.0 );
    if ( effModel_ != 0 ) {
        eff = effModel_->calcEfficiency( kinematics_ );
    }
    return eff;
}

Double_t LauIsobarDynamics::retrieveScfFraction( Int_t tagCat )
{
    Double_t scfFraction( 0.0 );

    // scf model and eff model are exactly the same, functionally
    // so we use an instance of LauAbsEffModel, and the method
    // calcEfficiency actually calculates the scf fraction
    if ( tagCat == -1 ) {
        if ( ! scfFractionModel_.empty() ) {
            scfFraction = scfFractionModel_[0]->calcEfficiency( kinematics_ );
        }
    } else {
        scfFraction = scfFractionModel_[tagCat]->calcEfficiency( kinematics_ );
    }
    return scfFraction;
}
