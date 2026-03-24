
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

/*! \file LauBkgndDPModel.cc
    \brief File containing implementation of LauBkgndDPModel class.
*/

#include "LauBkgndDPModel.hh"

#include "Lau2DHistDPPdf.hh"
#include "Lau2DSplineDPPdf.hh"
#include "LauDaughters.hh"
#include "LauFitDataTree.hh"
#include "LauKinematics.hh"
#include "LauRandom.hh"
#include "LauVetoes.hh"

#include "TRandom.h"
#include "TSystem.h"

#include <cstdlib>
#include <iostream>

LauBkgndDPModel::LauBkgndDPModel( LauDaughters* daughters, LauVetoes* vetoes ) :
    LauAbsBkgndDPModel( daughters, vetoes ),
    symmetricalDP_( kFALSE ),
    squareDP_( kFALSE ),
    bgHistDPPdf_( 0 ),
    curEvtHistVal_( 0.0 ),
    maxPdfHeight_( 1.0 ),
    pdfNorm_( 1.0 ),
    doneGenWarning_( kFALSE ),
    lowBinWarningIssued_( kFALSE )
{
    if ( daughters != 0 ) {
        symmetricalDP_ = daughters->gotSymmetricalDP();
    }
}

LauBkgndDPModel::~LauBkgndDPModel()
{
    if ( bgHistDPPdf_ != 0 ) {
        delete bgHistDPPdf_;
        bgHistDPPdf_ = 0;
    }
}

void LauBkgndDPModel::setBkgndHisto( const TH2* histo,
                                     Bool_t useInterpolation,
                                     Bool_t fluctuateBins,
                                     Bool_t useUpperHalfOnly,
                                     Bool_t squareDP )
{
    if ( ! histo ) {
        std::cerr << "WARNING in LauBkgndDPModel::setBkgndHisto : Supplied background histogram pointer is null, likelihood for this component will be flat in the Dalitz plot"
                  << std::endl;
    }

    Bool_t upperHalf = kFALSE;
    if ( symmetricalDP_ == kTRUE && useUpperHalfOnly == kTRUE ) {
        upperHalf = kTRUE;
    }
    std::cout << "INFO in LauBkgndDPModel::setBkgndHisto : Background histogram has upperHalf = "
              << static_cast<Int_t>( upperHalf ) << std::endl;

    squareDP_ = squareDP;

    LauKinematics* kinematics = this->getKinematics();
    const LauVetoes* vetoes   = this->getVetoes();
    bgHistDPPdf_              = new Lau2DHistDPPdf( histo,
                                       kinematics,
                                       vetoes,
                                       useInterpolation,
                                       fluctuateBins,
                                       upperHalf,
                                       squareDP );

    maxPdfHeight_ = bgHistDPPdf_->getMaxHeight();
    pdfNorm_      = bgHistDPPdf_->getHistNorm();
}

void LauBkgndDPModel::setBkgndSpline( const TH2* histo,
                                      Bool_t fluctuateBins,
                                      Bool_t useUpperHalfOnly,
                                      Bool_t squareDP )
{
    if ( ! histo ) {
        std::cerr << "WARNING in LauBkgndDPModel::setBkgndSpline : Supplied background histogram pointer is null, construction of spline will fail"
                  << std::endl;
    }

    Bool_t upperHalf = kFALSE;
    if ( symmetricalDP_ == kTRUE && useUpperHalfOnly == kTRUE ) {
        upperHalf = kTRUE;
    }
    std::cout << "INFO in LauBkgndDPModel::setBkgndSpline : Background histogram has upperHalf = "
              << static_cast<Int_t>( upperHalf ) << std::endl;

    squareDP_ = squareDP;

    LauKinematics* kinematics = this->getKinematics();
    const LauVetoes* vetoes   = this->getVetoes();
    bgHistDPPdf_ =
        new Lau2DSplineDPPdf( histo, kinematics, vetoes, fluctuateBins, upperHalf, squareDP );

    maxPdfHeight_ = bgHistDPPdf_->getMaxHeight();
    pdfNorm_      = bgHistDPPdf_->getHistNorm();
}

Double_t LauBkgndDPModel::calcHistValue( Double_t xVal, Double_t yVal )
{
    // Get the likelihood value of the background in the Dalitz plot.

    // Check that we have a valid histogram PDF
    if ( bgHistDPPdf_ == 0 ) {
        std::cerr << "WARNING in LauBkgndDPModel::calcHistValue : We don't have a histogram so assuming the likelihood is flat in the Dalitz plot."
                  << std::endl;
        this->setBkgndHisto( 0, kFALSE, kFALSE, kFALSE, kFALSE );
    }

    // Find out the un-normalised PDF value
    Double_t value = bgHistDPPdf_->interpolateXY( xVal, yVal );

    // Check that the value is greater than zero
    // If we're using a spline then negative values can be caused by adjacent bins that all contain a value of zero.
    // The spline requires the value, its first derivatives and the mixed second derivative to be continuous and to match the input histogram
    // at the bin centres. Derivatives are calculated using a finite difference approximation taking the difference between the neighbouring bins.
    // If two bins are zero but the third is not then the second bin will have a positive first derivative causing the spline to dip below zero
    // between the two zero bins to remain smooth. Such dips are unavoidable but are correctly removed here.
    if ( value < 0.0 ) {
        if ( ! lowBinWarningIssued_ ) {
            std::cerr << "WARNING in LauBkgndDPModel::calcHistValue : Value " << value
                      << " is less than 0 - setting to 0.  You may want to check your histogram!"
                      << std::endl
                      << "                                          : If you are using a spline then this could be caused by adjacent empty bins. Further warnings will be suppressed."
                      << std::endl;
            lowBinWarningIssued_ = kTRUE;
        }
        return 0.0;
    }

    LauKinematics* kinematics = this->getKinematics();

    // For square DP co-ordinates, need to divide by Jacobian
    if ( squareDP_ == kTRUE ) {
        // Make sure that square DP kinematics are correct, then get Jacobian
        kinematics->updateSqDPKinematics( xVal, yVal );
        Double_t jacobian  = kinematics->calcSqDPJacobian();
        value             /= jacobian;
    }

    return value;
}

void LauBkgndDPModel::initialise()
{
}

Bool_t LauBkgndDPModel::generate()
{
    // Routine to generate the background, using data provided by an
    // already defined histogram.

    LauKinematics* kinematics = this->getKinematics();

    Bool_t gotBG( kFALSE );

    while ( gotBG == kFALSE ) {

        if ( squareDP_ == kTRUE ) {
            // Generate a point in m', theta' space. By construction, this point
            // is already within the DP region.
            Double_t mPrime( 0.0 ), thetaPrime( 0.0 );
            kinematics->genFlatSqDP( mPrime, thetaPrime );

            // If we're in a symmetrical DP then we should only generate events in one half
            if ( symmetricalDP_ && thetaPrime > 0.5 ) {
                thetaPrime = 1.0 - thetaPrime;
            }

            // Calculate histogram height for DP point and
            // compare with the maximum height
            if ( bgHistDPPdf_ != 0 ) {
                Double_t bgContDP = bgHistDPPdf_->interpolateXY( mPrime, thetaPrime ) / maxPdfHeight_;
                if ( LauRandom::randomFun()->Rndm() < bgContDP ) {
                    kinematics->updateSqDPKinematics( mPrime, thetaPrime );
                    gotBG = kTRUE;
                }
            } else {
                if ( ! doneGenWarning_ ) {
                    std::cerr << "WARNING in LauBkgndDPModel::generate : We don't have a histogram so generating flat in the square DP, which won't be flat in the conventional DP!"
                              << std::endl;
                    std::cerr << "WARNING in LauBkgndDPModel::generate : This should never happen!! What have you done?!"
                              << std::endl;
                    doneGenWarning_ = kTRUE;
                }
                kinematics->updateSqDPKinematics( mPrime, thetaPrime );
                gotBG = kTRUE;
            }
        } else {
            // Generate a point in the Dalitz plot (phase-space).
            Double_t m13Sq( 0.0 ), m23Sq( 0.0 );
            kinematics->genFlatPhaseSpace( m13Sq, m23Sq );

            // If we're in a symmetrical DP then we should only generate events in one half
            if ( symmetricalDP_ && m13Sq > m23Sq ) {
                Double_t tmpSq = m13Sq;
                m13Sq          = m23Sq;
                m23Sq          = tmpSq;
            }

            // Calculate histogram height for DP point and
            // compare with the maximum height
            if ( bgHistDPPdf_ != 0 ) {
                Double_t bgContDP = bgHistDPPdf_->interpolateXY( m13Sq, m23Sq ) / maxPdfHeight_;
                if ( LauRandom::randomFun()->Rndm() < bgContDP ) {
                    kinematics->updateKinematics( m13Sq, m23Sq );
                    gotBG = kTRUE;
                }
            } else {
                if ( ! doneGenWarning_ ) {
                    std::cerr << "WARNING in LauBkgndDPModel::generate : We don't have a histogram so generating flat in the DP."
                              << std::endl;
                    doneGenWarning_ = kTRUE;
                }
                kinematics->updateKinematics( m13Sq, m23Sq );
                gotBG = kTRUE;
            }
        }
    }

    // Implement veto
    Bool_t vetoOK( kTRUE );
    const LauVetoes* vetoes = this->getVetoes();
    if ( vetoes ) {
        vetoOK = vetoes->passVeto( kinematics );
    }
    // Call this function recusively until we pass the veto.
    if ( vetoOK == kFALSE ) {
        this->generate();
    }
    return kTRUE;
}

void LauBkgndDPModel::fillDataTree( const LauFitDataTree& inputFitTree )
{
    // In LauFitDataTree, the first two variables should always be
    // m13^2 and m23^2. Other variables follow thus: charge.

    Int_t nBranches = inputFitTree.nBranches();
    if ( nBranches < 2 ) {
        std::cerr << "ERROR in LauBkgndDPModel::fillDataTree : Expecting at least 2 variables "
                  << "in input data tree, but there are " << nBranches << "!\n";
        std::cerr << "                                       : Make sure you have the right number of variables in your input data file!"
                  << std::endl;
        return;
    }

    const UInt_t nEvents { inputFitTree.nEvents() };

    // clear and resize the data vector
    bgData_.clear();
    bgData_.resize( nEvents );

    Double_t m13Sq { 0.0 }, m23Sq { 0.0 };

    for ( UInt_t iEvt = 0; iEvt < nEvents; iEvt++ ) {

        const LauFitData& dataValues = inputFitTree.getData( iEvt );
        m13Sq                        = dataValues.at( "m13Sq" );
        m23Sq                        = dataValues.at( "m23Sq" );

        bgData_[iEvt] = this->getUnNormValue( m13Sq, m23Sq );
    }
}

Double_t LauBkgndDPModel::getUnNormValue( const Double_t m13Sq, const Double_t m23Sq )
{
    // Update the kinematics. This will also update m' and theta' if squareDP = true
    LauKinematics* kinematics = this->getKinematics();
    kinematics->updateKinematics( m13Sq, m23Sq );

    if ( squareDP_ ) {
        const Double_t mPrime { kinematics->getmPrime() };
        const Double_t thetaPrime { kinematics->getThetaPrime() };
        curEvtHistVal_ = this->calcHistValue( mPrime, thetaPrime );
    } else {
        curEvtHistVal_ = this->calcHistValue( m13Sq, m23Sq );
    }

    return curEvtHistVal_;
}

Double_t LauBkgndDPModel::getLikelihood( const Double_t m13Sq, const Double_t m23Sq )
{
    this->getUnNormValue( m13Sq, m23Sq );
    return curEvtHistVal_ / this->getPdfNorm();
}

Double_t LauBkgndDPModel::getUnNormValue( UInt_t iEvt )
{
    // Retrieve the likelihood for the given event
    this->setDataEventNo( iEvt );
    return curEvtHistVal_;
}

Double_t LauBkgndDPModel::getLikelihood( UInt_t iEvt )
{
    // Retrieve the likelihood for the given event
    this->setDataEventNo( iEvt );
    Double_t llhd = curEvtHistVal_ / this->getPdfNorm();
    return llhd;
}

void LauBkgndDPModel::setDataEventNo( UInt_t iEvt )
{
    // Retrieve the data for event iEvt
    if ( bgData_.size() > iEvt ) {
        curEvtHistVal_ = bgData_[iEvt];
    } else {
        std::cerr << "ERROR in LauBkgndDPModel::setDataEventNo : Event index too large: " << iEvt
                  << " >= " << bgData_.size() << "." << std::endl;
    }
}
