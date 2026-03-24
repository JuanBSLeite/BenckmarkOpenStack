
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

/*! \file LauSimpleFitModel.cc
    \brief File containing implementation of LauSimpleFitModel class.
*/

#include "LauSimpleFitModel.hh"

#include "LauAbsBkgndDPModel.hh"
#include "LauAbsCoeffSet.hh"
#include "LauAbsPdf.hh"
#include "LauComplex.hh"
#include "LauConstants.hh"
#include "LauDaughters.hh"
#include "LauEffModel.hh"
#include "LauEmbeddedData.hh"
#include "LauFitNtuple.hh"
#include "LauGenNtuple.hh"
#include "LauIsobarDynamics.hh"
#include "LauKinematics.hh"
#include "LauPrint.hh"
#include "LauRandom.hh"
#include "LauScfMap.hh"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH2.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TVirtualFitter.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <typeinfo>

LauSimpleFitModel::LauSimpleFitModel( LauIsobarDynamics* sigModel ) :
    LauAbsFitModel(),
    sigDPModel_( sigModel ),
    kinematics_( sigModel ? sigModel->getKinematics() : 0 ),
    usingBkgnd_( kFALSE ),
    nSigComp_( 0 ),
    nSigDPPar_( 0 ),
    nExtraPdfPar_( 0 ),
    nNormPar_( 0 ),
    meanEff_( "meanEff", 0.0, 0.0, 1.0 ),
    dpRate_( "dpRate", 0.0, 0.0, 100.0 ),
    //signalEvents_("signalEvents",1.0,0.0,1.0),
    signalEvents_( 0 ),
    useSCF_( kFALSE ),
    useSCFHist_( kFALSE ),
    scfFrac_( "scfFrac", 0.0, 0.0, 1.0 ),
    scfFracHist_( 0 ),
    scfMap_( 0 ),
    compareFitData_( kFALSE ),
    signalTree_( 0 ),
    reuseSignal_( kFALSE ),
    useReweighting_( kFALSE ),
    sigDPLike_( 0.0 ),
    scfDPLike_( 0.0 ),
    sigExtraLike_( 0.0 ),
    scfExtraLike_( 0.0 ),
    sigTotalLike_( 0.0 ),
    scfTotalLike_( 0.0 )
{
}

LauSimpleFitModel::~LauSimpleFitModel()
{
    delete signalTree_;
    for ( LauBkgndEmbDataList::iterator iter = bkgndTree_.begin(); iter != bkgndTree_.end(); ++iter ) {
        delete ( *iter );
    }
    delete scfFracHist_;
}

void LauSimpleFitModel::setupBkgndVectors()
{
    UInt_t nBkgnds = this->nBkgndClasses();
    bkgndDPModels_.resize( nBkgnds );
    bkgndPdfs_.resize( nBkgnds );
    bkgndEvents_.resize( nBkgnds );
    bkgndTree_.resize( nBkgnds );
    reuseBkgnd_.resize( nBkgnds );
    bkgndDPLike_.resize( nBkgnds );
    bkgndExtraLike_.resize( nBkgnds );
    bkgndTotalLike_.resize( nBkgnds );
}

void LauSimpleFitModel::setNSigEvents( LauParameter* nSigEvents )
{
    if ( nSigEvents == 0 ) {
        std::cerr << "ERROR in LauSimpleFitModel::setNSigEvents : The signal yield LauParameter pointer is null."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }
    if ( signalEvents_ != 0 ) {
        std::cerr << "ERROR in LauSimpleFitModel::setNSigEvents : You are trying to overwrite the signal yield."
                  << std::endl;
        return;
    }

    signalEvents_ = nSigEvents;
    TString name  = signalEvents_->name();
    if ( ! name.Contains( "signalEvents" ) &&
         ! ( name.BeginsWith( "signal" ) && name.EndsWith( "Events" ) ) ) {
        signalEvents_->name( "signalEvents" );
    }
    Double_t value = nSigEvents->value();
    signalEvents_->range( -2.0 * ( TMath::Abs( value ) + 1.0 ), 2.0 * ( TMath::Abs( value ) + 1.0 ) );
}

void LauSimpleFitModel::setNBkgndEvents( LauAbsRValue* nBkgndEvents )
{
    if ( nBkgndEvents == 0 ) {
        std::cerr << "ERROR in LauSimpleFitModel::setNBkgndEvents : The background yield LauParameter pointer is null."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    if ( ! this->validBkgndClass( nBkgndEvents->name() ) ) {
        std::cerr << "ERROR in LauSimpleFitModel::setNBkgndEvents : Invalid background class \""
                  << nBkgndEvents->name() << "\"." << std::endl;
        std::cerr << "                                            : Background class names must be provided in \"setBkgndClassNames\" before any other background-related actions can be performed."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    UInt_t bkgndID = this->bkgndClassID( nBkgndEvents->name() );

    if ( bkgndEvents_[bkgndID] != 0 ) {
        std::cerr << "ERROR in LauSimpleFitModel::setNBkgndEvents : You are trying to overwrite the background yield."
                  << std::endl;
        return;
    }

    nBkgndEvents->name( nBkgndEvents->name() + "Events" );
    if ( nBkgndEvents->isLValue() ) {
        Double_t value      = nBkgndEvents->value();
        LauParameter* yield = dynamic_cast<LauParameter*>( nBkgndEvents );
        yield->range( -2.0 * ( TMath::Abs( value ) + 1.0 ), 2.0 * ( TMath::Abs( value ) + 1.0 ) );
    }
    bkgndEvents_[bkgndID] = nBkgndEvents;
}

void LauSimpleFitModel::splitSignalComponent( const TH2* dpHisto,
                                              const Bool_t upperHalf,
                                              const Bool_t fluctuateBins,
                                              LauScfMap* scfMap )
{
    if ( useSCF_ == kTRUE ) {
        std::cerr << "ERROR in LauSimpleFitModel::splitSignalComponent : Have already setup SCF."
                  << std::endl;
        return;
    }

    if ( dpHisto == 0 ) {
        std::cerr << "ERROR in LauSimpleFitModel::splitSignalComponent : The histogram pointer is null."
                  << std::endl;
        return;
    }

    const LauDaughters* daughters = sigDPModel_->getDaughters();
    scfFracHist_                  = new LauEffModel( daughters, 0 );
    scfFracHist_
        ->setEffHisto( dpHisto, kTRUE, fluctuateBins, 0.0, 0.0, upperHalf, daughters->squareDP() );

    scfMap_ = scfMap;

    useSCF_     = kTRUE;
    useSCFHist_ = kTRUE;
}

void LauSimpleFitModel::splitSignalComponent( const Double_t scfFrac, const Bool_t fixed )
{
    if ( useSCF_ == kTRUE ) {
        std::cerr << "ERROR in LauSimpleFitModel::splitSignalComponent : Have already setup SCF."
                  << std::endl;
        return;
    }

    scfFrac_.range( 0.0, 1.0 );
    scfFrac_.value( scfFrac );
    scfFrac_.initValue( scfFrac );
    scfFrac_.genValue( scfFrac );
    scfFrac_.fixed( fixed );

    useSCF_     = kTRUE;
    useSCFHist_ = kFALSE;
}

void LauSimpleFitModel::setBkgndDPModel( const TString& bkgndClass, LauAbsBkgndDPModel* bkgndDPModel )
{
    if ( bkgndDPModel == 0 ) {
        std::cerr << "ERROR in LauSimpleFitModel::setBkgndDPModel : The model pointer is null."
                  << std::endl;
        return;
    }

    // check that this background name is valid
    if ( ! this->validBkgndClass( bkgndClass ) ) {
        std::cerr << "ERROR in LauSimpleFitModel::setBkgndDPModel : Invalid background class \""
                  << bkgndClass << "\"." << std::endl;
        std::cerr << "                                            : Background class names must be provided in \"setBkgndClassNames\" before any other background-related actions can be performed."
                  << std::endl;
        return;
    }

    UInt_t bkgndID          = this->bkgndClassID( bkgndClass );
    bkgndDPModels_[bkgndID] = bkgndDPModel;

    usingBkgnd_ = kTRUE;
}

void LauSimpleFitModel::setSignalPdf( LauAbsPdf* pdf )
{
    if ( pdf == 0 ) {
        std::cerr << "ERROR in LauSimpleFitModel::setSignalPdf : The PDF pointer is null."
                  << std::endl;
        return;
    }
    signalPdfs_.push_back( pdf );
}

void LauSimpleFitModel::setSCFPdf( LauAbsPdf* pdf )
{
    if ( pdf == 0 ) {
        std::cerr << "ERROR in LauSimpleFitModel::setSCFPdf : The PDF pointer is null." << std::endl;
        return;
    }
    scfPdfs_.push_back( pdf );
}

void LauSimpleFitModel::setBkgndPdf( const TString& bkgndClass, LauAbsPdf* pdf )
{
    if ( pdf == 0 ) {
        std::cerr << "ERROR in LauSimpleFitModel::setBkgndPdf : The PDF pointer is null."
                  << std::endl;
        return;
    }

    // check that this background name is valid
    if ( ! this->validBkgndClass( bkgndClass ) ) {
        std::cerr << "ERROR in LauSimpleFitModel::setBkgndPdf : Invalid background class \""
                  << bkgndClass << "\"." << std::endl;
        std::cerr << "                                        : Background class names must be provided in \"setBkgndClassNames\" before any other background-related actions can be performed."
                  << std::endl;
        return;
    }

    UInt_t bkgndID = this->bkgndClassID( bkgndClass );
    bkgndPdfs_[bkgndID].push_back( pdf );

    usingBkgnd_ = kTRUE;
}

void LauSimpleFitModel::setAmpCoeffSet( LauAbsCoeffSet* coeffSet )
{
    // Resize the coeffPars vector if not already done
    if ( coeffPars_.empty() ) {
        coeffPars_.resize( sigDPModel_->getnTotAmp() );
        for ( std::vector<LauAbsCoeffSet*>::iterator iter = coeffPars_.begin();
              iter != coeffPars_.end();
              ++iter ) {
            ( *iter ) = 0;
        }
    }

    // Is there a component called compName in the signal model?
    TString compName( coeffSet->name() );
    const Int_t index = sigDPModel_->resonanceIndex( compName );
    if ( index < 0 ) {
        std::cerr << "ERROR in LauSimpleFitModel::setAmpCoeffSet : Signal DP model doesn't contain component \""
                  << compName << "\"." << std::endl;
        return;
    }

    // Do we already have it in our list of names?
    if ( coeffPars_[index] != 0 && coeffPars_[index]->name() == compName ) {
        std::cerr << "ERROR in LauSimpleFitModel::setAmpCoeffSet : Have already set coefficients for \""
                  << compName << "\"." << std::endl;
        return;
    }

    coeffSet->index( index );
    coeffPars_[index] = coeffSet;

    ++nSigComp_;

    std::cout << "INFO in LauSimpleFitModel::setAmpCoeffSet : Added coefficients for component A"
              << index << ": \"" << compName << "\" to the fit model." << std::endl;
    coeffSet->printParValues();
}

void LauSimpleFitModel::initialise()
{
    // Initialisation
    if ( ! this->useDP() && signalPdfs_.empty() ) {
        std::cerr << "ERROR in LauSimpleFitModel::initialise : Signal model doesn't exist for any variable."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    if ( this->useDP() ) {
        // Check that we have all the Dalitz-plot models
        if ( sigDPModel_ == 0 ) {
            std::cerr << "ERROR in LauSimpleFitModel::initialise : The pointer to the signal DP model is null.\n";
            std::cerr << "                                       : Removing the Dalitz Plot from the model."
                      << std::endl;
            this->useDP( kFALSE );
        }
        if ( usingBkgnd_ ) {
            for ( LauBkgndDPModelList::const_iterator dpmodel_iter = bkgndDPModels_.begin();
                  dpmodel_iter != bkgndDPModels_.end();
                  ++dpmodel_iter ) {
                if ( ( *dpmodel_iter ) == 0 ) {
                    std::cerr << "ERROR in LauSimpleFitModel::initialise : The pointer to one of the background DP models is null.\n";
                    std::cerr << "                                       : Removing the Dalitz Plot from the model."
                              << std::endl;
                    this->useDP( kFALSE );
                    break;
                }
            }
        }

        // Need to check that the number of components we have and that the dynamics has matches up
        UInt_t nAmp = sigDPModel_->getnTotAmp();
        if ( nAmp != nSigComp_ ) {
            std::cerr << "ERROR in LauSimpleFitModel::initialise : Number of signal DP components with magnitude and phase set not right."
                      << std::endl;
            gSystem->Exit( EXIT_FAILURE );
        }

        // From the initial parameter values calculate the coefficients
        // so they can be passed to the signal model
        this->updateCoeffs();

        // If all is well, go ahead and initialise them
        this->initialiseDPModels();
    }

    // Next check that, if a given component is being used we've got the
    // right number of PDFs for all the variables involved
    // TODO - should probably check variable names and so on as well

    UInt_t nsigpdfvars( 0 );
    for ( LauPdfPList::const_iterator pdf_iter = signalPdfs_.begin(); pdf_iter != signalPdfs_.end();
          ++pdf_iter ) {
        std::vector<TString> varNames = ( *pdf_iter )->varNames();
        for ( std::vector<TString>::const_iterator var_iter = varNames.begin();
              var_iter != varNames.end();
              ++var_iter ) {
            if ( ( *var_iter ) != "m13Sq" && ( *var_iter ) != "m23Sq" ) {
                ++nsigpdfvars;
            }
        }
    }
    if ( useSCF_ ) {
        UInt_t nscfpdfvars( 0 );
        for ( LauPdfPList::const_iterator pdf_iter = scfPdfs_.begin(); pdf_iter != scfPdfs_.end();
              ++pdf_iter ) {
            std::vector<TString> varNames = ( *pdf_iter )->varNames();
            for ( std::vector<TString>::const_iterator var_iter = varNames.begin();
                  var_iter != varNames.end();
                  ++var_iter ) {
                if ( ( *var_iter ) != "m13Sq" && ( *var_iter ) != "m23Sq" ) {
                    ++nscfpdfvars;
                }
            }
        }
        if ( nscfpdfvars != nsigpdfvars ) {
            std::cerr << "ERROR in LauSimpleFitModel::initialise : There are " << nsigpdfvars
                      << " TM signal PDF variables but " << nscfpdfvars
                      << " SCF signal PDF variables." << std::endl;
            gSystem->Exit( EXIT_FAILURE );
        }
    }
    if ( usingBkgnd_ ) {
        for ( LauBkgndPdfsList::const_iterator bgclass_iter = bkgndPdfs_.begin();
              bgclass_iter != bkgndPdfs_.end();
              ++bgclass_iter ) {
            UInt_t nbkgndpdfvars( 0 );
            const LauPdfPList& pdfList = ( *bgclass_iter );
            for ( LauPdfPList::const_iterator pdf_iter = pdfList.begin(); pdf_iter != pdfList.end();
                  ++pdf_iter ) {
                std::vector<TString> varNames = ( *pdf_iter )->varNames();
                for ( std::vector<TString>::const_iterator var_iter = varNames.begin();
                      var_iter != varNames.end();
                      ++var_iter ) {
                    if ( ( *var_iter ) != "m13Sq" && ( *var_iter ) != "m23Sq" ) {
                        ++nbkgndpdfvars;
                    }
                }
            }
            if ( nbkgndpdfvars != nsigpdfvars ) {
                std::cerr << "ERROR in LauSimpleFitModel::initialise : There are " << nsigpdfvars
                          << " signal PDF variables but " << nbkgndpdfvars
                          << " bkgnd PDF variables." << std::endl;
                gSystem->Exit( EXIT_FAILURE );
            }
        }
    }

    // Clear the vectors of parameter information so we can start from scratch
    this->clearFitParVectors();

    // Set the fit parameters for signal and background models
    this->setSignalDPParameters();

    // Set the fit parameters for the various extra PDFs
    this->setExtraPdfParameters();

    // Set the initial bg and signal events
    this->setFitNEvents();

    // Check that we have the expected number of fit variables
    const LauParameterPList& fitVars = this->fitPars();
    if ( fitVars.size() != ( nSigDPPar_ + nExtraPdfPar_ + nNormPar_ ) ) {
        std::cerr << "ERROR in LauSimpleFitModel::initialise : Number of fit parameters not of expected size."
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    this->setExtraNtupleVars();
}

void LauSimpleFitModel::initialiseDPModels()
{
    std::cout << "INFO in LauSimpleFitModel::initialiseDPModels : Initialising DP models"
              << std::endl;

    sigDPModel_->initialise( coeffs_ );

    if ( usingBkgnd_ == kTRUE ) {
        for ( LauBkgndDPModelList::iterator iter = bkgndDPModels_.begin();
              iter != bkgndDPModels_.end();
              ++iter ) {
            ( *iter )->initialise();
        }
    }
}

void LauSimpleFitModel::recalculateNormalisation()
{
    //std::cout << "INFO in LauSimpleFitModel::recalculateNormalizationInDPModels : Recalc Norm in DP model" << std::endl;
    sigDPModel_->recalculateNormalisation();
    sigDPModel_->modifyDataTree();
}

void LauSimpleFitModel::setSignalDPParameters()
{
    // Set the fit parameters for the signal model.

    nSigDPPar_ = 0;
    if ( ! this->useDP() ) {
        return;
    }

    std::cout << "INFO in LauSimpleFitModel::setSignalDPParameters : Setting the initial fit parameters for the signal DP model."
              << std::endl;

    // Place isobar coefficient parameters in vector of fit variables
    for ( UInt_t i = 0; i < nSigComp_; i++ ) {
        LauParameterPList pars  = coeffPars_[i]->getParameters();
        nSigDPPar_             += this->addFitParameters( pars, kTRUE );
    }

    // Obtain the resonance parameters and place them in the vector of fit variables and in a separate vector
    LauParameterPList& sigDPPars  = sigDPModel_->getFloatingParameters();
    nSigDPPar_                   += this->addResonanceParameters( sigDPPars );
}

void LauSimpleFitModel::setExtraPdfParameters()
{
    // Include all the parameters of the various PDFs in the fit
    // NB all of them are passed to the fit, even though some have been fixed through parameter.fixed(kTRUE)
    // With the new "cloned parameter" scheme only "original" parameters are passed to the fit.
    // Their clones are updated automatically when the originals are updated.

    nExtraPdfPar_ = 0;

    nExtraPdfPar_ += this->addFitParameters( signalPdfs_ );

    if ( useSCF_ == kTRUE ) {
        nExtraPdfPar_ += this->addFitParameters( scfPdfs_ );
    }

    if ( usingBkgnd_ == kTRUE ) {
        for ( LauBkgndPdfsList::iterator iter = bkgndPdfs_.begin(); iter != bkgndPdfs_.end();
              ++iter ) {
            nExtraPdfPar_ += this->addFitParameters( *iter );
        }
    }
}

void LauSimpleFitModel::setFitNEvents()
{
    if ( signalEvents_ == 0 ) {
        std::cerr << "ERROR in LauSimpleFitModel::setFitNEvents : Signal yield not defined."
                  << std::endl;
        return;
    }
    nNormPar_ = 0;

    // initialise the total number of events to be the sum of all the hypotheses
    Double_t nTotEvts = signalEvents_->value();
    for ( LauBkgndYieldList::const_iterator iter = bkgndEvents_.begin(); iter != bkgndEvents_.end();
          ++iter ) {
        nTotEvts += ( *iter )->value();
        if ( ( *iter ) == 0 ) {
            std::cerr << "ERROR in LauSimpleFitModel::setFitNEvents : Background yield not defined."
                      << std::endl;
            return;
        }
    }
    this->eventsPerExpt( TMath::FloorNint( nTotEvts ) );

    // if doing an extended ML fit add the number of signal events into the fit parameters
    if ( this->doEMLFit() ) {
        std::cout << "INFO in LauSimpleFitModel::setFitNEvents : Initialising number of events for signal and background components..."
                  << std::endl;
        // add the signal events to the list of fit parameters
        nNormPar_ += this->addFitParameters( signalEvents_ );
    } else {
        std::cout << "INFO in LauSimpleFitModel::setFitNEvents : Initialising number of events for background components (and hence signal)..."
                  << std::endl;
    }

    if ( useSCF_ && ! useSCFHist_ ) {
        nNormPar_ += this->addFitParameters( &scfFrac_ );
    }

    if ( usingBkgnd_ == kTRUE ) {
        nNormPar_ += this->addFitParameters( bkgndEvents_ );
    }
}

void LauSimpleFitModel::setExtraNtupleVars()
{
    // Set-up other parameters derived from the fit results, e.g. fit fractions.
    if ( this->useDP() != kTRUE ) {
        return;
    }

    // First clear the vectors so we start from scratch
    this->clearExtraVarVectors();

    LauParameterList& extraVars = this->extraPars();

    // Add a fit fraction for each signal component
    fitFrac_ = sigDPModel_->getFitFractions();
    if ( fitFrac_.size() != nSigComp_ ) {
        std::cerr << "ERROR in LauSimpleFitModel::setExtraNtupleVars : Initial Fit Fraction array of unexpected dimension: "
                  << fitFrac_.size() << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }
    for ( UInt_t i( 0 ); i < nSigComp_; ++i ) {
        if ( fitFrac_[i].size() != nSigComp_ ) {
            std::cerr << "ERROR in LauSimpleFitModel::setExtraNtupleVars : Initial Fit Fraction array of unexpected dimension: "
                      << fitFrac_[i].size() << std::endl;
            gSystem->Exit( EXIT_FAILURE );
        }
    }

    // Add the fit fraction that has not been corrected for the efficiency for each signal component
    fitFracEffUnCorr_ = sigDPModel_->getFitFractionsEfficiencyUncorrected();
    if ( fitFracEffUnCorr_.size() != nSigComp_ ) {
        std::cerr << "ERROR in LauSimpleFitModel::setExtraNtupleVars : Initial Fit Fraction array of unexpected dimension: "
                  << fitFracEffUnCorr_.size() << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }
    for ( UInt_t i( 0 ); i < nSigComp_; ++i ) {
        if ( fitFracEffUnCorr_[i].size() != nSigComp_ ) {
            std::cerr << "ERROR in LauSimpleFitModel::setExtraNtupleVars : Initial Fit Fraction array of unexpected dimension: "
                      << fitFracEffUnCorr_[i].size() << std::endl;
            gSystem->Exit( EXIT_FAILURE );
        }
    }

    for ( UInt_t i = 0; i < nSigComp_; i++ ) {
        for ( UInt_t j = i; j < nSigComp_; j++ ) {
            extraVars.push_back( fitFrac_[i][j] );
            extraVars.push_back( fitFracEffUnCorr_[i][j] );
        }
    }

    // Store any extra parameters/quantities from the DP model (e.g. K-matrix total fit fractions)
    std::vector<LauParameter> extraParams = sigDPModel_->getExtraParameters();
    std::vector<LauParameter>::iterator extraIter;
    for ( extraIter = extraParams.begin(); extraIter != extraParams.end(); ++extraIter ) {
        LauParameter extraParameter = ( *extraIter );
        extraVars.push_back( extraParameter );
    }

    // Now add in the DP efficiency value
    Double_t initMeanEff = sigDPModel_->getMeanEff().initValue();
    meanEff_.value( initMeanEff );
    meanEff_.initValue( initMeanEff );
    meanEff_.genValue( initMeanEff );
    extraVars.push_back( meanEff_ );

    // Also add in the DP rate
    Double_t initDPRate = sigDPModel_->getDPRate().initValue();
    dpRate_.value( initDPRate );
    dpRate_.initValue( initDPRate );
    dpRate_.genValue( initDPRate );
    extraVars.push_back( dpRate_ );
}

void LauSimpleFitModel::finaliseFitResults( const TString& tablePrefixName )
{
    // Retrieve parameters from the fit results for calculations and toy generation
    // and eventually store these in output root ntuples/text files

    // Now take the fit parameters and update them as necessary
    // e.g. to make mag > 0.0 and phase in the right range.
    // This function will also calculate any other values, such as the
    // fit fractions, using any errors provided by fitParErrors as appropriate.
    // Also obtain the pull values: (measured - generated)/(average error)

    if ( this->useDP() == kTRUE ) {
        for ( UInt_t i = 0; i < nSigComp_; i++ ) {
            // Check whether we have mag > 0.0, and phase in the right range
            coeffPars_[i]->finaliseValues();
        }
    }

    // update the pulls on the events
    if ( this->doEMLFit() ) {
        signalEvents_->updatePull();
    }
    if ( useSCF_ && ! useSCFHist_ ) {
        scfFrac_.updatePull();
    }
    if ( usingBkgnd_ == kTRUE ) {
        for ( LauBkgndYieldList::iterator iter = bkgndEvents_.begin(); iter != bkgndEvents_.end();
              ++iter ) {
            std::vector<LauParameter*> parameters = ( *iter )->getPars();
            for ( LauParameter* parameter : parameters ) {
                parameter->updatePull();
            }
        }
    }

    // Update the pulls on all the extra PDFs' parameters
    this->updateFitParameters( signalPdfs_ );
    if ( useSCF_ == kTRUE ) {
        this->updateFitParameters( scfPdfs_ );
    }
    if ( usingBkgnd_ == kTRUE ) {
        for ( LauBkgndPdfsList::iterator iter = bkgndPdfs_.begin(); iter != bkgndPdfs_.end();
              ++iter ) {
            this->updateFitParameters( *iter );
        }
    }

    // Fill the fit results to the ntuple for current experiment

    // update the coefficients and then calculate the fit fractions
    if ( this->useDP() == kTRUE ) {
        this->updateCoeffs();
        sigDPModel_->updateCoeffs( coeffs_ );
        sigDPModel_->calcExtraInfo();

        LauParArray fitFrac = sigDPModel_->getFitFractions();
        if ( fitFrac.size() != nSigComp_ ) {
            std::cerr << "ERROR in LauSimpleFitModel::finaliseFitResults : Fit Fraction array of unexpected dimension: "
                      << fitFrac.size() << std::endl;
            gSystem->Exit( EXIT_FAILURE );
        }
        for ( UInt_t i( 0 ); i < nSigComp_; ++i ) {
            if ( fitFrac[i].size() != nSigComp_ ) {
                std::cerr << "ERROR in LauSimpleFitModel::finaliseFitResults : Fit Fraction array of unexpected dimension: "
                          << fitFrac[i].size() << std::endl;
                gSystem->Exit( EXIT_FAILURE );
            }
        }

        LauParArray fitFracEffUnCorr = sigDPModel_->getFitFractionsEfficiencyUncorrected();
        if ( fitFracEffUnCorr.size() != nSigComp_ ) {
            std::cerr << "ERROR in LauSimpleFitModel::finaliseFitResults : Fit Fraction array of unexpected dimension: "
                      << fitFracEffUnCorr.size() << std::endl;
            gSystem->Exit( EXIT_FAILURE );
        }
        for ( UInt_t i( 0 ); i < nSigComp_; ++i ) {
            if ( fitFracEffUnCorr[i].size() != nSigComp_ ) {
                std::cerr << "ERROR in LauSimpleFitModel::finaliseFitResults : Fit Fraction array of unexpected dimension: "
                          << fitFracEffUnCorr[i].size() << std::endl;
                gSystem->Exit( EXIT_FAILURE );
            }
        }

        for ( UInt_t i( 0 ); i < nSigComp_; ++i ) {
            for ( UInt_t j( i ); j < nSigComp_; ++j ) {
                fitFrac_[i][j].value( fitFrac[i][j].value() );
                fitFracEffUnCorr_[i][j].value( fitFracEffUnCorr[i][j].value() );
            }
        }

        meanEff_.value( sigDPModel_->getMeanEff().value() );
        dpRate_.value( sigDPModel_->getDPRate().value() );

        this->clearExtraVarVectors();
        LauParameterList& extraVars = this->extraPars();

        // Then store the final fit parameters, and any extra parameters for
        // the signal model (e.g. fit fractions)
        for ( UInt_t i( 0 ); i < nSigComp_; ++i ) {
            for ( UInt_t j( i ); j < nSigComp_; ++j ) {
                extraVars.push_back( fitFrac_[i][j] );
                extraVars.push_back( fitFracEffUnCorr_[i][j] );
            }
        }

        // Store any extra parameters/quantities from the DP model (e.g. K-matrix total fit fractions)
        std::vector<LauParameter> extraParams = sigDPModel_->getExtraParameters();
        std::vector<LauParameter>::iterator extraIter;
        for ( extraIter = extraParams.begin(); extraIter != extraParams.end(); ++extraIter ) {
            LauParameter extraParameter = ( *extraIter );
            extraVars.push_back( extraParameter );
        }

        // Now add in the DP efficiency value
        extraVars.push_back( meanEff_ );

        // Also add in the DP rate
        extraVars.push_back( dpRate_ );

        this->printFitFractions( std::cout );
    }

    LauFitNtuple* ntuple = this->fitNtuple();
    ntuple->storeParsAndErrors( this->fitPars(), this->multiDimConstrainedPars(), this->extraPars() );

    // find out the correlation matrix for the parameters
    ntuple->storeCorrMatrix( this->iExpt(), this->fitStatus(), this->covarianceMatrix() );

    // Fill the data into ntuple
    ntuple->updateFitNtuple();

    // Print out the partial fit fractions, phases and the
    // averaged efficiency, reweighted by the dynamics (and anything else)
    if ( this->writeLatexTable() ) {
        TString sigOutFileName( tablePrefixName );
        sigOutFileName += "_";
        sigOutFileName += this->iExpt();
        sigOutFileName += "Expt.tex";
        this->writeOutTable( sigOutFileName );
    }
}

void LauSimpleFitModel::printFitFractions( std::ostream& output )
{
    // Print out Fit Fractions, total DP rate and mean efficiency
    for ( UInt_t i = 0; i < nSigComp_; i++ ) {
        output << "FitFraction for component " << i << " (" << coeffPars_[i]->name()
               << ") = " << fitFrac_[i][i] << std::endl;
    }
    output << "Overall DP rate (integral of matrix element squared) = " << dpRate_ << std::endl;
    output << "Average efficiency weighted by whole DP dynamics = " << meanEff_ << std::endl;
}

void LauSimpleFitModel::writeOutTable( const TString& outputFile )
{
    // Write out the results of the fit to a tex-readable table
    // TODO - need to include the yields in this table
    std::ofstream fout( outputFile );
    LauPrint print;

    std::cout << "INFO in LauSimpleFitModel::writeOutTable : Writing out results of the fit to the tex file "
              << outputFile << std::endl;

    if ( this->useDP() == kTRUE ) {
        // print the fit coefficients in one table
        coeffPars_.front()->printTableHeading( fout );
        for ( UInt_t i = 0; i < nSigComp_; i++ ) {
            coeffPars_[i]->printTableRow( fout );
        }
        fout << "\\hline" << std::endl;
        fout << "\\end{tabular}" << std::endl << std::endl;

        // print the fit fractions in another
        fout << "\\begin{tabular}{|l|c|}" << std::endl;
        fout << "\\hline" << std::endl;
        fout << "Component & FitFraction \\\\" << std::endl;
        fout << "\\hline" << std::endl;
        Double_t fitFracSum( 0.0 );
        for ( UInt_t i = 0; i < nSigComp_; i++ ) {
            TString resName = coeffPars_[i]->name();
            resName         = resName.ReplaceAll( "_", "\\_" );
            fout << resName << "  &  $";
            Double_t fitFrac  = fitFrac_[i][i].value();
            fitFracSum       += fitFrac;
            print.printFormat( fout, fitFrac );
            fout << "$ \\\\" << std::endl;
        }
        fout << "\\hline" << std::endl;

        // Also print out sum of fit fractions
        fout << "Fit Fraction Sum  &  $";
        print.printFormat( fout, fitFracSum );
        fout << "$ \\\\" << std::endl;
        fout << "\\hline" << std::endl;

        fout << "DP rate  &  $";
        print.printFormat( fout, dpRate_.value() );
        fout << "$ \\\\" << std::endl;
        fout << "\\hline" << std::endl;

        fout << "$< \\varepsilon >$  &  $";
        print.printFormat( fout, meanEff_.value() );
        fout << "$ \\\\" << std::endl;
        fout << "\\hline" << std::endl;
        fout << "\\end{tabular}" << std::endl << std::endl;
    }

    if ( ! signalPdfs_.empty() ) {
        fout << "\\begin{tabular}{|l|c|}" << std::endl;
        fout << "\\hline" << std::endl;
        if ( useSCF_ == kTRUE ) {
            fout << "\\Extra TM Signal PDFs' Parameters: & \\\\" << std::endl;
        } else {
            fout << "\\Extra Signal PDFs' Parameters: & \\\\" << std::endl;
        }
        this->printFitParameters( signalPdfs_, fout );
        if ( useSCF_ == kTRUE && ! scfPdfs_.empty() ) {
            fout << "\\hline" << std::endl;
            fout << "\\Extra SCF Signal PDFs' Parameters: & \\\\" << std::endl;
            this->printFitParameters( scfPdfs_, fout );
        }
        if ( usingBkgnd_ == kTRUE && ! bkgndPdfs_.empty() ) {
            fout << "\\hline" << std::endl;
            fout << "\\Extra Background PDFs' Parameters: & \\\\" << std::endl;
            for ( LauBkgndPdfsList::const_iterator iter = bkgndPdfs_.begin();
                  iter != bkgndPdfs_.end();
                  ++iter ) {
                this->printFitParameters( *iter, fout );
            }
        }
        fout << "\\hline \n\\end{tabular}" << std::endl << std::endl;
    }
}

void LauSimpleFitModel::checkInitFitParams()
{
    // Update the number of signal events to be total-sum(background events)
    this->updateSigEvents();

    // Check whether we want to have randomised initial fit parameters for the signal model
    if ( this->useRandomInitFitPars() == kTRUE ) {
        std::cout << "INFO in LauSimpleFitModel::checkInitFitParams : Setting random parameters for the signal DP model"
                  << std::endl;
        this->randomiseInitFitPars();
    }
}

void LauSimpleFitModel::randomiseInitFitPars()
{
    // Only randomise those parameters that are not fixed!
    std::cout << "INFO in LauSimpleFitModel::randomiseInitFitPars : Randomising the initial fit magnitudes and phases of the resonances..."
              << std::endl;

    for ( UInt_t i = 0; i < nSigComp_; i++ ) {
        coeffPars_[i]->randomiseInitValues();
    }
}

std::pair<LauSimpleFitModel::LauGenInfo, Bool_t> LauSimpleFitModel::eventsToGenerate()
{
    // Determine the number of events to generate for each hypothesis
    // If we're smearing then smear each one individually

    LauGenInfo nEvtsGen;

    // Keep track of whether any yield or asymmetry parameters are blinded
    Bool_t blind = kFALSE;

    // Signal
    if ( signalEvents_->blind() ) {
        blind = kTRUE;
    }

    Double_t evtWeight( 1.0 );
    Double_t nEvts = signalEvents_->genValue();
    if ( nEvts < 0 ) {
        evtWeight = -1.0;
        nEvts     = TMath::Abs( nEvts );
    }

    Int_t nEvtsToGen { static_cast<Int_t>( nEvts ) };
    if ( this->doPoissonSmearing() ) {
        nEvtsToGen = LauRandom::randomFun()->Poisson( nEvts );
    }

    nEvtsGen["signal"] = std::make_pair( nEvtsToGen, evtWeight );

    // Backgrounds
    const UInt_t nBkgnds = this->nBkgndClasses();
    for ( UInt_t bkgndID( 0 ); bkgndID < nBkgnds; ++bkgndID ) {
        const LauAbsRValue* evtsPar = bkgndEvents_[bkgndID];
        if ( evtsPar->blind() ) {
            blind = kTRUE;
        }

        evtWeight = 1.0;
        nEvts     = evtsPar->genValue();
        if ( nEvts < 0 ) {
            evtWeight = -1.0;
            nEvts     = TMath::Abs( nEvts );
        }

        nEvtsToGen = static_cast<Int_t>( nEvts );
        if ( this->doPoissonSmearing() ) {
            nEvtsToGen = LauRandom::randomFun()->Poisson( nEvts );
        }

        const TString& bkgndClass = this->bkgndClassName( bkgndID );
        nEvtsGen[bkgndClass]      = std::make_pair( nEvtsToGen, evtWeight );
    }

    return std::make_pair( nEvtsGen, blind );
}

Bool_t LauSimpleFitModel::genExpt()
{
    // Routine to generate toy Monte Carlo events according to the various models we have defined.

    // Determine the number of events to generate for each hypothesis
    std::pair<LauGenInfo, Bool_t> info = this->eventsToGenerate();
    LauGenInfo nEvts                   = info.first;
    const Bool_t blind                 = info.second;

    Bool_t genOK( kTRUE );
    Int_t evtNum( 0 );

    const UInt_t nBkgnds = this->nBkgndClasses();
    std::vector<TString> bkgndClassNames( nBkgnds );
    std::vector<TString> bkgndClassNamesGen( nBkgnds );
    for ( UInt_t iBkgnd( 0 ); iBkgnd < nBkgnds; ++iBkgnd ) {
        TString name( this->bkgndClassName( iBkgnd ) );
        bkgndClassNames[iBkgnd]    = name;
        bkgndClassNamesGen[iBkgnd] = "gen" + name;
    }

    const Bool_t storeSCFTruthInfo = ( useSCF_ || ( this->enableEmbedding() && signalTree_ != 0 &&
                                                    signalTree_->haveBranch( "mcMatch" ) ) );

    // Loop over the hypotheses and generate the requested number of events for each one
    for ( LauGenInfo::const_iterator iter = nEvts.begin(); iter != nEvts.end(); ++iter ) {

        const TString& type( iter->first );
        Int_t nEvtsGen( iter->second.first );
        Double_t evtWeight( iter->second.second );

        for ( Int_t iEvt( 0 ); iEvt < nEvtsGen; ++iEvt ) {

            this->setGenNtupleDoubleBranchValue( "evtWeight", evtWeight );
            // Add efficiency information
            this->setGenNtupleDoubleBranchValue( "efficiency", 1 );

            if ( type == "signal" ) {
                this->setGenNtupleIntegerBranchValue( "genSig", 1 );
                for ( UInt_t iBkgnd( 0 ); iBkgnd < nBkgnds; ++iBkgnd ) {
                    this->setGenNtupleIntegerBranchValue( bkgndClassNamesGen[iBkgnd], 0 );
                }
                genOK = this->generateSignalEvent();
                this->setGenNtupleDoubleBranchValue( "efficiency", sigDPModel_->getEvtEff() );
            } else {
                this->setGenNtupleIntegerBranchValue( "genSig", 0 );
                if ( storeSCFTruthInfo ) {
                    this->setGenNtupleIntegerBranchValue( "genTMSig", 0 );
                    this->setGenNtupleIntegerBranchValue( "genSCFSig", 0 );
                }
                UInt_t bkgndID( 0 );
                for ( UInt_t iBkgnd( 0 ); iBkgnd < nBkgnds; ++iBkgnd ) {
                    Int_t gen( 0 );
                    if ( bkgndClassNames[iBkgnd] == type ) {
                        gen     = 1;
                        bkgndID = iBkgnd;
                    }
                    this->setGenNtupleIntegerBranchValue( bkgndClassNamesGen[iBkgnd], gen );
                }
                genOK = this->generateBkgndEvent( bkgndID );
            }

            if ( ! genOK ) {
                // If there was a problem with the generation then break out and return.
                // The problem model will have adjusted itself so that all should be OK next time.
                break;
            }

            if ( this->useDP() == kTRUE ) {
                this->setDPBranchValues();
            }

            // Store the event number (within this experiment)
            this->setGenNtupleIntegerBranchValue( "iEvtWithinExpt", evtNum );
            // and then increment it
            ++evtNum;

            this->fillGenNtupleBranches();
            if ( ! blind && ( iEvt % 500 == 0 ) ) {
                std::cout << "INFO in LauSimpleFitModel::genExpt : Generated event number " << iEvt
                          << " out of " << nEvtsGen << " " << type << " events." << std::endl;
            }
        }

        if ( ! genOK ) {
            break;
        }
    }

    if ( this->useDP() && genOK ) {

        sigDPModel_->checkToyMC( kTRUE, kTRUE );

        // Get the fit fractions if they're to be written into the latex table
        if ( this->writeLatexTable() ) {
            LauParArray fitFrac = sigDPModel_->getFitFractions();
            if ( fitFrac.size() != nSigComp_ ) {
                std::cerr << "ERROR in LauSimpleFitModel::genExpt : Fit Fraction array of unexpected dimension: "
                          << fitFrac.size() << std::endl;
                gSystem->Exit( EXIT_FAILURE );
            }
            for ( UInt_t i( 0 ); i < nSigComp_; ++i ) {
                if ( fitFrac[i].size() != nSigComp_ ) {
                    std::cerr << "ERROR in LauSimpleFitModel::genExpt : Fit Fraction array of unexpected dimension: "
                              << fitFrac[i].size() << std::endl;
                    gSystem->Exit( EXIT_FAILURE );
                }
            }

            for ( UInt_t i( 0 ); i < nSigComp_; ++i ) {
                for ( UInt_t j = i; j < nSigComp_; j++ ) {
                    fitFrac_[i][j].value( fitFrac[i][j].value() );
                }
            }
            meanEff_.value( sigDPModel_->getMeanEff().value() );
            dpRate_.value( sigDPModel_->getDPRate().value() );
        }
    }

    // If we're reusing embedded events or if the generation is being
    // reset then clear the lists of used events
    if ( reuseSignal_ || ! genOK ) {
        if ( signalTree_ ) {
            signalTree_->clearUsedList();
        }
    }
    for ( UInt_t bkgndID( 0 ); bkgndID < nBkgnds; ++bkgndID ) {
        LauEmbeddedData* data = bkgndTree_[bkgndID];
        if ( reuseBkgnd_[bkgndID] || ! genOK ) {
            if ( data ) {
                data->clearUsedList();
            }
        }
    }

    return genOK;
}

Bool_t LauSimpleFitModel::generateSignalEvent()
{
    // Generate signal event
    Bool_t genOK( kTRUE );
    Bool_t genSCF( kFALSE );

    if ( this->useDP() ) {
        if ( this->enableEmbedding() && signalTree_ ) {
            if ( useReweighting_ ) {
                // Select a (random) event from the generated data. Then store the
                // reconstructed DP co-ords, together with other pdf information,
                // as the embedded data.
                genOK = signalTree_->getReweightedEvent( sigDPModel_ );
            } else {

                signalTree_->getEmbeddedEvent( kinematics_ );
            }
            if ( signalTree_->haveBranch( "mcMatch" ) ) {
                Int_t match = static_cast<Int_t>( signalTree_->getValue( "mcMatch" ) );
                if ( match ) {
                    this->setGenNtupleIntegerBranchValue( "genTMSig", 1 );
                    this->setGenNtupleIntegerBranchValue( "genSCFSig", 0 );
                    genSCF = kFALSE;
                } else {
                    this->setGenNtupleIntegerBranchValue( "genTMSig", 0 );
                    this->setGenNtupleIntegerBranchValue( "genSCFSig", 1 );
                    genSCF = kTRUE;
                }
            }
        } else {
            genOK = sigDPModel_->generate();
            if ( genOK && useSCF_ ) {
                Double_t frac( 0.0 );
                if ( useSCFHist_ ) {
                    frac = scfFracHist_->calcEfficiency( kinematics_ );
                } else {
                    frac = scfFrac_.genValue();
                }
                if ( frac < LauRandom::randomFun()->Rndm() ) {
                    this->setGenNtupleIntegerBranchValue( "genTMSig", 1 );
                    this->setGenNtupleIntegerBranchValue( "genSCFSig", 0 );
                    genSCF = kFALSE;
                } else {
                    this->setGenNtupleIntegerBranchValue( "genTMSig", 0 );
                    this->setGenNtupleIntegerBranchValue( "genSCFSig", 1 );
                    genSCF = kTRUE;

                    // Optionally smear the DP position
                    // of the SCF event
                    if ( scfMap_ != 0 ) {
                        Double_t xCoord( 0.0 ), yCoord( 0.0 );
                        if ( kinematics_->squareDP() ) {
                            xCoord = kinematics_->getmPrime();
                            yCoord = kinematics_->getThetaPrime();
                        } else {
                            xCoord = kinematics_->getm13Sq();
                            yCoord = kinematics_->getm23Sq();
                        }

                        // Find the bin number where this event is generated
                        Int_t binNo = scfMap_->binNumber( xCoord, yCoord );

                        // Retrieve the migration histogram
                        TH2* histo = scfMap_->trueHist( binNo );

                        const LauAbsEffModel* effModel = sigDPModel_->getEffModel();
                        do {
                            // Get a random point from the histogram
                            histo->GetRandom2( xCoord, yCoord );

                            // Update the kinematics
                            if ( kinematics_->squareDP() ) {
                                kinematics_->updateSqDPKinematics( xCoord, yCoord );
                            } else {
                                kinematics_->updateKinematics( xCoord, yCoord );
                            }
                        } while ( ! effModel->passVeto( kinematics_ ) );
                    }
                }
            }
        }
    } else {
        if ( this->enableEmbedding() && signalTree_ ) {
            signalTree_->getEmbeddedEvent( 0 );
            if ( signalTree_->haveBranch( "mcMatch" ) ) {
                Int_t match = static_cast<Int_t>( signalTree_->getValue( "mcMatch" ) );
                if ( match ) {
                    this->setGenNtupleIntegerBranchValue( "genTMSig", 1 );
                    this->setGenNtupleIntegerBranchValue( "genSCFSig", 0 );
                    genSCF = kFALSE;
                } else {
                    this->setGenNtupleIntegerBranchValue( "genTMSig", 0 );
                    this->setGenNtupleIntegerBranchValue( "genSCFSig", 1 );
                    genSCF = kTRUE;
                }
            }
        } else if ( useSCF_ ) {
            Double_t frac = scfFrac_.genValue();
            if ( frac < LauRandom::randomFun()->Rndm() ) {
                this->setGenNtupleIntegerBranchValue( "genTMSig", 1 );
                this->setGenNtupleIntegerBranchValue( "genSCFSig", 0 );
                genSCF = kFALSE;
            } else {
                this->setGenNtupleIntegerBranchValue( "genTMSig", 0 );
                this->setGenNtupleIntegerBranchValue( "genSCFSig", 1 );
                genSCF = kTRUE;
            }
        }
    }
    if ( genOK ) {
        if ( useSCF_ ) {
            if ( genSCF ) {
                this->generateExtraPdfValues( &scfPdfs_, signalTree_ );
            } else {
                this->generateExtraPdfValues( &signalPdfs_, signalTree_ );
            }
        } else {
            this->generateExtraPdfValues( &signalPdfs_, signalTree_ );
        }
    }
    // Check for problems with the embedding
    if ( this->enableEmbedding() && signalTree_ &&
         ( signalTree_->nEvents() == signalTree_->nUsedEvents() ) ) {
        std::cerr << "WARNING in LauSimpleFitModel::generateSignalEvent : Source of embedded signal events used up, clearing the list of used events."
                  << std::endl;
        signalTree_->clearUsedList();
    }
    return genOK;
}

Bool_t LauSimpleFitModel::generateBkgndEvent( UInt_t bkgndID )
{
    // Generate background event
    Bool_t genOK( kTRUE );

    LauAbsBkgndDPModel* model = bkgndDPModels_[bkgndID];
    LauEmbeddedData* embeddedData( 0 );
    if ( this->enableEmbedding() ) {
        embeddedData = bkgndTree_[bkgndID];
    }
    LauPdfPList* extraPdfs = &bkgndPdfs_[bkgndID];
    if ( this->useDP() ) {
        if ( embeddedData ) {
            embeddedData->getEmbeddedEvent( kinematics_ );
        } else {
            if ( model == 0 ) {
                std::cerr << "ERROR in LauSimpleFitModel::generateBkgndEvent : Can't find the DP model for background class \""
                          << this->bkgndClassName( bkgndID ) << "\"." << std::endl;
                gSystem->Exit( EXIT_FAILURE );
            }
            genOK = model->generate();
        }
    } else {
        if ( embeddedData ) {
            embeddedData->getEmbeddedEvent( 0 );
        }
    }
    if ( genOK ) {
        this->generateExtraPdfValues( extraPdfs, embeddedData );
    }
    // Check for problems with the embedding
    if ( embeddedData && ( embeddedData->nEvents() == embeddedData->nUsedEvents() ) ) {
        std::cerr << "WARNING in LauSimpleFitModel::generateBkgndEvent : Source of embedded "
                  << this->bkgndClassName( bkgndID )
                  << " events used up, clearing the list of used events." << std::endl;
        embeddedData->clearUsedList();
    }
    return genOK;
}

void LauSimpleFitModel::setupGenNtupleBranches()
{
    // Setup the required ntuple branches
    this->addGenNtupleDoubleBranch( "evtWeight" );
    this->addGenNtupleIntegerBranch( "genSig" );
    this->addGenNtupleDoubleBranch( "efficiency" );
    if ( useSCF_ ||
         ( this->enableEmbedding() && signalTree_ != 0 && signalTree_->haveBranch( "mcMatch" ) ) ) {
        this->addGenNtupleIntegerBranch( "genTMSig" );
        this->addGenNtupleIntegerBranch( "genSCFSig" );
    }
    const UInt_t nBkgnds = this->nBkgndClasses();
    for ( UInt_t iBkgnd( 0 ); iBkgnd < nBkgnds; ++iBkgnd ) {
        TString name( this->bkgndClassName( iBkgnd ) );
        name.Prepend( "gen" );
        this->addGenNtupleIntegerBranch( name );
    }
    if ( this->useDP() == kTRUE ) {
        this->addGenNtupleDoubleBranch( "m12" );
        this->addGenNtupleDoubleBranch( "m23" );
        this->addGenNtupleDoubleBranch( "m13" );
        this->addGenNtupleDoubleBranch( "m12Sq" );
        this->addGenNtupleDoubleBranch( "m23Sq" );
        this->addGenNtupleDoubleBranch( "m13Sq" );
        this->addGenNtupleDoubleBranch( "cosHel12" );
        this->addGenNtupleDoubleBranch( "cosHel23" );
        this->addGenNtupleDoubleBranch( "cosHel13" );
        if ( kinematics_->squareDP() ) {
            this->addGenNtupleDoubleBranch( "mPrime" );
            this->addGenNtupleDoubleBranch( "thPrime" );
        }
    }
    for ( LauPdfPList::const_iterator pdf_iter = signalPdfs_.begin(); pdf_iter != signalPdfs_.end();
          ++pdf_iter ) {
        std::vector<TString> varNames = ( *pdf_iter )->varNames();
        for ( std::vector<TString>::const_iterator var_iter = varNames.begin();
              var_iter != varNames.end();
              ++var_iter ) {
            if ( ( *var_iter ) != "m13Sq" && ( *var_iter ) != "m23Sq" ) {
                this->addGenNtupleDoubleBranch( ( *var_iter ) );
            }
        }
    }
}

void LauSimpleFitModel::setDPBranchValues()
{
    // Store all the DP information
    this->setGenNtupleDoubleBranchValue( "m12", kinematics_->getm12() );
    this->setGenNtupleDoubleBranchValue( "m23", kinematics_->getm23() );
    this->setGenNtupleDoubleBranchValue( "m13", kinematics_->getm13() );
    this->setGenNtupleDoubleBranchValue( "m12Sq", kinematics_->getm12Sq() );
    this->setGenNtupleDoubleBranchValue( "m23Sq", kinematics_->getm23Sq() );
    this->setGenNtupleDoubleBranchValue( "m13Sq", kinematics_->getm13Sq() );
    this->setGenNtupleDoubleBranchValue( "cosHel12", kinematics_->getc12() );
    this->setGenNtupleDoubleBranchValue( "cosHel23", kinematics_->getc23() );
    this->setGenNtupleDoubleBranchValue( "cosHel13", kinematics_->getc13() );
    if ( kinematics_->squareDP() ) {
        this->setGenNtupleDoubleBranchValue( "mPrime", kinematics_->getmPrime() );
        this->setGenNtupleDoubleBranchValue( "thPrime", kinematics_->getThetaPrime() );
    }
}

void LauSimpleFitModel::generateExtraPdfValues( LauPdfPList* extraPdfs, LauEmbeddedData* embeddedData )
{
    // Generate from the extra PDFs
    if ( extraPdfs ) {
        for ( LauPdfPList::iterator pdf_iter = extraPdfs->begin(); pdf_iter != extraPdfs->end();
              ++pdf_iter ) {
            LauFitData genValues;
            if ( embeddedData ) {
                genValues = embeddedData->getValues( ( *pdf_iter )->varNames() );
            } else {
                genValues = ( *pdf_iter )->generate( kinematics_ );
            }
            for ( LauFitData::const_iterator var_iter = genValues.begin();
                  var_iter != genValues.end();
                  ++var_iter ) {
                TString varName = var_iter->first;
                if ( varName != "m13Sq" && varName != "m23Sq" ) {
                    Double_t value = var_iter->second;
                    this->setGenNtupleDoubleBranchValue( varName, value );
                }
            }
        }
    }
}

void LauSimpleFitModel::propagateParUpdates()
{
    // Update the total normalisation for the signal likelihood
    if ( this->useDP() == kTRUE ) {
        this->updateCoeffs();
        sigDPModel_->updateCoeffs( coeffs_ );
    }

    // Update the signal events from the background events if not doing an extended fit
    if ( ! this->doEMLFit() && ! signalEvents_->fixed() ) {
        this->updateSigEvents();
    }
}

void LauSimpleFitModel::updateSigEvents()
{
    // The background parameters will have been set from Minuit.
    // We need to update the signal events using these.
    Double_t nTotEvts = this->eventsPerExpt();

    signalEvents_->range( -2.0 * nTotEvts, 2.0 * nTotEvts );
    for ( LauBkgndYieldList::iterator iter = bkgndEvents_.begin(); iter != bkgndEvents_.end();
          ++iter ) {
        LauAbsRValue* nBkgndEvents = ( *iter );
        if ( nBkgndEvents->isLValue() ) {
            LauParameter* yield = dynamic_cast<LauParameter*>( nBkgndEvents );
            yield->range( -2.0 * nTotEvts, 2.0 * nTotEvts );
        }
    }

    if ( signalEvents_->fixed() ) {
        return;
    }

    // Subtract background events (if any) from signal.
    Double_t signalEvents = nTotEvts;
    if ( usingBkgnd_ == kTRUE ) {
        for ( LauBkgndYieldList::const_iterator iter = bkgndEvents_.begin();
              iter != bkgndEvents_.end();
              ++iter ) {
            signalEvents -= ( *iter )->value();
        }
    }
    signalEvents_->value( signalEvents );
}

void LauSimpleFitModel::cacheInputFitVars()
{
    // Fill the internal data trees of the signal and backgrond models.

    LauFitDataTree* inputFitData = this->fitData();

    // First the Dalitz plot variables (m_ij^2)
    if ( this->useDP() == kTRUE ) {

        // need to append SCF smearing bins before caching DP amplitudes
        if ( scfMap_ != 0 ) {
            this->appendBinCentres( inputFitData );
        }
        sigDPModel_->fillDataTree( *inputFitData );

        if ( usingBkgnd_ == kTRUE ) {
            for ( LauBkgndDPModelList::iterator iter = bkgndDPModels_.begin();
                  iter != bkgndDPModels_.end();
                  ++iter ) {
                ( *iter )->fillDataTree( *inputFitData );
            }
        }
    }

    // ...and then the extra PDFs
    this->cacheInfo( signalPdfs_, *inputFitData );
    this->cacheInfo( scfPdfs_, *inputFitData );
    for ( LauBkgndPdfsList::iterator iter = bkgndPdfs_.begin(); iter != bkgndPdfs_.end(); ++iter ) {
        this->cacheInfo( ( *iter ), *inputFitData );
    }

    // and finally the SCF fractions and jacobians
    if ( useSCF_ && useSCFHist_ ) {
        if ( ! inputFitData->haveBranch( "m13Sq" ) || ! inputFitData->haveBranch( "m23Sq" ) ) {
            std::cerr << "ERROR in LauSimpleFitModel::cacheInputFitVars : Input data does not contain DP branches and so can't cache the SCF fraction."
                      << std::endl;
            gSystem->Exit( EXIT_FAILURE );
        }

        UInt_t nEvents = inputFitData->nEvents();
        recoSCFFracs_.clear();
        recoSCFFracs_.reserve( nEvents );
        if ( kinematics_->squareDP() ) {
            recoJacobians_.clear();
            recoJacobians_.reserve( nEvents );
        }
        for ( UInt_t iEvt = 0; iEvt < nEvents; iEvt++ ) {
            const LauFitData& dataValues        = inputFitData->getData( iEvt );
            LauFitData::const_iterator m13_iter = dataValues.find( "m13Sq" );
            LauFitData::const_iterator m23_iter = dataValues.find( "m23Sq" );
            kinematics_->updateKinematics( m13_iter->second, m23_iter->second );
            Double_t scfFrac = scfFracHist_->calcEfficiency( kinematics_ );
            recoSCFFracs_.push_back( scfFrac );
            if ( kinematics_->squareDP() ) {
                recoJacobians_.push_back( kinematics_->calcSqDPJacobian() );
            }
        }
    }
}

void LauSimpleFitModel::appendBinCentres( LauFitDataTree* inputData )
{
    // We'll be caching the DP amplitudes and efficiencies of the centres of the true bins.
    // To do so, we attach some fake points at the end of inputData, the number of the entry
    // minus the total number of events corresponding to the number of the histogram for that
    // given true bin in the LauScfMap object. (What this means is that when Laura is provided with
    // the LauScfMap object by the user, it's the latter who has to make sure that it contains the
    // right number of histograms and in exactly the right order!)

    // Get the x and y co-ordinates of the bin centres
    std::vector<Double_t> binCentresXCoords;
    std::vector<Double_t> binCentresYCoords;
    scfMap_->listBinCentres( binCentresXCoords, binCentresYCoords );

    // The SCF histograms could be in square Dalitz plot histograms.
    // The dynamics takes normal Dalitz plot coords, so we might have to convert them back.
    Bool_t sqDP  = kinematics_->squareDP();
    UInt_t nBins = binCentresXCoords.size();
    fakeSCFFracs_.clear();
    fakeSCFFracs_.reserve( nBins );
    if ( sqDP ) {
        fakeJacobians_.clear();
        fakeJacobians_.reserve( nBins );
    }

    for ( UInt_t iBin = 0; iBin < nBins; ++iBin ) {

        if ( sqDP ) {
            kinematics_->updateSqDPKinematics( binCentresXCoords[iBin], binCentresYCoords[iBin] );
            binCentresXCoords[iBin] = kinematics_->getm13Sq();
            binCentresYCoords[iBin] = kinematics_->getm23Sq();
            fakeJacobians_.push_back( kinematics_->calcSqDPJacobian() );
        } else {
            kinematics_->updateKinematics( binCentresXCoords[iBin], binCentresYCoords[iBin] );
        }

        fakeSCFFracs_.push_back( scfFracHist_->calcEfficiency( kinematics_ ) );
    }

    // Set up inputFitVars_ object to hold the fake events
    inputData->appendFakePoints( binCentresXCoords, binCentresYCoords );
}

Double_t LauSimpleFitModel::getTotEvtLikelihood( UInt_t iEvt )
{
    // Get the DP likelihood for signal and backgrounds
    this->getEvtDPLikelihood( iEvt );

    // Get the combined extra PDFs likelihood for signal and backgrounds
    this->getEvtExtraLikelihoods( iEvt );

    // If appropriate, combine the TM and SCF likelihoods
    Double_t sigLike = sigDPLike_ * sigExtraLike_;
    if ( useSCF_ ) {
        Double_t scfFrac( 0.0 );
        if ( useSCFHist_ ) {
            scfFrac = recoSCFFracs_[iEvt];
        } else {
            scfFrac = scfFrac_.unblindValue();
        }
        sigLike *= ( 1.0 - scfFrac );
        if ( ( scfMap_ != 0 ) && ( this->useDP() == kTRUE ) ) {
            // if we're smearing the SCF DP PDF then the SCF frac
            // is already included in the SCF DP likelihood
            sigLike += ( scfDPLike_ * scfExtraLike_ );
        } else {
            sigLike += ( scfFrac * scfDPLike_ * scfExtraLike_ );
        }
    }

    // Construct the total event likelihood
    Double_t likelihood  = signalEvents_->unblindValue() * sigLike;
    const UInt_t nBkgnds = this->nBkgndClasses();
    for ( UInt_t bkgndID( 0 ); bkgndID < nBkgnds; ++bkgndID ) {
        likelihood += ( bkgndEvents_[bkgndID]->unblindValue() * bkgndDPLike_[bkgndID] *
                        bkgndExtraLike_[bkgndID] );
    }

    return likelihood;
}

Double_t LauSimpleFitModel::getEventSum() const
{
    Double_t eventSum( 0.0 );
    eventSum += signalEvents_->unblindValue();
    for ( LauBkgndYieldList::const_iterator iter = bkgndEvents_.begin(); iter != bkgndEvents_.end();
          ++iter ) {
        eventSum += ( *iter )->unblindValue();
    }
    return eventSum;
}

void LauSimpleFitModel::getEvtDPLikelihood( UInt_t iEvt )
{
    // Function to return the signal and background likelihoods for the
    // Dalitz plot for the given event evtNo.

    if ( this->useDP() == kTRUE ) {
        sigDPModel_->calcLikelihoodInfo( iEvt );
        sigDPLike_ = sigDPModel_->getEvtLikelihood();

        if ( useSCF_ == kTRUE ) {
            if ( scfMap_ == 0 ) {
                // we're not smearing the SCF DP position
                // so the likelihood is the same as the TM
                scfDPLike_ = sigDPLike_;
            } else {
                // calculate the smeared SCF DP likelihood
                scfDPLike_ = this->getEvtSCFDPLikelihood( iEvt );
            }
        }

        const UInt_t nBkgnds = this->nBkgndClasses();
        for ( UInt_t bkgndID( 0 ); bkgndID < nBkgnds; ++bkgndID ) {
            if ( usingBkgnd_ == kTRUE ) {
                bkgndDPLike_[bkgndID] = bkgndDPModels_[bkgndID]->getLikelihood( iEvt );
            } else {
                bkgndDPLike_[bkgndID] = 0.0;
            }
        }
    } else {
        // There's always going to be a term in the likelihood for the
        // signal, so we'd better not zero it.
        sigDPLike_ = 1.0;
        scfDPLike_ = 1.0;

        const UInt_t nBkgnds = this->nBkgndClasses();
        for ( UInt_t bkgndID( 0 ); bkgndID < nBkgnds; ++bkgndID ) {
            if ( usingBkgnd_ == kTRUE ) {
                bkgndDPLike_[bkgndID] = 1.0;
            } else {
                bkgndDPLike_[bkgndID] = 0.0;
            }
        }
    }
}

Double_t LauSimpleFitModel::getEvtSCFDPLikelihood( UInt_t iEvt )
{
    Double_t scfDPLike( 0.0 );

    Double_t recoJacobian( 1.0 );
    Double_t xCoord( 0.0 );
    Double_t yCoord( 0.0 );
    Bool_t squareDP = kinematics_->squareDP();
    if ( squareDP ) {
        xCoord       = sigDPModel_->getEvtmPrime();
        yCoord       = sigDPModel_->getEvtthPrime();
        recoJacobian = recoJacobians_[iEvt];
    } else {
        xCoord = sigDPModel_->getEvtm13Sq();
        yCoord = sigDPModel_->getEvtm23Sq();
    }

    // Find the bin that our reco event falls in
    Int_t recoBin = scfMap_->binNumber( xCoord, yCoord );

    // Find out which true Bins contribute to the given reco bin
    const std::vector<Int_t>* trueBins = scfMap_->trueBins( recoBin );

    Int_t nDataEvents = this->eventsPerExpt();

    // Loop over the true bins
    for ( std::vector<Int_t>::const_iterator iter = trueBins->begin(); iter != trueBins->end();
          ++iter ) {
        Int_t trueBin = ( *iter );

        // prob of a true event in the given true bin migrating to the reco bin
        Double_t pRecoGivenTrue = scfMap_->prob( recoBin, trueBin );

        // We've cached the DP amplitudes and the efficiency for the
        // true bin centres, just after the data points
        sigDPModel_->calcLikelihoodInfo( nDataEvents + trueBin );
        Double_t pTrue = sigDPModel_->getEvtLikelihood();

        // Get the cached SCF fraction (and jacobian if we're using the square DP)
        Double_t scfFraction = fakeSCFFracs_[trueBin];
        Double_t jacobian( 1.0 );
        if ( squareDP ) {
            jacobian = fakeJacobians_[trueBin];
        }

        scfDPLike += pTrue * jacobian * scfFraction * pRecoGivenTrue;
    }

    // Divide by the reco jacobian
    scfDPLike /= recoJacobian;

    return scfDPLike;
}

void LauSimpleFitModel::getEvtExtraLikelihoods( UInt_t iEvt )
{
    // Function to return the signal and background likelihoods for the
    // extra variables for the given event evtNo.

    sigExtraLike_ = 1.0;

    for ( LauPdfPList::iterator iter = signalPdfs_.begin(); iter != signalPdfs_.end(); ++iter ) {
        ( *iter )->calcLikelihoodInfo( iEvt );
        sigExtraLike_ *= ( *iter )->getLikelihood();
    }
    if ( useSCF_ ) {
        scfExtraLike_ = 1.0;
        for ( LauPdfPList::iterator iter = scfPdfs_.begin(); iter != scfPdfs_.end(); ++iter ) {
            ( *iter )->calcLikelihoodInfo( iEvt );
            scfExtraLike_ *= ( *iter )->getLikelihood();
        }
    }

    const UInt_t nBkgnds = this->nBkgndClasses();
    for ( UInt_t bkgndID( 0 ); bkgndID < nBkgnds; ++bkgndID ) {
        if ( usingBkgnd_ ) {
            bkgndExtraLike_[bkgndID] = 1.0;
            LauPdfPList& pdfList     = bkgndPdfs_[bkgndID];
            for ( LauPdfPList::iterator pdf_iter = pdfList.begin(); pdf_iter != pdfList.end();
                  ++pdf_iter ) {
                ( *pdf_iter )->calcLikelihoodInfo( iEvt );
                bkgndExtraLike_[bkgndID] *= ( *pdf_iter )->getLikelihood();
            }
        } else {
            bkgndExtraLike_[bkgndID] = 0.0;
        }
    }
}

void LauSimpleFitModel::updateCoeffs()
{
    coeffs_.clear();
    coeffs_.reserve( nSigComp_ );
    for ( UInt_t i = 0; i < nSigComp_; i++ ) {
        coeffs_.push_back( coeffPars_[i]->particleCoeff() );
    }
}

void LauSimpleFitModel::setupSPlotNtupleBranches()
{
    // add branches for storing the experiment number and the number of
    // the event within the current experiment
    this->addSPlotNtupleIntegerBranch( "iExpt" );
    this->addSPlotNtupleIntegerBranch( "iEvtWithinExpt" );

    // Store the efficiency of the event (for inclusive BF calculations).
    if ( this->storeDPEff() ) {
        this->addSPlotNtupleDoubleBranch( "efficiency" );
        if ( sigDPModel_->usingScfModel() ) {
            this->addSPlotNtupleDoubleBranch( "scffraction" );
        }
    }

    // Store the total event likelihood for each species.
    if ( useSCF_ ) {
        this->addSPlotNtupleDoubleBranch( "sigTMTotalLike" );
        this->addSPlotNtupleDoubleBranch( "sigSCFTotalLike" );
        this->addSPlotNtupleDoubleBranch( "sigSCFFrac" );
    } else {
        this->addSPlotNtupleDoubleBranch( "sigTotalLike" );
    }
    if ( usingBkgnd_ ) {
        const UInt_t nBkgnds = this->nBkgndClasses();
        for ( UInt_t iBkgnd( 0 ); iBkgnd < nBkgnds; ++iBkgnd ) {
            TString name( this->bkgndClassName( iBkgnd ) );
            name += "TotalLike";
            this->addSPlotNtupleDoubleBranch( name );
        }
    }

    // Store the DP likelihoods
    if ( this->useDP() ) {
        if ( useSCF_ ) {
            this->addSPlotNtupleDoubleBranch( "sigTMDPLike" );
            this->addSPlotNtupleDoubleBranch( "sigSCFDPLike" );
        } else {
            this->addSPlotNtupleDoubleBranch( "sigDPLike" );
        }
        if ( usingBkgnd_ ) {
            const UInt_t nBkgnds = this->nBkgndClasses();
            for ( UInt_t iBkgnd( 0 ); iBkgnd < nBkgnds; ++iBkgnd ) {
                TString name( this->bkgndClassName( iBkgnd ) );
                name += "DPLike";
                this->addSPlotNtupleDoubleBranch( name );
            }
        }
    }

    // Store the likelihoods for each extra PDF
    if ( useSCF_ ) {
        this->addSPlotNtupleBranches( &signalPdfs_, "sigTM" );
        this->addSPlotNtupleBranches( &scfPdfs_, "sigSCF" );
    } else {
        this->addSPlotNtupleBranches( &signalPdfs_, "sig" );
    }
    if ( usingBkgnd_ ) {
        const UInt_t nBkgnds = this->nBkgndClasses();
        for ( UInt_t iBkgnd( 0 ); iBkgnd < nBkgnds; ++iBkgnd ) {
            const TString& bkgndClass  = this->bkgndClassName( iBkgnd );
            const LauPdfPList* pdfList = &( bkgndPdfs_[iBkgnd] );
            this->addSPlotNtupleBranches( pdfList, bkgndClass );
        }
    }
}

void LauSimpleFitModel::addSPlotNtupleBranches( const LauPdfPList* extraPdfs, const TString& prefix )
{
    if ( extraPdfs ) {
        // Loop through each of the PDFs
        for ( LauPdfPList::const_iterator pdf_iter = extraPdfs->begin(); pdf_iter != extraPdfs->end();
              ++pdf_iter ) {

            // Count the number of input variables that are not
            // DP variables (used in the case where there is DP
            // dependence for e.g. MVA)
            UInt_t nVars( 0 );
            std::vector<TString> varNames = ( *pdf_iter )->varNames();
            for ( std::vector<TString>::const_iterator var_iter = varNames.begin();
                  var_iter != varNames.end();
                  ++var_iter ) {
                if ( ( *var_iter ) != "m13Sq" && ( *var_iter ) != "m23Sq" ) {
                    ++nVars;
                }
            }

            if ( nVars == 1 ) {
                // If the PDF only has one variable then
                // simply add one branch for that variable
                TString varName = ( *pdf_iter )->varName();
                TString name( prefix );
                name += varName;
                name += "Like";
                this->addSPlotNtupleDoubleBranch( name );
            } else if ( nVars == 2 ) {
                // If the PDF has two variables then we
                // need a branch for them both together and
                // branches for each
                TString allVars( "" );
                for ( std::vector<TString>::const_iterator var_iter = varNames.begin();
                      var_iter != varNames.end();
                      ++var_iter ) {
                    allVars += ( *var_iter );
                    TString name( prefix );
                    name += ( *var_iter );
                    name += "Like";
                    this->addSPlotNtupleDoubleBranch( name );
                }
                TString name( prefix );
                name += allVars;
                name += "Like";
                this->addSPlotNtupleDoubleBranch( name );
            } else {
                std::cerr << "WARNING in LauSimpleFitModel::addSPlotNtupleBranches : Can't yet deal with 3D PDFs."
                          << std::endl;
            }
        }
    }
}

Double_t LauSimpleFitModel::setSPlotNtupleBranchValues( LauPdfPList* extraPdfs,
                                                        const TString& prefix,
                                                        UInt_t iEvt )
{
    // Store the PDF value for each variable in the list
    Double_t totalLike( 1.0 );
    Double_t extraLike( 0.0 );
    if ( extraPdfs ) {
        for ( LauPdfPList::iterator pdf_iter = extraPdfs->begin(); pdf_iter != extraPdfs->end();
              ++pdf_iter ) {

            // calculate the likelihood for this event
            ( *pdf_iter )->calcLikelihoodInfo( iEvt );
            extraLike  = ( *pdf_iter )->getLikelihood();
            totalLike *= extraLike;

            // Count the number of input variables that are not
            // DP variables (used in the case where there is DP
            // dependence for e.g. MVA)
            UInt_t nVars( 0 );
            std::vector<TString> varNames = ( *pdf_iter )->varNames();
            for ( std::vector<TString>::const_iterator var_iter = varNames.begin();
                  var_iter != varNames.end();
                  ++var_iter ) {
                if ( ( *var_iter ) != "m13Sq" && ( *var_iter ) != "m23Sq" ) {
                    ++nVars;
                }
            }

            if ( nVars == 1 ) {
                // If the PDF only has one variable then
                // simply store the value for that variable
                TString varName = ( *pdf_iter )->varName();
                TString name( prefix );
                name += varName;
                name += "Like";
                this->setSPlotNtupleDoubleBranchValue( name, extraLike );
            } else if ( nVars == 2 ) {
                // If the PDF has two variables then we
                // store the value for them both together
                // and for each on their own
                TString allVars( "" );
                for ( std::vector<TString>::const_iterator var_iter = varNames.begin();
                      var_iter != varNames.end();
                      ++var_iter ) {
                    allVars += ( *var_iter );
                    TString name( prefix );
                    name               += ( *var_iter );
                    name               += "Like";
                    Double_t indivLike  = ( *pdf_iter )->getLikelihood( ( *var_iter ) );
                    this->setSPlotNtupleDoubleBranchValue( name, indivLike );
                }
                TString name( prefix );
                name += allVars;
                name += "Like";
                this->setSPlotNtupleDoubleBranchValue( name, extraLike );
            } else {
                std::cerr << "WARNING in LauSimpleFitModel::setSPlotNtupleBranchValues : Can't yet deal with 3D PDFs."
                          << std::endl;
            }
        }
    }
    return totalLike;
}

LauSPlot::NameSet LauSimpleFitModel::variableNames() const
{
    LauSPlot::NameSet nameSet;
    if ( this->useDP() ) {
        nameSet.insert( "DP" );
    }
    // Loop through all the signal PDFs
    for ( LauPdfPList::const_iterator pdf_iter = signalPdfs_.begin(); pdf_iter != signalPdfs_.end();
          ++pdf_iter ) {
        // Loop over the variables involved in each PDF
        std::vector<TString> varNames = ( *pdf_iter )->varNames();
        for ( std::vector<TString>::const_iterator var_iter = varNames.begin();
              var_iter != varNames.end();
              ++var_iter ) {
            // If they are not DP coordinates then add them
            if ( ( *var_iter ) != "m13Sq" && ( *var_iter ) != "m23Sq" ) {
                nameSet.insert( ( *var_iter ) );
            }
        }
    }
    return nameSet;
}

LauSPlot::NumbMap LauSimpleFitModel::freeSpeciesNames() const
{
    LauSPlot::NumbMap numbMap;
    if ( ! signalEvents_->fixed() && this->doEMLFit() ) {
        numbMap["sig"] = signalEvents_->genValue();
    }
    if ( usingBkgnd_ ) {
        const UInt_t nBkgnds = this->nBkgndClasses();
        for ( UInt_t iBkgnd( 0 ); iBkgnd < nBkgnds; ++iBkgnd ) {
            const TString& bkgndClass = this->bkgndClassName( iBkgnd );
            const LauAbsRValue* par   = bkgndEvents_[iBkgnd];
            if ( ! par->fixed() ) {
                numbMap[bkgndClass] = par->genValue();
                if ( ! par->isLValue() ) {
                    std::cerr << "WARNING in LauSimpleFitModel::freeSpeciesNames : \"" << par->name()
                              << "\" is a LauFormulaPar, which implies it is perhaps not entirely free to float in the fit, so the sWeight calculation may not be reliable"
                              << std::endl;
                }
            }
        }
    }
    return numbMap;
}

LauSPlot::NumbMap LauSimpleFitModel::fixdSpeciesNames() const
{
    LauSPlot::NumbMap numbMap;
    if ( signalEvents_->fixed() && this->doEMLFit() ) {
        numbMap["sig"] = signalEvents_->genValue();
    }
    if ( usingBkgnd_ ) {
        const UInt_t nBkgnds = this->nBkgndClasses();
        for ( UInt_t iBkgnd( 0 ); iBkgnd < nBkgnds; ++iBkgnd ) {
            const TString& bkgndClass = this->bkgndClassName( iBkgnd );
            const LauAbsRValue* par   = bkgndEvents_[iBkgnd];
            if ( par->fixed() ) {
                numbMap[bkgndClass] = par->genValue();
            }
        }
    }
    return numbMap;
}

LauSPlot::TwoDMap LauSimpleFitModel::twodimPDFs() const
{
    LauSPlot::TwoDMap twodimMap;

    for ( LauPdfPList::const_iterator pdf_iter = signalPdfs_.begin(); pdf_iter != signalPdfs_.end();
          ++pdf_iter ) {
        // Count the number of input variables that are not DP variables
        UInt_t nVars( 0 );
        std::vector<TString> varNames = ( *pdf_iter )->varNames();
        for ( std::vector<TString>::const_iterator var_iter = varNames.begin();
              var_iter != varNames.end();
              ++var_iter ) {
            if ( ( *var_iter ) != "m13Sq" && ( *var_iter ) != "m23Sq" ) {
                ++nVars;
            }
        }
        if ( nVars == 2 ) {
            if ( useSCF_ ) {
                twodimMap.insert(
                    std::make_pair( "sigTM", std::make_pair( varNames[0], varNames[1] ) ) );
            } else {
                twodimMap.insert(
                    std::make_pair( "sig", std::make_pair( varNames[0], varNames[1] ) ) );
            }
        }
    }

    if ( useSCF_ ) {
        for ( LauPdfPList::const_iterator pdf_iter = scfPdfs_.begin(); pdf_iter != scfPdfs_.end();
              ++pdf_iter ) {
            // Count the number of input variables that are not DP variables
            UInt_t nVars( 0 );
            std::vector<TString> varNames = ( *pdf_iter )->varNames();
            for ( std::vector<TString>::const_iterator var_iter = varNames.begin();
                  var_iter != varNames.end();
                  ++var_iter ) {
                if ( ( *var_iter ) != "m13Sq" && ( *var_iter ) != "m23Sq" ) {
                    ++nVars;
                }
            }
            if ( nVars == 2 ) {
                twodimMap.insert(
                    std::make_pair( "sigSCF", std::make_pair( varNames[0], varNames[1] ) ) );
            }
        }
    }

    if ( usingBkgnd_ ) {
        const UInt_t nBkgnds = this->nBkgndClasses();
        for ( UInt_t iBkgnd( 0 ); iBkgnd < nBkgnds; ++iBkgnd ) {
            const TString& bkgndClass  = this->bkgndClassName( iBkgnd );
            const LauPdfPList& pdfList = bkgndPdfs_[iBkgnd];
            for ( LauPdfPList::const_iterator pdf_iter = pdfList.begin(); pdf_iter != pdfList.end();
                  ++pdf_iter ) {
                // Count the number of input variables that are not DP variables
                UInt_t nVars( 0 );
                std::vector<TString> varNames = ( *pdf_iter )->varNames();
                for ( std::vector<TString>::const_iterator var_iter = varNames.begin();
                      var_iter != varNames.end();
                      ++var_iter ) {
                    if ( ( *var_iter ) != "m13Sq" && ( *var_iter ) != "m23Sq" ) {
                        ++nVars;
                    }
                }
                if ( nVars == 2 ) {
                    twodimMap.insert(
                        std::make_pair( bkgndClass, std::make_pair( varNames[0], varNames[1] ) ) );
                }
            }
        }
    }

    return twodimMap;
}

void LauSimpleFitModel::storePerEvtLlhds()
{
    std::cout << "INFO in LauSimpleFitModel::storePerEvtLlhds : Storing per-event likelihood values..."
              << std::endl;

    // if we've not been using the DP model then we need to cache all
    // the info here so that we can get the efficiency from it

    LauFitDataTree* inputFitData = this->fitData();

    if ( ! this->useDP() && this->storeDPEff() ) {
        sigDPModel_->initialise( coeffs_ );
        sigDPModel_->fillDataTree( *inputFitData );
    }

    UInt_t evtsPerExpt( this->eventsPerExpt() );

    for ( UInt_t iEvt = 0; iEvt < evtsPerExpt; ++iEvt ) {

        this->setSPlotNtupleIntegerBranchValue( "iExpt", this->iExpt() );
        this->setSPlotNtupleIntegerBranchValue( "iEvtWithinExpt", iEvt );

        // the DP information
        this->getEvtDPLikelihood( iEvt );
        if ( this->storeDPEff() ) {
            if ( ! this->useDP() ) {
                sigDPModel_->calcLikelihoodInfo( iEvt );
            }
            this->setSPlotNtupleDoubleBranchValue( "efficiency", sigDPModel_->getEvtEff() );
            if ( sigDPModel_->usingScfModel() ) {
                this->setSPlotNtupleDoubleBranchValue( "scffraction",
                                                       sigDPModel_->getEvtScfFraction() );
            }
        }
        if ( this->useDP() ) {
            sigTotalLike_ = sigDPLike_;
            if ( useSCF_ ) {
                scfTotalLike_ = scfDPLike_;
                this->setSPlotNtupleDoubleBranchValue( "sigTMDPLike", sigDPLike_ );
                this->setSPlotNtupleDoubleBranchValue( "sigSCFDPLike", scfDPLike_ );
            } else {
                this->setSPlotNtupleDoubleBranchValue( "sigDPLike", sigDPLike_ );
            }
            if ( usingBkgnd_ ) {
                const UInt_t nBkgnds = this->nBkgndClasses();
                for ( UInt_t iBkgnd( 0 ); iBkgnd < nBkgnds; ++iBkgnd ) {
                    TString name  = this->bkgndClassName( iBkgnd );
                    name         += "DPLike";
                    this->setSPlotNtupleDoubleBranchValue( name, bkgndDPLike_[iBkgnd] );
                }
            }
        } else {
            sigTotalLike_ = 1.0;
            if ( useSCF_ ) {
                scfTotalLike_ = 1.0;
            }
            if ( usingBkgnd_ ) {
                const UInt_t nBkgnds = this->nBkgndClasses();
                for ( UInt_t iBkgnd( 0 ); iBkgnd < nBkgnds; ++iBkgnd ) {
                    bkgndTotalLike_[iBkgnd] = 1.0;
                }
            }
        }

        // the signal PDF values
        if ( useSCF_ ) {
            sigTotalLike_ *= this->setSPlotNtupleBranchValues( &signalPdfs_, "sigTM", iEvt );
            scfTotalLike_ *= this->setSPlotNtupleBranchValues( &scfPdfs_, "sigSCF", iEvt );
        } else {
            sigTotalLike_ *= this->setSPlotNtupleBranchValues( &signalPdfs_, "sig", iEvt );
        }

        // the background PDF values
        if ( usingBkgnd_ ) {
            const UInt_t nBkgnds = this->nBkgndClasses();
            for ( UInt_t iBkgnd( 0 ); iBkgnd < nBkgnds; ++iBkgnd ) {
                const TString& bkgndClass = this->bkgndClassName( iBkgnd );
                LauPdfPList& pdfs         = bkgndPdfs_[iBkgnd];
                bkgndTotalLike_[iBkgnd] *=
                    this->setSPlotNtupleBranchValues( &( pdfs ), bkgndClass, iEvt );
            }
        }

        // the total likelihoods
        if ( useSCF_ ) {
            Double_t scfFrac( 0.0 );
            if ( useSCFHist_ ) {
                scfFrac = recoSCFFracs_[iEvt];
            } else {
                scfFrac = scfFrac_.unblindValue();
            }
            this->setSPlotNtupleDoubleBranchValue( "sigSCFFrac", scfFrac );
            sigTotalLike_ *= ( 1.0 - scfFrac );
            if ( scfMap_ == 0 ) {
                scfTotalLike_ *= scfFrac;
            }
            this->setSPlotNtupleDoubleBranchValue( "sigTMTotalLike", sigTotalLike_ );
            this->setSPlotNtupleDoubleBranchValue( "sigSCFTotalLike", scfTotalLike_ );
        } else {
            this->setSPlotNtupleDoubleBranchValue( "sigTotalLike", sigTotalLike_ );
        }
        if ( usingBkgnd_ ) {
            const UInt_t nBkgnds = this->nBkgndClasses();
            for ( UInt_t iBkgnd( 0 ); iBkgnd < nBkgnds; ++iBkgnd ) {
                TString name  = this->bkgndClassName( iBkgnd );
                name         += "TotalLike";
                this->setSPlotNtupleDoubleBranchValue( name, bkgndTotalLike_[iBkgnd] );
            }
        }
        // fill the tree
        this->fillSPlotNtupleBranches();
    }
    std::cout << "INFO in LauSimpleFitModel::storePerEvtLlhds : Finished storing per-event likelihood values."
              << std::endl;
}

std::map<TString, LauComplex> LauSimpleFitModel::getDPAmps( const Double_t m13Sq,
                                                            const Double_t m23Sq )
{
    // Initialise the DP model, if not already done
    if ( coeffs_.empty() ) {
        this->updateCoeffs();
        this->initialiseDPModels();
    }

    if ( ! kinematics_->withinDPLimits( m13Sq, m23Sq ) ) {
        return {
            {"signal", {}}
        };
    }

    sigDPModel_->calcLikelihoodInfo( m13Sq, m23Sq );

    return {
        {"signal", sigDPModel_->getEvtDPAmp()}
    };
}

std::map<TString, Double_t> LauSimpleFitModel::getDPLikelihoods( const Double_t m13Sq,
                                                                 const Double_t m23Sq )
{
    // Initialise the DP model, if not already done
    if ( coeffs_.empty() ) {
        this->updateCoeffs();
        this->initialiseDPModels();
    }

    std::map<TString, Double_t> likelihoods;

    if ( ! kinematics_->withinDPLimits( m13Sq, m23Sq ) ) {
        likelihoods.emplace( "signal", 0.0 );
        if ( usingBkgnd_ ) {
            const UInt_t nBkgnds { this->nBkgndClasses() };
            for ( UInt_t bkgndID( 0 ); bkgndID < nBkgnds; ++bkgndID ) {
                likelihoods.emplace( this->bkgndClassName( bkgndID ), 0.0 );
            }
        }
        return likelihoods;
    }

    sigDPModel_->calcLikelihoodInfo( m13Sq, m23Sq );
    likelihoods.emplace( "signal", sigDPModel_->getEvtLikelihood() );

    // TODO - SCF signal
    static bool warningIssued { false };
    if ( useSCF_ && ! warningIssued ) {
        warningIssued = true;
        std::cerr << "WARNING in LauSimpleFitModel::getDPLikelihoods : calculation of SCF likelihoods not currently implemented in this function\n";
        std::cerr << "                                               : signal likelihood will just be the truth-matched value";
        std::cerr << std::endl;
    }

    if ( usingBkgnd_ ) {
        const UInt_t nBkgnds { this->nBkgndClasses() };
        for ( UInt_t bkgndID( 0 ); bkgndID < nBkgnds; ++bkgndID ) {
            likelihoods.emplace( this->bkgndClassName( bkgndID ),
                                 bkgndDPModels_[bkgndID]->getLikelihood( m13Sq, m23Sq ) );
        }
    }

    return likelihoods;
}

void LauSimpleFitModel::embedSignal( const TString& fileName,
                                     const TString& treeName,
                                     Bool_t reuseEventsWithinEnsemble,
                                     Bool_t reuseEventsWithinExperiment,
                                     Bool_t useReweighting )
{
    if ( signalTree_ ) {
        std::cerr << "ERROR in LauSimpleFitModel::embedSignal : Already embedding signal from a file."
                  << std::endl;
        return;
    }

    signalTree_ = new LauEmbeddedData( fileName, treeName, reuseEventsWithinExperiment );

    Bool_t dataOK = signalTree_->findBranches();
    if ( ! dataOK ) {
        delete signalTree_;
        signalTree_ = 0;
        std::cerr << "ERROR in LauSimpleFitModel::embedSignal : Problem creating data tree for embedding."
                  << std::endl;
        return;
    }

    reuseSignal_    = reuseEventsWithinEnsemble;
    useReweighting_ = useReweighting;

    if ( this->enableEmbedding() == kFALSE ) {
        this->enableEmbedding( kTRUE );
    }
}

void LauSimpleFitModel::embedBkgnd( const TString& bkgndClass,
                                    const TString& fileName,
                                    const TString& treeName,
                                    Bool_t reuseEventsWithinEnsemble,
                                    Bool_t reuseEventsWithinExperiment )
{
    if ( ! this->validBkgndClass( bkgndClass ) ) {
        std::cerr << "ERROR in LauSimpleFitModel::embedBkgnd : Invalid background class \""
                  << bkgndClass << "\"." << std::endl;
        std::cerr << "                                       : Background class names must be provided in \"setBkgndClassNames\" before any other background-related actions can be performed."
                  << std::endl;
        return;
    }

    UInt_t bkgndID = this->bkgndClassID( bkgndClass );

    if ( bkgndTree_[bkgndID] ) {
        std::cerr << "ERROR in LauSimpleFitModel::embedBkgnd : Already embedding background from a file."
                  << std::endl;
        return;
    }

    bkgndTree_[bkgndID] = new LauEmbeddedData( fileName, treeName, reuseEventsWithinExperiment );

    Bool_t dataOK = bkgndTree_[bkgndID]->findBranches();
    if ( ! dataOK ) {
        delete bkgndTree_[bkgndID];
        bkgndTree_[bkgndID] = 0;
        std::cerr << "ERROR in LauSimpleFitModel::embedBkgnd : Problem creating data tree for embedding."
                  << std::endl;
        return;
    }

    reuseBkgnd_[bkgndID] = reuseEventsWithinEnsemble;

    if ( this->enableEmbedding() == kFALSE ) {
        this->enableEmbedding( kTRUE );
    }
}

void LauSimpleFitModel::weightEvents( const TString& dataFileName, const TString& dataTreeName )
{
    // Routine to provide weights for events that are uniformly distributed
    // in the DP (or square DP) so as to reproduce the given DP model

    const Bool_t squareDP { kinematics_->squareDP() };
    if ( squareDP ) {
        std::cout << "INFO in LauSimpleFitModel::weightEvents : will store DP model weights and the square DP jacobian\n";
        std::cout << "                                        : the DP model weights can be used on their own to weight events that were generated flat in phase space\n";
        std::cout << "                                        : or they can be multiplied by the jacobian to weight events that were generated flat in the square DP\n";
        std::cout << "                                        : or they can be multiplied by max(1.0, jacobian) to weight events that were generated quasi-flat in the square DP"
                  << std::endl;
    } else {
        std::cout << "INFO in LauSimpleFitModel::weightEvents : will store DP model weights suitable for weighting events that were generated flat in phase space"
                  << std::endl;
    }

    // This reads in the given dataFile and creates an input
    // fit data tree that stores them for all events and experiments.
    const Bool_t dataOK { this->verifyFitData( dataFileName, dataTreeName ) };
    if ( ! dataOK ) {
        std::cerr << "ERROR in LauSimpleFitModel::weightEvents : Problem caching the data."
                  << std::endl;
        return;
    }

    LauFitDataTree* inputFitData { this->fitData() };

    if ( ! inputFitData->haveBranch( "m13Sq_MC" ) || ! inputFitData->haveBranch( "m23Sq_MC" ) ) {
        std::cerr << "ERROR in LauSimpleFitModel::weightEvents : Cannot find MC truth DP coordinate branches in supplied data, aborting."
                  << std::endl;
        return;
    }

    // Create the ntuple to hold the DP weights
    TString weightsFileName { dataFileName };
    const Ssiz_t index { weightsFileName.Last( '.' ) };
    weightsFileName.Insert( index, "_DPweights" );

    LauGenNtuple weightsTuple { weightsFileName, dataTreeName };
    weightsTuple.addIntegerBranch( "iExpt" );
    weightsTuple.addIntegerBranch( "iEvtWithinExpt" );
    weightsTuple.addDoubleBranch( "dpModelWeight" );
    if ( squareDP ) {
        weightsTuple.addDoubleBranch( "sqDPJacobian" );
    }

    const UInt_t nExpmt { this->nExpt() };
    const UInt_t firstExpmt { this->firstExpt() };
    for ( UInt_t iExpmt { firstExpmt }; iExpmt < ( firstExpmt + nExpmt ); ++iExpmt ) {

        inputFitData->readExperimentData( iExpmt );

        const UInt_t nEvents { inputFitData->nEvents() };
        if ( nEvents < 1 ) {
            std::cerr << "WARNING in LauSimpleFitModel::weightEvents : Zero events in experiment "
                      << iExpmt << ", skipping..." << std::endl;
            continue;
        }

        weightsTuple.setIntegerBranchValue( "iExpt", iExpmt );

        // Calculate and store the weights for the events in this experiment
        for ( UInt_t iEvent { 0 }; iEvent < nEvents; ++iEvent ) {

            weightsTuple.setIntegerBranchValue( "iEvtWithinExpt", iEvent );

            const LauFitData& evtData { inputFitData->getData( iEvent ) };

            const Double_t m13Sq_MC { evtData.at( "m13Sq_MC" ) };
            const Double_t m23Sq_MC { evtData.at( "m23Sq_MC" ) };

            Double_t dpModelWeight { 0.0 };
            Double_t jacobian { 1.0 };

            if ( kinematics_->withinDPLimits( m13Sq_MC, m23Sq_MC ) ) {

                kinematics_->updateKinematics( m13Sq_MC, m23Sq_MC );
                dpModelWeight = sigDPModel_->getEventWeight();

                if ( squareDP ) {
                    jacobian = kinematics_->calcSqDPJacobian();
                }

                dpModelWeight /= sigDPModel_->getDPNorm();
            }

            weightsTuple.setDoubleBranchValue( "dpModelWeight", dpModelWeight );
            if ( squareDP ) {
                weightsTuple.setDoubleBranchValue( "sqDPJacobian", jacobian );
            }
            weightsTuple.fillBranches();
        }
    }

    weightsTuple.buildIndex( "iExpt", "iEvtWithinExpt" );
    weightsTuple.addFriendTree( dataFileName, dataTreeName );
    weightsTuple.writeOutGenResults();
}

void LauSimpleFitModel::savePDFPlots( const TString& label )
{
    savePDFPlotsWave( label, 0 );
    savePDFPlotsWave( label, 1 );
    savePDFPlotsWave( label, 2 );

    std::cout << "LauCPFitModel::plot" << std::endl;
    //    ((LauIsobarDynamics*)sigDPModel_)->plot();

    //Double_t minm13 = sigDPModel_->getKinematics()->getm13Min();
    Double_t minm13 = 0.0;
    Double_t maxm13 = sigDPModel_->getKinematics()->getm13Max();
    //Double_t minm23 = sigDPModel_->getKinematics()->getm23Min();
    Double_t minm23 = 0.0;
    Double_t maxm23 = sigDPModel_->getKinematics()->getm23Max();

    Double_t mins13 = minm13 * minm13;
    Double_t maxs13 = maxm13 * maxm13;
    Double_t mins23 = minm23 * minm23;
    Double_t maxs23 = maxm23 * maxm23;

    Double_t s13, s23, chPdf;

    Int_t n13 = 200.00, n23 = 200.00;
    Double_t delta13, delta23;
    delta13     = ( maxs13 - mins13 ) / n13;
    delta23     = ( maxs23 - mins23 ) / n23;
    UInt_t nAmp = sigDPModel_->getnCohAmp();
    for ( UInt_t resID = 0; resID <= nAmp; ++resID ) {
        TGraph2D* dt    = new TGraph2D();
        TString resName = "TotalAmp";
        if ( resID != nAmp ) {
            TString tStrResID                = Form( "%d", resID );
            const LauIsobarDynamics* model   = sigDPModel_;
            const LauAbsResonance* resonance = model->getResonance( resID );
            resName                          = resonance->getResonanceName();
            std::cout << "resName = " << resName << std::endl;
        }

        resName.ReplaceAll( "(", "" );
        resName.ReplaceAll( ")", "" );
        resName.ReplaceAll( "*", "Star" );

        TCanvas* c = new TCanvas( "c" + resName + label, resName + " (" + label + ")", 0, 0, 600, 400 );
        dt->SetName( resName + label );
        dt->SetTitle( resName + " (" + label + ")" );
        Int_t count = 0;
        for ( Int_t i = 0; i < n13; i++ ) {
            s13 = mins13 + i * delta13;
            for ( Int_t j = 0; j < n23; j++ ) {
                s23 = mins23 + j * delta23;
                if ( sigDPModel_->getKinematics()->withinDPLimits2( s23, s13 ) ) {
                    //if (s13 > 4) continue;
                    sigDPModel_->calcLikelihoodInfo( s13, s23 );
                    LauComplex chAmp = sigDPModel_->getEvtDPAmp();
                    if ( resID != nAmp ) {
                        chAmp = sigDPModel_->getFullAmplitude( resID );
                    }
                    chPdf = chAmp.abs2();
                    //if ((z > 0.04)||(z < -0.03)) continue;
                    //z = TMath::Log(z);
                    if ( sigDPModel_->getDaughters()->gotSymmetricalDP() ) {
                        Double_t sLow  = s13;
                        Double_t sHigh = s23;
                        if ( sLow > sHigh ) {
                            continue;
                            //sLow = s23;
                            //sHigh = s13;
                        }

                        //if (sLow > 3.5) continue;
                        //if (i<10) std::cout << "SymmetricalDP" << std::endl;
                        //if (z>0.02) z = 0.02;
                        //if (z<-0.02) z = -0.02;

                        dt->SetPoint( count, sHigh, sLow, chPdf );
                        count++;
                    } else {
                        dt->SetPoint( count, s13, s23, chPdf );
                        count++;
                    }
                }
            }
        }
        gStyle->SetPalette( 1 );
        dt->GetXaxis()->SetTitle( "m_{K#pi}(low)" );
        dt->GetYaxis()->SetTitle( "m_{K#pi}(high)" );
        dt->Draw( "SURF1" );
        dt->GetXaxis()->SetTitle( "m_{K#pi}(low)" );
        dt->GetYaxis()->SetTitle( "m_{K#pi}(high)" );
        c->SaveAs( "plot_2D_" + resName + "_" + label + ".C" );
    }
}

void LauSimpleFitModel::savePDFPlotsWave( const TString& label, const Int_t& spin )
{

    std::cout << "label = " << label << ", spin = " << spin << std::endl;

    TString tStrResID = "S_Wave";
    if ( spin == 1 )
        tStrResID = "P_Wave";
    if ( spin == 2 )
        tStrResID = "D_Wave";

    std::cout << "LauSimpleFitModel::savePDFPlotsWave: " << tStrResID << std::endl;

    TCanvas* c = new TCanvas( "c" + tStrResID + label, tStrResID + " (" + label + ")", 0, 0, 600, 400 );

    Double_t minm13 = 0.0;
    Double_t maxm13 = sigDPModel_->getKinematics()->getm13Max();
    Double_t minm23 = 0.0;
    Double_t maxm23 = sigDPModel_->getKinematics()->getm23Max();

    Double_t mins13 = minm13 * minm13;
    Double_t maxs13 = maxm13 * maxm13;
    Double_t mins23 = minm23 * minm23;
    Double_t maxs23 = maxm23 * maxm23;

    Double_t s13, s23, chPdf;
    TGraph2D* dt = new TGraph2D();
    dt->SetName( tStrResID + label );
    dt->SetTitle( tStrResID + " (" + label + ")" );

    Int_t n13 = 200.00, n23 = 200.00;
    Double_t delta13, delta23;
    delta13     = ( maxs13 - mins13 ) / n13;
    delta23     = ( maxs23 - mins23 ) / n23;
    UInt_t nAmp = sigDPModel_->getnCohAmp();

    Int_t count = 0;
    for ( Int_t i = 0; i < n13; i++ ) {
        s13 = mins13 + i * delta13;
        for ( Int_t j = 0; j < n23; j++ ) {
            s23 = mins23 + j * delta23;
            if ( sigDPModel_->getKinematics()->withinDPLimits2( s23, s13 ) ) {
                //if (s13 > 4) continue;
                LauComplex chAmp( 0, 0 );
                Bool_t noWaveRes = kTRUE;
                sigDPModel_->calcLikelihoodInfo( s13, s23 );
                for ( UInt_t resID = 0; resID < nAmp; ++resID ) {
                    const LauIsobarDynamics* model   = sigDPModel_;
                    const LauAbsResonance* resonance = model->getResonance( resID );
                    Int_t spin_res                   = resonance->getSpin();
                    if ( spin != spin_res )
                        continue;
                    noWaveRes  = kFALSE;
                    chAmp     += sigDPModel_->getFullAmplitude( resID );
                }

                if ( noWaveRes )
                    return;

                chPdf = chAmp.abs2();
                if ( sigDPModel_->getDaughters()->gotSymmetricalDP() ) {
                    Double_t sLow  = s13;
                    Double_t sHigh = s23;
                    if ( sLow > sHigh ) {
                        continue;
                        //sLow = s23;
                        //sHigh = s13;
                    }
                    dt->SetPoint( count, sHigh, sLow, chPdf );
                    count++;
                } else {
                    dt->SetPoint( count, s13, s23, chPdf );
                    count++;
                }
            }
        }
    }
    gStyle->SetPalette( 1 );
    dt->GetXaxis()->SetTitle( "pipi" );
    dt->GetYaxis()->SetTitle( "pipi" );
    dt->Draw( "SURF1" );
    dt->GetXaxis()->SetTitle( "pipi" );
    dt->GetYaxis()->SetTitle( "pipi" );
    c->SaveAs( "plot_2D_" + tStrResID + "_" + label + ".C" );
}
