
/*
Copyright 2007 University of Warwick

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

/*! \file LauEmbeddedData.cc
    \brief File containing implementation of LauEmbeddedData class.
*/

#include <iostream>
#include <vector>
using std::cerr;
using std::endl;

#include "LauEmbeddedData.hh"
#include "LauIsobarDynamics.hh"
#include "LauKinematics.hh"
#include "LauRandom.hh"

#include "TRandom.h"

LauEmbeddedData::LauEmbeddedData( const TString& fileName,
                                  const TString& treeName,
                                  Bool_t allowReuseOfEvents ) :
    theDataTree_( new LauFitDataTree( fileName, treeName ) ),
    reuseEvents_( allowReuseOfEvents )
{
}

LauEmbeddedData::~LauEmbeddedData()
{
    delete theDataTree_;
}

Bool_t LauEmbeddedData::findBranches()
{
    if ( theDataTree_ == 0 ) {
        cerr << "ERROR in LauEmbeddedData::findBranches : Invalid pointer to the data tree object."
             << endl;
        return kFALSE;
    }

    Bool_t branchesReadOK = theDataTree_->findBranches();
    if ( ! branchesReadOK ) {
        cerr << "ERROR in LauEmbeddedData::findBranches : Problem finding branches." << endl;
        return kFALSE;
    }

    theDataTree_->readAllData();

    return kTRUE;
}

Bool_t LauEmbeddedData::getReweightedEvent( LauIsobarDynamics* dynamics )
{
    if ( ! theDataTree_ ) {
        cerr << "ERROR in LauEmbeddedData::getReweightedEvent : Invalid pointer to the data tree object."
             << endl;
        return kFALSE;
    }

    if ( ! dynamics ) {
        cerr << "ERROR in LauEmbeddedData::getReweightedEvent : Amplitude model is null." << endl;
        return kFALSE;
    }

    UInt_t numEvents = this->nEvents();
    UInt_t iEvt( 0 );

    // First select a random event within the embedded data sample
    if ( this->reuseEvents() ) {
        iEvt = LauRandom::randomFun()->Integer( numEvents );
    } else {
        if ( this->nUsedEvents() == numEvents ) {
            cerr << "ERROR in LauEmbeddedData::getReweightedEvent : Have already used all events in the tree."
                 << endl;
            return kFALSE;
        }
        Bool_t ok( kFALSE );
        while ( ! ok ) {
            iEvt = LauRandom::randomFun()->Integer( numEvents );
            ok   = usedEvents_.insert( iEvt ).second;
        }
    }

    // Retrieve the data for the selected event
    theData_ = theDataTree_->getData( iEvt );

    LauKinematics* kinematics = dynamics->getKinematics();

    if ( kinematics != 0 ) {

        // Get the true m13Sq, m23Sq variables.
        Double_t m13Sq_MC = this->getValue( "m13Sq_MC" );
        Double_t m23Sq_MC = this->getValue( "m23Sq_MC" );

        if ( kinematics->withinDPLimits( m13Sq_MC, m23Sq_MC ) ) {

            // Update the kinematics with the true m13Sq, m23Sq variables
            kinematics->updateKinematics( m13Sq_MC, m23Sq_MC );

            // Use this event or not, according to whether the amplitude model selects it.
            // The LauIsobarDynamics uses an accept/reject method based on the
            // ratio of the current to maximum ASq value. Other models may use alternative methods,
            // but they must provide a yes/no answer if the event is accepted.
            Bool_t gotReweightedEvent = dynamics->gotReweightedEvent();
            if ( gotReweightedEvent == kTRUE ) {

                // Event is accepted.
                // Update the kinematics to use the reco variables.
                Double_t m13Sq = this->getValue( "m13Sq" );
                Double_t m23Sq = this->getValue( "m23Sq" );
                kinematics->updateKinematics( m13Sq, m23Sq );

                if ( LauIsobarDynamics::GenOK != dynamics->checkToyMC( kFALSE, kFALSE ) ) {
                    return kFALSE;
                }

            } else {

                // Recursively run this function until we get an accepted event
                return this->getReweightedEvent( dynamics );
            }

        } else {

            cerr << "WARNING in LauEmbeddedData::getReweightedEvent : Skipping event " << iEvt
                 << ", which isn't within the DP boundary." << endl;
            return this->getReweightedEvent( dynamics );
        }
    }

    return kTRUE;
}

void LauEmbeddedData::getEmbeddedEvent( LauKinematics* kinematics )
{
    if ( ! theDataTree_ ) {
        cerr << "ERROR in LauEmbeddedData::getEmbeddedEvent : Invalid pointer to the data tree object."
             << endl;
        return;
    }
    UInt_t numEvents = this->nEvents();
    UInt_t iEvt( 0 );
    if ( this->reuseEvents() ) {
        iEvt = LauRandom::randomFun()->Integer( numEvents );
    } else {
        if ( this->nUsedEvents() == numEvents ) {
            cerr << "ERROR in LauEmbeddedData::getEmbeddedEvent : Have already used all events in the tree."
                 << endl;
            return;
        }
        Bool_t ok( kFALSE );
        while ( ! ok ) {
            iEvt = LauRandom::randomFun()->Integer( numEvents );
            ok   = usedEvents_.insert( iEvt ).second;
        }
    }
    theData_ = theDataTree_->getData( iEvt );

    if ( kinematics != 0 ) {
        Double_t m13Sq = this->getValue( "m13Sq" );
        Double_t m23Sq = this->getValue( "m23Sq" );
        if ( kinematics->withinDPLimits( m13Sq, m23Sq ) ) {
            kinematics->updateKinematics( m13Sq, m23Sq );
        } else {
            cerr << "WARNING in LauEmbeddedData::getEmbeddedEvent : Skipping event that isn't within the DP boundary."
                 << endl;
            this->getEmbeddedEvent( kinematics );
        }
    }
}

Double_t LauEmbeddedData::getValue( const TString& name ) const
{
    LauFitData::const_iterator iter = theData_.find( name );
    if ( iter == theData_.end() ) {
        cerr << "ERROR in LauEmbeddedData::getValue : Could not find branch \"" << name
             << "\" in embedded event." << endl;
        return 0.0;
    } else {
        return iter->second;
    }
}

LauFitData LauEmbeddedData::getValues( const std::vector<TString>& names ) const
{
    LauFitData retVal;

    for ( std::vector<TString>::const_iterator iter = names.begin(); iter != names.end(); ++iter ) {
        const TString& name                  = ( *iter );
        LauFitData::const_iterator data_iter = theData_.find( name );
        if ( data_iter == theData_.end() ) {
            cerr << "WARNING in LauEmbeddedData::getValues : Could not find branch \"" << name
                 << "\" in embedded event." << endl;
        } else {
            retVal.insert( *data_iter );
        }
    }

    return retVal;
}
