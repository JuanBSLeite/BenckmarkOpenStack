
/*
Copyright 2008 University of Warwick

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

/*! \file LauKMatrixPropFactory.cc
    \brief File containing implementation of LauKMatrixPropFactory class.
*/

// Class for storing K-matrix propagator objects
// using the factory method.

#include "LauKMatrixPropFactory.hh"

#include "LauKMatrixPropagator.hh"

#include <iostream>
using std::cout;
using std::endl;

// the singleton instance
LauKMatrixPropFactory* LauKMatrixPropFactory::theFactory_ = 0;

LauKMatrixPropFactory::LauKMatrixPropFactory()
{
    // Constructor
    map_.clear();
}

LauKMatrixPropFactory::~LauKMatrixPropFactory()
{
    // Destructor
    KMatrixPropMap::iterator iter;
    for ( iter = map_.begin(); iter != map_.end(); ++iter ) {
        LauKMatrixPropagator* thePropagator = iter->second;
        delete thePropagator;
    }
    map_.clear();
}

LauKMatrixPropFactory* LauKMatrixPropFactory::getInstance()
{
    if ( theFactory_ == 0 ) {
        theFactory_ = new LauKMatrixPropFactory();
    }
    return theFactory_;
}

LauKMatrixPropagator* LauKMatrixPropFactory::getPropagator( const TString& name,
                                                            const TString& paramFileName,
                                                            Int_t resPairAmpInt,
                                                            Int_t nChannels,
                                                            Int_t nPoles,
                                                            Int_t rowIndex )
{
    LauKMatrixPropagator* thePropagator( 0 );

    KMatrixPropMap::iterator iter = map_.find( name );

    if ( iter != map_.end() ) {
        // We have already made this propagator
        thePropagator = iter->second;
    } else {
        // The propagator does not exist. Create it and store it in the map.
        thePropagator =
            new LauKMatrixPropagator( name, paramFileName, resPairAmpInt, nChannels, nPoles, rowIndex );
        map_[name] = thePropagator;
    }

    return thePropagator;
}
