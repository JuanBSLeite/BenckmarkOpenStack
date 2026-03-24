
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

/*! \file LauKMatrixPropFactory.hh
    \brief File containing declaration of LauKMatrixPropFactory class.
*/

/*! \class LauKMatrixPropFactory
    \brief Factory class for the K-matrix propagators.

    Factory class for creating and storing K-matrix propagator objects.
*/

#ifndef LAU_KMATRIX_PROP_FACTORY
#define LAU_KMATRIX_PROP_FACTORY

class LauKMatrixPropagator;

#include "TString.h"

#include <map>

class LauKMatrixPropFactory {

  public:
    //! Destructor
    virtual ~LauKMatrixPropFactory();

    //! Get a static instance of this factory class. Only one is created per application.
    static LauKMatrixPropFactory* getInstance();

    //! Retrieve the propagator if it already exists, otherwise create one.
    /*!
        \param [in] name name of the propagator
        \param [in] paramFileName the parameter file
        \param [in] resPairAmpInt the number of the daughter not produced by the resonance
        \param [in] nChannels the number of channels
        \param [in] nPoles the number of poles
        \param [in] rowIndex this specifies which row of the propagator should be used when summing over the amplitude channels

        \return the propagator
    */
    LauKMatrixPropagator* getPropagator( const TString& name,
                                         const TString& paramFileName,
                                         Int_t resPairAmpInt,
                                         Int_t nChannels,
                                         Int_t nPoles,
                                         Int_t rowIndex );

  protected:
    //! A typedef to define a map used to keep track of which propagators have been created.
    typedef std::map<TString, LauKMatrixPropagator*> KMatrixPropMap;

  private:
    //! Private constructor (to ensure the singleton nature of this class)
    LauKMatrixPropFactory();

    //! Copy constructor (not implemented)
    LauKMatrixPropFactory( const LauKMatrixPropFactory& rhs );

    //! Copy assignment operator (not implemented)
    LauKMatrixPropFactory& operator=( const LauKMatrixPropFactory& rhs );

    //! The singleton instance
    static LauKMatrixPropFactory* theFactory_;

    //! The map used to store the propagator information
    KMatrixPropMap map_;

    ClassDef( LauKMatrixPropFactory, 0 )
};

#endif
