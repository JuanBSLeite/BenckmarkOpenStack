
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

/*! \file LauFitter.hh
    \brief File containing declaration of LauFitter class.
*/

#ifndef LAU_FITTER
#define LAU_FITTER

#include "LauPrint.hh"

#include "Rtypes.h"
#include "TString.h"

#include <memory>

class LauAbsFitter;

/*! \class LauFitter
    \brief Factory class for creating and providing access to the fitter.

    The fitter type and verbosity can be set before first access to determine which fitter is used.
*/

class LauFitter final {

  public:
    //! The types of fitter available
    enum class Type {
        Minuit /*!< the Minuit fitter */
    };

    //! Set the type of the fitter
    /*!
        \param [in] type the type of the fitter (default set to Minuit)
    */
    static void setFitterType( const Type type );

    //! Set the verbosity level of the fitter
    /*!
        \param [in] level the level of verbosity of the fitter (default set to Standard)
    */
    static void setFitterVerbosity( const LauOutputLevel level );

    //! Set the maximum number of parameters for the fitter
    /*!
        \param [in] maxPars the maximum number of parameters for the fitter (default set to 100)
    */
    static void setFitterMaxPars( const UInt_t maxPars );

    //! Method that provides access to the singleton fitter
    /*!
        \return a reference to a singleton LauAbsFitter object
    */
    static LauAbsFitter& fitter();

    //! Destroy the current fitter
    /*!
        A new fitter will be created on the next call to LauFitter::fitter
    */
    static void destroyFitter();

  private:
    //! Constructor
    LauFitter() = default;

    //! Destructor
    ~LauFitter() = default;

    //! Copy constructor (deleted)
    LauFitter( const LauFitter& ) = delete;

    //! Move constructor (deleted)
    LauFitter( LauFitter&& ) = delete;

    //! Copy assignment operator (deleted)
    LauFitter& operator=( const LauFitter& ) = delete;

    //! Move assignment operator (deleted)
    LauFitter& operator=( LauFitter&& ) = delete;

    //! Pointer to the singleton fitter instance
    static std::unique_ptr<LauAbsFitter> theInstance_;

    //! The fitter type
    static Type fitterType_;

    //! The fitter verbosity
    static LauOutputLevel fitterVerbosity_;

    //! The maximum number of parameters for the fitter
    static UInt_t fitterMaxPars_;

    ClassDef( LauFitter, 0 );
};

#endif
