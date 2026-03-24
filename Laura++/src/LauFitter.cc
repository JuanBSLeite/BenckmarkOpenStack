
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

/*! \file LauFitter.cc
    \brief File containing implementation of LauFitter methods.
*/

#include "LauFitter.hh"

#include "LauMinuit.hh"

#include <array>
#include <iostream>
#include <memory>

std::unique_ptr<LauAbsFitter> LauFitter::theInstance_ = nullptr;
LauFitter::Type LauFitter::fitterType_                = LauFitter::Type::Minuit;
LauOutputLevel LauFitter::fitterVerbosity_            = LauOutputLevel::Standard;
UInt_t LauFitter::fitterMaxPars_                      = 100;

void LauFitter::setFitterType( const Type type )
{
    if ( theInstance_ != nullptr ) {
        std::cerr << "ERROR in LauFitter::setFitterType : The fitter has already been created, cannot change the type now."
                  << std::endl;
        return;
    }

    fitterType_ = type;
}

void LauFitter::setFitterVerbosity( const LauOutputLevel level )
{
    if ( theInstance_ != nullptr ) {
        std::cerr << "ERROR in LauFitter::setFitterVerbosity : The fitter has already been created, cannot change the verbosity now."
                  << std::endl;
        return;
    }

    fitterVerbosity_ = level;
}

void LauFitter::setFitterMaxPars( const UInt_t maxPars )
{
    if ( theInstance_ != nullptr ) {
        std::cerr << "ERROR in LauFitter::setFitterMaxPars : The fitter has already been created, cannot change the maximum number of parameters now."
                  << std::endl;
        return;
    }

    fitterMaxPars_ = maxPars;
}

LauAbsFitter& LauFitter::fitter()
{
    // Returns a reference to a singleton LauAbsFitter object.
    // Creates the object the first time it is called.

    if ( theInstance_ == nullptr ) {
        if ( fitterType_ == Type::Minuit ) {
            // NB cannot use std::make_unique here since the LauMinuit constructor is private
            theInstance_.reset( new LauMinuit( fitterMaxPars_, fitterVerbosity_ ) );
        }
    }

    return *theInstance_;
}

void LauFitter::destroyFitter()
{
    // destroy the current fitter
    theInstance_.reset();

    // restore the default settings
    fitterType_      = LauFitter::Type::Minuit;
    fitterVerbosity_ = LauOutputLevel::Standard;
    fitterMaxPars_   = 100;
}
