
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

/*! \file LauRandom.hh
    \brief File containing LauRandom namespace.
*/

/*! \namespace LauRandom
    \brief Namespace for holding the random number generator functions.

    Holds one generator with a given (changable) seed and one where the seed is taken from the machine clock.
*/

#ifndef LAU_RANDOM
#define LAU_RANDOM

#include "Rtypes.h"

class TRandom;

namespace LauRandom {

    //! Access the singleton random number generator with a particular seed
    /*!
        The seed value can be set using LauRandom::setSeed

        \return pointer to the singleton random number generator
    */
    TRandom* randomFun();

    //! Access the singleton random number generator with seed set from machine clock time (within +-1 sec).
    /*!
        \return pointer to the singleton random number generator
    */
    TRandom* zeroSeedRandom();

    //! Set seed of random-number generator to a given value
    /*!
        \param [in] seed the seed to be used in the generation
    */
    void setSeed( UInt_t seed );

}    // namespace LauRandom

#endif
