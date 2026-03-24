
/*
Copyright 2006 University of Warwick

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

/*! \file LauParamFixed.hh
    \brief File containing declaration of LauParamFixed class.
*/

/*! \class LauParamFixed
    \brief Predicate to allow counting of the number of fixed parameters
*/

#ifndef LAU_PARAM_FIXED
#define LAU_PARAM_FIXED

#include "LauAbsRValue.hh"

class LauParamFixed {
  public:
    //! Predicate that check whether the parameter is fixed
    /*!
        \param [in] par the LauParameter that holds the information about the parameter
        \return whether or not the parameter is fixed
    */
    inline bool operator()( const LauAbsRValue* par ) const { return par->fixed(); }
};

#endif
