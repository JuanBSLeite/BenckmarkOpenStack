
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

/*! \file LauAbsBkgndDPModel.cc
    \brief File containing implementation of LauAbsBkgndDPModel class.
*/

#include "LauAbsBkgndDPModel.hh"

#include "LauDaughters.hh"

LauAbsBkgndDPModel::LauAbsBkgndDPModel( LauDaughters* daughters, LauVetoes* vetoes ) :
    daughters_( daughters ),
    kinematics_( daughters ? daughters->getKinematics() : 0 ),
    vetoes_( vetoes )
{
}
