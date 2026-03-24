
/*
Copyright 2022 University of Warwick

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

/*! \file LauASqMaxFinder.hh
    \brief File containing declaration of LauASqMaxFinder class.
*/

#ifndef LAU_ASQ_MAX_FINDER
#define LAU_ASQ_MAX_FINDER

#include "LauFitObject.hh"

#include <memory>
#include <vector>

class LauIsobarDynamics;
class LauParameter;

/*! \class LauASqMaxFinder
    \brief Class for locating the |A|^2 max for a given model

    Currently supports finding the |A|^2 max for a LauIsobarDynamics object
*/

class LauASqMaxFinder : public LauFitObject {

  public:
    //! Constructor
    /*!
        \param [in] iso the isobar dynamics for which we need to find the maximum |A|^2 value
    */
    explicit LauASqMaxFinder( LauIsobarDynamics& iso );

    //! Run the minimisation to locate the maximum |A|^2 value
    /*!
        \return the maximum |A|^2 value
    */
    Double_t find();

    //! This function sets the parameter values from Minuit
    /*!
        This function has to be public since it is called from the global FCN.
        It should not be called otherwise!

        \param [in] par an array storing the various parameter values
        \param [in] npar the number of free parameters
    */
    void setParsFromMinuit( Double_t* par, Int_t npar ) override;

    //! Calculate the new value of the negative log likelihood
    /*!
        This function has to be public since it is called from the global FCN.
        It should not be called otherwise!
    */
    Double_t getTotNegLogLikelihood() override;

  private:
    //! The isobar dynamics for which we need to find the maximum |A|^2 value
    LauIsobarDynamics& iso_;

    //! Flag to control printing of warnings if the fit wanders outside the DP
    Bool_t printMinimisationWarnings_ { kTRUE };

    //! The fit parameters (owned by us)
    std::vector<std::unique_ptr<LauParameter>> ownedParams_;

    //! The fit parameters in a form to be passed to the minimiser
    std::vector<LauParameter*> params_;

    ClassDefOverride( LauASqMaxFinder, 0 )
};

#endif
