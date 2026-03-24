
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

/*! \file LauFitNtuple.hh
    \brief File containing declaration of LauFitNtuple class.
*/

/*! \class LauFitNtuple
    \brief Class to store the results from the fit into an ntuple

    Class to store the results from the fit into an ntuple.
    The fitted values and their errors are stored.  In the case of toy MC, the true values are also stored, along with the pulls.
    Derived quantities can also be stored, such as the fit fractions.  What precisely is stored depends on the fit model.
*/

//****************************************************************************
// Class to store the results from the fit/toy MC into an ntuple
// -- CLASS DESCRIPTION [MISC] --
/// Class to store the results from the fit/toy MC into an ntuple

#ifndef LAU_FIT_NTUPLE
#define LAU_FIT_NTUPLE

#include "LauAbsFitter.hh"

#include "TMatrixDfwd.h"
#include "TString.h"

#include <set>
#include <vector>

class TFile;
class TTree;

class LauParameter;

class LauFitNtuple {

  public:
    //! Constructor
    /*!
        \param [in] fileName the name for the ntuple
        \param [in] storeAsymErrors whether or not to store the asymmetric error variables
    */
    LauFitNtuple( const TString& fileName, Bool_t storeAsymErrors );

    //! Destructor
    virtual ~LauFitNtuple();

    //! Store the correlation matrix and other fit information
    /*!
        \param [in] iExpt the experiment number
        \param [in] fitStatus the status of the fit
        \param [in] covMatrix the fit covariance matrix
    */
    void storeCorrMatrix( const UInt_t iExpt,
                          const LauAbsFitter::FitStatus& fitStatus,
                          const TMatrixD& covMatrix );

    //! Store parameters and their errors
    /*!
        \param [in] fitVars the fit parameters
        \param [in] constrainedVars the names of fit parameters that are included in multi-dimensional constraints
        \param [in] extraVars variables in addition to the fit parameters, e.g. derived quantities such as fit fractions
    */
    void storeParsAndErrors( const std::vector<LauParameter*>& fitVars,
                             const std::set<TString>& constrainedVars,
                             const std::vector<LauParameter>& extraVars );

    //! Update the fit ntuple
    void updateFitNtuple();

    //! Write out fit results
    void writeOutFitResults();

  private:
    //! Copy constructor (not implemented)
    LauFitNtuple( const LauFitNtuple& rhs );

    //! Copy assignment operator (not implemented)
    LauFitNtuple& operator=( const LauFitNtuple& rhs );

    //! Name of root file
    TString rootFileName_;
    //! Root file
    TFile* rootFile_;
    //! Fit results
    TTree* fitResults_;

    //! Fit variables
    std::vector<LauParameter*> fitVars_;
    //! Fit variables that are included in a multi-dimensional constraint
    std::set<TString> constrainedVars_;
    //! Extra variables
    std::vector<LauParameter> extraVars_;

    //! Correlation matrix
    std::vector<std::vector<Double_t>> corrMatrix_;

    //! Flags whether the fit tree has been defined
    Bool_t definedFitTree_;
    //! Flags whether or not to store the asymmetric error information
    Bool_t storeAsymErrors_;

    //! Status of fit
    LauAbsFitter::FitStatus fitStatus_;
    //! Number of fit parameters
    UInt_t nFitPars_;
    //! Number of free parameters
    UInt_t nFreePars_;
    //! Number of extra parameters
    UInt_t nExtraPars_;
    //! Experiment number
    Int_t iExpt_;

    ClassDef( LauFitNtuple, 0 )    // Fit/toyMC results ntuple
};

#endif
