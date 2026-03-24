
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

/*! \file LauParametricStepFuncPdf.hh
    \brief File containing declaration of LauParametricStepFuncPdf class.
*/

/*! \class LauParametricStepFuncPdf
    \brief Class for defining a Parametric Step Function PDF.

    Class that allows the user to define a Parametric Step Function PDF.
    Code mostly copied from RooFit/RooParametricStepFunction.
*/

/*****************************************************************************
 * Class based on RooFit/RooParametricStepFunction.                          *
 * Original copyright given below.                                           *
 *****************************************************************************
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef LAU_PAR_STEP_FUNC_PDF
#define LAU_PAR_STEP_FUNC_PDF

#include "LauAbsPdf.hh"

#include "TString.h"

#include <vector>

class LauParameter;

class LauParametricStepFuncPdf : public LauAbsPdf {

  public:
    //! Define the allowed options for the normalisation bin
    enum NormBin {
        First, /*!< The height of the first bin is determined from the rest to satisfy the normalisation */
        Last /*!< The height of the last bin is determined from the rest to satisfy the normalisation */
    };

    //! Constructor
    /*!
        \param [in] theVarName the name of the abscissa variable
        \param [in] params the bin contents of all the but the normalisation bin, so (N_bins - 1 entries)
        \param [in] limits the limits of the bins (lower limit of all bins and upper limit of the last), so (N_bins + 1) entries
        \param [in] normalisationBin the normalisation bin
    */
    LauParametricStepFuncPdf( const TString& theVarName,
                              const std::vector<LauAbsRValue*>& params,
                              const std::vector<Double_t>& limits,
                              NormBin normalisationBin = Last );

    //! Destructor
    virtual ~LauParametricStepFuncPdf();

    //! Calculate the likelihood (and intermediate info) for a given abscissa
    /*!
        \param [in] abscissas the values of the abscissa(s)
    */
    virtual void calcLikelihoodInfo( const LauAbscissas& abscissas );

    using LauAbsPdf::calcLikelihoodInfo;

    //! Calculate the normalisation
    virtual void calcNorm();

    //! Calculate the PDF height
    /*!
        \param [in] kinematics the current DP kinematics
    */
    virtual void calcPDFHeight( const LauKinematics* kinematics );

  protected:
    //! Number of bins
    /*!
        \return the number of bins
    */
    UInt_t nBins() const { return limits_.size() - 1; }

    //! Normalisation bin
    /*!
        \return the normalisation bin
    */
    NormBin normBin() const { return normBin_; }

  private:
    //! Copy constructor (not implemented)
    LauParametricStepFuncPdf( const LauParametricStepFuncPdf& rhs );

    //! Copy assignment operator (not implemented)
    LauParametricStepFuncPdf& operator=( const LauParametricStepFuncPdf& rhs );

    //! normalisation bin
    NormBin normBin_;
    //! limits of the bins
    std::vector<Double_t> limits_;

    ClassDef( LauParametricStepFuncPdf, 0 )    // Define the PSF PDF
};

#endif
