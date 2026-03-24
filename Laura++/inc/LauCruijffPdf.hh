
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

/*! \file LauCruijffPdf.hh
    \brief File containing declaration of LauCruijffPdf class.
*/

/*! \class LauCruijffPdf
    \brief Class for defining a Cruijff PDF.

    Class that allows the user to define a Cruijff PDF, a bifurcated Gaussian with asymmetric tails.
    The guts of the implementation have been copied from Wouter Hulsbergen's RooFit class.
*/

/*****************************************************************************
 * Class based on RooFit/RooCruijff.                                         *
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

#ifndef LAU_CRUIJFF_PDF
#define LAU_CRUIJFF_PDF

#include "LauAbsPdf.hh"
#include "LauParameter.hh"

#include "TRandom.h"
#include "TString.h"

#include <vector>

class LauCruijffPdf : public LauAbsPdf {

  public:
    //! Constructor
    /*!
        \param [in] theVarName the name of the abscissa variable
        \param [in] params the PDF parameters - mean, sigmaR, sigmaL, alphaR and alphaL
        \param [in] minAbscissa the minimum value of the abscissa
        \param [in] maxAbscissa the maximum value of the abscissa
    */
    LauCruijffPdf( const TString& theVarName,
                   const std::vector<LauAbsRValue*>& params,
                   Double_t minAbscissa,
                   Double_t maxAbscissa );

    //! Destructor
    virtual ~LauCruijffPdf();

    //! Calculate the likelihood (and intermediate info) for a given abscissa
    /*!
        \param [in] abscissas the values of the abscissa(s)
    */
    virtual void calcLikelihoodInfo( const LauAbscissas& abscissas );

    using LauAbsPdf::calcLikelihoodInfo;

    //! Calculate the PDF height
    /*!
        \param [in] kinematics the current DP kinematics
    */
    virtual void calcPDFHeight( const LauKinematics* kinematics );

  protected:

  private:
    //! Copy constructor (not implemented)
    LauCruijffPdf( const LauCruijffPdf& other );

    //! Copy assignment operator (not implemented)
    LauCruijffPdf& operator=( const LauCruijffPdf& other );

    //! Gaussian mean
    LauAbsRValue* mean_;
    //! Sigma of left Gaussian
    LauAbsRValue* sigmaL_;
    //! Sigma of right Gaussian
    LauAbsRValue* sigmaR_;
    //! Alpha of left Gaussian
    LauAbsRValue* alphaL_;
    //! Alpha of right Gaussian
    LauAbsRValue* alphaR_;

    ClassDef( LauCruijffPdf, 0 )    // Define the Cruijff PDF
};

#endif
