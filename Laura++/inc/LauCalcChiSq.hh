
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

#include "TH2Poly.h"
#include "TString.h"

#include <vector>

/*! \file LauCalcChiSq.hh
    \brief File containing declaration of LauCalcChiSq class.
*/

/*! \class LauCalcChiSq
    \brief Utility class to allow the calculation of the chisq of the fit to the Dalitz plot

    A utility class to allow the calculation of the chisq of the fit to the Dalitz plot.

    A text config file is provided that gives the datasets for the data and toy
    MC generated from the fit results. These can be in the traditional DP or
    the square DP.

    A sample config file is provided. The fields are:

    \<name of file containing data histogram\> \<name of data tree\> \<name of the x axis variable in data\> \<name of the y axis variable in data\><br>
    \<name of file containing toy MC histogram\> \<name of toy MC tree\> \<name of the x axis variable in toy MC\> \<name of the y axis variable in toy MC\><br>
    \<the minimum content for any bin\> \<the number of free parameter in the fit\> \<the scalefactor by which the toy MC bin content should be multiplied\> \<minimum x value\> \<maximum x value\> \<minimum y value\> \<maximum y value\>
*/

class LauCalcChiSq {

  public:
    //! Constructor
    /*!
        \param [in] inputFileName name of the config file
    */
    LauCalcChiSq( const TString& inputFileName = "chiSqInput.txt" );

    //! Destructor
    virtual ~LauCalcChiSq();

    //! Toggle verbose printout
    /*!
        \param [in] flag true to enable verbose printout, false to disable
    */
    inline void setVerbose( const Bool_t flag ) { verbose_ = flag; }

    //! Run the calculations
    void run();

  private:
    //! Read the config file, read the data and create histograms
    void initialiseHistos();

    //! Choose the binning scheme
    /*!
        \param [in] xs x coordinates of low statistics sample
        \param [in] ys y coordinates of low statistics sample
        \param [in] nEntries number of entries in low statistics sample
        \param [out] divisions resulting binning scheme
    */
    void pickBinning( const Double_t* xs,
                      const Double_t* ys,
                      const Int_t nEntries,
                      std::vector<Int_t>& divisions );

    //! Create the template histogram based on the binning scheme
    /*!
        This function is called recursively to perform subdivisions

        \param [in] xMin the minimum x coordinate of the region to be (sub)divided
        \param [in] xMax the maximum x coordinate of the region to be (sub)divided
        \param [in] yMin the minimum y coordinate of the region to be (sub)divided
        \param [in] yMax the maximum y coordinate of the region to be (sub)divided
        \param [in] xs x coordinates of low statistics sample
        \param [in] ys y coordinates of low statistics sample
        \param [in] nEntries number of entries in low statistics sample
        \param [in] divisions the binning scheme
        \param [in] iter indicates depth of the subdivisions
    */
    void getHisto( const Double_t xMin,
                   const Double_t xMax,
                   const Double_t yMin,
                   const Double_t yMax,
                   const Double_t* xs,
                   const Double_t* ys,
                   const Int_t nEntries,
                   const std::vector<Int_t>& divisions,
                   const UInt_t iter = 0 );

    //! Calculate the chisq from the data histograms
    void calculateChiSq();

    //! Create plots
    void makePlots();

    //! Name of the config file
    TString inputFileName_;

    //! Name of the low stats data file
    TString fileName1_;
    //! Name of the high stats data file
    TString fileName2_;
    //! Name of the low stats data tree
    TString treeName1_;
    //! Name of the high stats data tree
    TString treeName2_;

    //! Name of the x-coordinate branch in tree 1
    TString xName1_;
    //! Name of the x-coordinate branch in tree 2
    TString xName2_;
    //! Name of the y-coordinate branch in tree 1
    TString yName1_;
    //! Name of the y-coordinate branch in tree 2
    TString yName2_;

    //! The minimum bin content
    Float_t minContent_;

    //! Template histogram constructed from the binning scheme
    TH2Poly* theHisto_;
    //! Histogram (constructed from template) filled from tree 1
    TH2Poly* histo1_;
    //! Histogram (constructed from template) filled from tree 2
    TH2Poly* histo2_;
    //! Histogram (constructed from template) filled with pulls of tree1 vs tree2
    TH2Poly* pullHisto_;
    //! Histogram (constructed from template) filled with chisq of tree1 vs tree2
    TH2Poly* chiSqHisto_;
    //! Histogram (constructed from template) filled with signed chisq of tree1 vs tree2
    TH2Poly* chiSqSignedHisto_;

    //! Minimum x coordinate of histograms
    Float_t xMin_;
    //! Maximum x coordinate of histograms
    Float_t xMax_;
    //! Minimum y coordinate of histograms
    Float_t yMin_;
    //! Maximum y coordinate of histograms
    Float_t yMax_;

    //! Number of free parameters in fit (used for calculating the ndof)
    Int_t nParams_;

    //! Scalefactor between low and high stats data samples
    Float_t scaleFactor_;

    //! Verbose flag
    Bool_t verbose_;

    ClassDef( LauCalcChiSq, 0 )
};
