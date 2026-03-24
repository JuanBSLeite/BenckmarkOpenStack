
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

/*! \file LauCalcChiSq.cc
    \brief File containing implementation of LauCalcChiSq class.

    Code to produce an adaptive binning scheme and calculate the 2D chi-square
    between two datasets (e.g. low-stat data and high-stat toyMC).

    Note that the low and high stat histograms must have the same bin axes ranges
    and number of bins.

    It works by using the low stat (first) histogram to find a binning scheme such
    that the total number of entries in each bin is >= Min_bin_content. The number
    of entries in the histogram is divided by the desired minimum bin content to
    give a target number of bins. The largest number of bins that can be expressed
    as a product of powers of 4, 9, 25, 49 and 121 that does not exceed the target
    value is chosen. The histogram is the recursively subdivided in 2x2, 3x3, 5x5,
    7x7 or 11x11 bins. For each stage of the subdivision, each bin is first
    divided into equally populated bins in x then each of these is further divded
    into equally popiulated bins in y.

    The (Pearson) chi-squared is then the sum of the chi-squared contributions of
    all bins:
    (low_stat_number - high_stat_number)^2/(high_stat_number)
    The nDof = number of bins - number of free params - 1
*/

#include "LauCalcChiSq.hh"

#include "TAxis.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

LauCalcChiSq::LauCalcChiSq( const TString& inputFileName ) :
    inputFileName_( inputFileName ),
    fileName1_( "" ),
    fileName2_( "" ),
    treeName1_( "" ),
    treeName2_( "" ),
    xName1_( "" ),
    xName2_( "" ),
    yName1_( "" ),
    yName2_( "" ),
    minContent_( 10.0 ),
    histo1_( 0 ),
    histo2_( 0 ),
    chiSqHisto_( 0 ),
    chiSqSignedHisto_( 0 ),
    xMin_( 0.0 ),
    xMax_( 0.0 ),
    yMin_( 0.0 ),
    yMax_( 0.0 ),
    nParams_( 0 ),
    scaleFactor_( 1.0 ),
    verbose_( kFALSE )
{
}

LauCalcChiSq::~LauCalcChiSq()
{
}

void LauCalcChiSq::run()
{
    std::cout << "Running chi-squared algorithm" << std::endl;
    this->initialiseHistos();

    std::cout << "Calculating chi-squared" << std::endl;
    this->calculateChiSq();

    //make plots
    this->makePlots();

    // Output the various histograms
    std::cout << "Writing out histogram output" << std::endl;

    TFile* outFile = new TFile( "checkHistos.root", "recreate" );

    histo1_->SetDirectory( outFile );
    histo2_->SetDirectory( outFile );
    pullHisto_->SetDirectory( outFile );
    chiSqHisto_->SetDirectory( outFile );
    chiSqSignedHisto_->SetDirectory( outFile );

    outFile->Write();
    outFile->Close();

    std::cout << "Done" << std::endl;
}

void LauCalcChiSq::calculateChiSq()
{
    Float_t toyVal( 0. ), dataVal( 0. );
    Int_t minContent( 0 ), ndof( 0 );
    Float_t diff( 0. ), chiSq( 0. ), totalChiSq( 0. );

    minContent = histo1_->GetBinContent( 1 );

    for ( Int_t i = 1; i <= histo1_->GetNumberOfBins(); ++i ) {
        //keep track of actual minimum
        if ( histo1_->GetBinContent( i ) < minContent ) {
            minContent = histo1_->GetBinContent( i );
        }

        // Calculate Pearson chi-square for this bin, using the
        // second histogram for the expected distribution
        chiSq   = 0.;
        toyVal  = histo2_->GetBinContent( i );
        dataVal = histo1_->GetBinContent( i );
        diff    = dataVal - toyVal;
        if ( toyVal > 0 )
            chiSq = ( diff * diff ) / toyVal;
        totalChiSq += chiSq;

        chiSqHisto_->SetBinContent( i, chiSq );

        if ( diff > 0 ) {
            chiSqSignedHisto_->SetBinContent( i, chiSq );
            pullHisto_->SetBinContent( i, sqrt( chiSq ) );
        } else {
            chiSqSignedHisto_->SetBinContent( i, -chiSq );
            pullHisto_->SetBinContent( i, -sqrt( chiSq ) );
        }
    }

    ndof = histo1_->GetNumberOfBins() - nParams_ - 1;

    std::cout << "Total ChiSq/nDof = " << totalChiSq << "/" << ndof << " = " << totalChiSq / ndof
              << std::endl;
    std::cout << "Actual minimum entries per bin: " << minContent << std::endl;
}

void LauCalcChiSq::makePlots()
{
    gStyle->SetPalette( 1, 0 );
    const Int_t NRGBs     = 5;
    const Int_t NCont     = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable( NRGBs, stops, red, green, blue, NCont );
    gStyle->SetNumberContours( NCont );
    gStyle->SetOptStat( 0000 );

    TCanvas can;

    can.SetTopMargin( 0.05 );
    can.SetRightMargin( 0.17 );
    can.SetBottomMargin( 0.16 );
    can.SetLeftMargin( 0.14 );

    histo1_->SetLabelFont( 62, "x" );
    histo1_->SetLabelFont( 62, "y" );
    histo1_->SetTitleFont( 62, "x" );
    histo1_->SetTitleFont( 62, "y" );
    histo1_->SetTitleSize( 0.06, "x" );
    histo1_->SetTitleSize( 0.06, "y" );
    histo1_->SetLabelSize( 0.05, "x" );
    histo1_->SetLabelSize( 0.05, "y" );
    histo1_->SetXTitle( xName1_ );
    histo1_->SetYTitle( yName1_ );
    histo1_->Draw( "colz" );
    can.SaveAs( "data.pdf" );

    histo2_->SetLabelFont( 62, "x" );
    histo2_->SetLabelFont( 62, "y" );
    histo2_->SetTitleFont( 62, "x" );
    histo2_->SetTitleFont( 62, "y" );
    histo2_->SetTitleSize( 0.06, "x" );
    histo2_->SetTitleSize( 0.06, "y" );
    histo2_->SetLabelSize( 0.05, "x" );
    histo2_->SetLabelSize( 0.05, "y" );
    histo2_->SetXTitle( xName1_ );
    histo2_->SetYTitle( yName1_ );
    histo2_->Draw( "colz" );
    can.SaveAs( "toy.pdf" );

    if ( -1. * pullHisto_->GetMinimum() > pullHisto_->GetMaximum() )
        pullHisto_->SetMaximum( -1. * pullHisto_->GetMinimum() );
    else
        pullHisto_->SetMinimum( -1. * pullHisto_->GetMaximum() );

    pullHisto_->SetLabelFont( 62, "x" );
    pullHisto_->SetLabelFont( 62, "y" );
    pullHisto_->SetTitleFont( 62, "x" );
    pullHisto_->SetTitleFont( 62, "y" );
    pullHisto_->SetTitleSize( 0.06, "x" );
    pullHisto_->SetTitleSize( 0.06, "y" );
    pullHisto_->SetLabelSize( 0.05, "x" );
    pullHisto_->SetLabelSize( 0.05, "y" );
    pullHisto_->SetXTitle( xName1_ );
    pullHisto_->SetYTitle( yName1_ );
    pullHisto_->Draw( "colz" );
    can.SaveAs( "pull.pdf" );

    chiSqHisto_->SetLabelFont( 62, "x" );
    chiSqHisto_->SetLabelFont( 62, "y" );
    chiSqHisto_->SetTitleFont( 62, "x" );
    chiSqHisto_->SetTitleFont( 62, "y" );
    chiSqHisto_->SetTitleSize( 0.06, "x" );
    chiSqHisto_->SetTitleSize( 0.06, "y" );
    chiSqHisto_->SetLabelSize( 0.05, "x" );
    chiSqHisto_->SetLabelSize( 0.05, "y" );
    chiSqHisto_->SetXTitle( xName1_ );
    chiSqHisto_->SetYTitle( yName1_ );
    chiSqHisto_->Draw( "colz" );
    can.SaveAs( "chiSq.pdf" );

    chiSqSignedHisto_->SetLabelFont( 62, "x" );
    chiSqSignedHisto_->SetLabelFont( 62, "y" );
    chiSqSignedHisto_->SetTitleFont( 62, "x" );
    chiSqSignedHisto_->SetTitleFont( 62, "y" );
    chiSqSignedHisto_->SetTitleSize( 0.06, "x" );
    chiSqSignedHisto_->SetTitleSize( 0.06, "y" );
    chiSqSignedHisto_->SetLabelSize( 0.05, "x" );
    chiSqSignedHisto_->SetLabelSize( 0.05, "y" );
    chiSqSignedHisto_->SetXTitle( xName1_ );
    chiSqSignedHisto_->SetYTitle( yName1_ );
    chiSqSignedHisto_->Draw( "colz" );
    can.SaveAs( "chiSqSigned.pdf" );
}

void LauCalcChiSq::initialiseHistos()
{

    // Open the input control file:
    // Low_stat_file_name   Low_stat_tree_name   Low_stat_x_axis_name   Low_stat_y_axis_name
    // High_stat_file_name  High_stat_tree_name  High_stat_x_axis_name  High_stat_y_axis_name
    // Min_bin_content      N_free_params        Low/high_stat_histo_ratio  xMin  xMax  yMin  yMax
    std::ifstream getData( inputFileName_.Data() );

    // get the info on the low stat histo
    getData >> fileName1_ >> treeName1_ >> xName1_ >> yName1_;
    if ( ! getData.good() ) {
        std::cerr << "Error. Could not read first line of the input file " << inputFileName_
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // open the file that contains the low stat histogram
    TFile* file1 = TFile::Open( fileName1_.Data(), "read" );

    if ( file1 == 0 ) {
        gSystem->Exit( EXIT_FAILURE );
    }

    // retrieve the low stat histogram
    TTree* tree1 = dynamic_cast<TTree*>( file1->Get( treeName1_.Data() ) );

    if ( tree1 == 0 ) {
        std::cerr << "Error. Could not find the tree " << treeName1_ << " in the file "
                  << fileName1_ << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // get the info on the high stat histogram
    getData >> fileName2_ >> treeName2_ >> xName2_ >> yName2_;
    if ( ! getData.good() ) {
        std::cerr << "Error. Could not read the second line of the input file " << inputFileName_
                  << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // if both histograms are in the same file then just retrieve the
    // high stat histogram from the file we already have open
    TFile* file2( 0 );
    TTree* tree2( 0 );
    if ( fileName2_ == fileName1_ ) {
        tree2 = dynamic_cast<TTree*>( file1->Get( treeName2_.Data() ) );
    }
    // otherwise open the other file and retrieve the high stat histogram from there
    else {
        file2 = TFile::Open( fileName2_.Data(), "read" );

        if ( file2 == 0 ) {
            gSystem->Exit( EXIT_FAILURE );
        }

        tree2 = dynamic_cast<TTree*>( file2->Get( treeName2_.Data() ) );
    }

    if ( tree2 == 0 ) {
        std::cerr << "Error. Could not find the tree " << treeName2_ << " in the file "
                  << fileName2_ << std::endl;
        gSystem->Exit( EXIT_FAILURE );
    }

    // get the info on the minimum content, number of parameters and scalefactor
    Int_t nParameters( 0 );
    Float_t minContent( 0.0 ), scaleFactor( 1.0 ), xMin( 0.0 ), xMax( 1.0 ), yMin( 0.0 ), yMax( 1.0 );
    getData >> minContent >> nParameters >> scaleFactor >> xMin >> xMax >> yMin >> yMax;
    if ( getData.good() ) {
        minContent_  = minContent;
        nParams_     = nParameters;
        scaleFactor_ = scaleFactor;
        xMin_        = xMin;
        xMax_        = xMax;
        yMin_        = yMin;
        yMax_        = yMax;
    }

    // close the text file
    getData.close();

    std::cout << "Using the files and trees: " << fileName1_ << "; " << treeName1_ << " and "
              << fileName2_ << "; " << treeName2_ << std::endl;
    std::cout << "Relative scaling factor histo1/histo2 set to " << scaleFactor_ << std::endl;
    std::cout << "Minimum bin content is " << minContent_ << std::endl;
    std::cout << "Number of free parameters is " << nParams_ << std::endl;

    Double_t x;
    Double_t y;

    tree1->SetBranchAddress( xName1_.Data(), &x );
    tree1->SetBranchAddress( yName1_.Data(), &y );

    Int_t nEntries = tree1->GetEntries();

    Double_t* xs = new Double_t[nEntries];
    Double_t* ys = new Double_t[nEntries];

    for ( Int_t i = 0; i < nEntries; ++i ) {
        tree1->GetEntry( i );
        xs[i] = x;
        ys[i] = y;
    }

    theHisto_ = new TH2Poly( "theHisto_", "", xMin_, xMax_, yMin_, yMax_ );

    //select the number of divisions to get us closest to minContent entries per bin
    std::vector<Int_t> divisions;
    this->pickBinning( xs, ys, nEntries, divisions );

    //perform the adaptive bining based on  histo1_
    this->getHisto( xMin_, xMax_, yMin_, yMax_, xs, ys, nEntries, divisions );

    histo1_           = dynamic_cast<TH2Poly*>( theHisto_->Clone( "histo1_" ) );
    histo2_           = dynamic_cast<TH2Poly*>( theHisto_->Clone( "histo2_" ) );
    pullHisto_        = dynamic_cast<TH2Poly*>( theHisto_->Clone( "pullHisto_" ) );
    chiSqHisto_       = dynamic_cast<TH2Poly*>( theHisto_->Clone( "chiSqHisto_" ) );
    chiSqSignedHisto_ = dynamic_cast<TH2Poly*>( theHisto_->Clone( "chiSqSignedHisto_" ) );
    delete[] xs;
    delete[] ys;

    //fill the two histograms from the trees
    TString drawString1, drawString2, weightString2;
    drawString1   += yName1_;
    drawString1   += ":";
    drawString1   += xName1_;
    drawString1   += ">>histo1_";
    drawString2   += yName2_;
    drawString2   += ":";
    drawString2   += xName2_;
    drawString2   += ">>histo2_";
    weightString2 += scaleFactor_;

    tree1->Draw( drawString1 );
    tree2->Draw( drawString2, weightString2 );

    histo1_->SetDirectory( 0 );
    histo2_->SetDirectory( 0 );
    pullHisto_->SetDirectory( 0 );
    chiSqHisto_->SetDirectory( 0 );
    chiSqSignedHisto_->SetDirectory( 0 );

    // close the file(s) containing the trees
    if ( file1 != 0 ) {
        file1->Close();
    }
    delete file1;
    if ( file2 != 0 ) {
        file2->Close();
    }
    delete file2;
}

void LauCalcChiSq::pickBinning( const Double_t* xs,
                                const Double_t* ys,
                                const Int_t nEntries,
                                std::vector<Int_t>& divisions )
{
    //first check how many events we have within the histogram limits
    Int_t nIn( 0 );
    for ( Int_t i = 0; i < nEntries; ++i ) {
        if ( xs[i] < xMax_ && xs[i] >= xMin_ && ys[i] < yMax_ && ys[i] >= yMin_ ) {
            ++nIn;
        }
    }

    //aim to have exactly minContent events in each bin
    Int_t nBinsTarget = nIn / minContent_;

    std::cout << "Target is " << minContent_ << " entries per bin" << std::endl;
    std::cout << "Aiming to divide " << nIn << " entries between " << nBinsTarget << " bins"
              << std::endl;

    //we will iteratively sub-divide histogram bins into either 4, 9, 25, 49 or 121
    //here we figure out how many 4s, 9s, 25s, 49s and 121s to use to best match our target without exceeding it
    Int_t nDivisions( 0 ), nNines( 0 ), nTwentyFives( 0 ), nFortyNines( 0 ), nEleventyElevens( 0 ),
        nBins( 1 );
    Int_t nDivisionsBest( 0 ), nNinesBest( 0 ), nTwentyFivesBest( 0 ), nFortyNinesBest( 0 ),
        nEleventyElevensBest( 0 ), nBinsBest( 1 );

    do {
        ++nDivisions;
        for ( nNines = 0; nNines <= nDivisions; ++nNines ) {
            for ( nTwentyFives = 0; nTwentyFives <= nDivisions - nNines; ++nTwentyFives ) {
                for ( nFortyNines = 0; nFortyNines <= nDivisions - nNines - nTwentyFives;
                      ++nFortyNines ) {
                    for ( nEleventyElevens = 0;
                          nEleventyElevens <= nDivisions - nNines - nTwentyFives - nFortyNines;
                          ++nEleventyElevens ) {
                        nBins = TMath::Power( 4,
                                              nDivisions - nNines - nTwentyFives - nFortyNines -
                                                  nEleventyElevens ) *
                                TMath::Power( 9, nNines ) * TMath::Power( 25, nTwentyFives ) *
                                TMath::Power( 49, nFortyNines ) *
                                TMath::Power( 121, nEleventyElevens );
                        if ( nBins < nBinsTarget && nBins > nBinsBest ) {
                            //keep track of the best number of bins so far
                            nBinsBest            = nBins;
                            nDivisionsBest       = nDivisions;
                            nNinesBest           = nNines;
                            nTwentyFivesBest     = nTwentyFives;
                            nFortyNinesBest      = nFortyNines;
                            nEleventyElevensBest = nEleventyElevens;
                        }
                    }
                }
            }
        }
    } while ( TMath::Power( 4, nDivisions + 1 ) <
              nBinsTarget );    //if 4^n > target then we've gone far enough

    std::cout << "Using " << nBinsBest << " bins" << std::endl;

    //fill the vector with the divisions that we want to make
    for ( Int_t i = 0; i < nEleventyElevensBest; ++i ) {
        divisions.push_back( 11 );
    }
    for ( Int_t i = 0; i < nFortyNinesBest; ++i ) {
        divisions.push_back( 7 );
    }
    for ( Int_t i = 0; i < nTwentyFivesBest; ++i ) {
        divisions.push_back( 5 );
    }
    for ( Int_t i = 0; i < nNinesBest; ++i ) {
        divisions.push_back( 3 );
    }
    for ( Int_t i = 0; i < nDivisionsBest - nNinesBest - nTwentyFivesBest - nFortyNinesBest -
                               nEleventyElevensBest;
          ++i ) {
        divisions.push_back( 2 );
    }
}

void LauCalcChiSq::getHisto( const Double_t xMin,
                             const Double_t xMax,
                             const Double_t yMin,
                             const Double_t yMax,
                             const Double_t* xs,
                             const Double_t* ys,
                             const Int_t nEntries,
                             const std::vector<Int_t>& divisions,
                             const UInt_t iter )
{

    //If it's the last iteration create the bin and return
    if ( iter == divisions.size() ) {
        Double_t* x_new = new Double_t[5];
        Double_t* y_new = new Double_t[5];
        x_new[0]        = xMin;
        x_new[1]        = xMin;
        x_new[2]        = xMax;
        x_new[3]        = xMax;
        x_new[4]        = xMin;
        y_new[0]        = yMin;
        y_new[1]        = yMax;
        y_new[2]        = yMax;
        y_new[3]        = yMin;
        y_new[4]        = yMin;
        theHisto_->AddBin( 5, x_new, y_new );
        if ( verbose_ )
            std::cout << "Adding bin from (" << xMin << "," << yMin << ") to (" << xMax << ","
                      << yMax << ")" << std::endl;
        return;
    }

    //If not the last iteration then divide the bin
    Int_t n_divx = divisions[iter];
    Int_t n_divy = divisions[iter];

    if ( verbose_ )
        std::cout << "Dividing bin from (" << xMin << "," << yMin << ") to (" << xMax << "," << yMax
                  << ") into " << n_divx << " by " << n_divy << " subbins" << std::endl;

    Double_t* xIn = new Double_t[nEntries];
    Double_t* yIn = new Double_t[nEntries];
    Int_t* xIndex = new Int_t[nEntries + 2];
    Int_t* yIndex = new Int_t[nEntries + 2];

    Int_t xCountIn = 0;
    for ( Int_t i = 0; i < nEntries; ++i ) {
        if ( ( xs[i] < xMin ) || ( xs[i] > xMax ) || ( ys[i] < yMin ) || ( ys[i] > yMax ) )
            continue;
        xIn[xCountIn] = xs[i];
        ++xCountIn;
    }

    //find the delimitting x and y values for the sub bins
    Double_t xLimits[n_divx + 1];
    Double_t yLimits[n_divx][n_divy + 1];

    //first sort entries in x and divide bin into equally populated bins in x
    TMath::Sort( xCountIn, xIn, xIndex, false );

    xLimits[0]      = xMin;
    xLimits[n_divx] = xMax;
    for ( Int_t nDivx = 0; nDivx < n_divx; ++nDivx ) {
        if ( nDivx < ( n_divx - 1 ) ) {
            xLimits[nDivx + 1] = xIn[xIndex[xCountIn * ( nDivx + 1 ) / n_divx]];
        }

        //for each bin in x divide into equally populated bins in y
        yLimits[nDivx][0]      = yMin;
        yLimits[nDivx][n_divy] = yMax;
        Int_t yCountIn         = 0;

        for ( Int_t i = 0; i < nEntries; ++i ) {
            if ( ( xs[i] < xMin ) || ( xs[i] > xMax ) || ( ys[i] < yMin ) || ( ys[i] > yMax ) )
                continue;
            if ( ( xs[i] < xLimits[nDivx] ) || ( xs[i] >= xLimits[nDivx + 1] ) ||
                 ( ys[i] < yMin ) || ( ys[i] > yMax ) )
                continue;
            yIn[yCountIn] = ys[i];
            ++yCountIn;
        }

        TMath::Sort( yCountIn, yIn, yIndex, false );

        for ( Int_t nDivy = 1; nDivy < n_divy; ++nDivy ) {
            yLimits[nDivx][nDivy] = yIn[yIndex[yCountIn * nDivy / n_divy]];
        }
    }

    delete[] xIn;
    delete[] yIn;
    delete[] xIndex;
    delete[] yIndex;

    //call for each sub bin
    for ( Int_t nDivx = 0; nDivx < n_divx; ++nDivx ) {
        for ( Int_t nDivy = 0; nDivy < n_divy; ++nDivy ) {
            this->getHisto( xLimits[nDivx],
                            xLimits[nDivx + 1],
                            yLimits[nDivx][nDivy],
                            yLimits[nDivx][nDivy + 1],
                            xs,
                            ys,
                            nEntries,
                            divisions,
                            iter + 1 );
        }
    }
}
