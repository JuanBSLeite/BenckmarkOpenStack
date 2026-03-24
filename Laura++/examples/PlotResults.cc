
/*
Copyright 2009 University of Warwick

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

#include "TArrow.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TLeaf.h"
#include "TMath.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"

#include <cstdlib>
#include <fstream>
#include <iostream>

using std::cerr;
using std::cout;
using std::endl;

void plotDists( const TString& fileName, const TString& treeName, Int_t numExpts );
void plotErrors( const TString& fileName, const TString& treeName, Int_t numExpts );
void plotPulls( const TString& fileName, const TString& treeName, Int_t numExpts );
void plotCorrelations( const TString& fileName,
                       const TString& treeName,
                       Int_t numExpts,
                       Bool_t onlyYieldCorrs = kTRUE );

int main( const int argc, const char** argv )
{
    if ( argc < 2 ) {
        std::cerr << "Usage: " << argv[0]
                  << " <command = plotpulls, ploterrors, plotdists, plotcorrs, plotcorrsall> [nExpts] [inputFileName] [treeName]"
                  << std::endl;
    }

    gROOT->SetStyle( "Plain" );
    gStyle->SetOptStat( 1000000 );
    gStyle->SetOptFit( 0111 );

    Int_t nExpt( 1 );

    TString inputFileName( "output.root" );
    TString treeName( "fitResults" );

    Bool_t doPlotPulls( kFALSE );
    Bool_t doPlotErrors( kFALSE );
    Bool_t doPlotDists( kFALSE );
    Bool_t doPlotCorrs( kFALSE );

    Bool_t onlyYieldCorrs( kTRUE );

    if ( argc > 1 ) {
        TString command( argv[1] );
        if ( command.Contains( "plotpulls", TString::kIgnoreCase ) ) {
            doPlotPulls  = kTRUE;
            doPlotErrors = kFALSE;
            doPlotDists  = kFALSE;
            doPlotCorrs  = kFALSE;
        } else if ( command.Contains( "ploterrors", TString::kIgnoreCase ) ) {
            doPlotPulls  = kFALSE;
            doPlotErrors = kTRUE;
            doPlotDists  = kFALSE;
            doPlotCorrs  = kFALSE;
        } else if ( command.Contains( "plotdists", TString::kIgnoreCase ) ) {
            doPlotPulls  = kFALSE;
            doPlotErrors = kFALSE;
            doPlotDists  = kTRUE;
            doPlotCorrs  = kFALSE;
        } else if ( command.Contains( "plotcorrs", TString::kIgnoreCase ) ) {
            doPlotPulls  = kFALSE;
            doPlotErrors = kFALSE;
            doPlotDists  = kFALSE;
            doPlotCorrs  = kTRUE;
            if ( command.Contains( "plotcorrsall", TString::kIgnoreCase ) ) {
                onlyYieldCorrs = kFALSE;
            }
        } else {
            cerr << "Unknown command: \"" << command << "\", exiting..." << endl;
            return EXIT_FAILURE;
        }
    }
    if ( argc > 2 ) {
        nExpt = atoi( argv[2] );
    }
    if ( argc > 3 ) {
        inputFileName = argv[3];
    }
    if ( argc > 4 ) {
        treeName = argv[4];
    }

    if ( doPlotPulls ) {
        plotPulls( inputFileName, treeName, nExpt );
    }
    if ( doPlotErrors ) {
        plotErrors( inputFileName, treeName, nExpt );
    }
    if ( doPlotDists ) {
        plotDists( inputFileName, treeName, nExpt );
    }
    if ( doPlotCorrs ) {
        plotCorrelations( inputFileName, treeName, nExpt, onlyYieldCorrs );
    }

    return EXIT_SUCCESS;
}

void plotPulls( const TString& fileName, const TString& treeName, Int_t numExpts )
{

    TString texName( "pulls.tex" );
    cout << "Opening file \"" << texName << "\" for writing table." << endl;
    std::ofstream fout( texName );

    // write the table heading
    fout << "\\begin{tabular}{lccc}" << endl;
    fout << "\\hline" << endl;
    fout << "Fit Parameter & Entries & Pull Mean & Pull Width \\\\" << endl;
    fout << "\\hline" << endl;

    // open the file and get the tree
    TFile* file = TFile::Open( fileName );
    if ( ! file ) {
        cerr << "Can't open file: " << fileName << endl;
        return;
    }
    TTree* tree = (TTree*)file->Get( treeName );
    if ( ! tree ) {
        cerr << "Can't find tree: " << treeName << endl;
        return;
    }

    // get the array of leaves from the tree
    TObjArray* pArray = tree->GetListOfLeaves();
    if ( ! pArray ) {
        cerr << "Can't get leaf array from tree." << endl;
        return;
    }
    TObjArray& array = *pArray;

    // create a canvas, histogram and gaussian function
    TCanvas* c1 = new TCanvas( "c1", "c1" );
    Int_t nBins = static_cast<Int_t>( numExpts / 5.0 );
    TH1F* histo = new TH1F( "histo", "", nBins, -5.0, 5.0 );
    TF1 gaus( "gaus", "gaus(0)", -5.0, 5.0 );
    gaus.SetParName( 0, "N" );
    gaus.SetParName( 1, "#mu" );
    gaus.SetParName( 2, "#sigma" );

    Int_t arraySize( array.GetEntries() );
    cout << "Looping over array of size: " << arraySize << endl;
    for ( Int_t i( 0 ); i < arraySize; ++i ) {

        // reset the starting values of the Gaussian parameters
        gaus.SetParameter( 0, 5.0 );
        gaus.SetParameter( 1, 0.0 );
        gaus.SetParameter( 2, 1.0 );

        // grab the leaf from the array
        // check that it's valid
        TLeaf* leaf = (TLeaf*)array.At( i );
        if ( ! leaf ) {
            break;
        }

        // get the leaf name
        // check that it contains "Pull"
        TString name = leaf->GetName();
        if ( ! name.Contains( "Pull" ) ) {
            continue;
        }

        // create the plotting command from the leaf name
        TString plotCmd( name );
        plotCmd += ">>histo";

        // set the histogram axis title
        histo->SetXTitle( name );

        // draw the data into the histogram
        tree->Draw( plotCmd );

        // draw the histogram and fit it
        histo->Draw( "E1P" );
        histo->Fit( "gaus", "EQ" );
        histo->Draw( "E1P" );

        // save the plot to an eps file
        TString plotName( name );
        plotName += ".eps";
        c1->SaveAs( plotName );

        // write the pull data to the text file
        TString compName( name );
        compName.ReplaceAll( "_", "\\_" );
        fout << compName << " & " << histo->Integral() << " & "
             << "$" << gaus.GetParameter( 1 ) << " \\pm " << gaus.GetParError( 1 ) << "$"
             << " & "
             << "$" << gaus.GetParameter( 2 ) << " \\pm " << gaus.GetParError( 2 ) << "$"
             << " \\\\" << endl;
    }

    // write the table footer
    fout << "\\hline" << endl;
    fout << "\\end{tabular}" << endl;
    fout.close();

    delete histo;
    delete c1;
    delete tree;
    file->Close();
    delete file;
}

void plotErrors( const TString& fileName, const TString& treeName, Int_t numExpts )
{
    // open the file and get the tree
    TFile* file = TFile::Open( fileName );
    if ( ! file ) {
        cerr << "Can't open file: " << fileName << endl;
        return;
    }
    TTree* tree = (TTree*)file->Get( treeName );
    if ( ! tree ) {
        cerr << "Can't find tree: " << treeName << endl;
        return;
    }

    // get the array of leaves from the tree
    TObjArray* pArray = tree->GetListOfLeaves();
    if ( ! pArray ) {
        cerr << "Can't get leaf array from tree." << endl;
        return;
    }
    TObjArray& array = *pArray;

    // create a canvas, histogram pointer and gaussian function
    TCanvas* c1 = new TCanvas( "c1", "c1" );
    Int_t nBins = static_cast<Int_t>( numExpts / 5.0 );
    TH1F* histo( 0 );
    TF1 gaus( "gaus", "gaus(0)", -5.0, 5.0 );
    gaus.SetParName( 0, "N" );
    gaus.SetParName( 1, "#mu" );
    gaus.SetParName( 2, "#sigma" );

    Int_t arraySize( array.GetEntries() );
    cout << "Looping over array of size: " << arraySize << endl;
    for ( Int_t i( 0 ); i < arraySize; ++i ) {

        // grab the leaf from the array
        // check that it's valid
        TLeaf* leaf = (TLeaf*)array.At( i );
        if ( ! leaf ) {
            break;
        }

        // get the leaf name
        // check that it contains "Error"
        TString name = leaf->GetName();
        if ( ! name.Contains( "Error" ) ) {
            continue;
        }

        // create the plotting command from the leaf name
        TString plotCmd( name );
        plotCmd += ">>histo";

        // create the histogram from the minimum and maximum values of the leaf
        Double_t min = tree->GetMinimum( name );
        Double_t max = tree->GetMaximum( name );
        histo        = new TH1F( "histo", "", 50, min, max );
        tree->Draw( plotCmd );

        // now find the mean and RMS of the leaf
        Double_t mean = histo->GetMean();
        Double_t rms  = histo->GetRMS();

        // delete and recreate the histo as mean+/-5sigma
        min = mean - 5.0 * rms;
        max = mean + 5.0 * rms;
        delete histo;
        histo = new TH1F( "histo", "", nBins, min, max );
        tree->Draw( plotCmd );

        // set the histogram axis title
        histo->SetXTitle( name );

        // draw the data into the histogram
        tree->Draw( plotCmd );

        // reset the starting values of the Gaussian parameters
        gaus.SetRange( min, max );
        gaus.SetParameter( 0, histo->GetMaximum() * 50.0 / nBins );
        gaus.SetParameter( 1, histo->GetMean() );
        gaus.SetParameter( 2, histo->GetRMS() );

        // draw the histogram and fit it
        histo->Draw( "E1P" );
        histo->Fit( "gaus", "EQ" );
        histo->Draw( "E1P" );

        // save the plot to an eps file
        TString plotName( name );
        plotName += ".eps";
        c1->SaveAs( plotName );

        // delete the histogram
        delete histo;
    }

    delete c1;
    delete tree;
    file->Close();
    delete file;
}

void plotDists( const TString& fileName, const TString& treeName, Int_t numExpts )
{

    TString texName( "dists.tex" );
    cout << "Opening file \"" << texName << "\" for writing table." << endl;
    std::ofstream fout( texName );

    // write the table heading
    fout << "\\begin{tabular}{lcccccc}" << endl;
    fout << "\\hline" << endl;
    fout << "Fit Parameter & Entries & True Value & Distribution Mean & Distribution Width & Residual Mean & Residual Width \\\\"
         << endl;
    fout << "\\hline" << endl;

    TFile* file = TFile::Open( fileName );
    if ( ! file ) {
        cerr << "Can't open file: " << fileName << endl;
        return;
    }
    TTree* tree = (TTree*)file->Get( treeName );
    if ( ! tree ) {
        cerr << "Can't find tree: " << treeName << endl;
        return;
    }

    // get the array of leaves from the tree
    TObjArray* pArray = tree->GetListOfLeaves();
    if ( ! pArray ) {
        cerr << "Can't get leaf array from tree." << endl;
        return;
    }
    TObjArray& array = *pArray;

    // create a canvas and gaussian function
    TCanvas* c1 = new TCanvas( "c1", "c1" );
    TF1 gaus( "gaus", "gaus(0)", -5, 5 );
    gaus.SetParName( 0, "N" );
    gaus.SetParName( 1, "#mu" );
    gaus.SetParName( 2, "#sigma" );

    Int_t nBins = static_cast<Int_t>( numExpts / 5.0 );

    Int_t arraySize( array.GetEntries() );
    cout << "Looping over array of size: " << arraySize << endl;
    for ( Int_t i( 0 ); i < arraySize; ++i ) {

        // grab the leaf from the array
        // check that it's valid
        TLeaf* leaf = (TLeaf*)array.At( i );
        if ( ! leaf ) {
            break;
        }

        // get the leaf name
        // check that it is not a pull, error, truth value or correlation
        // also check that it isn't "iExpt" or "fitStatus"
        TString name = leaf->GetName();
        if ( name.EndsWith( "_Pull" ) || name.EndsWith( "_Error" ) || name.EndsWith( "_True" ) ||
             name.BeginsWith( "corr__" ) || name.EndsWith( "_GCC" ) || name == "iExpt" ||
             name == "fitStatus" ) {
            continue;
        }

        // only plot variables that are varying in the fit and fit fractions
        TString namePull( name );
        namePull += "_Pull";
        if ( ( ! name.Contains( "TotFF" ) ) && ( ! name.Contains( "FitFrac" ) ) &&
             ( ! tree->GetBranch( namePull ) ) ) {
            continue;
        }

        // create the plotting command from the leaf name
        TString plotCmd( name );
        plotCmd += ">>histo";

        // create the histogram from the minimum and maximum values of the leaf
        Double_t min = tree->GetMinimum( name );
        Double_t max = tree->GetMaximum( name );
        TH1F* histo  = new TH1F( "histo", "", 50, min, max );
        tree->Draw( plotCmd );

        // now find the mean and RMS of the leaf
        // the recreate the histo as mean+/-5sigma
        Double_t mean = histo->GetMean();
        Double_t rms  = histo->GetRMS();
        Double_t xMin = mean - 5.0 * rms;
        Double_t xMax = mean + 5.0 * rms;
        if ( ( name.Contains( "TotFF" ) || ( name.Contains( "FitFrac" ) && name.Contains( "Sq" ) ) ) &&
             xMin < 0.0 ) {
            xMin = 0.0;
        }
        delete histo;
        histo = new TH1F( "histo", "", nBins, xMin, xMax );
        tree->Draw( plotCmd );

        // set the Gaussian parameters to the mean and RMS
        gaus.SetParameter( 0, 5.0 );
        gaus.SetParameter( 1, histo->GetMean() );
        gaus.SetParameter( 2, histo->GetRMS() );

        // set the histogram axis title
        histo->SetXTitle( name );

        // draw the histogram and fit it
        histo->Draw( "E1P" );
        histo->Fit( "gaus", "EQ" );
        histo->Draw( "E1P" );

        // plot the true value as an arrow
        TString trueName( name );
        trueName += "_True";

        Double_t maxTrueVal = tree->GetMaximum( trueName );
        Double_t minTrueVal = tree->GetMinimum( trueName );
        Bool_t goodTrueVal  = kTRUE;
        if ( TMath::Abs( maxTrueVal - minTrueVal ) > 1e-6 ) {
            goodTrueVal = kFALSE;
        }
        Double_t trueVal     = maxTrueVal;
        Double_t arrowHeight = histo->GetMaximum() * 0.2;
        TArrow arrow( trueVal, 0.0, trueVal, arrowHeight, 0.05, "<" );
        arrow.SetLineColor( kRed );
        if ( goodTrueVal ) {
            arrow.Draw();
        }

        // save the plot to an eps file
        TString plotName( name );
        plotName += ".eps";
        c1->SaveAs( plotName );

        // grab the values
        Double_t entries      = histo->Integral();
        Double_t distMean     = gaus.GetParameter( 1 );
        Double_t distMeanErr  = gaus.GetParError( 1 );
        Double_t distWidth    = gaus.GetParameter( 2 );
        Double_t distWidthErr = gaus.GetParError( 2 );

        // plot the residual
        delete histo;
        histo    = new TH1F( "histo", "", nBins, xMin - trueVal, xMax - trueVal );
        plotCmd  = name;
        plotCmd += "-";
        plotCmd += trueName;
        plotCmd += ">>histo";
        tree->Draw( plotCmd );

        TString xTitle  = name;
        xTitle         += "_Residual";
        histo->SetXTitle( xTitle );

        Double_t trueValue( 0.0 );
        Double_t residualMean( 0.0 );
        Double_t residualMeanErr( 0.0 );
        Double_t residualWidth( 0.0 );
        Double_t residualWidthErr( 0.0 );

        // set the Gaussian parameters to the mean and RMS
        gaus.SetParameter( 0, 5.0 );
        gaus.SetParameter( 1, histo->GetMean() );
        gaus.SetParameter( 2, histo->GetRMS() );

        // draw the histogram and fit it
        histo->Draw( "E1P" );
        histo->Fit( "gaus", "EQ" );
        histo->Draw( "E1P" );

        plotName  = name;
        plotName += "_Residual.eps";
        c1->SaveAs( plotName );

        // grab the values
        trueValue        = tree->GetMinimum( trueName );
        residualMean     = gaus.GetParameter( 1 );
        residualMeanErr  = gaus.GetParError( 1 );
        residualWidth    = gaus.GetParameter( 2 );
        residualWidthErr = gaus.GetParError( 2 );

        // write the distribution data to the text file
        TString compName( name );
        compName.ReplaceAll( "_", "\\_" );
        fout << compName << " & " << entries << " & " << trueValue << " & ";
        fout << "$" << distMean << " \\pm " << distMeanErr << "$"
             << " & ";
        fout << "$" << distWidth << " \\pm " << distWidthErr << "$"
             << " & ";
        fout << "$" << residualMean << " \\pm " << residualMeanErr << "$"
             << " & ";
        fout << "$" << residualWidth << " \\pm " << residualWidthErr << "$"
             << " & ";
        fout << " \\\\" << endl;

        // delete the histogram
        delete histo;
    }

    // write the table footer
    fout << "\\hline" << endl;
    fout << "\\end{tabular}" << endl;
    fout.close();

    delete c1;
    delete tree;
    file->Close();
    delete file;
}

void plotCorrelations( const TString& fileName,
                       const TString& treeName,
                       Int_t numExpts,
                       Bool_t onlyYieldCorrs )
{
    TFile* file = TFile::Open( fileName );
    if ( ! file ) {
        cerr << "Can't open file: " << fileName << endl;
        return;
    }
    TTree* tree = (TTree*)file->Get( treeName );
    if ( ! tree ) {
        cerr << "Can't find tree: " << treeName << endl;
        return;
    }

    // get the array of leaves from the tree
    TObjArray* pArray = tree->GetListOfLeaves();
    if ( ! pArray ) {
        cerr << "Can't get leaf array from tree." << endl;
        return;
    }
    TObjArray& array = *pArray;

    // create a canvas and gaussian function
    TCanvas* c1 = new TCanvas( "c1", "c1" );
    TF1 gaus( "gaus", "gaus(0)", -1, 1 );
    gaus.SetParName( 0, "N" );
    gaus.SetParName( 1, "#mu" );
    gaus.SetParName( 2, "#sigma" );

    Int_t nBins = static_cast<Int_t>( numExpts / 5.0 );

    Int_t arraySize( array.GetEntries() );
    cout << "Looping over array of size: " << arraySize << endl;
    for ( Int_t i( 0 ); i < arraySize; ++i ) {

        // grab the leaf from the array
        // check that it's valid
        TLeaf* leaf = (TLeaf*)array.At( i );
        if ( ! leaf ) {
            break;
        }

        // get the leaf name and check that it is a correlation
        TString name = leaf->GetName();
        if ( ! name.BeginsWith( "corr__" ) && ! name.EndsWith( "GCC" ) ) {
            continue;
        }

        // check that if it is not a global correlation coeff then
        // it is a correlation between two yield parameters (if requested)
        if ( onlyYieldCorrs && ! name.EndsWith( "GCC" ) ) {
            TString namePart( "" );
            Ssiz_t a( 0 );
            name.Tokenize( namePart, a, "_" );
            name.Tokenize( namePart, a, "_" );
            if ( ! namePart.Contains( "Events" ) ) {
                continue;
            }
            name.Tokenize( namePart, a, "_" );
            if ( ! namePart.Contains( "Events" ) ) {
                continue;
            }
        }

        // create the plotting command from the leaf name
        TString plotCmd( name );
        plotCmd += ">>histo";

        // create the histogram from the minimum and maximum values of the leaf
        Double_t min = tree->GetMinimum( name );
        Double_t max = tree->GetMaximum( name );
        TH1F* histo  = new TH1F( "histo", "", 50, min, max );
        tree->Draw( plotCmd );

        // now find the mean and RMS of the leaf
        // the recreate the histo as mean+/-5sigma
        Double_t mean = histo->GetMean();
        Double_t rms  = histo->GetRMS();
        delete histo;
        histo = new TH1F( "histo", "", nBins, mean - 5.0 * rms, mean + 5.0 * rms );
        tree->Draw( plotCmd );

        // set the Gaussian parameters to the mean and RMS
        gaus.SetParameter( 0, 5.0 );
        gaus.SetParameter( 1, histo->GetMean() );
        gaus.SetParameter( 2, histo->GetRMS() );

        // set the histogram axis title
        histo->SetXTitle( name );

        // draw the histogram and fit it
        histo->Draw( "E1P" );
        histo->Fit( "gaus", "EQ" );
        histo->Draw( "E1P" );

        // save the plot to an eps file
        TString plotName( name );
        plotName += ".eps";
        c1->SaveAs( plotName );

        // delete the histogram
        delete histo;
    }

    delete c1;
    delete tree;
    file->Close();
    delete file;
}
