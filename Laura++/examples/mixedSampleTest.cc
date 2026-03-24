
/*
Copyright 2015 University of Warwick

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

// Code to perform the mixed sample unbinned goodness of fit test on the data and simulation provided.
//
// If the data and simulation are derived from the same underlying PDF then they will mix optimally.
// This is used to define the test statistic T = (1 / n_k (n_d + n_s) ) SUM_i^(n_d+n_s) SUM_k^(n_k) I(i,k)
// where I(i,k) is 1 if the ith event and its kth nearest neighbour are from the same sample and zero otherwise,
// n_d is the number of events in the data sample, n_s is the number of events in the simulation sample and
// n_k is the number of nearest neighbours to consider.
//
// For optimal mixing T should be Gaussian distributed with mean and width given by:
// mu_T      = n_d(n_d - 1) + n_s(n_s - 1) / (n_d + n_s)(n_d + n_s - 1)
// sigma_T^2 = ( 1 / n_k(n_d+n_s) ) ( (n_d n_s)/(n_d + n_s)^2 + 4(n_d^2 n_s^2)/(n_d + n_s)^4 )
//
// This script uses 500 bootstrapped samples taken from the data and simulation provided,
// calculates T from each sample, performs a Gaussian fit to (T-mu_T)/sigma_T and reports the
// chi^2 of that fit.

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TH1D.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TTree.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>

TRandom3* r( 0 );

//Class to store a single event
class Entry {
  public:
    Entry( Bool_t mc, Double_t mPrime, Double_t thPrime ) :
        _mc( mc ),
        _mPrime( mPrime ),
        _thPrime( thPrime )
    {
    }

    Bool_t isMC() { return _mc; }

    Double_t getDistSq( Entry* e )
    {
        return ( this->_mPrime - e->_mPrime ) * ( this->_mPrime - e->_mPrime ) +
               ( this->_thPrime - e->_thPrime ) * ( this->_thPrime - e->_thPrime );
    }

  private:
    Bool_t _mc;
    Double_t _mPrime;
    Double_t _thPrime;
};

//Function to calculate the test statistic, T, for a single bootstrapped sample
Double_t calcT( TString file0Name,
                TString tree0Name,
                TString file1Name,
                TString tree1Name,
                Int_t nd,
                Int_t nmc,
                Int_t nk )
{

    //Open the data and simulation files
    TFile *f0( 0 ), *f1( 0 );

    f0 = TFile::Open( file0Name );
    f1 = TFile::Open( file1Name );

    TTree* t0 = (TTree*)f0->Get( tree0Name );
    TTree* t1 = (TTree*)f1->Get( tree1Name );

    Double_t mprime( 0. ), thprime( 0. );

    t0->SetBranchAddress( "mPrime", &mprime );
    t0->SetBranchAddress( "thPrime", &thprime );
    t1->SetBranchAddress( "mPrime", &mprime );
    t1->SetBranchAddress( "thPrime", &thprime );

    Int_t n0 = t0->GetEntries();
    Int_t n1 = t1->GetEntries();

    std::vector<Int_t> entries0;
    std::vector<Int_t> entries1;

    //Get randomised bootstrapped samples of data and simulation
    for ( Int_t i = 0; i < nmc; ++i ) {
        entries0.push_back( (Int_t)( r->Rndm() * n0 ) );
    }
    for ( Int_t i = 0; i < nd; ++i ) {
        entries1.push_back( (Int_t)( r->Rndm() * n1 ) );
    }

    //Sort entries to speed up extraction
    std::sort( entries0.begin(), entries0.end() );
    std::sort( entries1.begin(), entries1.end() );

    //Extract required entries from TTrees and store
    std::vector<Entry*> entries;

    for ( Int_t i = 0; i < nmc; ++i ) {
        t0->GetEntry( entries0[i] );
        entries.push_back( new Entry( true, mprime, thprime ) );
    }
    for ( Int_t i = 0; i < nd; ++i ) {
        t1->GetEntry( entries1[i] );
        entries.push_back( new Entry( false, mprime, thprime ) );
    }

    f0->Close();
    f1->Close();

    entries0.clear();
    entries1.clear();

    //Construct a matrix of the distances between all selected points
    TMatrixD distances( nmc + nd, nmc + nd );

    for ( Int_t i = 0; i < nmc + nd; ++i ) {
        for ( Int_t j = i; j < nmc + nd; ++j ) {
            distances( i, j ) = entries[i]->getDistSq( entries[j] );
        }
    }

    //Now calculate T
    Double_t T( 0. );

    for ( Int_t i = 0; i < nmc + nd; ++i ) {
        //For each entry, first sort all other entries by their distance
        std::multimap<Double_t, Int_t> distsSorted;
        for ( Int_t j = 0; j < i; ++j ) {
            distsSorted.insert( std::pair<Double_t, Int_t>( distances( j, i ), j ) );
        }
        for ( Int_t j = i + 1; j < nmc + nd; ++j ) {
            distsSorted.insert( std::pair<Double_t, Int_t>( distances( i, j ), j ) );
        }
        std::multimap<Double_t, Int_t>::iterator it = distsSorted.begin();
        //Use k nearest neighbours to calculate T
        for ( Int_t k = 0; k < nk; ++k ) {
            if ( entries[( *it ).second]->isMC() == entries[i]->isMC() )
                ++T;
            ++it;
        }
    }

    //Apply the normalisation factor
    T /= nk * ( nd + nmc );

    std::vector<Entry*>::iterator it  = entries.begin();
    std::vector<Entry*>::iterator end = entries.end();

    //Clean up
    for ( ; it != end; ++it ) {
        delete ( *it );
        *it = 0;
    }
    entries.clear();

    return T;
}

int main( int argc, char** argv )
{

    gStyle->SetOptStat( 0 );

    Int_t nd = 100;    // the number of data entries to use per sample
    Int_t nmcOvernd =
        10;    // the number of simulation entries to use divided by the number of data entries
    Int_t nk = 10;                               // the number of nearest neighbours to consider
    TString outFile( "mixedSampleTest.pdf" );    // the name of the output file for the histogram
    TString mcFile( "" );                        // the file containing the simulation
    TString mcTree( "" );                        // the tree containing the simulation
    TString dataFile( "" );                      // the file containing the data
    TString dataTree( "" );                      // the tree containing the data
    Int_t seed                = 0;               // seed for the bootstrapping
    Int_t floatGaussianParams = 0;    // whether to float the mean and width of the Gaussian or not

    if ( argc > 4 ) {
        mcFile   = argv[1];
        mcTree   = argv[2];
        dataFile = argv[3];
        dataTree = argv[4];
        if ( argc > 5 ) {
            seed = atoi( argv[5] );
            if ( argc > 6 ) {
                nd = atoi( argv[6] );
                if ( argc > 7 ) {
                    nmcOvernd = atoi( argv[7] );
                    if ( argc > 8 ) {
                        nk = atoi( argv[8] );
                        if ( argc > 9 ) {
                            floatGaussianParams = atoi( argv[9] );
                            if ( argc > 10 ) {
                                outFile = argv[10];
                            }
                        }
                    }
                }
            }
        }
    } else {
        std::cout << "Usage:\n";
        std::cout << argv[0]
                  << " <mcFile> <mcTree> <dataFile> <dataTree> [seed = 0] [nd = 100] [nmcOvernd = 10] [nk = 10] [floatGaussianParams = 0] [outFile = \"mixedSampleTest.pdf\"]"
                  << std::endl;
        return 0;
    }

    Int_t nmc = nmcOvernd * nd;

    //Calculate expected values of mu_T and sigma_T
    Double_t na( nd ), nb( nmc ), n( nd + nmc );
    Double_t muT  = ( na * ( na - 1 ) + nb * ( nb - 1 ) ) / ( n * ( n - 1 ) );
    Double_t sigT = TMath::Sqrt(
        ( na * nb / ( n * n ) + 4 * na * na * nb * nb / ( n * n * n * n ) ) / ( (n)*nk ) );

    TH1D* h = new TH1D( "h", "", 50, -5, 5 );
    r       = new TRandom3( seed );

    //Calculate the test statistic for 500 bootstrapped samples and fill the histogram with (T-mu_T)/sigma_T
    for ( Int_t i = 0; i < 500; ++i ) {
        if ( i % 10 == 0 )
            std::cout << "Processing sample " << i + 1 << " of 500..." << std::endl;
        Double_t val = ( calcT( mcFile, mcTree, dataFile, dataTree, nd, nmc, nk ) - muT ) / sigT;
        std::cout << val << std::endl;
        h->Fill( val );
    }

    //Perform a Gaussian fit to (T-mu_T)/sigma_T
    TCanvas can;
    TF1 gaus( "gaus", "gaus(0)", -5.0, 5.0 );
    gaus.SetParameter( 0, 100. );
    if ( floatGaussianParams == 1 ) {
        gaus.SetParameter( 0, 500.0 );
        gaus.SetParameter( 1, 0.0 );
        gaus.SetParameter( 2, 1.0 );
    } else {
        gaus.SetParameter( 0, 500.0 );
        gaus.FixParameter( 1, 0.0 );
        gaus.FixParameter( 2, 1.0 );
    }
    TFitResultPtr p = h->Fit( "gaus", "S" );
    h->GetFunction( "gaus" )->SetLineColor( kRed );
    h->GetXaxis()->SetTitle( "(T-#mu_{T})/#sigma_{T}" );
    h->Draw( "E" );
    can.SaveAs( outFile );

    //Report the chi^2 of the Gaussian fit
    std::cout << "chi2/ndf=" << p->Chi2() << "/" << p->Ndf() << std::endl;

    delete r;
    delete h;
}
