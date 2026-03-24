
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

// Code to perform the point to point similarity unbinned goodness of fit test on a subsample of the data and simulation provided.
//
// A test statistic is constructed, T = (1/n_d^2) SUM_(i,j>i)^(n_d) psi(|x_i^d - x_j^d|) - 1/(n_d * n_s) SUM_(i,j)^(n_d,n_s) psi(|x_i^d - x_j^s|),
// where x_i^d is the position of the ith data entry in the phase space and x_i^s is the position of the ith simulated entry,
// n_d and n_s are the number of data and simulation entries, respectively,
// and psi is a function that drops off for larger distances. This script uses exp(-(m'^2+theat'^2)/(2sigma^2)) where sigma should be set to roughly the
// size of the narrowest structure in the square Dalitz plot.
//
// The distribution of this test statistic is unknown so the permutation test is used.
// The data and simulation samples are divided into n subsamples and a p-value is determined for subsample i.
// For a good match between data and simulation, the distribution of p-values for the n samples should be uniform between 0 and 1.
// For a bad match, the p-values will peak at zero.

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
#include <set>
#include <vector>

TRandom3* r( 0 );

//object to store a single entry
class Entry {
  public:
    Entry( Bool_t data, Double_t mPrime, Double_t thPrime, Int_t nPerms ) :
        _data( data ),
        _permutedData( nPerms, 0 ),
        _mPrime( mPrime ),
        _thPrime( thPrime )
    {
    }

    Bool_t isData() { return _data; }

    Bool_t isPermutedData( Int_t i ) { return _permutedData[i]; }

    void setPermutedData( Int_t i, Bool_t permutedData ) { _permutedData[i] = permutedData; }

    //get squared distance from entry e
    Double_t getDistSq( Entry* e )
    {
        return ( this->_mPrime - e->_mPrime ) * ( this->_mPrime - e->_mPrime ) +
               ( this->_thPrime - e->_thPrime ) * ( this->_thPrime - e->_thPrime );
    }

  private:
    Bool_t _data;                         //whether entry is data
    std::vector<Bool_t> _permutedData;    //whether entry is considered data in each permuted dataset
    Double_t _mPrime;                     //SDP co-ordinate
    Double_t _thPrime;                    //SDP co-ordinate
};

//function to load entries from data and MC files for sample i of n
Int_t loadEntries( std::vector<Entry*>& entries,
                   TString file0Name,
                   TString tree0Name,
                   TString file1Name,
                   TString tree1Name,
                   Int_t sample,
                   Int_t nSamples,
                   Int_t nPerms )
{

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

    Int_t n0 = t0->GetEntries() / nSamples;
    Int_t n1 = t1->GetEntries() / nSamples;

    entries.reserve( n0 + n1 );

    //load independent entries for each sample
    for ( Int_t i = 0; i < n0; ++i ) {
        t0->GetEntry( i * nSamples + sample );
        entries.push_back( new Entry( false, mprime, thprime, nPerms ) );
    }
    for ( Int_t i = 0; i < n1; ++i ) {
        t1->GetEntry( i * nSamples + sample );
        entries.push_back( new Entry( true, mprime, thprime, nPerms ) );
    }

    f0->Close();
    f1->Close();

    return n1;
}

//method to setup nPerms permuted datasets
void makePermutations( std::vector<Entry*>& entries, Int_t nd, Int_t nPerms )
{
    Int_t n = entries.size();

    for ( Int_t i = 0; i < nPerms; ++i ) {
        std::set<Int_t> dataSet;
        Int_t entry( 0 );
        //pick nd random entries
        for ( Int_t j = 0; j < nd; ++j ) {
            //throw random numbers until we find one that hasn't been inserted already
            do {
                entry = r->Integer( n );
            } while ( ! dataSet.insert( entry ).second );
        }

        std::set<Int_t>::iterator it  = dataSet.begin();
        std::set<Int_t>::iterator end = dataSet.end();

        //set the flags
        for ( ; it != end; ++it ) {
            entries[( *it )]->setPermutedData( i, true );
        }
    }
}

//method to calculate the p-value for this sample
Double_t calcPVal( std::vector<Entry*>& entries, Int_t nd, Int_t nPerms, Double_t sigma )
{

    Int_t n   = entries.size();
    Int_t nmc = n - nd;

    Double_t distSq( 0. );
    Double_t psi( 0. );
    Double_t twoSigSq( 2 * sigma * sigma );

    Double_t T1( 0. ), T2( 0. );
    std::vector<Double_t> T1Perms( nPerms, 0. );
    std::vector<Double_t> T2Perms( nPerms, 0. );

    //loop over all pairs of entries
    for ( Int_t i = 0; i < n; ++i ) {
        std::cout << i << " of " << n << std::endl;
        for ( Int_t j = i; j < n; ++j ) {
            distSq = entries[i]->getDistSq( entries[j] );
            psi    = TMath::Exp( -distSq / twoSigSq );
            //calculate contribution from this pair to the test statistic for data
            if ( entries[i]->isData() ) {
                if ( entries[j]->isData() ) {
                    T1 += psi;
                } else {
                    T2 += psi;
                }
            } else if ( entries[j]->isData() ) {
                T2 += psi;
            }
            //calculate contributions to test statistics for the permutations of the data
            for ( Int_t k = 0; k < nPerms; ++k ) {
                if ( entries[i]->isPermutedData( k ) ) {
                    if ( entries[j]->isPermutedData( k ) ) {
                        T1Perms[k] += psi;
                    } else {
                        T2Perms[k] += psi;
                    }
                } else if ( entries[j]->isPermutedData( k ) ) {
                    T2Perms[k] += psi;
                }
            }
        }
    }

    Double_t T = T1 / ( nd * nd ) - T2 / ( nd * nmc );
    Double_t pVal( 0. );
    for ( Int_t i = 0; i < nPerms; ++i ) {
        //count how many permutations give a worse test statistic than the data
        if ( T1Perms[i] / ( nd * nd ) - T2Perms[i] / ( nd * nmc ) > T )
            ++pVal;
    }
    pVal /= nPerms;

    return pVal;
}

int main( int argc, char** argv )
{

    gStyle->SetOptStat( 0 );

    TString mcFile( "" );      // the file to load simulated events from
    TString mcTree( "" );      // the tree containing the simulation
    TString dataFile( "" );    // the file to load data events from
    TString dataTree( "" );    // the tree containing the data
    Int_t seed     = 0;        // seed the random number used to produce permutations
    Int_t iSample  = 0;        // the index of this sample
    Int_t nSamples = 1;        // the number of samples to split the events into
    Int_t nPerms   = 100;      // the number of permutations to test against
    Double_t sigma( 0.01 );    // the width of the Gaussian decay used to weight events by distance

    if ( argc > 2 ) {
        mcFile   = argv[1];
        mcTree   = argv[2];
        dataFile = argv[3];
        dataTree = argv[4];
        if ( argc > 5 ) {
            seed = atoi( argv[5] );
            if ( argc > 6 ) {
                iSample = atoi( argv[6] );
                if ( argc > 7 ) {
                    nSamples = atoi( argv[7] );
                    if ( argc > 8 ) {
                        nPerms = atoi( argv[8] );
                        if ( argc > 9 ) {
                            sigma = atof( argv[9] );
                        }
                    }
                }
            }
        }
    } else {
        std::cout << "Usage:\n";
        std::cout << argv[0]
                  << "<mcFile> <mcTree> <dataFile> <dataTree> [seed = 0] [iSample = 0] [nSamples = 1] [nPerms = 100] [sigma = 0.01]"
                  << std::endl;
        return 0;
    }

    r = new TRandom3( seed );

    std::vector<Entry*> entries;

    //load entries from MC and from data
    Int_t nd = loadEntries( entries, mcFile, mcTree, dataFile, dataTree, iSample, nSamples, nPerms );

    //make the permuted-data sets
    makePermutations( entries, nd, nPerms );

    //calculate the p-value
    std::cout << calcPVal( entries, nd, nPerms, sigma );

    delete r;
}
