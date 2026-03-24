
/*
Copyright 2010 University of Warwick

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

void plotDataIsobars( TString fileName = "test.root" )
{
    gStyle->SetOptStat( 0 );

    TFile* _file0     = TFile::Open( fileName );
    TTree* fitResults = dynamic_cast<TTree*>( _file0->Get( "fitResults" ) );
    fitResults->GetEntry( 0 );

    TCanvas* c1 = new TCanvas( "c1", "c1", 500, 500 );
    TH2F histo( "histo", "histo", 150, -1.5, 1.5, 150, -1.5, 1.5 );
    histo.Draw();

    // draw grid
    gStyle->SetLineWidth( 1 );
    gStyle->SetLineStyle( 3 );
    TEllipse ellipse2( 0.0, 0.0, 1.0, 1.0 );
    ellipse2.Draw();
    TEllipse ellipse1( 0.0, 0.0, 0.5, 0.5 );
    ellipse1.Draw();
    TLine lh( -1.5, 0.0, 1.5, 0.0 );
    lh.Draw();
    TLine lv( 0.0, -1.5, 0.0, 1.5 );
    lv.Draw();

    // draw results
    TObjArray* listOfLeaves = fitResults->GetListOfLeaves();

    const UInt_t maxNVar( 100 );
    TString varNames[100];

    const UInt_t nLeaves = listOfLeaves->GetEntries();
    UInt_t nVar( 0 );
    for ( UInt_t iLeaf = 0; iLeaf < nLeaves; ++iLeaf ) {
        TLeaf* leaf  = dynamic_cast<TLeaf*>( ( *listOfLeaves )[iLeaf] );
        TString name = leaf->GetName();
        if ( name.BeginsWith( "A" ) && name.EndsWith( "_X" ) ) {
            TString shortName( name( 0, name.Length() - 2 ) );
            varNames[nVar] = shortName;
            ++nVar;

            if ( maxNVar < nVar ) {
                std::cerr << "Too many variables!!!" << std::endl;
                break;
            }
        }
    }

    for ( UInt_t iVar = 0; iVar < nVar; ++iVar ) {
        fitResults->SetMarkerColor( iVar + 2 );
        fitResults->SetMarkerStyle( 2 );

        TString drawCmd = varNames[iVar] + "_Y:" + varNames[iVar] + "_X";
        fitResults->Draw( drawCmd, "", "same" );
    }

    c1->SaveAs( "data-isobars.eps" );
}
