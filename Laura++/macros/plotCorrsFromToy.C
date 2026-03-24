
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

void plotCorrsFromToy( TString fileName = "test.root", const Bool_t plotIndividualCorrs = kFALSE )
{
    gStyle->SetOptStat( 0 );
    gStyle->SetPalette( 1 );
    gStyle->SetPadTopMargin( 0.15 );
    gStyle->SetPadBottomMargin( 0.15 );
    gStyle->SetPadLeftMargin( 0.15 );
    gStyle->SetPadRightMargin( 0.15 );

    TCanvas* c1 = new TCanvas( "c1", "c1" );

    TFile* _file0     = TFile::Open( fileName );
    TTree* fitResults = dynamic_cast<TTree*>( _file0->Get( "fitResults" ) );
    fitResults->GetEntry( 0 );

    if ( fitResults->GetEntries() == 1 ) {
        std::cout << "Only one entry in fitResults tree ...\n"
                  << "\t--> you probably want to use plotCorrs.C instead" << std::endl;
        return;
    }

    TObjArray* listOfLeaves = fitResults->GetListOfLeaves();

    const UInt_t maxNVar( 100 );
    TString varNames[100];

    UInt_t nVar( 0 );
    for ( UInt_t iLeaf = 0; iLeaf < listOfLeaves->GetEntries(); ++iLeaf ) {
        TLeaf* leaf  = dynamic_cast<TLeaf*>( ( *listOfLeaves )[iLeaf] );
        TString name = leaf->GetName();
        if ( name.EndsWith( "_Error" ) ) {
            TString shortName( name( 0, name.Length() - 6 ) );
            //      std::cout << name << " " << shortName << std::endl;
            varNames[nVar] = shortName;
            ++nVar;

            if ( maxNVar < nVar ) {
                std::cerr << "Too many variables!!!" << std::endl;
                break;
            }
        }
    }

    TH2I corrHisto( "corrHisto", "correlations", nVar, -0.5, nVar - 0.5, nVar, -0.5, nVar - 0.5 );

    for ( UInt_t iVar = 0; iVar < nVar; ++iVar ) {
        corrHisto.GetXaxis()->SetBinLabel( iVar + 1, varNames[iVar] );
        corrHisto.GetYaxis()->SetBinLabel( iVar + 1, varNames[iVar] );

        for ( UInt_t jVar = 0; jVar < nVar; ++jVar ) {

            if ( jVar == iVar ) {
                corrHisto.SetBinContent( iVar + 1, jVar + 1, 100 );
                continue;
            }

            TString corrName = "corr__" + varNames[iVar] + "__" + varNames[jVar];

            if ( plotIndividualCorrs ) {
                fitResults->Draw( corrName );
                TString filename = corrName + ".eps";
                c1->SaveAs( filename );
            }

            // get correlation factor from the distribution of Y vs X
            TString drawCmd = varNames[jVar] + ":" + varNames[iVar] + ">>bob";
            fitResults->Draw( drawCmd );
            TH2* histo    = (TH2*)gDirectory->Get( "bob" );
            Double_t corr = histo->GetCorrelationFactor();
            corrHisto.SetBinContent( iVar + 1, jVar + 1, TMath::Nint( corr * 100 ) );
            delete histo;

            //         // get average correlation factor from toy fit results (MINUIT)
            //         Double_t avValue(0.0);
            //         for (UInt_t iResult=0; iResult<fitResults->GetEntries(); ++iResult) {
            //         fitResults->GetEntry(iResult);
            //         TLeaf *corrLeaf = dynamic_cast<TLeaf *>(fitResults->GetListOfLeaves()->FindObject(corrName.Data()));
            //         avValue += corrLeaf->GetValue();

            //         };
            //         avValue /= fitResults->GetEntries();
            //         corrHisto.SetBinContent( iVar+1, jVar+1, TMath::Nint(avValue*100) );
        }
    }

    corrHisto.SetMinimum( -100.0 );
    corrHisto.SetMaximum( 100.0 );

    corrHisto.GetZaxis()->SetTitle( "Linear correlation coefficient (%)" );
    corrHisto.GetZaxis()->SetTitleSize( 0.033 );
    corrHisto.GetZaxis()->SetLabelSize( 0.033 );
    corrHisto.GetXaxis()->SetLabelSize( 0.033 );
    corrHisto.GetYaxis()->SetLabelSize( 0.033 );
    corrHisto.GetXaxis()->LabelsOption( "v" );

    corrHisto.Draw( "colz" );
    corrHisto.Draw( "text same" );

    c1->SaveAs( "corrHisto.eps" );
}
