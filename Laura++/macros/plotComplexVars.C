
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

void plotComplexVars( TString fileName = "test.root" )
{
    TFile* _file0     = TFile::Open( fileName );
    TTree* fitResults = dynamic_cast<TTree*>( _file0->Get( "fitResults" ) );
    fitResults->GetEntry( 0 );

    TCanvas* c1 = new TCanvas( "c1", "c1", 500, 500 );
    TH2F histo( "histo", "histo", 120, -1.2, 1.2, 120, -1.2, 1.2 );
    histo.Draw();

    // draw grid
    gStyle->SetLineWidth( 1 );
    gStyle->SetLineStyle( 3 );
    TEllipse ellipse2( 0.0, 0.0, 1.0, 1.0 );
    ellipse2.Draw();
    TEllipse ellipse1( 0.0, 0.0, 0.5, 0.5 );
    ellipse1.Draw();
    TLine lh( -1.2, 0.0, 1.2, 0.0 );
    lh.Draw();
    TLine lv( 0.0, -1.2, 0.0, 1.2 );
    lv.Draw();

    // draw results
    fitResults->SetMarkerColor( 2 );
    fitResults->SetMarkerStyle( 1 );
    fitResults->Draw( "A1_Y:A1_X", "", "same" );
    fitResults->SetMarkerColor( 1 );
    fitResults->SetMarkerStyle( 2 );
    fitResults->Draw( "A1_Y_True:A1_X_True", "", "same" );

    fitResults->SetMarkerColor( 3 );
    fitResults->SetMarkerStyle( 1 );
    fitResults->Draw( "A2_Y:A2_X", "", "same" );
    fitResults->SetMarkerColor( 1 );
    fitResults->SetMarkerStyle( 2 );
    fitResults->Draw( "A2_Y_True:A2_X_True", "", "same" );

    fitResults->SetMarkerColor( 4 );
    fitResults->SetMarkerStyle( 1 );
    fitResults->Draw( "A3_Y:A3_X", "", "same" );
    fitResults->SetMarkerColor( 1 );
    fitResults->SetMarkerStyle( 2 );
    fitResults->Draw( "A3_Y_True:A3_X_True", "", "same" );

    c1->SaveAs( "complexVars.eps" );
}
