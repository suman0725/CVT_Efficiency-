/* 
 * File:   DrawEffPlots.cc
 * Author: rafopar
 *
 * Created on February 13, 2023, 11:45 AM
 */

#include <cstdlib>
#include <string>

using namespace std;

/*
 * 
 */
void DrawEffPlots( int run, double IBeam, std::string UsrComment ) {

    gStyle->SetOptStat(0);

    //int run = 4899; // 80 nA
    //int run = 4900; // 65 nA
    //int run = 6616;
    //int run = 6769;
    //int run = 16290;
    //int run = 615250;
    //int run = 615275;

    const int nCVT_Layers = 6;
    const int nSVT_Layers = 6;
    const int nBMT_Layers = 6;
    const int nBMT_Sec = 3;

    std::map<int, double> I_bam_Run;
    I_bam_Run[4899] = 80.;
    I_bam_Run[4900] = 60.;
    I_bam_Run[6616] = 5.;
    I_bam_Run[6618] = 10.;
    I_bam_Run[6769] = 50.;
    I_bam_Run[16290] = 4.;
    I_bam_Run[61525] = 5.;
    I_bam_Run[615250] = 50.;
    I_bam_Run[615275] = 75.;
    
    std::map<int, std::string> m_comment_Run;
    m_comment_Run[4899] = "F18 RG-A";
    m_comment_Run[4900] = "F18 RG-A";
    m_comment_Run[6616] = "S19 RG-A";
    m_comment_Run[6618] = "S19 RG-A";
    m_comment_Run[6769] = "S19 RG-A";
    m_comment_Run[16290] = "RG-C 12C";
    m_comment_Run[61525] = "RG-B";
    m_comment_Run[615250] = "RG-B";
    m_comment_Run[615275] = "RG-B";


    const double z_MaxBMT[nBMT_Layers] = {  17.,  18.,  22.,  22. ,  22. ,  22.};
    const double z_MinBMT[nBMT_Layers] = { -13., -14., -15., -15.5, -16.5, -17 };
    const double phi_max[nBMT_Sec] = { 220., 100., 340.};
    const double phi_min[nBMT_Sec] = { 140., 20., 260. };
    
    
    std::string BMT_TTileNames[nCVT_Layers] = {"4C(1)", "4Z(2)", "5Z(3)", "5C(4)", "6Z(5)", "6C(6)"};

    TLine *line1 = new TLine();
    line1->SetLineColor(2);
    line1->SetLineWidth(3);

    TLatex *lat1 = new TLatex();
    lat1->SetNDC();

    TFile *file_in = new TFile(Form("CVT_Eff_%d.root", run));
    if( !file_in->IsOpen() ){
        cout<<Form("The file CVT_Eff_%d.root can not be opened.", run)<<endl;
        cout<<"Exiting"<<endl;
        exit(1);
    }

    std::map<int, int> m_canv_pad;
    m_canv_pad[0] = 1;
    m_canv_pad[1] = 4;
    m_canv_pad[2] = 5;
    m_canv_pad[3] = 2;
    m_canv_pad[4] = 6;
    m_canv_pad[5] = 3;

    double z_Gap_SVT[nSVT_Layers] = {-2.5, -2.3, 1.35, 1.65, 5.5, 5.5}; // Those indicate z coordinates of the gap between sub-sensors.
    int skip_zBins_SVT = 3;

    TH2D * h_z_phi_layers3_[nCVT_Layers];
    TH2D * h_z_phi_layers4_[nCVT_Layers];
    TH2D * h_z_phi_layers_Eff_[nCVT_Layers];
    TH1D * h_z_layers3_[nCVT_Layers];
    TH1D * h_z_layers4_[nCVT_Layers];
    TH1D * h_phi_layers3_[nCVT_Layers];
    TH1D * h_phi_layers4_[nCVT_Layers];
    TH1D * h_phi_Dep_layer_EffBMT_[nCVT_Layers];
    TH1D * h_z_Dep_layer_EffBMT_[nCVT_Layers];

    TH1D * h_BMT_Residual_[nCVT_Layers];
    TH1D * h_BMT_Residual_OnTrk_[nCVT_Layers];

    TCanvas *c1 = new TCanvas("c1", "", 1800, 1200);
    c1->Divide(3, 2, 0., 0.);
    TCanvas *c2 = new TCanvas("c2", "", 1800, 1200);
    c2->Divide(3, 2, 0., 0.);
    TCanvas *c3 = new TCanvas("c3", "", 1800, 1200);
    c3->Divide(3, 2, 0., 0.);
    TCanvas *c4 = new TCanvas("c4", "", 1800, 1200);
    c4->Divide(3, 2, 0.0001, 0.01);

    
    ofstream out_BMTTileEff(Form("BMT_Tile_Eff_%d.dat", run));
    for (int i = 0; i < nCVT_Layers; i++) {
        h_z_phi_layers3_[i] = (TH2D*) file_in->Get(Form("h_z_phi_layers3_%d", i));
        h_z_phi_layers3_[i]->SetTitle("; # phi [deg]; z [cm]");
        h_z_layers3_[i] = (TH1D*) h_z_phi_layers3_[i]->ProjectionY(Form("h_z_layers3_%d", i), 1, h_z_phi_layers3_[i]->GetNbinsX());
        
        int zBin_Max = h_z_phi_layers3_[i]->GetYaxis()->FindBin(z_MaxBMT[i]);
        int zBin_Min = h_z_phi_layers3_[i]->GetYaxis()->FindBin(z_MinBMT[i]);

        h_z_phi_layers4_[i] = (TH2D*) file_in->Get(Form("h_z_phi_layers4_%d", i));
        h_z_phi_layers4_[i]->SetTitle("; # phi [deg]; z [cm]");
        h_z_layers4_[i] = (TH1D*) h_z_phi_layers4_[i]->ProjectionY(Form("h_z_layers4_%d", i), 1, h_z_phi_layers4_[i]->GetNbinsX());

        double eff[nBMT_Sec];
        for( int iSec = 0; iSec < nBMT_Sec; iSec++ ){
            int phiBinMax = h_z_phi_layers3_[i]->GetXaxis()->FindBin(phi_max[iSec]) ;
            int phiBinMin = h_z_phi_layers3_[i]->GetXaxis()->FindBin(phi_min[iSec]) ;
            
            double counts_HasHit = h_z_phi_layers3_[i]->Integral(phiBinMin, phiBinMax, zBin_Min, zBin_Max);
            double counts_All = h_z_phi_layers4_[i]->Integral(phiBinMin, phiBinMax, zBin_Min, zBin_Max);
            eff[iSec] = counts_HasHit/counts_All;
            out_BMTTileEff<<i+1<<setw(5)<<iSec + 1<<setw(15)<<eff[iSec]<<endl;
        }
        int zMaxBin = i >= nCVT_Layers - 2 ? h_z_phi_layers4_[i]->GetYaxis()->FindBin(23.) : h_z_phi_layers4_[i]->GetNbinsY();

        h_phi_layers3_[i] = h_z_phi_layers3_[i]->ProjectionX(Form("h_phi_layers3__%d", i), 1, zMaxBin);
        h_phi_layers4_[i] = h_z_phi_layers4_[i]->ProjectionX(Form("h_phi_layers4__%d", i), 1, zMaxBin);

        h_phi_Dep_layer_EffBMT_[i] = (TH1D*) h_phi_layers3_[i]->Clone(Form("h_phi_Dep_layer_EffBMT_%d", i));
        h_phi_Dep_layer_EffBMT_[i]->Divide(h_phi_layers3_[i], h_phi_layers4_[i], 1, 1, "B");

        h_z_phi_layers_Eff_[i] = (TH2D*) h_z_phi_layers3_[i]->Clone(Form("h_z_phi_layers_Eff_%d", i));
        h_z_phi_layers_Eff_[i]->Divide(h_z_phi_layers4_[i]);

        c1->cd(i + 1);
        //h_z_phi_layers_Eff_[i]->SetMinimum(0.8);
        h_z_phi_layers_Eff_[i]->SetAxisRange(-25, 25, "Y");
        h_z_phi_layers_Eff_[i]->Draw("colz");

        c2->cd(i + 1)->SetGridy();
        h_z_Dep_layer_EffBMT_[i] = (TH1D*) h_z_layers3_[i]->Clone(Form("h_z_Dep_layer_EffBMT_%d", i));
        //h_z_Dep_layer_EffBMT_[i]->Divide(h_z_layers4_[i]);
        h_z_Dep_layer_EffBMT_[i]->Divide(h_z_layers3_[i], h_z_layers4_[i], 1, 1, "B");
        h_z_Dep_layer_EffBMT_[i]->SetAxisRange(0., 1.05, "Y");
        h_z_Dep_layer_EffBMT_[i]->GetYaxis()->SetNdivisions(511);
        h_z_Dep_layer_EffBMT_[i]->Draw("E");

        c3->cd(i + 1)->SetGridy();
        h_phi_Dep_layer_EffBMT_[i]->SetAxisRange(0., 1.05, "Y");
        h_phi_Dep_layer_EffBMT_[i]->GetYaxis()->SetNdivisions(511);
        h_phi_Dep_layer_EffBMT_[i]->SetLineWidth(3);
        h_phi_Dep_layer_EffBMT_[i]->Draw("E Hist");
        lat1->DrawLatex(0.15, 0.5, Form("%s", BMT_TTileNames[i].c_str()));
        lat1->DrawLatex(0.12, 0.3, "Sec. A(2)");
        lat1->DrawLatex(0.42, 0.3, "Sec. B(1)");
        lat1->DrawLatex(0.72, 0.3, "Sec. C(3)");
        for( int iSec = 0; iSec < nBMT_Sec; iSec++ ){
            line1->DrawLine(phi_min[iSec], eff[iSec], phi_max[iSec], eff[iSec] );
        }
        
        line1->DrawLine(120., 0., 120., 1.);
        line1->DrawLine(240., 0., 240., 1.);

        c4->cd(m_canv_pad[i])->SetLogy();
        h_BMT_Residual_[i] = (TH1D*) file_in->Get(Form("h_BMT_Residual_%d", i));
        std::string title = m_canv_pad[i] < 4 ? "; #Delta z [cm]" : "; #Delta #phi [deg]";
        h_BMT_Residual_[i]->SetTitle(title.c_str());
        h_BMT_Residual_[i]->SetLineColor(2);
        h_BMT_Residual_OnTrk_[i] = (TH1D*) file_in->Get(Form("h_BMT_Residual_OnTrk_%d", i));
        h_BMT_Residual_[i]->Draw();
        h_BMT_Residual_OnTrk_[i]->Draw("Same");
    }

    c1->cd();
    lat1->DrawLatex(0.35, 0.96, Form("Run %d, I = %1.1f nA  %s", run, IBeam, UsrComment.c_str()));
    c2->cd();
    lat1->DrawLatex(0.35, 0.96, Form("Run %d, I = %1.1f nA  %s", run, IBeam, UsrComment.c_str()));
    c3->cd();
    lat1->DrawLatex(0.35, 0.96, Form("Run %d, I = %1.1f nA  %s", run, IBeam, UsrComment.c_str()));

    c1->Print(Form("Figs/BMT_2D_HitEff_%d.pdf", run));
    c1->Print(Form("Figs/BMT_2D_HitEff_%d.png", run));
    c1->Print(Form("Figs/BMT_2D_HitEff_%d.root", run));

    c2->Print(Form("Figs/BMT_HitEff_zDep_%d.pdf", run));
    c2->Print(Form("Figs/BMT_HitEff_zDep_%d.png", run));
    c2->Print(Form("Figs/BMT_HitEff_zDep_%d.root", run));

    c3->Print(Form("Figs/BMT_HitEff_phi_%d.pdf", run));
    c3->Print(Form("Figs/BMT_HitEff_phi_%d.png", run));
    c3->Print(Form("Figs/BMT_HitEff_phi_%d.root", run));

    c4->Print(Form("Figs/BMT_Residuals_%d.pdf", run));
    c4->Print(Form("Figs/BMT_Residuals_%d.png", run));
    c4->Print(Form("Figs/BMT_Residuals_%d.root", run));



    TH2D * h_BST_z_phi_layers3_[nSVT_Layers];
    TH2D * h_BST_z_phi_layers4_[nSVT_Layers];
    TH2D * h_BST_z_phi_Eff_layers_[nSVT_Layers];
    TH1D * h_BST_phi_Eff_layers_[nSVT_Layers];

    TH1D * h_SVT_z_layers3_[nSVT_Layers];
    TH1D * h_SVT_z_layers4_[nSVT_Layers];

    TCanvas *c5 = new TCanvas("c5", "", 1800, 1200);
    c5->Divide(3, 2, 0., 0.);
    TCanvas *c6 = new TCanvas("c6", "", 1800, 1200);
    c6->Divide(3, 2, 0., 0.);

    for (int i = 0; i < nSVT_Layers; i++) {

        h_BST_z_phi_layers3_[i] = (TH2D*) file_in->Get(Form("h_BST_z_phi_layers3_%d", i));
        h_BST_z_phi_layers4_[i] = (TH2D*) file_in->Get(Form("h_BST_z_phi_layers4_%d", i));
        h_BST_z_phi_Eff_layers_[i] = (TH2D*) h_BST_z_phi_layers3_[i]->Clone(Form("h_BST_z_phi_Eff_layers_%d", i));

        h_BST_z_phi_Eff_layers_[i]->Divide(h_BST_z_phi_layers4_[i]);

        c5->cd(i + 1);

        h_BST_z_phi_Eff_layers_[i]->Draw("colz");

        int ybin1 = h_BST_z_phi_layers3_[i]->GetYaxis()->FindBin(z_Gap_SVT[i]) - skip_zBins_SVT;
        int ybin2 = h_BST_z_phi_layers3_[i]->GetYaxis()->FindBin(z_Gap_SVT[i]) + skip_zBins_SVT;
        TH1D *h_tmp1 = h_BST_z_phi_layers3_[i]->ProjectionX("h_tmp1", 1, ybin1);
        TH1D *h_tmp2 = h_BST_z_phi_layers3_[i]->ProjectionX("h_tmp2", ybin2, h_BST_z_phi_layers3_[i]->GetNbinsY() );
        
        h_SVT_z_layers3_[i] = (TH1D*)h_tmp1->Clone(Form("h_SVT_z_layers3_%d", i));
        h_SVT_z_layers3_[i]->Add(h_tmp2);
        
        delete h_tmp1;
        delete h_tmp2;
        
        h_tmp1 = h_BST_z_phi_layers4_[i]->ProjectionX("h_tmp1", 1, ybin1);
        h_tmp2 = h_BST_z_phi_layers4_[i]->ProjectionX("h_tmp2", ybin2, h_BST_z_phi_layers4_[i]->GetNbinsY() );
        
        h_SVT_z_layers4_[i] = (TH1D*)h_tmp1->Clone(Form("h_SVT_z_layers4_%d", i));
        h_SVT_z_layers4_[i]->Add(h_tmp2);
        
        delete h_tmp1;
        delete h_tmp2;
        
        h_BST_phi_Eff_layers_[i] = (TH1D*)h_SVT_z_layers3_[i]->Clone(Form("h_BST_phi_Eff_layers_%d", i));
        h_BST_phi_Eff_layers_[i]->Divide(h_SVT_z_layers3_[i], h_SVT_z_layers4_[i], 1, 1, "B");
        
        c6->cd(i+1)->SetGridy();
        h_BST_phi_Eff_layers_[i]->SetAxisRange(0., 1.05, "Y");
        h_BST_phi_Eff_layers_[i]->GetYaxis()->SetNdivisions(511);
        h_BST_phi_Eff_layers_[i]->SetTitle("; #phi [deg]; Efficiency");
        h_BST_phi_Eff_layers_[i]->SetLineWidth(3);
        h_BST_phi_Eff_layers_[i]->Draw("E hist");
        lat1->DrawLatex(0.45, 0.15, Form("Layer %d", i+1));
    }

    c5->cd();
    lat1->DrawLatex(0.35, 0.96, Form("Run %d, I = %1.1f nA   %s", run, IBeam, UsrComment.c_str()));
    c6->cd();
    lat1->DrawLatex(0.35, 0.96, Form("Run %d, I = %1.1f nA   %s", run, IBeam, UsrComment.c_str()));
    
    c5->Print(Form("Figs/SVT_2D_HitEff_%d.pdf", run));
    c5->Print(Form("Figs/SVT_2D_HitEff_%d.png", run));
    c5->Print(Form("Figs/SVT_2D_HitEff_%d.root", run));
    
    c6->Print(Form("Figs/SVT_HitEff_phi_%d.pdf", run));
    c6->Print(Form("Figs/SVT_HitEff_phi_%d.png", run));
    c6->Print(Form("Figs/SVT_HitEff_phi_%d.root", run));

    return;
}
