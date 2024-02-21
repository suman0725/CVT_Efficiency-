/* 
 * File:   CVT_DetEff.cc
 * Author: rafopar
 *
 * Created on February 1, 2023, 7:29 PM
 */

#include <cstdlib>
#include <iostream>

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TMath.h>

// ===== Hipo headers =====
#include <reader.h>
#include <writer.h>
#include <dictionary.h>
#include <vector>

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

    std::cout << " reading file example program (HIPO) " << __cplusplus << std::endl;

    char inputFile[256];
    int run = -9999;

    if (argc > 1) {
        run = atoi(argv[1]);
        sprintf(inputFile, "Data/rec_clas_%d.hipo", run);
    } else {
        std::cout << " *** please provide the run number..." << std::endl;
        cout << "Exiting ..." << endl;
        exit(0);
    }


    const int detBMT = 1;
    const int detCVT = 5;
    const int detCTOF = 4;
    const double zCTOF_max = 25.;
    const double zCTOF_min = -22.;
    const int nBMTLayers = 6;
    const int nBSTLayers = 6;
    const int nBMTSec = 3;

    const double trk_dist_cuts_[nBMTLayers] = {0.5, 2., 2., 0.5, 2., 0.5};

    const double chi2NDFCut = 5.;

    TFile *file_out = new TFile(Form("CVT_Eff_%d.root", run), "Recreate");

    TH1D h_nCVTTr("h_nCVTTr", "", 11, -0.5, 10.5);
    TH1D h_nCVTTraj("h_nCVTTraj", "", 11, -0.5, 10.5);

    TH1D h_Chi2NDF1("h_Chi2NDF1", "", 200, 0., 40.);

    TH1D h_trkID_1("h_trkID_1", "", 11, -0.5, 10.5);
    TH1D h_trkID_2("h_trkID_2", "", 11, -0.5, 10.5);

    TH1D h_nBSTCl1("h_nBSTCl1", "", 21, -0.5, 20.5);
    TH2D h_BST_Layer_phi1("h_BST_Layer_phi1", "", 200, 0., 360., 7, -0.5, 6.5);
    TH2D h_BST_Layer_phi2("h_BST_Layer_phi2", "", 200, 0., 360., 7, -0.5, 6.5);
    TH2D h_BST_Layer_phi3("h_BST_Layer_phi3", "", 200, 0., 360., 7, -0.5, 6.5);
    TH2D h_BST_Layer_z1("h_BST_Layer_z1", "", 200, -30, 20., 7, -0.5, 6.5);
    TH2D h_BST_Layer_z2("h_BST_Layer_z2", "", 200, -30, 20., 7, -0.5, 6.5);
    TH1D h_nBSTClOnTrk1("h_nBSTClOnTrk1", "", 21, -0.5, 20.5);
    TH1D h_BSTLayers1("h_BSTLayers1", "", 65, -0.5, 64.5);
    TH2D h_nBMT_Crs_Cl1("h_nBMT_Crs_Cl1", "", 46, -0.5, 45.5, 46, -0.5, 45.5);
    TH1D h_nBMTClOnTrk1("h_nBMTClOnTrk1", "", 21, -0.5, 20.5);
    TH1D h_BMTLayers1("h_BMTLayers1", "", 65, -0.5, 64.5);
    TH1D h_BMTLayers2("h_BMTLayers2", "", 65, -0.5, 64.5);
    TH2D h_nBMTBST_Cls1("h_nBMTBST_Cls1", "", 11, -0.5, 10.5, 11, -0.5, 10.5);
    TH2D h_BMTvsBST_layers("h_BMTvsBST_layers", "", 65, -0.5, 64.5, 65, -0.5, 64.5);

    TH2D h_BMT_z2_z1("h_BMT_z2_z1", "", 200, -30., 30., 200, -30., 30.);
    TH2D h_BMT_z2_layer("h_BMT_z2_layer", "", 200, -30., 30., 7, -0.5, 6.5);
    TH2D h_BMT_z1_layer("h_BMT_z1_layer", "", 200, -30., 30., 7, -0.5, 6.5);

    TH2D h_BMT_yx2("h_BMT_yx2", "", 200, -25., 25., 200, -25., 25);
    TH2D h_BMT_yx1("h_BMT_yx1", "", 200, -25., 25., 200, -25., 25);
    TH2D h_ph2_phi1("h_ph2_phi1", "", 200, 0., 360, 200, 0., 360);

    TH1D h_z_CTOF1("h_z_CTOF1", "", 200, -30., 50.);
    TH2D h_N_BMT_Cl_AllVsOnTrk1_[3][nBMTLayers];

    for (int isec = 0; isec < 3; isec++) {

        for (int il = 0; il < nBMTLayers; il++) {
            h_N_BMT_Cl_AllVsOnTrk1_[isec][il] = TH2D(Form("h_N_BMT_Cl_AllVsOnTrk1_%d_%d", isec, il), "", 5, -0.5, 4.5, 5, -0.5, 4.5);
        }
    }

    TH1D h_phi_layers1_[nBMTLayers];
    TH1D h_phi_layers2_[nBMTLayers];
    TH2D h_z_phi_layers1_[nBMTLayers];
    TH2D h_z_phi_layers2_[nBMTLayers];
    TH2D h_z_phi_layers3_[nBMTLayers];
    TH2D h_z_phi_layers4_[nBMTLayers];
    TH1D h_BMT_Residual_[nBMTLayers];
    TH1D h_BMT_Residual_OnTrk_[nBMTLayers];
    TH2D h_BST_z_phi_layers3_[nBSTLayers];
    TH2D h_BST_z_phi_layers4_[nBSTLayers];

    for (int il = 0; il < nBMTLayers; il++) {
        h_phi_layers1_[il] = TH1D(Form("h_phi_layers1_%d", il), "", 60, 0., 360.);
        h_phi_layers2_[il] = TH1D(Form("h_phi_layers2_%d", il), "", 60, 0., 360.);

        h_z_phi_layers1_[il] = TH2D(Form("h_z_phi_layers1_%d", il), "", 30, 0., 360., 30, -40., 40.);
        h_z_phi_layers2_[il] = TH2D(Form("h_z_phi_layers2_%d", il), "", 30, 0., 360., 30, -40., 40.);
        h_z_phi_layers3_[il] = TH2D(Form("h_z_phi_layers3_%d", il), "", 30, 0., 360., 30, -40., 40.);
        h_z_phi_layers4_[il] = TH2D(Form("h_z_phi_layers4_%d", il), "", 30, 0., 360., 30, -40., 40.);
        h_BMT_Residual_[il] = TH1D(Form("h_BMT_Residual_%d", il), "", 2000, -50., 50.);
        h_BMT_Residual_OnTrk_[il] = TH1D(Form("h_BMT_Residual_OnTrk_%d", il), "", 2000, -50., 50.);
    }

    for (int il = 0; il < nBSTLayers; il++) {
        h_BST_z_phi_layers3_[il] = TH2D(Form("h_BST_z_phi_layers3_%d", il), "", 100, 0., 360., 100, -15., 15.);
        h_BST_z_phi_layers4_[il] = TH2D(Form("h_BST_z_phi_layers4_%d", il), "", 100, 0., 360., 100, -15., 15.);
    }

    hipo::reader reader;
    reader.open(inputFile);

    hipo::dictionary factory;

    reader.readDictionary(factory);

    factory.show();

    hipo::event event;
    int evCounter = 0;

    hipo::bank bCVTRecTr(factory.getSchema("CVTRec::Tracks"));
    hipo::bank bCVTRecTraj(factory.getSchema("CVTRec::Trajectory"));
    hipo::bank bBMTRecCrosses(factory.getSchema("BMTRec::Crosses"));
    hipo::bank bBMTRecClusters(factory.getSchema("BMTRec::Clusters"));
    hipo::bank bBSTRecCrosses(factory.getSchema("BSTRec::Crosses"));
    hipo::bank bBSTRecClusters(factory.getSchema("BSTRec::Clusters"));

    const int nMaxEv = 100000;

    try {

        while (reader.next() == true) {
            reader.read(event);

            evCounter = evCounter + 1;

            if (evCounter % 10000 == 0) {
                cout.flush() << "Processed " << evCounter << " events \r";
            }

//            if (evCounter > nMaxEv) {
//                break;
//            }

            event.getStructure(bCVTRecTr);
            event.getStructure(bBMTRecCrosses);
            event.getStructure(bBMTRecClusters);
            event.getStructure(bBSTRecCrosses);
            event.getStructure(bBSTRecClusters);
            event.getStructure(bCVTRecTraj);

            int nCVTTr = bCVTRecTr.getRows();
            h_nCVTTr.Fill(nCVTTr);


            double chi2NDF = 9999.;
            int trk_ID = -9999;
            for (int itr = 0; itr < nCVTTr; itr++) {
                trk_ID = bCVTRecTr.getInt("ID", itr);
                double chi2 = double(bCVTRecTr.getFloat("chi2", itr));
                int ndf = bCVTRecTr.getInt("ndf", itr);
                chi2NDF = chi2 / double(ndf);

                h_Chi2NDF1.Fill(chi2NDF);

                h_trkID_1.Fill(trk_ID);
            }

            // Use only clean events: when there is only 1 track in the event
            if (nCVTTr != 1) {
                continue;
            }

            h_trkID_2.Fill(trk_ID);

            if (chi2NDF > chi2NDFCut) {
                continue;
            }


            /*
            Loop over BST clusters and make sure there are 6 BMT clusters corresponding to the track
             */
            int nBSTCl = bBSTRecClusters.getRows();
            h_nBSTCl1.Fill(nBSTCl);

            int nBSTClInTrk1 = 0;
            int BST_Cl_Layers = 0;
            int bst_sec[3] = {0};
            for (int iBSTCl = 0; iBSTCl < nBSTCl; iBSTCl++) {
                int trkID = bBSTRecClusters.getInt("trkID", iBSTCl);
                int layer = bBSTRecClusters.getInt("layer", iBSTCl);
                int region = (layer - 1) / 2; // 0, 1 and 2
                int sec = bBSTRecClusters.getInt("sector", iBSTCl);
                float xCl1 = bBSTRecClusters.getFloat("x1", iBSTCl);
                float yCl1 = bBSTRecClusters.getFloat("y1", iBSTCl);
                float zCl1 = bBSTRecClusters.getFloat("z1", iBSTCl);
                float phi_cl1 = atan2(yCl1, xCl1) * TMath::RadToDeg() - 30.;
                phi_cl1 = phi_cl1 < 0. ? phi_cl1 + 360 : phi_cl1;
                h_BST_Layer_phi1.Fill(phi_cl1, layer);
                float xCl2 = bBSTRecClusters.getFloat("x2", iBSTCl);
                float yCl2 = bBSTRecClusters.getFloat("y2", iBSTCl);
                float zCl2 = bBSTRecClusters.getFloat("z2", iBSTCl);
                float phi_cl2 = atan2(yCl2, xCl2) * TMath::RadToDeg() - 30.;
                phi_cl2 = phi_cl2 < 0. ? phi_cl2 + 360 : phi_cl2;
                h_BST_Layer_phi2.Fill(phi_cl2, layer);

                h_BST_Layer_z1.Fill(zCl1, layer);
                h_BST_Layer_z2.Fill(zCl2, layer);

                float xClAvg = (xCl1 + xCl2) / 2.;
                float yClAvg = (yCl1 + yCl2) / 2.;

                float phi_clavg = atan2(yClAvg, xClAvg) * TMath::RadToDeg() - 30.;
                phi_clavg = phi_clavg < 0. ? phi_clavg + 360 : phi_clavg;

                h_BST_Layer_phi3.Fill(phi_clavg, layer);

                if (layer < 1) {
                    cout << "Wrong Layer!!!!!   " << layer << endl;
                    exit(1);
                }

                if (trkID == 1) {
                    nBSTClInTrk1 = nBSTClInTrk1 + 1;
                    BST_Cl_Layers = BST_Cl_Layers | 1 << (layer - 1);
                    bst_sec[region] = sec;
                }
            }

            h_BSTLayers1.Fill(BST_Cl_Layers);
            h_nBSTClOnTrk1.Fill(nBSTClInTrk1);

            int nBMTCross = bBMTRecCrosses.getRows();
            int nBMTCl = bBMTRecClusters.getRows();
            int nCVTRecTraj = bCVTRecTraj.getRows();

            h_nCVTTraj.Fill(nCVTRecTraj);

            h_nBMT_Crs_Cl1.Fill(nBMTCl, nBMTCross);

            int sector = -1;
            int nBMTClInTrk1 = 0;
            int BMT_Cl_Layers_OnTrk = 0;
            int BMT_Cl_Layers_CloseToTrk = 0;
            int nBMTClAll_Layers[3][nBMTLayers] = {0};
            int nBMTClInTrk_inLayers[3][nBMTLayers] = {0};

            vector<double> BMT_Cl_Coord_[nBMTLayers]; // For C layers, the coordinate is the "z coordinate" of the strip, and for Z layers the coordinate is the phi angle of the strip
            vector<double> BMT_Cl_CoordOnTrk_[nBMTLayers]; // For C layers, the coordinate is the "z coordinate" of the strip, and for Z layers the coordinate is the phi angle of the strip

            for (int iCl = 0; iCl < nBMTCl; iCl++) {
                int ID = bBMTRecClusters.getInt("ID", iCl);
                int trkID = bBMTRecClusters.getInt("trkID", iCl);
                int layer = bBMTRecClusters.getInt("layer", iCl);

                float x1 = bBMTRecClusters.getFloat("x1", iCl);
                float x2 = bBMTRecClusters.getFloat("x2", iCl);
                float y1 = bBMTRecClusters.getFloat("y1", iCl);
                float y2 = bBMTRecClusters.getFloat("y2", iCl);
                float z1 = bBMTRecClusters.getFloat("z1", iCl);
                float z2 = bBMTRecClusters.getFloat("z2", iCl);

                double phi1 = atan2(y1, x1) * TMath::RadToDeg() - 30.;
                phi1 = phi1 < 0 ? phi1 + 360. : phi1;

                double phi2 = atan2(y2, x2) * TMath::RadToDeg() - 30.;
                phi2 = phi2 < 0 ? phi2 + 360. : phi2;

                double strip_coord = (layer == 1 || layer == 4 || layer == 6) ? (z1 + z2) / 2. : (phi1 + phi2) / 2.;
                BMT_Cl_Coord_[layer - 1].push_back(strip_coord);

                h_ph2_phi1.Fill(phi1, phi2);

                h_BMT_yx1.Fill(x1, y1);
                h_BMT_yx2.Fill(x2, y2);
                h_BMT_z2_z1.Fill(z1, z2);
                h_BMT_z2_layer.Fill(z2, layer);
                h_BMT_z1_layer.Fill(z1, layer);

                int inter_sector = bBMTRecClusters.getInt("sector", iCl);

                nBMTClAll_Layers[inter_sector - 1][layer - 1] = nBMTClAll_Layers[inter_sector - 1][layer - 1] + 1;

                if (trkID == 1) {
                    sector = inter_sector;
                    nBMTClInTrk1 = nBMTClInTrk1 + 1;
                    BMT_Cl_Layers_OnTrk = BMT_Cl_Layers_OnTrk | 1 << (layer - 1);
                    nBMTClInTrk_inLayers[sector - 1][layer - 1] = nBMTClInTrk_inLayers[sector - 1][layer - 1] + 1;
                    BMT_Cl_CoordOnTrk_[layer - 1].push_back(strip_coord);
                }
            }



            h_nBMTClOnTrk1.Fill(nBMTClInTrk1);
            h_BMTLayers1.Fill(BMT_Cl_Layers_OnTrk);
            h_nBMTBST_Cls1.Fill(nBSTClInTrk1, nBMTClInTrk1);


            h_BMTvsBST_layers.Fill(BST_Cl_Layers, BMT_Cl_Layers_OnTrk);

            if (nBMTClInTrk1 < 4) {
                continue;
            }

            float zCTOF = -999.;
            float z_layers_[nBMTLayers] = {-999.};
            float phi_layers_[nBMTLayers] = {-999.};

            float z_layers_BST_[nBMTLayers] = {-999.};
            float phi_layers_BST_[nBMTLayers] = {-999.};

            for (int iTraj = 0; iTraj < nCVTRecTraj; iTraj++) {
                int det = bCVTRecTraj.getInt("detector", iTraj);
                int trkId = bCVTRecTraj.getInt("id", iTraj);
                int layer_SVT = bCVTRecTraj.getInt("layer", iTraj) - 1; // layers 1 to 6 correspond to SVT, and 7 to 12 correspond to BMT

                int layer_BMT = layer_SVT - 6;

                float z = bCVTRecTraj.getFloat("z", iTraj);
                float y = bCVTRecTraj.getFloat("y", iTraj);
                float x = bCVTRecTraj.getFloat("x", iTraj);

                float phi = atan2(y, x) * TMath::RadToDeg() - 30.;
                phi = phi < 0. ? phi + 360. : phi;

                if (det == detCTOF && trkId == 1) {
                    h_z_CTOF1.Fill(z);
                    zCTOF = z;
                } else if (det == detCVT && trkId == 1 && layer_BMT >= 0) {
                    z_layers_[layer_BMT] = z;
                    phi_layers_[layer_BMT] = phi;
                    h_phi_layers1_[layer_BMT].Fill(phi);

                    double resid_coord = (layer_BMT == 0 || layer_BMT == 3 || layer_BMT == 5) ? z : phi;

                    for (int icl = 0; icl < BMT_Cl_Coord_[layer_BMT].size(); icl++) {
                        h_BMT_Residual_[layer_BMT].Fill(BMT_Cl_Coord_[layer_BMT].at(icl) - resid_coord);

                        if (TMath::Abs(BMT_Cl_Coord_[layer_BMT].at(icl) - resid_coord) < trk_dist_cuts_[layer_BMT]) {

                            BMT_Cl_Layers_CloseToTrk = BMT_Cl_Layers_CloseToTrk | 1 << layer_BMT;
                        }
                    }

                    for (int icl = 0; icl < BMT_Cl_CoordOnTrk_[layer_BMT].size(); icl++) {

                        h_BMT_Residual_OnTrk_[layer_BMT].Fill(BMT_Cl_CoordOnTrk_[layer_BMT].at(icl) - resid_coord);
                    }

                } else if (det == detCVT && trkId == 1 && layer_SVT < nBSTLayers) {
                    z_layers_BST_[layer_SVT] = z;
                    phi_layers_BST_[layer_SVT] = phi;
                }
            }

            for (int il = 0; il < nBMTLayers; il++) {
                h_N_BMT_Cl_AllVsOnTrk1_[sector - 1][il].Fill(nBMTClAll_Layers[sector - 1][il], nBMTClInTrk_inLayers[sector - 1][il]);
                h_phi_layers2_[il].Fill(phi_layers_[il]);
                h_z_phi_layers2_[il].Fill(phi_layers_[il], z_layers_[il]);
            }



            /**
             *                        BMT Efficiency
             */
            if (BST_Cl_Layers == 63 && nBSTClInTrk1 == 6 && nBMTClInTrk1 >= 4) {
                h_BMTLayers2.Fill(BMT_Cl_Layers_OnTrk);

                for (int il = 0; il < nBMTLayers; il++) {
                    //if (BMT_Cl_Layers_OnTrk & (1 << il)) {
                    if (BMT_Cl_Layers_CloseToTrk & (1 << il)) {
                        h_z_phi_layers3_[il].Fill(phi_layers_[il], z_layers_[il]);
                    }

                    h_z_phi_layers4_[il].Fill(phi_layers_[il], z_layers_[il]);
                }
            }

            /**
             *                        BST Efficiency
             */

            if ( /*BMT_Cl_Layers_OnTrk == 63 &&*/ nBMTClInTrk1 >= 4 && ((BST_Cl_Layers | 3 << 0 == 63) || (BST_Cl_Layers | 3 << 2 == 63) || (BST_Cl_Layers | 3 << 4 == 63))) {

                for (int il = 0; il < nBSTLayers; il++) {
                    if (BST_Cl_Layers & 1 << il) {
                        h_BST_z_phi_layers3_[il].Fill(phi_layers_BST_[il], z_layers_BST_[il]);
                    }
                    h_BST_z_phi_layers4_[il].Fill(phi_layers_BST_[il], z_layers_BST_[il]);
                }

            }

        }
    } catch (const char *msg) {
        cerr << msg << endl;
    }

    gDirectory->Write();

    return 0;
}