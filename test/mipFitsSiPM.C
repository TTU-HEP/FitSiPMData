//#include "TROOT.h"
//#include "TObject.h"
//#include "TH1.h"
//#include "TCanvas.h"
//#include "TF1.h"
//#include "TLegend.h"
//#include "TMath.h"
//#include "TFile.h"
//#include "TDirectory.h"
//#include "TList.h"
//#include "TLatex.h"
#include "TChain.h"
//#include <cmath>
//#include <cstdio>
//#include <vector>
#include "SPEfunc.h"
#include "../include/NTupleReader.h"

struct TreeVars {
    float ped;
    float meanPE;
    float gain;
    float ctProb;
    std::string cName;
};

//This function fits the SiPM MIP distribution
void fitSPEMIP(TH1* hfit, int run, TreeVars& vars, TTree* tree)
{
    TCanvas c1("c1","c1",800,800);
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.12);
    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.14);
    //c1.SetLogx();
    c1.SetLogy();
    //hfit->GetXaxis()->SetRangeUser(10, 400);
    hfit->GetXaxis()->SetRangeUser(1, 1500);
    hfit->SetMinimum(0.1);
    hfit->SetMaximum(70000);
    
    std::string ffname;
    double nPeaks = 6;
    
    PPEFunc background1(3,false);
    background1.h = hfit;    
    PPEFunc background(nPeaks,false);
    background.h = hfit;
    PPEFunc background_MIP(nPeaks,true);
    background_MIP.h = hfit;
    hfit->Draw("E hist");


                                                       ////////////////////////////////////////////////////////////
    double min0  =  4500; double max0  =    5300;      //p[0]  : pedestal amplitude-------------//               //
    double min1  =   2.0; double max1  =    10.0;      //p[1]  : pedestal width-----------------//               //
    double min2  =     0; double max2  =     200;      //p[2]  : Overall Shift------------------//(avg. pedestal)//
    double min3  =     0; double max3  =    2000;      //p[3]  : background amplitude-----------//               //
    double min4  =     0; double max4  =      10;      //p[4]  : background width---------------//               //

    double min5  =  3000; double max5  =    6000;      //p[5]  : overall PE amplitude ----------//               //
    double min6  =  0.05; double max6  =     0.3;      //p[6]  : poisson mean number of PE------//               //
    double min7  =   140; double max7  =     150;      //p[7]  : PE peak spacing----------------//(Gain)         //
    double min8  =    12; double max8  =      18;      //p[8]  : PE peak width------------------//               //
    double min9  =  0.02; double max9  =    0.06;      //p[9]  : pixel cross-talk probability---//               //

    double min10 =     0; double max10 =       1;      //p[10] : Total MIP peak area------------//               //
    double min11 =     0; double max11 =       1;      //p[11] : Most probable value of Landau--//(MPV)          //
    double min12 =     1; double max12 =     150;      //p[12] : Landau width parameter---------//               //
    double min13 =     1; double max13 =     200;      //p[13] : gaussian width-----------------//               //
                                                       ////////////////////////////////////////////////////////////

    double set0  =  4900;
    double set1  =     8;
    double set2  =   140;
    double set3  =  2000;
    double set4  =     5;

    double set5  =  5000;
    double set6  =  0.15;
    double set7  =   140;
    double set8  =  15.0;
    double set9  =  0.05;

    double set10 =     0;
    double set11 =     0;
    double set12 =    50;
    double set13 =   100;
    

    //Fit Background Peak
    ffname = "sff" + std::string(hfit->GetName());
    TF1* fit0 = new TF1(ffname.c_str(), background1, 0.0, 1000.0, 10);
    fit0->SetParLimits(0,  min0, max0); 
    fit0->SetParLimits(1,  min1, max1); 
    fit0->SetParLimits(2,  min2, max2); 
    fit0->SetParLimits(3,  min3, max3); 
    fit0->SetParLimits(4,  min4, max4); 
    fit0->SetParLimits(5,  min5, max5); 
    fit0->SetParLimits(6,  min6, max6); 
    fit0->SetParLimits(7,  min7, max7); 
    fit0->SetParLimits(8,  min8, max8); 
    fit0->SetParLimits(9,  min9, max9); 
    fit0->FixParameter(0,  0.01);
    fit0->FixParameter(1,  set1);
    fit0->FixParameter(2,  set2);
    fit0->SetParameter(3,  set3);
    fit0->SetParameter(4,  set4);
    fit0->FixParameter(5,  0.01);
    fit0->FixParameter(6,  set6);
    fit0->FixParameter(7,  set7);
    fit0->FixParameter(8,  set8);
    fit0->FixParameter(9,  set9);
    hfit->Fit(fit0, "RQML", "", 0, 1000);
    fit0->SetLineColor(kGreen);
    fit0->SetLineWidth(2);
    // fit0->Draw("same");

    //Fit Pedestal Peak
    ffname = "sff" + std::string(hfit->GetName());
    TF1* fit1 = new TF1(ffname.c_str(), background1, 130.0, 160.0, 10);
    fit1->SetParLimits(0,  min0, max0); 
    fit1->SetParLimits(1,  min1, max1); 
    fit1->SetParLimits(2,  min2, max2); 
    fit1->SetParLimits(3,  min3, max3); 
    fit1->SetParLimits(4,  min4, max4); 
    fit1->SetParLimits(5,  min5, max5); 
    fit1->SetParLimits(6,  min6, max6); 
    fit1->SetParLimits(7,  min7, max7); 
    fit1->SetParLimits(8,  min8, max8); 
    fit1->SetParLimits(9,  min9, max9); 
    fit1->SetParameter(0,  set0);
    fit1->SetParameter(1,  set1);
    fit1->SetParameter(2,  set2);
    fit1->FixParameter(3,  fit0->GetParameter(3));
    fit1->FixParameter(4,  fit0->GetParameter(4));
    fit1->FixParameter(5,  0.01);
    fit1->FixParameter(6,  set6);
    fit1->FixParameter(7,  set7);
    fit1->FixParameter(8,  set8);
    fit1->FixParameter(9,  set9);
    hfit->Fit(fit1, "RQML", "", 130, 160);
    fit1->SetLineColor(kRed);
    fit1->SetLineWidth(2);
    // fit1->Draw("same");
    
    //fit PE Peaks 
    ffname = "sff2" + std::string(hfit->GetName());
    TF1* fit2 = new TF1(ffname.c_str(), background1, 140.0, 500.0, 10);
    fit2->SetParLimits(0,  min0, max0); 
    fit2->SetParLimits(1,  min1, max1); 
    fit2->SetParLimits(2,  min2, max2); 
    fit2->SetParLimits(3,  min3, max3); 
    fit2->SetParLimits(4,  min4, max4); 
    fit2->SetParLimits(5,  min5, max5); 
    fit2->SetParLimits(6,  min6, max6); 
    fit2->SetParLimits(7,  min7, max7); 
    fit2->SetParLimits(8,  min8, max8); 
    fit2->SetParLimits(9,  min9, max9); 
    fit2->FixParameter(0,  fit1->GetParameter(0));
    fit2->FixParameter(1,  fit1->GetParameter(1));
    fit2->FixParameter(2,  fit1->GetParameter(2));
    fit2->SetParameter(3,  fit1->GetParameter(3));
    fit2->SetParameter(4,  fit1->GetParameter(4));
    fit2->SetParameter(5,  set5);
    fit2->SetParameter(6,  set6);
    fit2->SetParameter(7,  set7);
    fit2->SetParameter(8,  set8);
    fit2->SetParameter(9,  set9);
    hfit->Fit(fit2, "RQML", "", fit1->GetParameter(2), 450);
    fit2->SetLineColor(kBlack);
    fit2->SetLineWidth(2);
    // fit2->Draw("same");

    //Fine tune PE peaks
    ffname = "sff3" + std::string(hfit->GetName());
    TF1* fit3 = new TF1(ffname.c_str(), background, 300.0, 1000.0, 10);
    double down1 = 0.5;
    double up1   = 1.5;
    fit3->SetParLimits(0,  down1*fit2->GetParameter(0), up1*fit2->GetParameter(0)); 
    fit3->SetParLimits(1,  down1*fit2->GetParameter(1), up1*fit2->GetParameter(1)); 
    fit3->SetParLimits(2,  down1*fit2->GetParameter(2), up1*fit2->GetParameter(2)); 
    fit3->SetParLimits(3,  down1*fit2->GetParameter(3),     fit2->GetParameter(3)); 
    fit3->SetParLimits(4,  down1*fit2->GetParameter(4),     fit2->GetParameter(4)); 
    fit3->SetParLimits(5,  down1*fit2->GetParameter(5), up1*fit2->GetParameter(5)); 
    fit3->SetParLimits(6,  down1*fit2->GetParameter(6), up1*fit2->GetParameter(6)); 
    fit3->SetParLimits(7,   0.85*fit2->GetParameter(7), up1*fit2->GetParameter(7)); 
    fit3->SetParLimits(8,  down1*fit2->GetParameter(8), up1*fit2->GetParameter(8)); 
    fit3->SetParLimits(9,  down1*fit2->GetParameter(9), up1*fit2->GetParameter(9)); 
    fit3->FixParameter(0,  fit2->GetParameter(0));
    fit3->FixParameter(1,  fit2->GetParameter(1));
    fit3->FixParameter(2,  fit2->GetParameter(2));
    fit3->FixParameter(3,  fit2->GetParameter(3));
    fit3->FixParameter(4,  fit2->GetParameter(4));
    fit3->FixParameter(5,  fit2->GetParameter(5));
    fit3->SetParameter(6,  fit2->GetParameter(6));
    fit3->SetParameter(7,  fit2->GetParameter(7));
    fit3->SetParameter(8,  fit2->GetParameter(8));
    fit3->SetParameter(9,  fit2->GetParameter(9));
    hfit->Fit(fit3, "RQML", "", 450, 1000);
    fit3->SetLineColor(kOrange);
    fit3->SetLineWidth(2);
    // fit3->Draw("same");
    
    //Fit MIP peak
    ffname = "sff4" + std::string(hfit->GetName());
    TF1* fit4 = new TF1(ffname.c_str(), background_MIP, 130.0, 1000.0, 14);
    fit4->SetParLimits(0,  min0,  max0); 
    fit4->SetParLimits(1,  min1,  max1); 
    fit4->SetParLimits(2,  min2,  max2); 
    fit4->SetParLimits(3,  min3,  max3); 
    fit4->SetParLimits(4,  min4,  max4); 
    fit4->SetParLimits(5,  min5,  max5); 
    fit4->SetParLimits(6,  min6,  max6); 
    fit4->SetParLimits(7,  min7,  max7); 
    fit4->SetParLimits(8,  min8,  max8); 
    fit4->SetParLimits(9,  min9,  max9);
    fit4->SetParLimits(10, min10, max10); 
    fit4->SetParLimits(11, min11, max11); 
    fit4->SetParLimits(12, min12, max12); 
    fit4->SetParLimits(13, min13, max13); 
    fit4->FixParameter(0,  fit3->GetParameter(0));
    fit4->FixParameter(1,  fit3->GetParameter(1));
    fit4->FixParameter(2,  fit3->GetParameter(2));
    fit4->FixParameter(3,  fit3->GetParameter(3));
    fit4->FixParameter(4,  fit3->GetParameter(4));
    fit4->FixParameter(5,  fit3->GetParameter(5));
    fit4->FixParameter(6,  fit3->GetParameter(6));
    fit4->FixParameter(7,  fit3->GetParameter(7));
    fit4->FixParameter(8,  fit3->GetParameter(8));
    fit4->FixParameter(9,  fit3->GetParameter(9));
    fit4->FixParameter(10, set10);
    fit4->FixParameter(11, set11);
    fit4->FixParameter(12, set12);
    fit4->FixParameter(13, set13);
    hfit->Fit(fit4, "RQML", "", fit3->GetParameter(2), 1000);
    fit4->SetLineColor(kOrange);
    fit4->SetLineWidth(2);
    //fit4->Draw("same");
    
    //Fine Tunning All Fit
    ffname = "sff5" + std::string(hfit->GetName());
    TF1* fit5 = new TF1(ffname.c_str(), background_MIP, 350.0, 1000.0, 14);
    double down = 0.95;
    double up   = 1.05;
    fit5->SetParLimits(0, down*fit4->GetParameter(0),  up*fit4->GetParameter(0));
    fit5->SetParLimits(1, down*fit4->GetParameter(1),  up*fit4->GetParameter(1));
    fit5->SetParLimits(2, down*fit4->GetParameter(2),  up*fit4->GetParameter(2));
    fit5->SetParLimits(3, down*fit4->GetParameter(3),     fit4->GetParameter(3));
    fit5->SetParLimits(4, down*fit4->GetParameter(4),     fit4->GetParameter(4));
    fit5->SetParLimits(5, down*fit4->GetParameter(5),  up*fit4->GetParameter(5));
    fit5->SetParLimits(6, down*fit4->GetParameter(6),  up*fit4->GetParameter(6));
    fit5->SetParLimits(7, down*fit4->GetParameter(7),  up*fit4->GetParameter(7));
    fit5->SetParLimits(8, down*fit4->GetParameter(8),     fit4->GetParameter(8));
    fit5->SetParLimits(9, down*fit4->GetParameter(9),  up*fit4->GetParameter(9));
    fit5->SetParLimits(10,down*fit4->GetParameter(10), up*fit4->GetParameter(10));
    fit5->SetParLimits(11,down*fit4->GetParameter(11), up*fit4->GetParameter(11));
    fit5->SetParLimits(12,down*fit4->GetParameter(12), up*fit4->GetParameter(12));
    fit5->SetParLimits(13,down*fit4->GetParameter(13), up*fit4->GetParameter(13));
    fit5->SetParameter(0,  fit4->GetParameter(0));
    fit5->SetParameter(1,  fit4->GetParameter(1));
    fit5->SetParameter(2,  fit4->GetParameter(2));
    fit5->SetParameter(3,  fit4->GetParameter(3));
    fit5->SetParameter(4,  fit4->GetParameter(4));
    fit5->SetParameter(5,  fit4->GetParameter(5));
    fit5->SetParameter(6,  fit4->GetParameter(6));
    fit5->SetParameter(7,  fit4->GetParameter(7));
    fit5->SetParameter(8,  fit4->GetParameter(8));
    fit5->SetParameter(9,  fit4->GetParameter(9));
    fit5->FixParameter(10, fit4->GetParameter(10));
    fit5->FixParameter(11, fit4->GetParameter(11)); 
    fit5->FixParameter(12, fit4->GetParameter(12));
    fit5->FixParameter(13, fit4->GetParameter(13));
    hfit->Fit(fit5, "RQML", "", fit4->GetParameter(2), 1500);
    fit5->SetLineColor(kBlack);
    fit5->SetLineWidth(2);
    //fit5->Draw("same");
        
    vars.ped  = fit5->GetParameter(2);
    vars.meanPE = fit5->GetParameter(6);
    vars.gain = fit5->GetParameter(7);
    vars.ctProb = fit5->GetParameter(9);
    vars.cName = hfit->GetName();
    tree->Fill();

    //printf("Chi^2:%10.4f Gain:%10.4f Unc:%10.4f MPV:%10.4f Unc:%10.4f", fit5->GetChisquare(), fit5->GetParameter(7), fit5->GetParError(7), fit5->GetParameter(11), fit5->GetParError(11));
    //printf("  Run: %i Hname: %s\n",run, hfit->GetName());
    printf("Chi^2:%10.4f, P0:%10.4f, P1:%10.4f", fit5->GetChisquare(), fit5->GetParameter(0), fit5->GetParameter(1));
    printf(" P2:%10.4f, P3:%10.4f, P4:%10.4f, P5:%10.4f", fit5->GetParameter(2), fit5->GetParameter(3), fit5->GetParameter(4), fit5->GetParameter(5));
    printf(" P6:%10.4f, P7:%10.4f, P8:%10.4f, P9:%10.4f", fit5->GetParameter(6), fit5->GetParameter(7), fit5->GetParameter(8), fit5->GetParameter(9));
    printf(" P10:%10.4f, P11:%10.4f, P12:%10.4f, P13:%10.4f\n", fit5->GetParameter(10), fit5->GetParameter(11), fit5->GetParameter(12), fit5->GetParameter(13));
    hfit->GetYaxis()->SetTitle("Events / ADC bin");
    hfit->GetXaxis()->SetTitle("ADC Counts");
    hfit->SetTitleOffset(1,"X");
    hfit->SetTitleOffset(1.2,"Y");
    hfit->SetTitleSize(0.05,"X");
    hfit->SetTitleSize(0.05,"Y");
    
    TF1* draw5PE = new TF1("BackGround_MIP", background_MIP, 50.0, 2000.0, 14);
    draw5PE->FixParameter(0,  fit5->GetParameter(0));
    draw5PE->FixParameter(1,  fit5->GetParameter(1));
    draw5PE->FixParameter(2,  fit5->GetParameter(2));
    draw5PE->FixParameter(3,  fit5->GetParameter(3));
    draw5PE->FixParameter(4,  fit5->GetParameter(4));
    draw5PE->FixParameter(5,  fit5->GetParameter(5));
    draw5PE->FixParameter(6,  fit5->GetParameter(6));
    draw5PE->FixParameter(7,  fit5->GetParameter(7));
    draw5PE->FixParameter(8,  fit5->GetParameter(8));
    draw5PE->FixParameter(9,  fit5->GetParameter(9));
    draw5PE->FixParameter(10, fit5->GetParameter(10));
    draw5PE->FixParameter(11, fit5->GetParameter(11));
    draw5PE->FixParameter(12, fit5->GetParameter(12));
    draw5PE->FixParameter(13, fit5->GetParameter(13));
    draw5PE->SetLineWidth(2);
    draw5PE->SetLineColor(kBlue);
    draw5PE->Draw("same");
    
    TF1* draw5Mip = new TF1("Lang", langaufun, 50.0, 2000.0, 4);
    draw5Mip->FixParameter(0, fit5->GetParameter(12));
    draw5Mip->FixParameter(1, fit5->GetParameter(11)+fit5->GetParameter(2));
    draw5Mip->FixParameter(2, fit5->GetParameter(10));
    draw5Mip->FixParameter(3, fit5->GetParameter(13));
    draw5Mip->SetLineWidth(2);
    draw5Mip->SetLineColor(kGreen+2);
    // draw5Mip->Draw("same");
    
    //Make Plots Pretty
    hfit->SetStats(false);
    hfit->SetTitle("");
    
    TLatex* CMSPrelim1 = new TLatex(0.14, 0.91, "3mm SiPM");
    //TLatex* CMSPrelim1 = new TLatex(0.14, 0.91, "Cherenkov Fiber");
    CMSPrelim1->SetNDC();
    CMSPrelim1->SetTextFont(42);
    
    TLatex* testbeam = new TLatex(0.95, 0.91, "HG-DREAM 2025");
    testbeam->SetNDC();
    testbeam->SetTextFont(42);
    testbeam->SetTextAlign(31);
    
    TLatex* SiPMTitle = new TLatex(0.93, 0.86, "SiPM (Silicon Photomultiplier)");
    SiPMTitle->SetNDC();
    SiPMTitle->SetTextFont(42);
    SiPMTitle->SetTextAlign(32);
    
    TLatex* Muon = new TLatex(0.93, 0.87, "Cosmics");
    Muon->SetNDC();
    Muon->SetTextFont(42);
    Muon->SetTextAlign(32);
    
    // char chan [100];
    // int iEta, iPhi, iDepth;
    // //sscanf (hfit->GetName(),"beam_adc_%d_%d_%d", &iEta, &iPhi, &iDepth);
    // //sscanf (hfit->GetName(),"adc_nosub_%d_%d_%d", &iEta, &iPhi, &iDepth);
    // //sscanf (hfit->GetName(),"adc_nosub_constBin_%d_%d_%d", &iEta, &iPhi, &iDepth);
    // sscanf (hfit->GetName(),"adc_nosub_binChris_%d_%d_%d", &iEta, &iPhi, &iDepth);
    // //sscanf (hfit->GetName(),"ped_adc_%d_%d_%d", &iEta, &iPhi, &iDepth);
    // sprintf (chan, "Channel: %d,%d,%d", iEta, iPhi, iDepth);
    
    // TLatex* channel = new TLatex(0.93, 0.86, chan);
    // channel->SetNDC();
    // channel->SetTextFont(42);
    // channel->SetTextAlign(32);
    
    char mpv [100];
    float intmpv = fit5->GetParameter(11);
    sprintf (mpv,"MPV: %0.3f", intmpv);
    
    TLatex* MPV = new TLatex(0.93, 0.8, mpv);
    MPV->SetNDC();
    MPV->SetTextFont(42);
    MPV->SetTextAlign(31);
    
    char gain [100];
    float intgain = fit5->GetParameter(7);
    sprintf (gain,"SiPM Gain: %0.3f", intgain);
    
    TLatex* Gain = new TLatex(0.93, 0.75, gain);
    Gain->SetNDC();
    Gain->SetTextFont(42);
    Gain->SetTextAlign(31);

    char ped [100];
    float intped = fit5->GetParameter(2);
    sprintf (ped,"Pedestal: %0.3f", intped);
    
    TLatex* Ped = new TLatex(0.93, 0.7, ped);
    Ped->SetNDC();
    Ped->SetTextFont(42);
    Ped->SetTextAlign(31);    

    //CMSPrelim1->Draw();
    testbeam->Draw();
    //SiPMTitle->Draw();
    Muon->Draw();
    //channel->Draw();
    MPV->Draw();
    Gain->Draw();
    Ped->Draw();
    
    std::string oname = std::string(hfit->GetName()) + ".pdf";
    //sprintf(oname, "%s_SiPMRuns_3030to3475.pdf", hfit->GetName());
    //sprintf(oname, "%s_HBRuns_3526to3534.pdf", hfit->GetName());
    c1.Print(oname.c_str());
}

int main()
{
    //char baseFile[]         = "/Users/mad24679/Documents/TTU-Research/CaloX/PulseShapeProcesses/run0583_small.root";
    //char baseFile[]         = "/Users/mad24679/Documents/TTU-Research/CaloX/PulseShapeProcesses/run0595_250610144350.root";
    char baseFile[]         = "/Users/mad24679/Documents/TTU-Research/CaloX/PulseShapeProcesses/run1020_250703180532.root";
    char treeName[]         = "EventTree";
    TChain *chBase = new TChain(treeName);
    chBase->Add(baseFile);

    // Create output file, tree, and branches
    TFile* file = new TFile("output.root", "RECREATE");
    TTree* tree = new TTree("myTree", "A simple tree");

    TreeVars vars;
    tree->Branch("ped", &vars.ped, "ped/F");
    tree->Branch("meanPE", &vars.meanPE, "meanPE/F");
    tree->Branch("gain", &vars.gain, "gain/F");
    tree->Branch("ctProb", &vars.ctProb, "ctProb/F");
    tree->Branch("cName", &vars.cName);

    try
    {
        NTupleReader tr(chBase, {"FERS_Board0_energyHG"});

        std::vector<double> binEdges = {
            //0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190,
            0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 
            100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195,
            200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265, 270, 275, 280, 285, 290, 295, 300, 305, 310, 315, 320, 325, 330, 335, 340, 345, 350, 355, 360, 365, 370, 375, 380, 385, 390, 395,

            400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 510, 520, 530, 540, 550, 560, 570, 580, 590, 
            600, 610, 620, 630, 640, 650, 660, 670, 680, 690, 700, 710, 720, 730, 740, 750, 760, 770, 780, 790, 
            800, 810, 820, 830, 840, 850, 860, 870, 880, 890, 900, 910, 920, 930, 940, 950, 960, 970, 980, 990, 
            1000, 1010, 1020, 1030, 1040, 1050, 1060, 1070, 1080, 1090, 1100, 1110, 1120, 1130, 1140, 1150, 1160, 1170, 1180, 1190, 1200, 1210, 1220, 1230, 1240, 1250, 1260, 1270, 1280, 1290, 1300, 1310, 1320, 1330, 1340, 1350, 1360, 1370, 1380, 1390, 1400, 1410, 1420, 1430, 1440, 1450, 1460, 1470, 1480, 1490,
            1500, 1520, 1540, 1560, 1580, 1600, 1620, 1640, 1660, 1680, 1700, 1720, 1740, 1760, 1780, 1800, 1820, 1840, 1860, 1880, 1900, 1920, 1940, 1960, 1980,
            2000, 2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100, 2110, 2120, 2130, 2140, 2150, 2160, 2170, 2180, 2190, 2200, 2210, 2220, 2230, 2240, 2250, 2260, 2270, 2280, 2290, 2300, 2310, 2320, 2330, 2340, 2350, 2360, 2370, 2380, 2390, 2400, 2410, 2420, 2430, 2440, 2450, 2460, 2470, 2480, 2490, 2500, 2510, 2520, 2530, 2540, 2550, 2560, 2570, 2580, 2590, 2600, 2610, 2620, 2630, 2640, 2650, 2660, 2670, 2680, 2690, 2700, 2710, 2720, 2730, 2740, 2750, 2760, 2770, 2780, 2790, 2800, 2810, 2820, 2830, 2840, 2850, 2860, 2870, 2880, 2890, 2900, 2910, 2920, 2930, 2940, 2950, 2960, 2970, 2980, 2990, 
            3000, 3010, 3020, 3030, 3040, 3050, 3060, 3070, 3080, 3090, 3100, 3110, 3120, 3130, 3140, 3150, 3160, 3170, 3180, 3190, 3200, 3210, 3220, 3230, 3240, 3250, 3260, 3270, 3280, 3290, 3300, 3310, 3320, 3330, 3340, 3350, 3360, 3370, 3380, 3390, 3400, 3410, 3420, 3430, 3440, 3450, 3460, 3470, 3480, 3490, 3500, 3510, 3520, 3530, 3540, 3550, 3560, 3570, 3580, 3590, 3600, 3610, 3620, 3630, 3640, 3650, 3660, 3670, 3680, 3690, 3700, 3710, 3720, 3730, 3740, 3750, 3760, 3770, 3780, 3790, 3800, 3810, 3820, 3830, 3840, 3850, 3860, 3870, 3880, 3890, 3900, 3910, 3920, 3930, 3940, 3950, 3960, 3970, 3980, 3990, 
            4000, 4010, 4020, 4030, 4040, 4050, 4060, 4070, 4080, 4090, 4100, 4110, 4120, 4130, 4140, 4150, 4160, 4170, 4180, 4190, 4200, 4210, 4220, 4230, 4240, 4250, 4260, 4270, 4280, 4290, 4300, 4310, 4320, 4330, 4340, 4350, 4360, 4370, 4380, 4390, 4400, 4410, 4420, 4430, 4440, 4450, 4460, 4470, 4480, 4490, 4500, 4510, 4520, 4530, 4540, 4550, 4560, 4570, 4580, 4590, 4600, 4610, 4620, 4630, 4640, 4650, 4660, 4670, 4680, 4690, 4700, 4710, 4720, 4730, 4740, 4750, 4760, 4770, 4780, 4790, 4800, 4810, 4820, 4830, 4840, 4850, 4860, 4870, 4880, 4890, 4900, 4910, 4920, 4930, 4940, 4950, 4960, 4970, 4980, 4990, 5000
            //400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 820, 840, 860, 880, 900, 920, 940, 960, 980,
            //1000, 1025, 1050, 1075, 1100, 1125, 1150, 1175, 1200, 1225, 1250, 1275, 1300, 1325, 1350, 1375, 1400, 1425, 1450, 1475, 1500, 1525, 1550, 1575, 1600, 1625, 1650, 1675, 1700, 1725, 1750, 1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000

        };

        std::vector<std::pair<std::string, std::vector<std::shared_ptr<TH1D>>>> hVec;
        hVec.push_back(std::make_pair("FERS_Board3",  std::vector<std::shared_ptr<TH1D>>()));
        hVec.push_back(std::make_pair("FERS_Board11", std::vector<std::shared_ptr<TH1D>>()));

        // Loop over the events in the tree
        while(tr.getNextEvent())
        {
            for(auto& h : hVec)
            {
                const auto& hg = tr.getVec<unsigned short>(h.first + "_energyHG");
                if(h.second.empty()) // Initialize histograms on first event
                {
                    for(unsigned int i = 0; i < hg.size(); i++)
                    {
                        std::string name = h.first + "_Channel" + std::to_string(i);
                        h.second.push_back(std::make_shared<TH1D>(name.c_str(), name.c_str(), 1000, 0, 2000));
                        //h.second.push_back(std::make_shared<TH1D>(name3.c_str(), name3.c_str(),  binEdges.size()-1, binEdges.data()));
                    }
                }
                for(unsigned int i = 0; i < hg.size(); i++)
                {
                    if(hg[i] > 20 && hg[i] < 5000) h.second[i]->Fill(hg[i]);
                }
            }
        }

        // Run the fit for each histogram
        for(auto& h : hVec)
        {
            for(unsigned int i = 0; i < h.second.size(); i++)
            {
                auto nEvents = h.second[i]->Integral();
                for(unsigned int j = 0; j < h.second[i]->GetNbinsX()+1; j++)
                {
                    auto binWidth = h.second[i]->GetBinWidth(j);
                    auto binVal = h.second[i]->GetBinContent(j);
                    auto binError = h.second[i]->GetBinError(j);
                    h.second[i]->SetBinContent(j, binVal/binWidth);
                    h.second[i]->SetBinError(j, binError/binWidth);
                }
                h.second[i]->Scale(nEvents/h.second[i]->Integral());

                fitSPEMIP(h.second[i].get(), i, vars, tree);
            }
        }

        // Save the histograms to the file
        tree->Write();
        file->Close();
        delete file;
        delete chBase;
    }
    catch(const NTRException& e)
    {
        e.print();
    }

    return 0;

}
