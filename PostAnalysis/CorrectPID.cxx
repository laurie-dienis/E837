#ifndef CorrectPID_cxx
#define CorrectPID_cxx

#include "ActColors.h"
#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActPIDCorrector.h"
#include "ActSilMatrix.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDataFrame.hxx"
#include "Rtypes.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TROOT.h"
#include "TString.h"

#include <fstream>
#include <iostream>
#include <string>

#include "./HistConfig.h"
#include "./Utils.cxx"

void CorrectPID(const std::string& light = "4He")
{
    ROOT::EnableImplicitMT();
    // Manager of  cuts
    ActRoot::CutsManager<std::string> cuts;
    // Read data
    ActRoot::DataManager datman {"/home/eactar/Analysis_e837/Analysis/configs/data_post.conf",
                                 ActRoot::ModeType::EMerge};
    auto chain {datman.GetJoinedData()};
    ROOT::RDataFrame d {*chain};
    // Apply basic cut on multiplicity
    auto df {d.Filter("fSilLayers.size() == 1")}; // && fThetaLight<20")};
    df = df.Define("ESil", "(double)fSilEs.front()");

    // Book histograms
    auto hPID {df.Histo2D(HistConfig::PID, "ESil", "fQave")};
    // df.Foreach([](const ActRoot::MergerData& d) { d.fQave; }, {"MergerData"});


    // Read PID
    cuts.ReadCut(light, TString::Format("./Cuts/LightPID/pid_%s.root", light.c_str()));

    // Init PID corrector
    auto* hModel {new TH2D {"HModel", "PID model", 200, 0, 350, 600, 0, 2000}};
    ActPhysics::PIDCorrector corr {"corr", cuts.GetListOfKeys(), hModel};
    double minE {8};
    double maxE {9};
    // Fill corrector
    df.Foreach(
        [&](const ActRoot::MergerData& d)
        {
            if(cuts.IsInside(light, d.fSilEs.front(), d.fQave))
            {
                corr.FillHisto(light, d.fSP.Z(), d.fQave, d.fSilEs.front(), minE, maxE);
            }
        },
        {"MergerData"});

    // SP per silicon
    // int nsil {20};
    // std::vector<TH1D*> h1ds;
    // std::vector<TH2D*> hs, hsq;
    // for(int s = 0; s < nsil; s++)
    // {
    //     hs.push_back(new TH2D {TString::Format("hSP%d", s), "SP colored", 200, 0, 300, 200, 0, 500});
    // }
    // // Fill colored histo
    // vetoed.Foreach([&](const ActRoot::MergerData& d) { hs[d.fSilNs.front()]->Fill(d.fSP.Y(), d.fSP.Z()); },
    //                {"MergerData"});

    // TString pidfile {};
    // pidfile = TString::Format("./Cuts/LightPID/pid_%s.root", light.c_str());
    // cuts.ReadCut(light, pidfile);
    // cuts.ReadCut("debug", "./Cuts/Debug/pid_central.root");
    // std::cout << BOLDCYAN << "Reading light PID in : " << pidfile << RESET << '\n';


    // plotting
    auto* c10 {new TCanvas("c10", "Pipe1 canvas 0")};
    hPID->DrawClone("colz");
    cuts.DrawAll();

    // Fit and plot
    corr.GetProfiles();
    corr.FitProfiles(150, 350);
    corr.Draw();
}
#endif
