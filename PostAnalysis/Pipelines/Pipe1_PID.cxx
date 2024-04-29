#ifndef Pipe1_PID_cxx
#define Pipe1_PID_cxx

#include "ActColors.h"
#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActSilMatrix.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDataFrame.hxx"
#include "Rtypes.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TROOT.h"
#include "TString.h"

#include <fstream>
#include <iostream>
#include <string>

#include "../HistConfig.h"
#include "../Utils.cxx"

void Pipe1_PID(const std::string& beam, const std::string& target, const std::string& light, double e_beam_i,
               int pressure)
{
    ROOT::EnableImplicitMT();
    // Manager of  cuts
    ActRoot::CutsManager<std::string> cuts;
    // Read data
    ActRoot::DataManager datman {"/home/dienis/Analysis_e837/configs/data_post.conf",
                                 ActRoot::ModeType::EMerge};
    auto chain {datman.GetJoinedData()};
    auto chainTpc {datman.GetJoinedData(ActRoot::ModeType::EReadTPC)};
    chain->AddFriend(chainTpc.get());
    ROOT::RDataFrame d {*chain};
    // Apply basic cut on multiplicity
    auto df {d.Filter("fSilLayers.size() == 1")}; // && fThetaLight<20")};

    // Preliminary cut to gate 8He
    TH1::AddDirectory(false);
    auto clean {[&](const ActRoot::TPCData& d)
                {
                    //d.Print();
                    // 1-> Fill Qproj X
                    TH1D p {"pQx", "Charge projection onto X;X [pad];Counts", 128, 0, 128};
                    // Fill from clusters
                    for(auto& cl : d.fClusters)
                    {
                        for(auto& v : cl.GetVoxels())
                            p.Fill(v.GetPosition().X(), v.GetCharge());
                    }
                    // 2-> Find maximum
                    auto maxBin = p.GetMaximumBin();
                    auto xMax = p.GetBinCenter(maxBin);
                    auto yMax = p.GetBinContent(maxBin);
                    auto yMax_afterpeak = p.GetBinContent(maxBin + 20.);
                    if(yMax_afterpeak > 10.)
                    {
                        return false;
                        //streamer << d.fRun << " " << d.fEntry << '\n';
                    }
                    // std::cout << "xmax=" << xMax << "\n";
                    // std::cout << "ymax=" << yMax << "\n";
                    // double fitWidth {3};
                    // double minPad {15};
                    // double maxPad {113};
                    // double sigmaThresh {3};
                    // bool fit {minPad <= xMax && xMax <= maxPad};
                    // if(fit)
                    // {
                    //     auto fitMin = xMax - fitWidth;
                    //     auto fitMax = xMax + fitWidth;
                    //     p.Fit("gaus", "0QR", "", fitMin, fitMax);
                    //     // Get function
                    //     auto* func {p.GetFunction("gaus")};
                    //     if(func)
                    //     {
                    //         double sigma {func->GetParameter("Sigma")};
                    //         if(sigma < sigmaThresh)
                    //             return false;
                    //     }
                    // }
                    return true;
                }};

    // Filter of 8He
    ROOT::RDF::RNode vetoed {df};
    // Spped up executation : ony apply this heavy cut for p = 700 mbar
    // For p = 900 mbar the separation is way better by itself
    if(pressure == 800)
        vetoed = df.Filter(clean, {"TPCData"});
    // Define ESil as alias of .front of energy vector
    vetoed = vetoed.Define("ESil", "(double)fSilEs.front()");

    // Book histograms
    auto hPID {vetoed.Histo2D(HistConfig::PID, "ESil", "fQave")};
    auto hSP {vetoed.Histo2D(HistConfig::SP, "fSP.fCoordinates.fY", "fSP.fCoordinates.fZ")};
    cuts.ReadCut("circular", "./Cuts/Debug/circular_beam.root");
    auto noHe8 = vetoed.Filter(
        [&](const ActRoot::MergerData& d)
        {
            // if(cuts.IsInside("circular", d.fSP.Y(), d.fSP.Z()))
            //     return false;
            // else
            //     return true;
            auto N {d.fSilNs.front()};
            if(N == 0 || N == 2 || N == 6 || N == 8 || N == 9 || N == 10 || N == 11)
            {
                return true;
            }
            else
                return false;
        },
        {"MergerData"});
    auto hPID_selected {noHe8.Histo2D(HistConfig::PID_selected, "ESil", "fQave")};

    // SP per silicon
    int nsil {20};
    std::vector<TH1D*> h1ds;
    std::vector<TH2D*> hs, hsq;
    for(int s = 0; s < nsil; s++)
    {
        hs.push_back(new TH2D {TString::Format("hSP%d", s), "SP colored", 200, 0, 300, 200, 0, 500});
    }
    // Fill colored histo
    vetoed.Foreach([&](const ActRoot::MergerData& d) { hs[d.fSilNs.front()]->Fill(d.fSP.Y(), d.fSP.Z()); },
                   {"MergerData"});

    TString pidfile {};
    pidfile = TString::Format("./Cuts/LightPID/pid_%s_%dmbar.root", light.c_str(), pressure);
    cuts.ReadCut(light, pidfile);
    cuts.ReadCut("debug", "./Cuts/Debug/carbon.root");
    std::cout << BOLDCYAN << "Reading light PID in : " << pidfile << RESET << '\n';

    std::ofstream streamer {"./carbon.dat"};
    vetoed.Foreach(
        [&](const ActRoot::MergerData& d)
        {
            if(cuts.IsInside("debug", d.fSilEs.front(), d.fQave))
                // if(cut.IsInside("debug", d.fSP.Y(), d.fSP.Z()))
                streamer << d.fRun << " " << d.fEntry << '\n';
        },
        {"MergerData"});
    streamer.close();
    if(cuts.GetCut(light))
    {
        // Filter
        auto pid {vetoed.Filter([&](const ActRoot::MergerData& d, double ESil)
                                { return cuts.IsInside(light, ESil, d.fQave); },
                                {"MergerData", "ESil"})};
        auto filename {E837Utils::GetFileName(1, pressure, beam, target, light)};
        std::cout << BOLDCYAN << "Saving PID_Tree in file : " << filename << '\n';
        std::cout << "  with Nentries : " << pid.Count().GetValue() << '\n';
        pid.Snapshot("PID_Tree", filename, {"MergerData", "ESil"});
    }

    // plotting
    auto* c10 {new TCanvas("c10", "Pipe1 canvas 0")};
    c10->DivideSquare(4);
    c10->cd(1);
    hPID->DrawClone("colz");
    cuts.SetLineAttributes(light, kMagenta, 2);
    cuts.DrawAll();
    c10->cd(2);
    hSP->DrawClone("colz");
    cuts.DrawAll();
    c10->cd(3);
    std::vector<int> colors {kYellow, kOrange,       kRed,     kMagenta - 10, kMagenta,   kViolet, kBlue + 2,
                             kAzure,  kCyan,         kTeal,    kSpring,       kGreen + 2, kYellow, kOrange,
                             kRed,    kMagenta - 10, kMagenta, kViolet,       kBlue + 2,  kAzure};
    for(int i = 0; i < 20; i++)
    {
        hs[i]->SetMarkerStyle(6);
        hs[i]->SetMarkerColor(colors[i]);
        if(i == 0)
            hs[i]->Draw("scat");
        else
            hs[i]->Draw("scat same");
    }
    c10->cd(4);
    hPID_selected->DrawClone("colz");
    // hTPC->DrawClone("colz");
}
#endif
