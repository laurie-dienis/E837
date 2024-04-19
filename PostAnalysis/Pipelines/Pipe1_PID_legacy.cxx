#ifndef Pipe1_PID_cxx
#define Pipe1_PID_cxx

#include "ActColors.h"
#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActSilMatrix.h"
#include "ActTypes.h"

#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDataFrame.hxx"
#include "Rtypes.h"

#include "TCanvas.h"
#include "TROOT.h"
#include "TString.h"

#include <fstream>
#include <iostream>
#include <string>

#include "../HistConfig.h"
#include "../Utils.cxx"

void Pipe1_PID(const std::string& beam, const std::string& target, const std::string& light)
{
    // ROOT::EnableImplicitMT();
    // Read data
    ActRoot::DataManager datman {"/home/eactar/Analysis_e837/Analysis/configs/data_post.conf",
                                 ActRoot::ModeType::EMerge};
    auto chain {datman.GetJoinedData()};
    ROOT::RDataFrame df {*chain};

    // Apply cuts
    auto vetoed {df.Filter("fSilLayers.size() == 1")};
    // vetoed = vetoed.Define("ESil", "fSilEs.front()* (9.5 / 5860)");
    vetoed = vetoed.Define("ESil", "(double)fSilEs.front()");

    // Book histograms
    auto hPID {vetoed.Histo2D(HistConfig::PID, "ESil", "fQave")};
    auto hSP {vetoed.Histo2D(HistConfig::SP, "fSP.fCoordinates.fY", "fSP.fCoordinates.fZ")};

    // Read PID cut
    ActRoot::CutsManager<std::string> cut;
    TString pidfile {};
    pidfile = TString::Format("./Cuts/LightPID/pid_%s.root", light.c_str());
    cut.ReadCut(light, pidfile);
    cut.ReadCut("debug", "./Cuts/Debug/pid_high.root");
    std::cout << BOLDCYAN << "Reading light PID in : " << pidfile << RESET << '\n';

    std::ofstream streamer {"./debug_pid.dat"};
    vetoed.Foreach(
        [&](const ActRoot::MergerData& d)
        {
            if(cut.IsInside("debug", d.fSilEs.front(), d.fQave))
                streamer << d.fRun << " " << d.fEntry << '\n';
        },
        {"MergerData"});
    streamer.close();
    if(cut.GetCut(light))
    {
        // Filter
        auto pid {vetoed.Filter([&](const ActRoot::MergerData& d, double ESil)
                                { return cut.IsInside(light, ESil, d.fQave); },
                                {"MergerData", "ESil"})};
        auto filename {E837Utils::GetFileName(1, beam, target, light)};
        std::cout << BOLDCYAN << "Saving PID_Tree in file : " << filename << '\n';
        std::cout << "  with Nentries : " << pid.Count().GetValue() << '\n';
        pid.Snapshot("PID_Tree", filename);

        // // Write
        // std::ofstream streamer {"./Hes_veto.dat"};
        // pid.Foreach([&](const ActRoot::MergerData& d) { streamer << d.fRun << " " << d.fEntry << '\n'; },
        //             {"MergerData"});
        // streamer.close();
    }

    // plotting
    auto* c10 {new TCanvas("c10", "Pipe1 canvas 0")};
    c10->DivideSquare(2);
    c10->cd(1);
    hPID->DrawClone("colz");
    cut.SetLineAttributes(light, kMagenta, 2);
    cut.DrawAll();
    c10->cd(2);
    hSP->DrawClone("colz");
}
#endif
