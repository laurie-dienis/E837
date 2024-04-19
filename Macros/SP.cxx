#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TChain.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2D.h"
#include "THStack.h"
#include "TMultiGraph.h"
#include "TString.h"

#include <fstream>
#include <vector>
void SP()
{
    // Get data
    ActRoot::DataManager datman {"../configs/data_post.conf", ActRoot::ModeType::EMerge};
    auto chain {datman.GetJoinedData()};

    // Set number of silicons
    const int nsil {12};
    std::vector<TH1D*> hs;
    std::vector<TGraph*> gs;
    std::vector<TH2D*> hsp;
    for(int s = 0; s < nsil; s++)
    {
        hs.push_back(new TH1D {TString::Format("hSil%d", s), TString::Format("Sil %d raw;Channel", s), 300, 0, 30});
        gs.push_back(new TGraph);
        gs.back()->SetTitle(TString::Format("SP for %d", s));
        hsp.push_back(new TH2D {TString::Format("hSP%d", s), TString::Format("SP for %d;Y [mm];Z [mm]", s), 120, 0, 250,
                                120, 0, 400});
    }

    // Debug cuts
    ActRoot::CutsManager<int> cuts;
    cuts.ReadCut(0, "./Cuts/debug_sp.root");
    std::ofstream streamer {"./debug_sp.dat"};
    // Fill
    ActRoot::MergerData* data {new ActRoot::MergerData};
    chain->SetBranchAddress("MergerData", &data);
    for(auto i = 0; i < chain->GetEntries(); i++)
    {
        chain->GetEntry(i);
        // Condition on size
        if(data->fSilNs.size() < 1)
            continue;
        // Fill silicon energy for first hit
        auto N {data->fSilNs.front()};
        auto E {data->fSilEs.front()};
        hs[N]->Fill(E);
        // Fill SP
        gs[N]->SetPoint(gs[N]->GetN(), data->fSP.Y(), data->fSP.Z());
        hsp[N]->Fill(data->fSP.Y(), data->fSP.Z());

        // Print to file
        if(cuts.IsInside(0, data->fSP.Y(), data->fSP.Z()))
        {
            streamer << data->fRun << " " << data->fEntry << '\n';
        }
    }
    streamer.close();


    // Plot
    auto* c0 {new TCanvas {"c0", "Raw sil data"}};
    c0->DivideSquare(hs.size());
    for(int s = 0; s < hs.size(); s++)
    {
        c0->cd(s + 1);
        hs[s]->Draw();
    }

    // Multigraph for graphs
    auto* mg {new TMultiGraph};
    for(auto& g : gs)
    {
        g->SetMarkerStyle(6);
        mg->Add(g);
    }

    auto* c1 {new TCanvas {"c1", "SP per Sil"}};
    mg->Draw("ap plc pmc");

    // HStack for histograms
    auto* stack {new THStack};
    stack->SetTitle("SP;Y [mm];Z [mm]");
    for(auto& h : hsp)
    {
        stack->Add(h);
    }

    // Plot
    auto* c2 {new TCanvas {"c2", "SP per sil histo"}};
    stack->Draw("nostack plc pmc");
    cuts.DrawAll();
}
