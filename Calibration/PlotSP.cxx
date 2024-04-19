#include "ActMergerData.h"

#include "TCanvas.h"
#include "TChain.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH1D.h"
#include "TMultiGraph.h"
#include "TString.h"

#include <vector>
void PlotSP()
{
    // Get data
    auto* chain {new TChain {"ACTAR_Merged"}};
    chain->Add("/home/eactar/Analysis_e837/Analysis/RootFiles/Merger/Merged_Run_0031.root");

    // Set number of silicons
    const int nsil {12};
    std::vector<TH1D*> hs;
    std::vector<TGraph*> gs;
    for(int s = 0; s < nsil; s++)
    {
        hs.push_back(new TH1D {TString::Format("hSil%d", s), TString::Format("Sil %d raw;Channel", s), 2000, 0, 8192});
        gs.push_back(new TGraph);
        gs.back()->SetTitle(TString::Format("SP for %d", s));
    }

    // Fill
    ActRoot::MergerData* data {new ActRoot::MergerData};
    chain->SetBranchAddress("MergerData", &data);
    for(auto i = 0; i < chain->GetEntries(); i++)
    {
        chain->GetEntry(i);
        // Fill silicon energy for first hit
        auto N {data->fSilNs.front()};
        auto E {data->fSilEs.front()};
        hs[N]->Fill(E);
        // Fill SP
        gs[N]->SetPoint(gs[N]->GetN(), data->fSP.Y(), data->fSP.Z());
    }


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
}
