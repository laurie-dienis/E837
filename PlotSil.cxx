#include "ActDataManager.h"
#include "ActModularData.h"
#include "ActSilData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TString.h"

#include <vector>

void PlotSil()
{
    ActRoot::DataManager datman {"./configs/data.conf", ActRoot::ModeType::EReadSilMod};
    auto chain {datman.GetJoinedData()};

    // Init DF
    ROOT::RDataFrame df {*chain};

    // Gatconf just to check
    auto hGATCONF {df.Define("GATCONF", [](ActRoot::ModularData& d) { return d.fLeaves["GATCONF"]; }, {"ModularData"})
                       .Histo1D("GATCONF")};

    // Init histograms
    std::string layer {"f0"};
    int nsil {12};
    std::vector<TH1D*> hs;
    for(int s = 0; s < nsil; s++)
    {
        hs.push_back(new TH1D(TString::Format("hSil%d", s), TString::Format("Raw sil %d;Channel", s), 4000, 0, 8192));
    }

    // Do!
    df.Foreach(
        [&](ActRoot::SilData& d)
        {
            for(int i = 0; i < d.fSiN[layer].size(); i++)
            {
                auto N {d.fSiN[layer][i]};
                auto E {d.fSiE[layer][i]};
                hs[N]->Fill(E);
            }
        },
        {"SilData"});


    // Plot
    auto* c0 {new TCanvas {"c0", "Raw sil data"}};
    c0->DivideSquare(hs.size());
    for(int i = 0; i < hs.size(); i++)
    {
        c0->cd(i + 1);
        hs[i]->Draw();
    }

    auto c1 {new TCanvas {"c1", "GATCONF"}};
    hGATCONF->DrawClone();
}
