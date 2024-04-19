#include "ActDataManager.h"
#include "ActMergerData.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TROOT.h"
void QaveX()
{
    // Read data
    ActRoot::DataManager datman {"../configs/data_post.conf", ActRoot::ModeType::EMerge};
    auto chain {datman.GetJoinedData()};

    auto* h {new TH2D {"hRP", "Qave against X;X [mm];Qave [mm^{-1}]", 150, 0, 260, 1000, 400, 1000}};
    // Process
    auto* m {new ActRoot::MergerData};
    chain->SetBranchAddress("MergerData", &m);
    for(auto i = 0; i < chain->GetEntries(); i++)
    {
        chain->GetEntry(i);
        auto rp {m->fRP};
        if(rp.X() == -1)
            continue;
        h->Fill(rp.X(), m->fQave);
    }


    // Plot
    auto* c0 {new TCanvas {"c0", "Qave dependence on X"}};
    h->DrawClone("colz");
}
