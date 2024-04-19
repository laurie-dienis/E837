#include "ActDataManager.h"
#include "ActMergerData.h"


#include "TCanvas.h"
#include "TH2.h"
#include "TROOT.h"

#include "Math/Point3D.h"

void RP()
{
    // Read data
    // Get data
    ActRoot::DataManager datman {"../configs/data_post.conf", ActRoot::ModeType::EMerge};
    auto chain {datman.GetJoinedData()};

    auto* hRP {new TH2D {"hRP", "RP;X [mm];Y [mm]", 150, 0, 260, 150, 0, 260}};
    // Process
    auto* m {new ActRoot::MergerData};
    chain->SetBranchAddress("MergerData", &m);
    for(auto i = 0; i < chain->GetEntries(); i++)
    {
        chain->GetEntry(i);
        auto rp {m->fRP};
        if(rp.X() == -1)
            continue;
        hRP->Fill(rp.X(), rp.Y());
    }

    auto* c0 {new TCanvas {"c0", "rp canvas"}};
    hRP->Draw("colz");
}
