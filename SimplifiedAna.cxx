#include "ActDataManager.h"
#include "ActModularData.h"
#include "ActSilData.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"
#include "Rtypes.h"

#include "TCanvas.h"
#include "THStack.h"
#include "TString.h"

#include "Math/Point3Dfwd.h"

#include <algorithm>
#include <iterator>
#include <numeric>
#include <string>

struct SimplifiedData
{
    ROOT::Math::XYZPointF fSP {-1, -1, -1};
    double fQtotal {-1};
    double fSilE {-1};
    int fSilN {-1};
};

void SimplifiedAna()
{
    ActRoot::DataManager datman {"./configs/data.conf", ActRoot::ModeType::EFilter};
    auto chain {datman.GetJoinedData()};
    auto silChain {datman.GetJoinedData(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(silChain.get());

    // DF
    ROOT::RDataFrame df {*chain};

    // Define placement in pad units
    double xactar {128};
    double xsil {35};
    double driftFactor {8.33};
    std::string layer {"f0"};
    // Histogram parameters
    bool converted {(driftFactor != 1)};
    double xyfactor {(converted) ? 2. : 1.};
    int nbinsy {250};
    double ymax {(converted) ? 260. : 128.};
    int nbinsz {250};
    double zmax {(converted) ? 650. : 200.};

    // Define
    auto def {df.Define("SimplifiedData",
                        [&](ActRoot::TPCData& dtpc, ActRoot::SilData& dsil)
                        {
                            SimplifiedData res;
                            if(dtpc.fClusters.size() == 1)
                            {
                                // Reference to cluster
                                auto& track {dtpc.fClusters.front()};
                                // 1-> Compute SP
                                const auto& line {track.GetLine()};
                                auto sp {line.MoveToX(xactar + xsil)};
                                res.fSP = {sp.X() * xyfactor, sp.Y() * xyfactor, sp.Z() * driftFactor};
                                // 2-> Extract SilN and E
                                // Find maximum
                                auto maxE {std::max_element(dsil.fSiE[layer].begin(), dsil.fSiE[layer].end())};
                                auto pos {std::distance(dsil.fSiE[layer].begin(), maxE)};
                                if(maxE != dsil.fSiE[layer].end())
                                {
                                    res.fSilN = dsil.fSiN[layer][pos];
                                    res.fSilE = *maxE;
                                }
                                // 3-> Compute totalQ
                                auto qTotal {std::accumulate(track.GetVoxels().begin(), track.GetVoxels().end(), 0.,
                                                             [](double sum, const ActRoot::Voxel& v)
                                                             { return sum + v.GetCharge(); })};
                                res.fQtotal = qTotal;
                            }
                            return res;
                        },
                        {"TPCData", "SilData"})};
    // Book histograms
    auto hSP {
        def.Filter([](const SimplifiedData& d) { return d.fSilN == 4; }, {"SimplifiedData"})
            .Define("SPy", [](SimplifiedData& d) { return d.fSP.Y(); }, {"SimplifiedData"})
            .Define("SPz", [](SimplifiedData& d) { return d.fSP.Z(); }, {"SimplifiedData"})
            .Histo2D({"hSP", "Silicon point N = 4;Y [pad];Z [pad]", nbinsy, 0, ymax, nbinsz, 0, zmax}, "SPy", "SPz")};

    // Filling
    int nsil {12};
    std::vector<TH2D*> hs, hsq;
    for(int s = 0; s < nsil; s++)
    {
        hs.push_back(
            new TH2D {TString::Format("hSP%d", s), TString::Format("SP for %d", s), nbinsy, 0, ymax, nbinsz, 0, zmax});
        hsq.push_back(new TH2D {TString::Format("hQ%d", s), TString::Format("Q vs ESil for %d", s), 500, 0, 4000, 1000,
                                0, 500000});
    }
    def.Foreach(
        [&](const SimplifiedData& d)
        {
            if(d.fSilN != -1)
            {
                hs[d.fSilN]->Fill(d.fSP.Y(), d.fSP.Z());
                hsq[d.fSilN]->Fill(d.fSilE, d.fQtotal);
            }
        },
        {"SimplifiedData"});

    // Plot
    auto* c0 {new TCanvas {"c0", "SP canvas"}};
    c0->DivideSquare(2);
    c0->cd(1);
    hSP->DrawClone("colz");
    c0->cd(2);
    std::vector<int> colors {kYellow,   kOrange, kRed,  kMagenta - 10, kMagenta, kViolet,
                             kBlue + 2, kAzure,  kCyan, kTeal,         kSpring,  kGreen + 3};
    for(int i = 0; i < nsil; i++)
    {
        hs[i]->SetMarkerStyle(6);
        hs[i]->SetMarkerColor(colors[i]);
        if(i == 0)
            hs[i]->Draw("scat");
        else
            hs[i]->Draw("scat same");
    }

    auto* c1 {new TCanvas {"c1", "Q vs E sil canvas"}};
    // Stack
    auto* qstack {new THStack};
    for(auto& h : hsq)
        qstack->Add(h, "colz");

    qstack->Draw("pads");
    // hs[0]->SetMarkerStyle(8);
    // hs[0]->SetMarkerSize(0.5);
    // hs[0]->SetMarkerColor(kRed);
    // hs[0]->Draw("scat");
}
