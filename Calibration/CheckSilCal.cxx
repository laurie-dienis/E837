#include "ActSilData.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TH1.h"
#include "TString.h"

#include <string>
#include <vector>
void CheckSilCal()
{
    // Read the data
    ROOT::RDataFrame df {"ACTAR_Data", "/scratch/dienis/RootFiles/Data/Data_Run_0094_Calibrated.root"};
    ROOT::RDataFrame raw {"ACTAR_Data", "/scratch/dienis/RootFiles/Data/Data_Run_0094_Uncalibrated.root"};

    // Set parameters
    const int nsil {12};
    // Init histograms
    std::vector<TH1D*> hs, hraw;
    for(int s = 0; s < nsil; s++)
    {
        hs.push_back(
            new TH1D {TString::Format("hSil%d", s), TString::Format("Calibrated Sil %d;E_{Sil} [MeV]", s), 500, 0, 10});
        hraw.push_back(
            new TH1D {TString::Format("hRaw%d", s), TString::Format("Raw Sil %d;E_{Sil} [MeV]", s), 800, 0, 8164});
    }

    // Fill the histograms
    const std::string layer {"f0"};
    df.Foreach(
        [&](const ActRoot::SilData& d)
        {
            if(d.fSiE.size() < 1)
                return;
            for(int i = 0, size = d.fSiE.at(layer).size(); i < size; i++)
            {
                auto N {d.fSiN.at(layer).at(i)};
                auto E {d.fSiE.at(layer).at(i)};
                hs[N]->Fill(E);
            }
        },
        {"SilData"});

    raw.Foreach(
        [&](const ActRoot::SilData& d)
        {
            for(int i = 0, size = d.fSiE.at(layer).size(); i < size; i++)
            {
                auto N {d.fSiN.at(layer).at(i)};
                auto E {d.fSiE.at(layer).at(i)};
                hraw[N]->Fill(E);
            }
        },
        {"SilData"});


    // Plot
    auto* c0 {new TCanvas {"c0", "Check sil canvas"}};
    c0->DivideSquare(hs.size());
    for(int s = 0; s < hs.size(); s++)
    {
        c0->cd(s + 1);
        hs[s]->Draw();
    }
    auto* c1 {new TCanvas {"c1", "Raw sil canvas"}};
    c1->DivideSquare(hraw.size());
    for(int s = 0; s < hraw.size(); s++)
    {
        c1->cd(s + 1);
        hraw[s]->Draw();
    }
}
