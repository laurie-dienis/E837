#include "TCanvas.h"
#include "ActSilData.h"
#include "ROOT/RDataFrame.hxx"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1.h"
#include "TUUID.h"

#include "CalibrationRunner.h"
#include "CalibrationSource.h"

#include <iostream>
#include <vector>

#include "CalibrationRunner.cxx"
#include "CalibrationSource.cxx"

std::vector<TH1D*> ReadData()
{
   // Read the data
    ROOT::RDataFrame df {"ACTAR_Data", "/scratch/dienis/RootFiles/Data/Data_Run_0094_Uncalibrated.root"};
   // Set parameters
    const int nsil {12};
    // Init histograms
    std::vector<TH1D*> hs;
    for(int s = 0; s < nsil; s++)
    {
        hs.push_back(
            new TH1D {TString::Format("hSil%d", s), TString::Format("Calibrated Sil %d;E_{Sil} [MeV]", s), 500, 0, 10});
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
    return hs;
}

void e837()
{
    auto hs {ReadData()};

    // Create source
    Calibration::Source source;
    source.Print();
    std::cout << "number = " << hs.size();

    auto* h2rebin {(TH1D*)hs[2]->Clone("h2Rebin")};
    h2rebin->Rebin(2);
    auto* h7rebin {(TH1D*)hs[7]->Clone("h7rebin")};
    h7rebin->Rebin(2);
    // Si_3
    Calibration::Runner runner2 {&source, h2rebin, hs[2],false};
    // Si_4
    Calibration::Runner runner7 {&source, h7rebin, hs[7],false};
    //Run for all
    for(auto& runner : {&runner2})
    {
        runner->SetRange(4000, 7000);
        runner->DoIt();
        //runner->Draw(new TCanvas);
        //runner->PrintRes();
    }


    // Plot
    //auto* c0 {new TCanvas {"c0", "Silicon canvas"}};
    //c0->DivideSquare(hs.size());
    //for(int i = 0; i < hs.size(); i++)
    //{
    //    c0->cd(i + 1);
    //    hs[i]->Draw();
    //}
}
