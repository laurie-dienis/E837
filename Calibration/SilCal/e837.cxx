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
            new TH1D {TString::Format("hSil%d", s), TString::Format("Calibrated Sil %d;E_{Sil} [MeV]", s), 800, 0, 8154});
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
            	// std::cout << "size = " << d.fSiE.at(layer).size()<<"\n";
                auto N {d.fSiN.at(layer).at(i)};
                // std::cout << "N = " << N<<"\n";;
                auto E {d.fSiE.at(layer).at(i)};
                // std::cout << "E = " << E<<"\n";;
                            // Skip invalid indices
            if (N == 12) 
            {
                std::cerr << "Warning: Invalid index N = " << N << ", skipping...\n";
                continue;
            }
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

std::vector<Calibration::Runner> runners;

for (int i = 0; i <= 11; ++i) {
    if (i == 4) continue; // Exclude hs[4]

    auto* hRebin = (TH1D*)hs[i]->Clone(("h" + std::to_string(i) + "Rebin").c_str());
    hRebin->Rebin(1);

    // Create Calibration::Runner for the current histogram
    runners.emplace_back(&source, hRebin, hs[i], false);
}

// Run for all created runners
for (int i = 0; i < runners.size(); ++i) {
    std::cout << "\n--- Detector " << i << " ---" << std::endl; // Add space and detector number
    runners[i].SetRange(4700, 6700);
    runners[i].DoIt(); // Perform calibration (this should print parameters)
    runners[i].Draw(new TCanvas);
    runners[i].PrintRes();
}

    //Plot
    auto* c0 {new TCanvas {"c0", "Silicon canvas"}};
    c0->DivideSquare(hs.size());
    for(int i = 0; i < hs.size(); i++)
    {
       c0->cd(i + 1);
       hs[i]->Draw();
    }
}
