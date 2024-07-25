#include "ActCalibrationManager.h"
#include "ActTPCDetector.h"
#include "ActTPCLegacyData.h"

#include "TCanvas.h"
#include "TChain.h"
#include "TH2D.h"
#include "TString.h"

#include <set>

void FillHistogram(ActRoot::CalibrationManager* calman, ActRoot::TPCParameters* pars, MEventReduced* evt, TH2D* h,
                   bool isMatched)
{
    // iterate over hits
    for(int it = 0, size = evt->CoboAsad.size(); it < size; it++)
    {
        int co = evt->CoboAsad[it].globalchannelid >> 11;
        int as = (evt->CoboAsad[it].globalchannelid - (co << 11)) >> 9;
        int ag = (evt->CoboAsad[it].globalchannelid - (co << 11) - (as << 9)) >> 7;
        int ch = evt->CoboAsad[it].globalchannelid - (co << 11) - (as << 9) - (ag << 7);
        int where = co * pars->GetNBASAD() * pars->GetNBAGET() * pars->GetNBCHANNEL() +
                    as * pars->GetNBAGET() * pars->GetNBCHANNEL() + ag * pars->GetNBCHANNEL() + ch;

        if((co != 31) && (co != 16))
        {
            auto xval {calman->ApplyLookUp(where, 4)};
            auto yval {calman->ApplyLookUp(where, 5)};
            for(int hit = 0, otherSize = evt->CoboAsad[it].peakheight.size(); hit < otherSize; hit++)
            {
                if(yval != -1)
                {
                    double z_position {evt->CoboAsad[it].peaktime[hit]};
                    if(z_position > 0.)
                    {
                        auto Qiaux {evt->CoboAsad[it].peakheight[hit]};
                        // Fill histogram
                        if(isMatched)
                            Qiaux = calman->ApplyPadAlignment(where, Qiaux);
                        h->Fill(where, Qiaux);
                    }
                }
            }
        }
    }
}

void ReadGainMatching(bool isMatched = false)
{
    // Get the data
    auto* chain {new TChain("ACTAR_TTree")};
    std::set<int> runs {4, 5, 6, 7, 8, 9, 10, 11, 12};
    for(const auto& run : runs)
        chain->Add(TString::Format("../RootFiles/Raw/Tree_Run_%04d_Merged.root", run));

    // Set the parameters
    ActRoot::TPCParameters tpc {"Actar"};
    // Calibration manager
    ActRoot::CalibrationManager calman {};
    calman.ReadLookUpTable("../Actar/LT.dat");
    if(isMatched)
        calman.ReadPadAlign("../Calibration/Outputs/gain_matching_v0.dat");

    // Create histogram
    auto* h {new TH2D {"h", (isMatched) ? "Matched pads;Channel;Q" : "Unmatched pads;Channel;Q", 17408, 0, 17408, 800,
                       0, 5000}};
    // Set MEventReduced
    MEventReduced* evt {new MEventReduced};
    chain->SetBranchAddress("data", &evt);
    // Run
    for(long int entry = 0, maxEntry = chain->GetEntries(); entry < maxEntry; entry++)
    {
        std::cout << "\r"
                  << "At entry : " << entry << std::flush;
        chain->GetEntry(entry);
        FillHistogram(&calman, &tpc, evt, h, isMatched);
    }

    // Plot
    auto* c0 {new TCanvas {"c0", "Gain matching canvas"}};
    h->Draw("colz");

    // Save histogram
    if(!isMatched)
        h->SaveAs("./Inputs/gain.root");
}
