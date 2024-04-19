#include "TCanvas.h"
#include "TH1.h"

#include "CalibrationRunner.h"
#include "CalibrationSource.h"

#include <fstream>

#include "CalibrationRunner.cxx"
#include "CalibrationSource.cxx"

TH1D* GetData()
{
    auto* ret {new TH1D {"hRaw", "Raw data;Channel;Counts", 8192, 0, 8192}};
    // std::ifstream streamer {"/media/Data/SilCal/Inputs/Bea_triple_alfa_1h_1Marzo.txt"};
    std::ifstream streamer {"/media/Data/SilCal/Inputs/Ivan_triple_alfa_3h_20m_29Feb.txt"};
    int bin {1};
    int val {};
    while(streamer >> val)
    {
        ret->AddBinContent(bin, val);
        bin++;
    }
    streamer.close();
    return ret;
}

void actar()
{
    // Init source
    Calibration::Source source {"3AlphaUSC"};
    source.Print();

    // Get data
    auto* data {GetData()};
    // Rebin it for initial fit
    auto* rebinned {(TH1D*)data->Clone()};
    rebinned->SetTitle("Rebinned raw data");
    rebinned->Rebin(4);

    // Run!
    Calibration::Runner runner {&source, rebinned, data, false};
    runner.SetRange(3200, 4200);
    runner.SetGaussPreWidth(100);
    runner.DoIt();
    runner.Draw(new TCanvas);
    runner.PrintRes();

    // // Plotting
    // auto* c0 {new TCanvas {"c0", "Cal test"}};
    // data->Draw("hist");
    // for(auto* o : *data->GetListOfFunctions())
    //     if(o)
    //         o->Draw("same");
}
