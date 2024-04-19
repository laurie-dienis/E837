#include "TCanvas.h"
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

std::vector<TH1D*> ReadData(const std::string& file)
{
    auto* f {new TFile {file.c_str()}};
    // f->ls();
    auto updir {f->Get<TDirectory>("Raw")};
    // updir->ls();
    auto lowdir {updir->Get<TDirectory>("SI")};
    // lowdir->ls();
    // Read
    std::vector<TH1D*> ret;
    auto keys {lowdir->GetListOfKeys()};
    for(auto* key : *keys)
    {
        ret.push_back((TH1D*)lowdir->Get<TH1I>(key->GetName()));
    }
    return ret;
}
void ganil()
{
    auto hs {ReadData("./Inputs/Ganil/USC_Si_test_3624-5-2_3_200V_1p2uA.root")};

    // Create source
    Calibration::Source source;
    source.Print();

    auto* h3rebin {(TH1D*)hs[3]->Clone("h3Rebin")};
    h3rebin->Rebin(2);
    auto* h4rebin {(TH1D*)hs[4]->Clone("h4rebin")};
    h4rebin->RebinX(2);
    // Si_3
    Calibration::Runner runner3 {&source, h3rebin, hs[3]};
    // Si_4
    Calibration::Runner runner4 {&source, h4rebin, hs[4]};
    // Run for all
    for(auto& runner : {&runner3, &runner4})
    {
        runner->SetRange(5000, 7000);
        runner->DoIt();
        runner->Draw(new TCanvas);
        runner->PrintRes();
    }


    // Plot
    auto* c0 {new TCanvas {"c0", "Silicon canvas"}};
    c0->DivideSquare(hs.size());
    for(int i = 0; i < hs.size(); i++)
    {
        c0->cd(i + 1);
        hs[i]->Draw();
    }
}
