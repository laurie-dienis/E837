#include "ActDataManager.h"
#include "ActModularData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TROOT.h"

#include <atomic>

void Pipe0_Beam()
{
    ROOT::EnableImplicitMT();
    // Read data
    ActRoot::DataManager datman {"../../configs/data.conf", ActRoot::ModeType::EReadSilMod}; 
    auto chain {datman.GetJoinedData()};
    ROOT::RDataFrame df {*chain};
    // df.Describe().Print();

    // Get GATCONF values
    auto def {df.Define("GATCONF", [](ActRoot::ModularData& mod) { return mod.fLeaves["GATCONF"]; }, {"ModularData"})};

    // Book histograms
    auto hGATCONF {def.Histo1D("GATCONF")};

    // And cound CFA triggers
    std::atomic<unsigned long int> cfa {};
    def.Foreach(
        [&](ActRoot::ModularData& mod)
        {
            auto gatconf {static_cast<int>(mod.fLeaves["GATCONF"])};
            if(gatconf == 1) //change the value depending on the GATCONF
                cfa++;
        },
        {"ModularData"});

    // Draw
    auto* c0 {new TCanvas {"c00", "Pipe 0 canvas 0"}};
    hGATCONF->DrawClone();

    // Print report
    std::cout << "===== GATCONF report =====" << '\n';
    std::cout << "-> CFA/div = " << cfa << '\n';
    std::cout << "==========================" << '\n';
}
