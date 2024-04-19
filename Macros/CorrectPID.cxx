#include "ActCutsManager.h"
#include "ActDataManager.h"

#include "ROOT/RDataFrame.hxx"

#include "TROOT.h"

#include "/home/eactar/Analysis_e837/Analysis/PostAnalysis/HistConfig.h"
void CorrectPID()
{
    ROOT::EnableImplicitMT();
    // Manager of  cuts
    ActRoot::CutsManager<std::string> cuts;
    // Read data
    ActRoot::DataManager datman {"/home/eactar/Analysis_e837/Analysis/configs/data_post.conf",
                                 ActRoot::ModeType::EMerge};
    auto chain {datman.GetJoinedData()};
    ROOT::RDataFrame d {*chain};
    auto df {d.Filter("fSilLayers.size() == 1")}; // && fThetaLight<20")};

    // Plot PID
    auto hPID {df.Histo2D(HistConfig::PID, "ESil", "fQave")};
}
