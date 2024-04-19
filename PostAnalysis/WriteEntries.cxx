#ifndef WriteEntries_cxx
#define WriteEntries_cxx

#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActMergerData.h"

#include "ROOT/RDataFrame.hxx"

#include <fstream>
#include <string>

#include "./Utils.cxx"

void WriteEntries(const std::string& beam, const std::string& target, const std::string& light)
{
    // Read data
    // ActRoot::DataManager datman {"/home/eactar/Analysis_e837/Analysis/configs/data_post.conf",
    //                              ActRoot::ModeType::EMerge};
    // auto chain {datman.GetJoinedData()};
    // ROOT::RDataFrame d {*chain};
    ROOT::RDataFrame d {"PID_Tree", E837Utils::GetFileName(1, beam, target, light)};

    // Apply any filter function
    // auto df {d};
    // Read cut
    ActRoot::CutsManager<int> cut;
    cut.ReadCut(0, "./Cuts/Debug/low_angle.root");
    auto df {d.Filter(
        [&](const ActRoot::MergerData& m)
        {
            if(m.fSilEs.size() > 0)
                if(cut.IsInside(0, m.fThetaLight, m.fThetaHeavy))
                    return true;
            return false;
        },
        {"MergerData"})};
    // auto df {d.Filter("fThetaBeam > 1.5")}; // in deg

    // Write to file
    std::ofstream streamer {"./debug_angle.dat"};
    df.Foreach([&](const ActRoot::MergerData& data) { streamer << data.fRun << " " << data.fEntry << '\n'; },
               {"MergerData"});
    streamer.close();


    std::cout << "Written " << df.Count().GetValue() << " entries to file!" << '\n';
}
#endif
