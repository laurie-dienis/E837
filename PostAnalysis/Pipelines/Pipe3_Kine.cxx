#ifndef Pipe3_Kine_cxx
#define Pipe3_Kine_cxx

#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActParticle.h"
#include "ActSRIM.h"

#include "ROOT/RDataFrame.hxx"
#include "Rtypes.h"

#include "TCanvas.h"
#include "TROOT.h"
#include "TString.h"

#include <string>

#include "../HistConfig.h"
#include "../Utils.cxx"

void Pipe3_Kine(const std::string& beam, const std::string& target, const std::string& light, double e_beam_ini,
                int pressure)
{
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"PID_Tree", E837Utils::GetFileName(1, beam, target, light)};

    // Read srim tables
    auto* srim {new ActPhysics::SRIM};
    srim->ReadTable("light",
                    TString::Format("/home/dienis/Analysis_e837/Analysis/Inputs/SRIM/%s_He97_butane_%dmbar.txt",
                                    light.c_str(), pressure)
                        .Data());
    auto vetoed {df.Filter("fSilLayers.size() == 1")};
    // vetoed = vetoed.Define("ESil", "fSilEs.front()* (9.5 / 5860)");
    vetoed = vetoed.Define("EVertex",
                           [&](const ActRoot::MergerData& d, double esil)
                           { return srim->EvalInitialEnergy("light", esil, d.fTrackLength); },
                           {"MergerData", "ESil"});

    // Book histograms
    auto hKine {vetoed.Histo2D(HistConfig::Kine, "EVertex", "fThetaLight")};
    auto hAngle {vetoed.Histo2D(HistConfig::Angle, "fThetaLight", "fThetaHeavy")};

    // plotting
    auto* c30 {new TCanvas("c30", "Pipe3 canvas 0")};
    c30->DivideSquare(2);
    c30->cd(1);
    hKine->DrawClone("colz");
    c30->cd(2);
    ActPhysics::Particle pb {beam};
    ActPhysics::Kinematics kin {beam, target, light, e_beam_ini * pb.GetAMU()};
    auto* g {kin.GetTheta3vs4Line()};
    hAngle->DrawClone("colz");
    g->Draw("l");
}
#endif
