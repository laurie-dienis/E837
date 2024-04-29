#ifndef Pipe3_Kine_cxx
#define Pipe3_Kine_cxx

#include "ActColors.h"
#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActSilMatrix.h"
#include "ActTPCData.h"
#include "ActTypes.h"
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

void set_plot_style()
{
    const Int_t NRGBs = 5;
    const Int_t NCont = 455;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}

void Pipe3_Kine(const std::string& beam, const std::string& target, const std::string& light, double e_beam_ini,
                int pressure)
{
    ROOT::EnableImplicitMT();
    //ROOT::RDataFrame df {"PID_Tree", E837Utils::GetFileName(1, pressure, beam, target, light)};

    // Read data
    ActRoot::DataManager datman {"/home/dienis/Analysis_e837/configs/data_post.conf",
                                 ActRoot::ModeType::EMerge};
    auto chain {datman.GetJoinedData()};
    auto chainTpc {datman.GetJoinedData(ActRoot::ModeType::EReadTPC)};
    chain->AddFriend(chainTpc.get());
    ROOT::RDataFrame d {*chain};
    
    // Read srim tables
    auto* srim {new ActPhysics::SRIM};
    srim->ReadTable("light",
                    TString::Format("/home/dienis/Analysis_e837/Inputs/SRIM/%s_He97_butane_%dmbar.txt",
                                    light.c_str(), pressure)
                        .Data());
    auto vetoed {d.Filter("fSilLayers.size() == 1")};
    // vetoed = vetoed.Define("ESil", "fSilEs.front()* (9.5 / 5860)");
    // Define ESil as alias of .front of energy vector
    vetoed = vetoed.Define("ESil", "(double)fSilEs.front()");
    vetoed = vetoed.Define("EVertex",
                           [&](const ActRoot::MergerData& d, double esil)
                           { return srim->EvalInitialEnergy("light", esil, d.fTrackLength); },
                           {"MergerData", "ESil"});

    
    // Manager of  cuts
    ActRoot::CutsManager<std::string> cuts;
    TString pidfile {};
    pidfile = TString::Format("./Cuts/LightPID/pid_%s_%dmbar.root", light.c_str(), pressure);
    cuts.ReadCut(light, pidfile);
    TString carbonfile {};
    carbonfile = TString::Format("./Cuts/Debug/carbon_%dmbar.root", pressure);
    cuts.ReadCut("carbon", carbonfile);
    std::cout << BOLDCYAN << "Reading light PID in : " << pidfile << RESET << '\n';

    auto vetoed_light {vetoed.Filter([&](const ActRoot::MergerData& d)
            { return cuts.IsInside(light, d.fSilEs.front(), d.fQave); },
            {"MergerData"})};

    auto vetoed_carbon {vetoed.Filter([&](const ActRoot::MergerData& d)
            { return cuts.IsInside("carbon", d.fSilEs.front(), d.fQave); },
            {"MergerData"})};

    // Book histograms
    auto hKine {vetoed_carbon.Histo2D(HistConfig::Kine, "fThetaLight", "EVertex")};
    auto hAngle {vetoed_carbon.Histo2D(HistConfig::Angle, "fThetaHeavy", "fThetaLight")};
    auto hPID {vetoed_carbon.Histo2D(HistConfig::PID, "ESil", "fQave")};

    // plotting
    set_plot_style();
    auto* c30 {new TCanvas("c30", "Pipe3 canvas 0")};
    c30->DivideSquare(3);

    c30->cd(1);
    hKine->DrawClone("colz");
    ActPhysics::Particle pb1 {beam};
    ActPhysics::Kinematics kin1 {beam, target, light, e_beam_ini * pb1.GetAMU()};
    auto* g1 {kin1.GetKinematicLine3()};
    g1->SetLineColor(kOrange+6);
    g1->Draw("l");
    ActPhysics::Kinematics kin2 {beam,"12C","8He", e_beam_ini * pb1.GetAMU()};
    auto* g2 {kin2.GetKinematicLine3()};
    g2->SetLineColor(kViolet-4);
    g2->Draw("l");
    ActPhysics::Kinematics kin3 {beam,"4He","6He", e_beam_ini * pb1.GetAMU()};
    auto* g5 {kin3.GetKinematicLine3()};
    g5->SetLineColor(kAzure+10);
    g5->Draw("l");

    TLegend *legend = new TLegend(0.1,0.75,0.3,0.9);
    legend->AddEntry(g1, "4He","l");
    legend->AddEntry(g2, "8He","l");
    legend->Draw("l");

    c30->cd(2);
    ActPhysics::Particle pb {beam};
    ActPhysics::Kinematics kin {beam, target, light, e_beam_ini * pb.GetAMU()};
    auto* g3 {kin.GetTheta3vs4Line()};
    std::cout<<"type of g"<<typeid(g3).name();
    g3->SetLineColor(kOrange+6);
    auto* g4 {kin2.GetTheta3vs4Line()};
    g4->SetLineColor(kViolet-4);
    auto* g6 {kin3.GetTheta3vs4Line()};
    g6->SetLineColor(kAzure+10);
    hAngle->DrawClone("colz");
    g3->Draw("l");
    g4->Draw("l");
    g6->Draw("l");

    c30->cd(3);
    hPID->DrawClone();
}
#endif
