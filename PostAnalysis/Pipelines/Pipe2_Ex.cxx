#ifndef Pipe2_Ex_cxx
#define Pipe2_Ex_cxx

#include "ActColors.h"
#include "ActCutsManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActParticle.h"
#include "ActSRIM.h"

#include "ROOT/RDataFrame.hxx"
#include "Rtypes.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TMath.h"
#include "TROOT.h"
#include "TString.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "../HistConfig.h"
#include "../Utils.cxx"

TH1D* GetProjectionX(TH2D* h, double xmin, double xmax, TString name = "projX_Ex")
{
    // Find bins
    int bin_low = h->GetYaxis()->FindBin(xmin);
    int bin_up = h->GetYaxis()->FindBin(xmax);

    auto p = h->ProjectionX(name, bin_low, bin_low);
    p->SetTitle(TString::Format("Proj #Theta [%.2f,%.2f]", xmin, xmax));
    return p;
}

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

double ComputeExcitationEnergy(double evertex, double ebeam, double theta, double mlight, double mheavy){
    // std::cout << "-> Evertex : " << evertex << '\n';
    // std::cout << "-> Ebeam : " << ebeam << '\n';
    // std::cout << "-> theta : " << theta << '\n';
    // std::cout << "-> mlight : " << mlight << '\n';
    // std::cout << "-> mheavy : " << mheavy << '\n';
    return evertex*((2*TMath::Cos(theta*TMath::DegToRad())*TMath::Sqrt((ebeam*mlight)/(evertex*mheavy)))-((mlight/mheavy)+1));
}

void Pipe2_Ex(const std::string& beam, const std::string& target, const std::string& light, double ebeam_i,
              int pressure)
{
    ROOT::EnableImplicitMT();
    // Read data
    ROOT::RDataFrame d {"PID_Tree", E837Utils::GetFileName(1, pressure, beam, target, light)};
    auto df = d.Filter("ESil>0");

    // Settings
    const double iter_threshold = 1; // keV
    double kin_particle_threshold {};
    // if(light == "4He")
    // {
    kin_particle_threshold = 8.956; // MeV
    // }
    // if(light == "6He")
    // {
    //     kin_particle_threshold = 8.956; // MeV
    // }
    // else
    // {
    //     throw std::runtime_error("Pipe 2 : cannot set threshold for particle " + light);
    // }
    // df.Describe().Print();

    // Build all the kinematics
    std::vector<ActPhysics::Kinematics> vkins;
    for(int s = 0; s < df.GetNSlots(); s++)
        vkins.push_back({beam, target, light, -1, 0});

    // Read srim tables
    auto* srim {new ActPhysics::SRIM};
    srim->ReadTable("light",
                    TString::Format("/home/dienis/Analysis_e837/Inputs/SRIM/%s_He97_butane_%dmbar.txt",
                                    light.c_str(), pressure)
                        .Data());
    srim->ReadTable("beam", TString::Format("/home/dienis/Analysis_e837/Inputs/SRIM/%s_He97_butane_%dmbar.txt",
                                            beam.c_str(), pressure)
                                .Data());

    ActPhysics::Particle pb {beam};
    ActPhysics::Particle pt {target};
    ActPhysics::Particle pl {light};
    double mass_beam = pb.GetAMU();
    double mass_target = pt.GetAMU();
    double mass_light = pl.GetAMU();
    double ebeam_i_MeV = ebeam_i * mass_beam;
    double range_beam_i = srim->EvalDirect("beam", ebeam_i_MeV);
    std::cout << "ebeam_i = " << ebeam_i_MeV << "MeV"
              << "\n";
    std::cout << "range_i = " << range_beam_i << "mm"
              << "\n";

    auto def = df.Define("EVertex",
                         [&](const ActRoot::MergerData& d, double esil)
                         {
                             // auto corr {30. / TMath::Cos(d.fThetaLight * TMath::DegToRad())};
                             // std::cout << "elight_ = " << srim->EvalInitialEnergy("light", esil, d.fTrackLength) <<
                             // "MeV\n";
                             return srim->EvalInitialEnergy("light", esil, d.fTrackLength);
                         },
                         {"MergerData", "ESil"});

    // Iterative process for determining Ereac and Elight
    def = def.Define("iter",
                     [&](const ActRoot::MergerData& merger, double esil, double esil_guess)
                     {
                         std::pair<double, double> ret;
                         double ereac, range_beam_f, dist_vertex, dist_sil, eloss, eloss_previous;
                         double elight = esil_guess;
                         // std::cout << "esilguess = " << elight << "MeV"
                         //           << "\n";
                         double theta = merger.fThetaLight * TMath::DegToRad();
                         double range_ini_sil = srim->EvalDirect("light", esil);
                         // std::cout << "------------------------------" << '\n';
                         // std::cout << "-> Theta : " << theta * TMath::RadToDeg() << '\n';
                         // std::cout << "-> ESil  : " << esil << '\n';
                         for(int i = 0; i < 15; i++)
                         {
                             // std::cout << "elight_" << i << " = " << elight << "MeV"
                             //         << "\n";
                             ereac = mass_beam * std::pow(((mass_light / mass_beam) + 1), 2) * elight /
                                     (4 * mass_light * std::pow(std::cos(theta), 2));
                             // std::cout << "ereac_" << i << " = " << ereac << "MeV"
                             //   << "\n";
                             // if (ereac>ebeam_i_MeV) {
                             //   elight = -11.;
                             //   ereac = -11;
                             //   break;
                             // }

                             //     ereac = ebeam_i_MeV;
                             // }
                             range_beam_f = srim->EvalDirect("beam", ereac);
                             // std::cout << "range_f_" << i << " = " << range_beam_f << "mm"
                             //   << "\n";
                             dist_vertex = range_beam_i - range_beam_f;
                             // if(dist_target > 350.)
                             // {
                             //     dist_target = 250.;
                             // }
                             dist_sil = (merger.fSP.X() - dist_vertex) / std::cos(theta);
                             eloss_previous = eloss;
                             auto range_vertex_sil = range_ini_sil + dist_sil;
                             elight = srim->EvalInverse("light", range_vertex_sil);
                             // Absolute value?
                             eloss = esil - elight;
                             // eloss = esil - srim->EvalInverse("light", dist_sil);
                             // std::cout << "eloss_" << i << " = " << eloss << "MeV"
                             //      << "\n";
                             // elight = esil + eloss;
                             // std::cout << "elight_" << i << " = " << elight << "MeV"
                             //          << "\n";
                             // std::cout << "elossdiff" << i << " = " << (eloss - eloss_previous) * 1000 << "MeV"
                             //        << "\n\n";
                             //  Print everything together
                             // std::cout << "-> Iter    : " << i << '\n';
                             // std::cout << "   Ereac   : " << ereac << '\n';
                             // std::cout << "   RBeamF  : " << range_beam_f << '\n';
                             // std::cout << "   DistVer : " << dist_vertex << '\n';
                             // std::cout << "   DistSil : " << dist_sil << '\n';
                             // std::cout << "   ESil    : " << esil << '\n';
                             // std::cout << "   ELoss   : " << eloss << '\n';
                             // std::cout << "   ELight  : " << elight << '\n';
                             // std::cout << '\n';

                             if(abs(eloss - eloss_previous) * 1000 < iter_threshold)
                                 break;
                         }
                         ret.first = elight;
                         ret.second = ereac;
                         return ret;
                     },
                     {"MergerData", "ESil", "EVertex"});

    // Ereac, energy of the beam at the reaction point
    def = def.Define("Ereac", [](const std::pair<double, double>& pair) { return pair.second; }, {"iter"});

    // Elight, energy of the light particle at the reaction point
    def = def.Define("Elight", [](const std::pair<double, double>& pair) { return pair.first; }, {"iter"});

    // Ex of compound nucleus
    def = def.Define(
        "Ex", [&](double ereac) { return (ereac * (mass_light / (mass_light + mass_beam))) + kin_particle_threshold; },
        {"Ereac"});

    // comparison with other method for reconstruction Ex

    def = def.Define("EBeam_range",
                     [&](const ActRoot::MergerData& d)
                     { return srim->Slow("beam", ebeam_i_MeV, d.fRP.X(), d.fThetaBeam * TMath::DegToRad()); },
                     {"MergerData"});

    def = def.DefineSlot("EBeam_Si",
        [&](unsigned int slot, double evertex, float theta)
        {
            return vkins[slot].ReconstructBeamEnergyFromLabKinematics(evertex, theta * TMath::DegToRad());
        },
        {"EVertex", "fThetaLight"});


    // ThetaCM
    def = def.DefineSlot("ThetaCM_range",
                         [&](unsigned int slot, double ereac, double evertex, float theta)
                         {
                             // return -1.;
                             if(std::isnan(ereac) || ereac < vkins[slot].GetT1Thresh())
                             {
                                 return -11.;
                             }
                             vkins[slot].SetBeamEnergy(ereac);
                             // std::cout << "ThetaCM = " << 180 - (vkins[slot].ReconstructTheta3CMFromLab(evertex,
                             // theta*TMath::DegToRad())*TMath::RadToDeg())<< "\n";
                             return 180 - (vkins[slot].ReconstructTheta3CMFromLab(evertex, theta * TMath::DegToRad()) *
                                           TMath::RadToDeg());
                         },
                         {"EBeam_range", "EVertex", "fThetaLight"});
    
    // ThetaCM
    def = def.DefineSlot("ThetaCM_Si",
                         [&](unsigned int slot, double ereac, double evertex, float theta)
                         {
                             // return -1.;
                             if(std::isnan(ereac) || ereac < vkins[slot].GetT1Thresh())
                             {
                                 return -11.;
                             }
                             vkins[slot].SetBeamEnergy(ereac);
                             // std::cout << "ThetaCM = " << 180 - (vkins[slot].ReconstructTheta3CMFromLab(evertex,
                             // theta*TMath::DegToRad())*TMath::RadToDeg())<< "\n";
                             return 180 - (vkins[slot].ReconstructTheta3CMFromLab(evertex, theta * TMath::DegToRad()) *
                                           TMath::RadToDeg());
                         },
                         {"EBeam_Si", "EVertex", "fThetaLight"});

    // Ex of 8He
    def = def.DefineSlot("Eex_range",
                         [&](unsigned int slot, const ActRoot::MergerData& d, double ebeam, double evertex)
                         {
                             if(std::isnan(ebeam) || ebeam < vkins[slot].GetT1Thresh())
                             {
                                 return -11.;
                             }
                             // std::cout << "E_4He = " << evertex << "\n";
                             // std::cout << "Theta_4He = " << d.fThetaLight << "\n";
                             // std::cout << "Ebeam from range = " << ebeam << "\n";
                             double Ex2 = ComputeExcitationEnergy(evertex, ebeam, d.fThetaLight, mass_light, mass_beam);
                             vkins[slot].SetBeamEnergy(ebeam);
                             // std::cout << "Ex from range = " << 
                             // vkins[slot].ReconstructExcitationEnergy(evertex, d.fThetaLight * TMath::DegToRad()) << "\n";
                             // std::cout << "Ex from range 2 = " << Ex2 << "\n";
                             //return vkins[slot].ReconstructExcitationEnergy(evertex, d.fThetaLight * TMath::DegToRad());
                             return Ex2;
                         },
                         {"MergerData", "EBeam_range", "EVertex"});


    // Ex of 8He
    def = def.DefineSlot("Eex_Si",
                         [&](unsigned int slot, const ActRoot::MergerData& d, double ebeam, double evertex)
                         {
                             if(std::isnan(ebeam) || ebeam < vkins[slot].GetT1Thresh())
                             {
                                 return -11.;
                             }
                             // std::cout << "Ebeam from Si = " << ebeam << "\n";
                             double Ex2 = ComputeExcitationEnergy(evertex, ebeam, d.fThetaLight, mass_light, mass_beam);
                             vkins[slot].SetBeamEnergy(ebeam);
                             // std::cout << "Ex from Si = " << 
                             // vkins[slot].ReconstructExcitationEnergy(evertex, d.fThetaLight * TMath::DegToRad()) << "\n";
                             // std::cout << "Ex from Si 2 = " << Ex2 << "\n";
                             // std::cout << '\n';
                             //return vkins[slot].ReconstructExcitationEnergy(evertex, d.fThetaLight * TMath::DegToRad());
                             return Ex2;
                         },
                         {"MergerData", "EBeam_Si", "EVertex"});

    //Ereac non interative
    def = def.Define(
        "Ereac_check_range",
        [&](double ebeam)
        {
            // std::cout << "ereac = " << vkins[slot].ReconstructBeamEnergyFromLabKinematics(evertex, theta *
            // TMath::DegToRad()) *
            //             (mass_target / (mass_target + mass_beam)) +
            //       kin_particle_threshold << " MeV" << "\n";
            return (ebeam * (mass_target / (mass_target + mass_beam)) + kin_particle_threshold);
        },
        {"EBeam_range"});

    //Ereac non interative
    def = def.Define(
        "Ereac_check_Si",
        [&](double ebeam)
        {
            // std::cout << "ereac = " << vkins[slot].ReconstructBeamEnergyFromLabKinematics(evertex, theta *
            // TMath::DegToRad()) *
            //             (mass_target / (mass_target + mass_beam)) +
            //       kin_particle_threshold << " MeV" << "\n";
            return (ebeam * (mass_target / (mass_target + mass_beam)) + kin_particle_threshold);
        },
        {"EBeam_Si"});

    //Ereac filtered
    // def = def.Define(
    //     "Ereac_check_filtered", 
    //     [&](double ereac, double Ex)
    //     {
    //         if (Ex<1 ) {
    //         return ereac;
    //         }
    //     },
    //     {"Ereac_check_Si","Eex_Si"});

    //debug
    // std::ofstream streamer {"./debug_Ex.dat"};
    // def.Foreach(
    //     [&](const ActRoot::MergerData& d, double Ex)
    //     {
    //         if(Ex>1)
    //             streamer << d.fRun << " " << d.fEntry << '\n';
    //     },
    //     {"MergerData", "Eex_Si"});
    // streamer.close();

    // // Write cuts to file
    // ActRoot::CutsManager<std::string> cuts;
    // cuts.ReadCut("debug", "./Cuts/Debug/kinematics.root");
    // std::ofstream streamer {"./debug_kinematics.dat"};
    // def.Foreach(
    //     [&](const ActRoot::MergerData& d, double evertex)
    //     {
    //         if(cuts.IsInside("debug", d.fThetaLight, evertex))
    //             streamer << d.fRun << " " << d.fEntry << '\n';
    //     },
    //     {"MergerData", "EVertex"});
    // streamer.close();
    // diff between the two methods
    // def = def.Define("Ereac_diff",
    //[&](double ereac, double ereac_check)
    //{ return std::abs(ereac-ereac_check); },
    //{"Ereac", "Ereac_check"});

    ActRoot::CutsManager<std::string> cuts;
    TString cutfile {};
    cutfile = TString::Format("/home/dienis/Analysis_e837/PostAnalysis/Cuts/Debug/CorrectEx_%dmbar.root", pressure);
    cuts.ReadCut("Ebeamcorrect", cutfile);
    
    cuts.ReadCut("weird", "./Cuts/Debug/false_resonance.root");
    std::ofstream streamer {"./debug_false_resonances.dat"};
    def.Foreach(
        [&](const ActRoot::MergerData& d, double Ereac_check_Si, double ThetaCM_Si)
        {
            if(cuts.IsInside("weird", Ereac_check_Si, ThetaCM_Si))
                // if(cut.IsInside("debug", d.fSP.Y(), d.fSP.Z()))
                streamer << d.fRun << " " << d.fEntry << '\n';
        },
        {"MergerData","Ereac_check_Si","ThetaCM_Si"});
    streamer.close();
    
    // def = def.Define(
    //     "Ereac_check_filtered", 
    //     [&](double ereac, double EBeam_range, double EBeam_Si)
    //     {
    //         if (cuts.IsInside("Ebeamcorrect", EBeam_range, EBeam_Si)) {
    //         return ereac;
    //         }
    //     },
    //     {"Ereac_check_range","EBeam_range","EBeam_Si"});
    // Filter
    auto vetoed {def.Filter([&](double EBeam_range, double EBeam_Si)
            { return cuts.IsInside("Ebeamcorrect", EBeam_range, EBeam_Si); },
            {"EBeam_range", "EBeam_Si"})};
    
    // Book histograms
    auto hEreac {def.Histo1D(HistConfig::Ereac, "EBeam_Si")};
    auto hTheta {d.Histo1D("fThetaLight")};
    auto hElight {def.Histo1D(HistConfig::Elight, "Elight")};
    auto hEex_range {def.Histo1D(HistConfig::Ex, "Eex_range")};
    auto hEex_Si {def.Histo1D(HistConfig::Ex2, "Eex_Si")};
    // auto hEdiff{def.Histo1D(HistConfig::Ediff, "Ereac_diff")};
    auto hExlab {def.Histo2D(HistConfig::Ex21Nalab, "fThetaLight", "Ex")};
    auto hEx {def.Histo2D(HistConfig::Ex21Na, "ThetaCM_Si", "Ex")};
    auto hEx2_range {def.Histo2D(HistConfig::Ex21Na2, "Ereac_check_range", "ThetaCM_range")};
    auto hEx2_Si {vetoed.Histo2D(HistConfig::Ex21Na2, "Ereac_check_Si", "ThetaCM_Si")};
    auto hProj = GetProjectionX(hEx2_Si.GetPtr(), 167, 176);
    hProj->Rebin(3);
    auto hKin {def.Histo2D(HistConfig::KinEl, "fThetaLight", "EVertex")};
    auto hEex_Si_vs_range {def.Histo2D(HistConfig::EexSiVsRange, "Eex_range","Eex_Si")};
    auto hEbeam_Si_vs_range {def.Histo2D(HistConfig::EbeamSiVsRange, "EBeam_range","EBeam_Si")};

    // plotting
    //gStyle->SetPalette(kCool);
    set_plot_style();
    auto* c20 {new TCanvas("c20", "Pipe2 canvas 0")};
    c20->DivideSquare(6);
    c20->cd(1);
    hEex_range->DrawClone("colz");
    c20->cd(2);
    hEx2_range->DrawClone("colz");
    c20->cd(3);
    hEex_Si_vs_range->DrawClone("colz");
    c20->cd(4);
    hEex_Si->DrawClone("colz");
    c20->cd(5);
    hEx2_Si->DrawClone("colz");
    cuts.DrawAll();
    c20->cd(6);
    hEbeam_Si_vs_range->DrawClone("colz");
    cuts.DrawAll();
    TLine *line = new TLine(5.,5.,13.,13.);
    line->SetLineColor(kOrange);
    line->SetLineWidth(2);
    line->Draw("same");


    // c20->cd(1);
    // hEreac->DrawClone("colz");
    // c20->cd(2);
    // hEex->DrawClone("colz");
    // hProj->DrawClone("colz");
    // c20->cd(3);
    // hEx2->DrawClone("colz");
    // c20->cd(4);
    // hKin->DrawClone("colz");
    // hEx->DrawClone("colz");
    // Draw
    // TGraph* g {};
    // if(light != "1H")
    // {
    //     vkins[0].SetBeamEnergyAndEx(ebeam_i_MeV, 0);
    //     g = vkins[0].GetKinematicLine3();
    //     g->Draw("l");
    // }
    // // cuts.DrawAll();
    // c20->cd(5);
    // hProj->DrawClone("colz");
}
#endif
