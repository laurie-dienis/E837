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

TH1D *GetProjectionX(TH2D *h, double xmin, double xmax,
                     TString name = "projX_Ex") {
  // Find bins
  int bin_low = h->GetYaxis()->FindBin(xmin);
  int bin_up = h->GetYaxis()->FindBin(xmax);

  auto p = h->ProjectionX(name, bin_low, bin_low);
  p->SetTitle(TString::Format("Proj #Theta [%.2f,%.2f]", xmin, xmax));
  return p;
}

void set_plot_style() {
  const Int_t NRGBs = 5;
  const Int_t NCont = 455;

  Double_t stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
  Double_t red[NRGBs] = {0.00, 0.00, 0.87, 1.00, 0.51};
  Double_t green[NRGBs] = {0.00, 0.81, 1.00, 0.20, 0.00};
  Double_t blue[NRGBs] = {0.51, 1.00, 0.12, 0.00, 0.00};
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}

double ComputeExcitationEnergy_elastic(double evertex, double ebeam,
                                       double theta, double mlight,
                                       double mheavy) {
  // std::cout << "-> Evertex : " << evertex << '\n';
  // std::cout << "-> Ebeam : " << ebeam << '\n';
  // std::cout << "-> theta : " << theta << '\n';
  // std::cout << "-> mlight : " << mlight << '\n';
  // std::cout << "-> mheavy : " << mheavy << '\n';
  return evertex * ((2 * TMath::Cos(theta * TMath::DegToRad()) *
                     TMath::Sqrt((ebeam * mlight) / (evertex * mheavy))) -
                    ((mlight / mheavy) + 1));
}

double ComputeExcitationEnergy(double evertex, double ebeam, double theta,
                               double mlight, double mbeam, double mtarget) {
  // std::cout << "-> Evertex : " << evertex << '\n';
  // std::cout << "-> Ebeam : " << ebeam << '\n';
  // std::cout << "-> theta : " << theta << '\n';
  // std::cout << "-> mlight : " << mlight << '\n';
  // std::cout << "-> mbeam : " << mbeam << '\n';
  double delta_m = 2 * mlight - mbeam - mtarget;
  // std::cout << "-> delta_m : " << delta_m << '\n';
  double Ex = ebeam * (1 - (mbeam / mlight)) +
              evertex * ((2 * TMath::Cos(theta * TMath::DegToRad()) *
                          TMath::Sqrt((ebeam * mbeam) / (evertex * mlight))) -
                         ((mlight / mlight) + 1)) +
              delta_m;
  // std::cout << "-> Ex : " << Ex << '\n';
  return Ex;
}

void Pipe2_Ex(const std::string &beam, const std::string &target,
              const std::string &light, double ebeam_i, int pressure) {
  // ROOT::EnableImplicitMT();
  //  Read data
  //  Open the main data file and create a RDataFrame
  ROOT::RDataFrame d{"PID_Tree",
                     E837Utils::GetFileName(1, pressure, beam, target, light)};

  // Apply an initial filter
  auto df = d.Filter("ESil>0");

  // Step 1: Load the labels from classification.root
  ROOT::RDataFrame label_df("Class_Tree", "./Input/classification_900.root");

  // Extract the labels into a vector
  std::vector<long long> labels = *label_df.Take<long long>("label");

  // Filter the indices where label == 1
  std::vector<int> selected_indices;
  for (int i = 0; i < labels.size(); ++i) {
    if (labels[i] == 1 || labels[i] == 0) {
      selected_indices.push_back(i);
    }
  }

  // Step 2: Add an entry index column to the dataframe
  auto df_with_index =
      df.Define("entryIndex", [](int entry) { return entry; }, {"fEntry"});

  // Step 3: Define a lambda function for filtering
  auto filter_func = [selected_indices](int entryIndex) {
    return std::find(selected_indices.begin(), selected_indices.end(),
                     entryIndex) != selected_indices.end();
  };

  // Apply the filter using the entry index
  auto filtered_df = df_with_index.Filter(filter_func, {"entryIndex"});

  std::ofstream streamer1{"./6He.dat"};
  df.Foreach(
      [&](const ActRoot::MergerData &d) {
        streamer1 << d.fRun << " " << d.fEntry << '\n';
      },
      {"MergerData"});
  streamer1.close();

  // Appy normalization with simulation ?
  bool normalization{1};
  // Settings
  double angle_min{134.};
  double angle_max{144.};
  const double iter_threshold = 1; // keV
  double kin_particle_threshold{};
  kin_particle_threshold = 8.956; // MeV
  // df.Describe().Print();

  // Build all the kinematics
  std::vector<ActPhysics::Kinematics> vkins;
  for (int s = 0; s < df.GetNSlots(); s++)
    vkins.push_back({beam, target, light, -1, 0});

  // Read srim tables
  auto *srim{new ActPhysics::SRIM};
  srim->ReadTable(
      "light",
      TString::Format(
          "/home/laurie/Analysis_e837/Inputs/SRIM/%s_He97_butane_%dmbar.txt",
          light.c_str(), pressure)
          .Data());
  srim->ReadTable(
      "beam",
      TString::Format(
          "/home/laurie/Analysis_e837/Inputs/SRIM/%s_He97_butane_%dmbar.txt",
          beam.c_str(), pressure)
          .Data());

  ActPhysics::Particle pb{beam};
  ActPhysics::Particle pt{target};
  ActPhysics::Particle pl{light};
  double mass_beam = pb.GetAMU();
  double mass_target = pt.GetAMU();
  double mass_light = pl.GetAMU();
  double ebeam_i_MeV = ebeam_i * mass_beam;
  double range_beam_i = srim->EvalDirect("beam", ebeam_i_MeV);
  std::cout << "ebeam_i = " << ebeam_i_MeV << "MeV"
            << "\n";
  std::cout << "range_i = " << range_beam_i << "mm"
            << "\n";

  auto def =
      df.Define("EVertex",
                [&](const ActRoot::MergerData &d, double esil) {
                  // auto corr {30. / TMath::Cos(d.fThetaLight *
                  // TMath::DegToRad())}; std::cout << "elight_ = " <<
                  // srim->EvalInitialEnergy("light", esil, d.fTrackLength) <<
                  // "MeV\n";
                  // std::cout << "-> esil : " << esil << '\n';
                  // std::cout << "-> evertex : " <<
                  // srim->EvalInitialEnergy("light", esil, d.fTrackLength) <<
                  // '\n';
                  return srim->EvalInitialEnergy("light", esil, d.fTrackLength);
                },
                {"MergerData", "ESil"});

  // Iterative process for determining Ereac and Elight
  // def = def.Define(
  //     "iter",
  //     [&](const ActRoot::MergerData &merger, double esil, double esil_guess)
  //     {
  //       std::pair<double, double> ret;
  //       double ereac, range_beam_f, dist_vertex, dist_sil, eloss,
  //           eloss_previous;
  //       double elight = esil_guess;
  //       std::cout << "esilguess = " << elight << "MeV"
  //                 << "\n";
  //       double theta = merger.fThetaLight * TMath::DegToRad();
  //       double range_ini_sil = srim->EvalDirect("light", esil);
  //       // std::cout << "------------------------------" << '\n';
  //       // std::cout << "-> Theta : " << theta * TMath::RadToDeg() << '\n';
  //       // std::cout << "-> ESil  : " << esil << '\n';
  //       for (int i = 0; i < 10; i++) {
  //          std::cout << "elight_" << i << " = " << elight << "MeV"
  //                  << "\n";
  //         ereac = mass_beam * std::pow(((mass_light / mass_beam) + 1), 2) *
  //                 elight / (4 * mass_light * std::pow(std::cos(theta), 2));
  //         // std::cout << "ereac_" << i << " = " << ereac << "MeV"
  //         //   << "\n";
  //         if (ereac>ebeam_i_MeV) {
  //           elight = 1.;
  //           ereac = 1;
  //           break;
  //         }
  //         //     ereac = ebeam_i_MeV;
  //         // }
  //         range_beam_f = srim->EvalDirect("beam", ereac);
  //         // std::cout << "range_f_" << i << " = " << range_beam_f << "mm"
  //         //   << "\n";
  //         dist_vertex = range_beam_i - range_beam_f;
  //         if (dist_vertex<-350){
  //           dist_vertex =250.;
  //         }
  //         if(dist_vertex > 350.)
  //         {
  //             dist_vertex = 250.;
  //          }
  //         dist_sil = (merger.fSP.X() - dist_vertex) / std::cos(theta);
  //         eloss_previous = eloss;
  //         auto range_vertex_sil = range_ini_sil + dist_sil;
  //         double elight_new = srim->EvalInverse("light", range_vertex_sil);
  //         // Absolute value?
  //         eloss = esil - elight;
  //         double alpha = 0.05;
  //         elight = alpha * elight_new + (1 - alpha) * elight;
  //         // eloss = esil - srim->EvalInverse("light", dist_sil);
  //         // std::cout << "eloss_" << i << " = " << eloss << "MeV"
  //         //      << "\n";
  //         // elight = esil + eloss;
  //         // std::cout << "elight_" << i << " = " << elight << "MeV"
  //         //          << "\n";
  //         // std::cout << "elossdiff" << i << " = " << (eloss -
  //         eloss_previous)
  //         // * 1000 << "MeV"
  //         //        << "\n\n";
  //         //  Print everything together
  //         std::cout << "\n" << "-> Iter    : " << i << '\n';
  //         std::cout << "   Ereac   : " << ereac << '\n';
  //         std::cout << "   RBeamF  : " << range_beam_f << '\n';
  //         std::cout << "   DistVer : " << dist_vertex << '\n';
  //         std::cout << "   DistSil : " << dist_sil << '\n';
  //         std::cout << "   ESil    : " << esil << '\n';
  //         std::cout << "   ELoss   : " << eloss << '\n';
  //         std::cout << "   ELight  : " << elight << '\n';
  //         std::cout << '\n';

  //         if (abs(eloss - eloss_previous) * 1000 < iter_threshold)
  //           break;
  //       }
  //       ret.first = elight;
  //       ret.second = ereac;
  //       return ret;
  //     },
  //     {"MergerData", "ESil", "EVertex"});

  // Iterative process for determining Ereac and Elight
  def = def.Define(
      "iter",
      [&](const ActRoot::MergerData &merger, double esil, double esil_guess) {
        std::pair<double, double> ret;
        double range_beam_i = srim->EvalDirect("beam", ebeam_i_MeV);
        std::cout << "range_i = " << range_beam_i << "mm"
                  << "\n";
        double ereac, range_beam_f, dist_target, dist_sil, eloss,
            eloss_previous;
        double elight = esil_guess;
        double theta = std::abs(merger.fThetaLight) * TMath::DegToRad();
        double range_ini_silicon = srim->EvalDirect("light", esil);
        std::cout << "theta = " << theta * TMath::RadToDeg() << "deg"
                  << "\n";
        for (int i = 0; i < 200; i++) {
          std::cout << "elight_" << i << " = " << elight << "MeV"
                    << "\n";
          ereac = mass_beam * std::pow(((mass_light / mass_beam) + 1), 2) *
                  elight / (4 * mass_light * std::pow(std::cos(theta), 2));
          std::cout << "ereac_" << i << " = " << ereac << "MeV"
                    << "\n";
          if (ereac < 3.457) {
            ereac = 3.5;
            break;
          }
          range_beam_f = srim->EvalDirect("beam", ereac);
          std::cout << "range_f_" << i << " = " << range_beam_f << "mm"
                    << "\n";
          dist_target = range_beam_i - range_beam_f;
          std::cout << "dist_target_" << i << " = " << dist_target << "mm"
                    << "\n";
          if (dist_target < -350) {
            dist_target = 250.;
          }
          if (dist_target > 350.) {
            dist_target = 250.;
          }
          dist_sil = (merger.fSP.X() - dist_target) / std::cos(theta);
          eloss_previous = eloss;
          // eloss = esil - srim->EvalInverse("light", dist_sil);
          std::cout << "eloss_" << i << " = " << eloss << "MeV"
                    << "\n";
          double range_vertex_silicon = range_ini_silicon + dist_sil;
          double elight_new = srim->EvalInverse("light", range_vertex_silicon);
          eloss = esil - elight;
          double alpha = 0.05;
          elight = alpha * elight_new + (1 - alpha) * elight;
          std::cout << "elight_" << i << " = " << elight << "MeV"
                    << "\n";
          std::cout << "elossdiff" << i << " = "
                    << (eloss - eloss_previous) * 1000 << "MeV"
                    << "\n\n";
          if (abs(eloss - eloss_previous) * 1000 < iter_threshold)
            break;
        }
        std::cout << "\n\n\n";
        ret.first = elight;
        ret.second = ereac;
        return ret;
      },
      {"MergerData", "ESil", "EVertex"});

  // Ereac, energy of the beam at the reaction point
  def = def.Define(
      "Ereac",
      [](const std::pair<double, double> &pair) { return pair.second; },
      {"iter"});

  // Elight, energy of the light particle at the reaction point
  def = def.Define(
      "Elight",
      [](const std::pair<double, double> &pair) { return pair.first; },
      {"iter"});

  // Ex of compound nucleus
  def = def.Define("Ex",
                   [&](double ereac) {
                     // std::cout << "ex = " << (ereac * (mass_light /
                     // (mass_light + mass_beam))) +
                     //       kin_particle_threshold << "MeV";
                     return (ereac * (mass_light / (mass_light + mass_beam))) +
                            kin_particle_threshold;
                   },
                   {"Ereac"});

  // comparison with other method for reconstruction Ex

  def = def.Define("EBeam_range",
                   [&](const ActRoot::MergerData &d) {
                     //  std::cout << "-> ebeam : " << srim->Slow("beam",
                     //  ebeam_i_MeV, d.fRP.X(),
                     //                    d.fThetaBeam * TMath::DegToRad())  <<
                     //                    '\n';
                     return srim->Slow("beam", ebeam_i_MeV, d.fRP.X(),
                                       d.fThetaBeam * TMath::DegToRad());
                   },
                   {"MergerData"});

  def = def.Define("Angle_heavy",
                   [&](const ActRoot::MergerData &d) {
                     //  std::cout << "-> ebeam : " << srim->Slow("beam",
                     //  ebeam_i_MeV, d.fRP.X(),
                     //                    d.fThetaBeam * TMath::DegToRad())  <<
                     //                    '\n';
                     return d.fThetaHeavy;
                   },
                   {"MergerData"});

  def = def.Define("Angle_light",
                   [&](const ActRoot::MergerData &d) {
                     //  std::cout << "-> ebeam : " << srim->Slow("beam",
                     //  ebeam_i_MeV, d.fRP.X(),
                     //                    d.fThetaBeam * TMath::DegToRad())  <<
                     //                    '\n';
                     return d.fThetaLight;
                   },
                   {"MergerData"});

  def = def.DefineSlot(
      "EVertex_calc",
      [&](unsigned int slot, double ebeam, float theta) {
        return ebeam * 2 * 4 * TMath::Cos(theta * TMath::DegToRad()) *
               TMath::Cos(theta * TMath::DegToRad()) /
               ((1. + mass_beam / mass_light) * (1. + mass_beam / mass_light));
      },
      {"EBeam_range", "fThetaLight"});

  def = def.DefineSlot(
      "EBeam_Si",
      [&](unsigned int slot, double evertex, float theta) {
        return vkins[slot].ReconstructBeamEnergyFromLabKinematics(
            evertex, theta * TMath::DegToRad());
      },
      {"EVertex", "fThetaLight"});

  // ThetaCM
  def = def.DefineSlot(
      "ThetaCM",
      [&](unsigned int slot, double ereac, double evertex, float theta) {
        if (std::isnan(ereac)) {
          return -11.;
        }
        vkins[slot].SetBeamEnergy(ereac);
        return (vkins[slot].ReconstructTheta3CMFromLab(
                    evertex, theta * TMath::DegToRad()) *
                TMath::RadToDeg());
      },
      {"Ereac", "Elight", "fThetaLight"});

  // // ThetaCM
  def = def.DefineSlot(
      "ThetaCM_range",
      [&](unsigned int slot, double ereac, double evertex, float theta) {
        // return -1.;
        if (std::isnan(ereac) || ereac < vkins[slot].GetT1Thresh()) {
          return -11.;
        }
        vkins[slot].SetBeamEnergy(ereac);
        // std::cout << "ThetaCM = " << 180 -
        // (vkins[slot].ReconstructTheta3CMFromLab(evertex,
        // theta*TMath::DegToRad())*TMath::RadToDeg())<< "\n";
        return (vkins[slot].ReconstructTheta3CMFromLab(
                    evertex, theta * TMath::DegToRad()) *
                TMath::RadToDeg());
      },
      {"EBeam_range", "EVertex", "fThetaLight"});

  // ThetaCM
  def = def.DefineSlot(
      "ThetaCM_Si",
      [&](unsigned int slot, double ereac, double evertex, float theta) {
        // return -1.;
        if (std::isnan(ereac) || ereac < vkins[slot].GetT1Thresh()) {
          return -11.;
        }
        vkins[slot].SetBeamEnergy(ereac);
        // std::cout << "ThetaCM = " << 180 -
        // (vkins[slot].ReconstructTheta3CMFromLab(evertex,
        // theta*TMath::DegToRad())*TMath::RadToDeg())<< "\n";
        return (vkins[slot].ReconstructTheta3CMFromLab(
                    evertex, theta * TMath::DegToRad()) *
                TMath::RadToDeg());
      },
      {"EBeam_Si", "EVertex", "fThetaLight"});

  // Ex of 8He
  def = def.DefineSlot(
      "Eex_range",
      [&](unsigned int slot, const ActRoot::MergerData &d, double ebeam,
          double evertex) {
        if (std::isnan(ebeam)) {
          return -11.;
        }
        double Ex2;
        // std::cout << "Theta_4He = " << d.fThetaLight << "\n";
        //  if (light.c_str()=="4He")
        if (strncmp(light.c_str(), "4He", strlen("4He")) == 0) {
          // std::cout << "Evertex = " << evertex << "\n";
          Ex2 = ComputeExcitationEnergy_elastic(evertex, ebeam, d.fThetaLight,
                                                mass_light, mass_beam);
        }
        // std::cout << "Ex = " << Ex2 << "\n";
        //  if (light.c_str()=="6He"){
        if (strncmp(light.c_str(), "6He", strlen("6He")) == 0) {
          // Ex2 = ComputeExcitationEnergy_elastic(evertex, ebeam,
          // d.fThetaLight,
          //                                      mass_light, mass_beam);
          Ex2 = ComputeExcitationEnergy(evertex, ebeam, d.fThetaLight,
                                        mass_light, mass_beam, mass_target);
        }
        // vkins[slot].SetBeamEnergy(ebeam);
        //  std::cout << "Ex from range = " <<
        //  vkins[slot].ReconstructExcitationEnergy(evertex, d.fThetaLight *
        //  TMath::DegToRad()) << "\n"; std::cout << "Ex from range 2 = " << Ex2
        //  << "\n";￼tionEnergy(evertex,
        //  d.fThetaLight
        //  * TMath::DegToRad());
        return Ex2;
      },
      {"MergerData", "EBeam_range", "EVertex"});

  // Ex of 8He
  def = def.DefineSlot(
      "Eex_Si",
      [&](unsigned int slot, const ActRoot::MergerData &d, double ebeam,
          double evertex) {
        if (std::isnan(ebeam) || ebeam < vkins[slot].GetT1Thresh()) {
          return -11.;
        }
        // std::cout << "Ebeam from Si = " << ebeam << "\n";
        double Ex2 = ComputeExcitationEnergy_elastic(
            evertex, ebeam, d.fThetaLight, mass_light, mass_beam);
        vkins[slot].SetBeamEnergy(ebeam);
        // std::cout << "Ex from Si = " <<
        // vkins[slot].ReconstructExcitationEnergy(evertex, d.fThetaLight *
        // TMath::DegToRad()) << "\n"; std::cout << "Ex from Si 2 = " << Ex2 <<
        // "\n"; std::cout << '\n';
        // return vkins[slot].ReconstructExcitationEnergy(evertex, d.fThetaLight
        // * TMath::DegToRad());
        return Ex2;
      },
      {"MergerData", "EBeam_Si", "EVertex"});

  // Ereac non interative
  def = def.Define("Ereac_check_range",
                   [&](double ebeam, double Ex) {
                     // std::cout << "ereac = " <<
                     // vkins[slot].ReconstructBeamEnergyFromLabKinematics(evertex,
                     // theta * TMath::DegToRad()) *
                     //             (mass_target / (mass_target + mass_beam)) +
                     //       kin_particle_threshold << " MeV" << "\n";
                     if (Ex>-2){
                     return (ebeam * (mass_target / (mass_target + mass_beam)) +
                             kin_particle_threshold);
                     }
                     else {
                    return -10.;
                     }
                   },
                   {"EBeam_range","Eex_range"});

  // Ereac non interative
  def = def.Define("Ereac_check_Si",
                   [&](double ebeam) {
                     // std::cout << "ereac = " <<
                     // vkins[slot].ReconstructBeamEnergyFromLabKinematics(evertex,
                     // theta * TMath::DegToRad()) *
                     //             (mass_target / (mass_target + mass_beam)) +
                     //       kin_particle_threshold << " MeV" << "\n";
                     return (ebeam * (mass_target / (mass_target + mass_beam)) +
                             kin_particle_threshold);
                   },
                   {"EBeam_Si"});

  def = def.Define("Ex_proj",
                   [&](double theta3CM, double ex) {
                     if (theta3CM >= angle_min && theta3CM <= angle_max) {
                       return ex;
                     } else {
                       return 0.;
                     }
                   },
                   {"ThetaCM_range", "Ereac_check_range"});

  // normalization of the excitation energy distribution
  TFile *norm = TFile::Open(
      TString::Format("./Input/Norm_%s_%dmbar_%.0f-%.0fdeg.root", light.c_str(),
                      pressure, angle_min, angle_max),
      "READ");
  if (!norm || norm->IsZombie()) {
    // Handle error if the file cannot be opened
    std::cerr << "Error: Could not open file" << std::endl;
    return;
  }
  // Retrieve the TH2D histogram named "hnorm" from the file
  TH2D *hnorm = dynamic_cast<TH2D *>(norm->Get("hnorm"));
  if (!hnorm) {
    // Handle error if the histogram cannot be found or if it's of the wrong
    // type
    std::cerr << "Error: Histogram 'hnorm' not found or is not a TH2D"
              << std::endl;
    // Close the file before exiting
    norm->Close();
    return;
  }
  // Retrieve the TH1D histogram named "hnorm_proj" from the file
  TH1D *hnorm_proj = dynamic_cast<TH1D *>(norm->Get("hprojection"));
  if (!hnorm_proj) {
    // Handle error if the histogram cannot be found or if it's of the wrong
    // type
    std::cerr << "Error: Histogram 'hnorm_proj' not found or is not a TH1D"
              << std::endl;
    // Close the file before exiting
    norm->Close();
    return;
  }

  ActRoot::CutsManager<std::string> cuts;
  TString cutfile{};
  cutfile = TString::Format("Cuts/"
                            "Debug/CorrectEx_%dmbar.root",
                            pressure);
  cuts.ReadCut("Ebeamcorrect", cutfile);

  TString cutfile2{};
  cutfile2 = TString::Format("Cuts/"
                             "Debug/Secondcleaning2_%dmbar.root",
                             pressure);
  cuts.ReadCut("Cleaning2", cutfile2);

  cuts.ReadCut("weird", "./Cuts/Debug/false_resonance.root");
  std::ofstream streamer{"./debug_false_resonances.dat"};
  def.Foreach(
      [&](const ActRoot::MergerData &d, double Ereac_check_Si,
          double ThetaCM_Si) {
        if (cuts.IsInside("weird", Ereac_check_Si, ThetaCM_Si))
          // if(cut.IsInside("debug", d.fSP.Y(), d.fSP.Z()))
          streamer << d.fRun << " " << d.fEntry << '\n';
      },
      {"MergerData", "Ereac_check_Si", "ThetaCM_Si"});
  streamer.close();

  // Filter
  auto vetoed{def.Filter(
      [&](double EVertex_calc, double EVertex) {
        return cuts.IsInside("Ebeamcorrect", EVertex_calc, EVertex);
      },
      {"EVertex_calc", "EVertex"})};

  auto vetoed2{vetoed.Filter(
      [&](float Angle_light, float Angle_heavy) {
        return cuts.IsInside("Cleaning2", Angle_light, Angle_heavy);
      },
      {"Angle_heavy", "Angle_light"})};


  // Book histograms
  auto hEreac{def.Histo1D(HistConfig::Ereac, "EBeam_Si")};
  auto hTheta{d.Histo1D("fThetaLight")};
  // auto hElight{def.Histo1D(HistConfig::Elight, "Elight")};
  auto hEx_proj{vetoed2.Histo1D(HistConfig::Exproj, "Ex_proj")}; // vetoed-def
  // auto hEex_range{def.Histo1D(HistConfig::Ex, "Eex_range")}; //vetoed-def
  auto hEex_range{vetoed2.Histo1D(HistConfig::Ex, "Eex_range")}; // vetoed-def
  auto hEex_Si{def.Histo1D(HistConfig::Ex2, "Eex_Si")};
  // auto hEdiff{def.Histo1D(HistConfig::Ediff, "Ereac_diff")};
  // auto hExlab{def.Histo2D(HistConfig::Ex21Nalab, "fThetaLight", "Ex")};
  // auto hEx{def.Histo2D(HistConfig::Ex21Na, "ThetaCM_Si", "Ex")};
  auto hEx2_range{vetoed2.Histo2D(HistConfig::Ex21Na3, "Ereac_check_range",
                                  "ThetaCM_range")}; // vetoed-def
  auto hEx2_Si{
      vetoed.Histo2D(HistConfig::Ex21Na2, "Angle_heavy", "Angle_light")};
  auto hProj = GetProjectionX(hEx2_range.GetPtr(), angle_min, angle_max);
  // hProj->Rebin(2);
  auto hKin{def.Histo2D(HistConfig::KinEl, "fThetaLight", "EVertex")};
  auto hEex_Si_vs_range{
      def.Histo2D(HistConfig::EexSiVsRange, "Eex_range", "Eex_Si")};
  auto hEbeam_Si_vs_range{
      def.Histo2D(HistConfig::EbeamSiVsRange, "EVertex_calc", "EVertex")};

  // Check compatibility
  if (hEx2_Si->GetNbinsX() != hnorm->GetNbinsX() or
      hEx2_Si->GetNbinsY() != hnorm->GetNbinsY()) {
    std::cerr << "Histograms have incompatible binning!" << std::endl;
    return;
  }

  // Determine the normalization factor
  double normalization_factor = hnorm->Integral() / hEx2_range->Integral();
  std::cout << "normalization factor = " << normalization_factor << "\n";
  std::cout << "integral = " << hEx2_range->Integral() << "\n";
  std::cout << "integral norm = " << hnorm->Integral() << "\n";

  int bin_min = hEx2_range->GetYaxis()->FindBin(angle_min);
  int bin_max = hEx2_range->GetYaxis()->FindBin(angle_max);
  TH2F *hEx2_range_original = (TH2F *)hEx2_range->Clone("hEx2_range_original");
  auto hProj2_nonorm = hEx2_range_original->ProjectionX("hProjY_nonorm", bin_min, bin_max);

  // Apply normalization
  for (int ix = 1; ix <= hEx2_range->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= hEx2_range->GetNbinsY(); ++iy) {
      double bin_content_phy = hEx2_range->GetBinContent(ix, iy);
      double bin_content_norm = hnorm->GetBinContent(ix, iy);
      if (bin_content_norm != 0 && bin_content_phy != 0) {
        hEx2_range->SetBinContent(ix, iy, (bin_content_phy / bin_content_norm)*50);
      } else {
        hEx2_range->SetBinContent(ix, iy, 0);
      }
    }
  }

  auto hProj2 = hEx2_range->ProjectionX("hProjY_norm", bin_min, bin_max);
  // Calcul des erreurs normalisées
  for (int i = 1; i <= hProj2->GetNbinsX(); ++i) {
      double bin_content_nonorm = hProj2_nonorm->GetBinContent(i); // Contenu
      //non normalisé 
      double bin_content_norm = hProj2->GetBinContent(i); //
      //Contenu normalisé

      // Vérification pour éviter division par zéro
      if (bin_content_nonorm > 0 && bin_content_norm > 0) {
          double norm_factor = bin_content_nonorm / bin_content_norm;
          std::cout << "bin content no norm = " << bin_content_nonorm << "\n";
          std::cout << "norm factor = " << norm_factor << "\n";
          double error = TMath::Sqrt(bin_content_nonorm) / norm_factor;
          std::cout << "error no norm = " << TMath::Sqrt(bin_content_nonorm)  << "\n";
          std::cout << "error norm = " << error  << "\n\n";
          hProj2->SetBinError(i, error);
      } else {
          hProj2->SetBinError(i, 0); // Pas d'erreur si le contenu est nul
      }
  }

  // Check compatibility
  // if (hProj->GetNbinsX() != hnorm_proj->GetNbinsX()) {
  //   std::cerr << "Histograms have incompatible binning!" << std::endl;
  //   return;
  // }
  // Check compatibility
  if (hEx_proj->GetNbinsX() != hnorm_proj->GetNbinsX()) {
    std::cerr << "Histograms have incompatible binning!" << std::endl;
    return;
  }

  // Cloning the histogram
  TH1D *hEx_proj_nonorm =
      dynamic_cast<TH1D *>(hEx_proj->Clone("hEx_proj_nonorm"));

  if (normalization == 1) {
    hProj->Divide(hnorm_proj);

    // Apply normalization
    for (int ix = 1; ix <= hEx_proj->GetNbinsX(); ++ix) {
      double bin_content_1 = hEx_proj->GetBinContent(ix);
      double bin_content_2 = hnorm_proj->GetBinContent(ix);
      if (bin_content_2 != 0 && bin_content_1 != 0) {
        hEx_proj->SetBinContent(ix, bin_content_1 / bin_content_2);
      } else {
        hEx_proj->SetBinContent(ix, 0);
      }
    }
  }
  
  gErrorIgnoreLevel = kError; 
  // plotting
  set_plot_style();
  auto *c20{new TCanvas("c20", "Pipe2 canvas 0")};
  c20->DivideSquare(6);
  c20->cd(1);
  hEex_range->DrawClone("colz");
  c20->cd(2);
  hEx2_range->DrawClone("colz");
  c20->cd(3);
  // hProj->DrawClone("");
  // hEx_proj_nonorm->SetLineColor(kAzure-2); // Couleur des lignes et barres
  // d'erreurs
  hEx_proj_nonorm->DrawClone("colz");
  // hEx_proj_nonorm->DrawClone("E1 SAME");
  c20->cd(4);
  hEbeam_Si_vs_range->DrawClone("colz");
  cuts.DrawAll();
  // hEex_Si->DrawClone("colz");
  c20->cd(5);
  hEx2_Si->DrawClone("colz");
  ActPhysics::Kinematics kin1{beam, target, light, ebeam_i_MeV};
  auto *g1{kin1.GetTheta3vs4Line()};
  g1->SetLineColor(kOrange + 6);
  g1->Draw("l");
  ActPhysics::Kinematics kin2{beam, "12C", "8He", ebeam_i_MeV};
  auto *g2{kin2.GetTheta3vs4Line()};
  g2->SetLineColor(kViolet - 4);
  g2->Draw("l");
  cuts.DrawAll();
  c20->cd(6);
  // hEbeam_Si_vs_range->DrawClone("colz");
  hProj2->SetTitle(
      TString::Format("Projection of Ex(12Be) between %.2f and %.2f #circ",
                      angle_min, angle_max));
// Créer le graphe avec erreurs
int nBins = hProj2->GetNbinsX();
TGraphErrors *graph = new TGraphErrors(nBins);

for (int i = 1; i <= nBins; ++i) {
    double x = hProj2->GetBinCenter(i);          // Coordonnée X (centre de la bin)
    double y = hProj2->GetBinContent(i);         // Coordonnée Y (contenu de la bin)
    double ex = 0;                               // Pas d'erreur sur X
    double ey = hProj2->GetBinError(i);          // Erreur sur Y
    graph->SetPoint(i - 1, x, y);                // Ajouter le point
    graph->SetPointError(i - 1, ex, ey);         // Ajouter les barres d'erreurs
}

// Personnalisation des barres d'erreurs et des points
graph->SetMarkerStyle(20);  // Style des points (cercles pleins)
graph->SetMarkerSize(0.5);  // Taille des points
graph->SetMarkerColor(kBlack); // Couleur des points et des barres d'erreurs
graph->SetLineColor(kBlack);   // Couleur des barres d'erreurs
graph->SetLineWidth(2);
graph->SetTitle(
      TString::Format("Proj of Ex between %.2f and %.2f#circ",
                      angle_min, angle_max));
graph->GetXaxis()->SetTitle("Ex_12Be (MeV)");

// Dessiner uniquement les barres d'erreurs et les points
graph->Draw("AP"); // "A" pour axes, "P" pour points et barres d'erreurs


// Utiliser TGraphSmooth pour ajouter une courbe lissée
TGraphSmooth *graphSmooth = new TGraphSmooth("smooth");
TGraph *smoothGraph = graphSmooth->SmoothKern(graph, "normal", 0.1); // méthode Kernel

// Personnalisation de la courbe lissée
smoothGraph->SetLineColor(kRed);  // Couleur de la courbe (bleu)
smoothGraph->SetLineWidth(2);      // Épaisseur de la courbe

// Dessiner la courbe lissée par-dessus
smoothGraph->Draw("C SAME"); // "C" pour courbe lisse, "SAME" pour superposer

// // Ajouter une courbe spline par-dessus
// TGraph *curve = new TGraph(nBins);
// for (int i = 1; i <= nBins; ++i) {
//     double x = hProj2->GetBinCenter(i);
//     double y = hProj2->GetBinContent(i);
//     curve->SetPoint(i - 1, x, y); // Ajouter les points pour la courbe
// }

// // Personnalisation de la courbe
// curve->SetLineColor(kRed);  // Couleur de la courbe (rouge)
// curve->SetLineWidth(2);     // Épaisseur de la courbe

// // Dessiner la courbe spline par-dessus
// curve->Draw("C SAME"); // "C" pour courbe lisse, "SAME" pour superposer


  // Get the y-axis range of the histogram
  double y_min = hEx_proj->GetMinimum();
  double y_max = hEx_proj->GetMaximum();

  // Create TLine objects for each vertical line
  TLine *line1 = new TLine(11.7, y_min, 11.7, y_max);
  TLine *line2 = new TLine(11.85, y_min, 11.85, y_max);
  TLine *line3 = new TLine(12.1, y_min, 12.1, y_max);
  TLine *line4 = new TLine(12.2, y_min, 12.2, y_max);
  TLine *line5 = new TLine(12.4, y_min, 12.4, y_max);
  TLine *line6 = new TLine(12.7, y_min, 12.7, y_max);

  // Set line properties
  line1->SetLineColor(kRed);
  line2->SetLineColor(kOrange + 7);
  line3->SetLineColor(kYellow);
  line4->SetLineColor(kGreen - 3);
  line5->SetLineColor(kAzure + 7);
  line6->SetLineColor(kViolet - 1);

  line1->SetLineWidth(2);
  line2->SetLineWidth(2);
  line3->SetLineWidth(2);
  line4->SetLineWidth(2);
  line5->SetLineWidth(2);
  line6->SetLineWidth(2);

  line2->SetLineStyle(2);
  line4->SetLineStyle(2);
  line5->SetLineStyle(2);
  line6->SetLineStyle(2);

  // Draw the lines on the canvas
  // line1->Draw("same");
  // line2->Draw("same");
  // line3->Draw("same");
  // line4->Draw("same");
  // line5->Draw("same");
  // line6->Draw("same");

  // Create a legend and add entries for each line
  TLegend *legend =
      new TLegend(0.65, 0.65, 0.9, 0.9); // Adjust position as needed
  legend->AddEntry(line1, "x = 11.7", "l");
  legend->AddEntry(line2, "x = 11.85 ?", "l");
  legend->AddEntry(line3, "x = 12.1", "l");
  legend->AddEntry(line4, "x = 12.2 ?", "l");
  legend->AddEntry(line5, "x = 12.4 ?", "l");
  legend->AddEntry(line6, "x = 12.7 ?", "l");

  // Set legend properties 
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
