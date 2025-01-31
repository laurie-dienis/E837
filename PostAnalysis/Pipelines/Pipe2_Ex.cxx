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

// Function to export histogram data to a file
void ExportToAZURE_2D(TH2D *hEx2_range, const std::string &outputFilename) {
  // Check if the histogram pointer is null
  if (!hEx2_range) {
    std::cerr << "Error: The provided histogram is null." << std::endl;
    return;
  }

  // Open the output file
  std::ofstream outFile(outputFilename);
  if (!outFile.is_open()) {
    std::cerr << "Error: Unable to create the output file." << std::endl;
    return;
  }

  // Set the output format to scientific notation with precision
  outFile << std::scientific << std::setprecision(3);

  // Loop over the bins of the histogram
  for (int ix = 1; ix <= hEx2_range->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= hEx2_range->GetNbinsY(); ++iy) {
      double energy = hEx2_range->GetXaxis()->GetBinCenter(
          ix); // Bin center for X-axis (Energy)
      double angle = hEx2_range->GetYaxis()->GetBinCenter(
          iy); // Bin center for Y-axis (Angle)
      double crossSection =
          hEx2_range->GetBinContent(ix, iy)*1e-3; // Bin content (Cross-section)
      double uncertainty =
          hEx2_range->GetBinError(ix, iy)*1e-3; // Bin error (Uncertainty)

      // Optional: Only write bins with non-zero content
      if (crossSection != 0) {
        outFile << energy << "\t" << angle << "\t" << crossSection << "\t"
                << uncertainty << std::endl;
      }
    }
  }

  // Close the file
  outFile.close();
  std::cout << "File " << outputFilename << " created successfully."
            << std::endl;
}

// Function to export histogram data to a file
void ExportToAZURE_1D_EProj(TH1D *hEx2_range, const std::string &outputFilename,
                            const double &angle_min) {
  // Check if the histogram pointer is null
  if (!hEx2_range) {
    std::cerr << "Error: The provided histogram is null." << std::endl;
    return;
  }

  // Open the output file
  std::ofstream outFile(outputFilename);
  if (!outFile.is_open()) {
    std::cerr << "Error: Unable to create the output file." << std::endl;
    return;
  }

  // Set the output format to scientific notation with precision
  outFile << std::scientific << std::setprecision(3);

  // Loop over the bins of the histogram
  for (int ix = 1; ix <= hEx2_range->GetNbinsX(); ++ix) {
    double energy = hEx2_range->GetXaxis()->GetBinCenter(ix); // Bin center for X-axis (Energy)
    double crossSection = hEx2_range->GetBinContent(ix)*1e-3; // Bin content (Cross-section)
    double uncertainty = hEx2_range->GetBinError(ix)*1e-3; // Bin error (Uncertainty)


    // Optional: Only write bins with non-zero content
    if (crossSection != 0) {
      outFile << energy << "\t" << angle_min << "\t" << crossSection << "\t"
              << uncertainty << std::endl;
    }
}

// Close the file
outFile.close();
std::cout << "File " << outputFilename << " created successfully." << std::endl;
}

// Function to export histogram data to a file
void ExportToAZURE_1D_ThetaProj(TH1D *hEx2_range, const std::string &outputFilename,
                            const double &energy_min) {
  // Check if the histogram pointer is null
  if (!hEx2_range) {
    std::cerr << "Error: The provided histogram is null." << std::endl;
    return;
  }

  // Open the output file
  std::ofstream outFile(outputFilename);
  if (!outFile.is_open()) {
    std::cerr << "Error: Unable to create the output file." << std::endl;
    return;
  }

  // Set the output format to scientific notation with precision
  outFile << std::scientific << std::setprecision(3);

  // Loop over the bins of the histogram
  for (int ix = 1; ix <= hEx2_range->GetNbinsX(); ++ix) {
    double angle = hEx2_range->GetXaxis()->GetBinCenter(
        ix); // Bin center for X-axis (Energy)
    double crossSection =
        hEx2_range->GetBinContent(ix)*1e-3; // Bin content (Cross-section)
    double uncertainty = hEx2_range->GetBinError(ix)*1e-3; // Bin error (Uncertainty)

    // Optional: Only write bins with non-zero content
    if (crossSection != 0) {
      outFile << energy_min << "\t" << angle<< "\t" << crossSection << "\t"
              << uncertainty << std::endl;
    }
}

// Close the file
outFile.close();
std::cout << "File " << outputFilename << " created successfully." << std::endl;
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

double ComputeExcitationEnergy_angle(double ebeam, double theta_heavy,
                                     double theta_light, double mlight,
                                     double mheavy) {
  double evertex = ebeam * (1 - (mheavy / (mlight + mheavy))) +
                   0.5 * mheavy * TMath::Cos(theta_heavy * TMath::DegToRad());
  return evertex * ((2 * TMath::Cos(theta_light * TMath::DegToRad()) *
                     TMath::Sqrt((ebeam * mlight) / (evertex * mheavy))) -
                    ((mlight / mheavy) + 1));
}

double ComputeExcitationEnergy_elastic(double evertex, double ebeam,
                                       double theta, double mlight,
                                       double mheavy) {
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
  // Setting
  double angle_min{130.};
  double angle_max{140.};

  double angle_min_azure{20.7}; //19.3 ou 26
  double angle_max_azure{27.3}; //23.3 ou 30
  
  double energy_min{11.61}; // 11.61 or 12.22
  double energy_max{12.17}; // 12.17 or 12.78

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
      "Eex_angle",
      [&](unsigned int slot, const ActRoot::MergerData &d, double ebeam) {
        if (std::isnan(ebeam)) {
          return -11.;
        }
        double Ex2;
        Ex2 = ComputeExcitationEnergy_angle(ebeam, d.fThetaHeavy, d.fThetaLight,
                                            mass_light, mass_beam);
        return Ex2;
      },
      {"MergerData", "EBeam_range"});

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
  def =
      def.Define("Ereac_check_range",
                 [&](double ebeam, double Ex) {
                   // std::cout << "ereac = " <<
                   // vkins[slot].ReconstructBeamEnergyFromLabKinematics(evertex,
                   // theta * TMath::DegToRad()) *
                   //             (mass_target / (mass_target + mass_beam)) +
                   //       kin_particle_threshold << " MeV" << "\n";
                   if (Ex > -0.2 && Ex < 0.2) {
                     return (ebeam * (mass_target / (mass_target + mass_beam)) +
                             kin_particle_threshold);
                   } else {
                     return -10.;
                   }
                 },
                 {"EBeam_range", "Eex_range"});

  // Elab alpha direct kinematic for AZURE2
  def = def.Define("Elab",
                   [&](double ebeam) {
                     double Ecm =
                         ebeam * (mass_target / (mass_target + mass_beam));
                     double Elab =
                         Ecm * ((mass_target + mass_beam) / mass_beam);
                     return Elab;
                   },
                   {"EBeam_range"});

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

  // Retrieve the TH2D histogram named "hnorm_azure" from the file
  TH2D *hnorm_azure = dynamic_cast<TH2D *>(norm->Get("hnorm_azure"));
  if (!hnorm_azure) {
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
                             "Debug/Secondcleaning3_%dmbar.root",
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
  auto hEex_Si{vetoed2.Histo1D(HistConfig::Ex2, "Eex_Si")};
  // auto hEdiff{def.Histo1D(HistConfig::Ediff, "Ereac_diff")};
  // auto hExlab{def.Histo2D(HistConfig::Ex21Nalab, "fThetaLight", "Ex")};
  // auto hEx{def.Histo2D(HistConfig::Ex21Na, "ThetaCM_Si", "Ex")};
  auto hEx2_range{vetoed2.Histo2D(HistConfig::Ex21Na3, "Ereac_check_Si",
                                  "ThetaCM_range")}; // vetoed-def
  auto hEx2_azure{vetoed2.Histo2D(HistConfig::Ex21Na4, "Elab",
                                  "Angle_light")}; // vetoed-def
  auto hEx2_Si{def.Histo2D(HistConfig::Ex21Na2, "Angle_heavy", "Angle_light")};
  auto hProj = GetProjectionX(hEx2_range.GetPtr(), angle_min, angle_max);
  // hProj->Rebin(2);
  auto hKin{def.Histo2D(HistConfig::KinEl, "fThetaLight", "EVertex")};
  auto hEex_Si_vs_range{
      def.Histo2D(HistConfig::EexSiVsRange, "Eex_range", "Eex_Si")};
  auto hEbeam_Si_vs_range{
      vetoed2.Histo2D(HistConfig::EbeamSiVsRange, "EVertex_calc", "EVertex")};

  //////////////////////////////////////////////////////////////////////  
  /////////////Normalize and pass to cross_section//////////////////////
  ////////////////////////////////////////////////////////////////////
  // Check compatibility
  if (hEx2_Si->GetNbinsX() != hnorm->GetNbinsX() or
      hEx2_Si->GetNbinsY() != hnorm->GetNbinsY()) {
    std::cerr << "Histograms have incompatible binning!" << std::endl;
    return;
  }

  int bin_min = hEx2_range->GetYaxis()->FindBin(angle_min);
  int bin_max = hEx2_range->GetYaxis()->FindBin(angle_max);
  TH2F *hEx2_range_original = (TH2F *)hEx2_range->Clone("hEx2_range_original");
  auto hProj2_nonorm =
      hEx2_range_original->ProjectionX("hProjY_nonorm", bin_min, bin_max);

  int bin_min_theta = hEx2_range->GetXaxis()->FindBin(energy_min);
  int bin_max_theta = hEx2_range->GetXaxis()->FindBin(energy_max);
  auto hProjTheta_nonorm = hEx2_range_original->ProjectionY(
      "hProjX_nonorm", bin_min_theta, bin_max_theta);

  // Apply normalization
  for (int ix = 1; ix <= hEx2_range->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= hEx2_range->GetNbinsY(); ++iy) {
      double bin_content_phy = hEx2_range->GetBinContent(ix, iy);
      double bin_content_norm = hnorm->GetBinContent(ix, iy);
      double angle_deg = hEx2_range->GetYaxis()->GetBinCenter(iy);
      double cross_section_factor_2;
      if (pressure == 700) {
        cross_section_factor_2 =
            1e24 * 1e3 /
            (2 * TMath::Pi() * TMath::Sin((angle_deg)*TMath::DegToRad()) *
             0.0002390941 * 2.58e9 * 1.29e19);
      }
      if (pressure == 900) {
        cross_section_factor_2 =
            1e24 * 1e3 /
            (2 * TMath::Pi() * TMath::Sin((angle_deg)*TMath::DegToRad()) *
             0.0002390941 * 1.09e9 * 2.21e19);
      }
      if (bin_content_norm != 0 && bin_content_phy != 0) {
        // std::cout<<"norm = "<<bin_content_norm<<"\n";
        // std::cout<<"counts = "<<bin_content_phy<<"\n";
        // std::cout<<"angle = "<<angle_deg<<"\n";
        // std::cout<<"cross_section_factor = "<<cross_section_factor_2<<"\n";
        // std::cout<<"\n";
        hEx2_range->SetBinContent(ix, iy,
                                  cross_section_factor_2 *
                                      (bin_content_phy / bin_content_norm));
      } else {
        hEx2_range->SetBinContent(ix, iy, 0);
      }
    }
  }

  //////////////////////////////////////////////////
  // Normalization and cross section for AZURE2 ////
  //////////////////////////////////////////////

  int bin_min_AZURE = hEx2_azure->GetYaxis()->FindBin(angle_min_azure);
  int bin_max_AZURE = hEx2_azure->GetYaxis()->FindBin(angle_max_azure);
  TH2F *hEx2_azure_original = (TH2F *)hEx2_range->Clone("hEx2_azure_original");
  auto hProj_azure_nonorm = hEx2_azure_original->ProjectionX(
      "hProjY_nonorm_azure", bin_min_AZURE, bin_max_AZURE);

  int bin_min_theta_AZURE = hEx2_range->GetXaxis()->FindBin((energy_min-8.96)*1.5);
  int bin_max_theta_AZURE = hEx2_range->GetXaxis()->FindBin((energy_max-8.96)*1.5);
  auto hProj_azure_Theta_nonorm = hEx2_azure_original->ProjectionY(
      "hProjX_nonorm_azure", bin_min_theta_AZURE, bin_max_theta_AZURE);

  /// Apply normalization
  for (int ix = 1; ix <= hEx2_azure->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= hEx2_azure->GetNbinsY(); ++iy) {
      double bin_content_phy = hEx2_azure->GetBinContent(ix, iy);
      double bin_content_norm = hnorm_azure->GetBinContent(ix, iy);
      double angle_deg = hEx2_azure->GetYaxis()->GetBinCenter(iy);
      double cross_section_factor_2;
      if (pressure == 700) {
        cross_section_factor_2 =
            1e24 * 1e3 /
            (2 * TMath::Pi() * TMath::Sin((angle_deg)*TMath::DegToRad()) *
             0.0002390941 * 2.58e9 * 1.29e19);
      }
      if (pressure == 900) {
        cross_section_factor_2 =
            1e24 * 1e3 /
            (2 * TMath::Pi() * TMath::Sin((angle_deg)*TMath::DegToRad()) *
             0.0002390941 * 1.09e9 * 2.21e19);
      }
      if (bin_content_norm != 0 && bin_content_phy != 0) {
        //std::cout<<"norm = "<<bin_content_norm<<"\n";
        // std::cout<<"counts = "<<bin_content_phy<<"\n";
        // std::cout<<"angle = "<<angle_deg<<"\n";
        // std::cout<<"cross_section_factor = "<<cross_section_factor_2<<"\n";
        // std::cout<<"\n";
        hEx2_azure->SetBinContent(ix, iy,
                                  cross_section_factor_2 *
                                      (bin_content_phy / bin_content_norm));
      } else {
        hEx2_azure->SetBinContent(ix, iy, 0);
      }
    }
  }

  // Set the output filename dynamically
  std::string outputFilename = std::string(
      TString::Format("ACTAR_differential_%d.dat", pressure).Data());

  // Get a pointer to the histogram
  TH2D *hEx2_azure_ptr = hEx2_azure.GetPtr();

  ExportToAZURE_2D(hEx2_azure_ptr, outputFilename);

  //////////////////////////////////////////
  //////// Projection Ex ///////////////////
  //////////////////////////////////////////
  auto hProj2 = hEx2_range->ProjectionX("hProjY_norm", bin_min, bin_max);
  // Calcul des erreurs normalisées
  for (int i = 1; i <= hProj2->GetNbinsX(); ++i) {
    double bin_content_nonorm = hProj2_nonorm->GetBinContent(i); // Contenu
    // non normalisé
    double bin_content_norm = hProj2->GetBinContent(i) / 10.; //
    // Contenu normalisé

    // Vérification pour éviter division par zéro
    if (bin_content_nonorm > 0 && bin_content_norm > 0) {
      hProj2->SetBinContent(i, bin_content_norm);
      double norm_factor = bin_content_norm / bin_content_nonorm;
      // std::cout << "bin content no norm = " << bin_content_nonorm << "\n";
      // std::cout << "norm factor = " << norm_factor << "\n";
      double error = TMath::Sqrt(bin_content_nonorm) * norm_factor;
      // std::cout << "error no norm = " << TMath::Sqrt(bin_content_nonorm)  <<
      // "\n"; std::cout << "error norm = " << error  << "\n\n";
      hProj2->SetBinError(i, error);
    } else {
      hProj2->SetBinError(i, 0); // Pas d'erreur si le contenu est nul
    }
  }

  //////////////////////////////////////////
  //////Projection Elab for azure///////////
  /////////////////////////////////////////
// Définition des paramètres de transformation
double a = kin_particle_threshold;  // Remplace par ta valeur
double b = (mass_beam+mass_target)/mass_beam;  // Remplace par ta valeur

// Création d'un nouvel histogramme pour stocker la transformation
TH1D* hProj2_transformed = (TH1D*)hProj2->Clone("hProj2_transformed");
hProj2_transformed->Reset(); // On vide l'histogramme pour le remplir correctement

// Application de la transformation sur chaque bin
for (int i = 1; i <= hProj2->GetNbinsX(); ++i) {
    double x = hProj2->GetBinCenter(i);  // Coordonnée x du centre du bin
    double new_x = (x - a) * b;          // Transformation
    std::cout<<"energie = "<<new_x<<"\n";

    double content = hProj2->GetBinContent(i); // Contenu du bin d'origine
    double error = hProj2->GetBinError(i);     // Erreur associée

    // Remplissage du nouvel histogramme
    int new_bin = hProj2_transformed->FindBin(new_x);
    hProj2_transformed->SetBinContent(new_bin, content);
    hProj2_transformed->SetBinError(new_bin, error);
} 


  auto hProj_azure_elab = hEx2_azure->ProjectionX("hProjY_azure_elab",
                                                  bin_min_AZURE, bin_max_AZURE);
  // Calcul des erreurs normalisées
  for (int i = 1; i <= hProj_azure_nonorm->GetNbinsX(); ++i) {
    double bin_content_nonorm = hProj_azure_nonorm->GetBinContent(i); // Contenu
    // non normalisé
    double bin_content_norm = hProj_azure_elab->GetBinContent(i) / 10.; //
    // Contenu normalisé

    // Vérification pour éviter division par zéro
    if (bin_content_nonorm > 0 && bin_content_norm > 0) {
      hProj_azure_elab->SetBinContent(i, bin_content_norm);
      double norm_factor = bin_content_norm / bin_content_nonorm;
      // std::cout << "bin content no norm = " << bin_content_nonorm << "\n";
      // std::cout << "norm factor = " << norm_factor << "\n";
      double error = TMath::Sqrt(bin_content_nonorm) * norm_factor;
      // std::cout << "error no norm = " << TMath::Sqrt(bin_content_nonorm)  <<
      // "\n"; std::cout << "error norm = " << error  << "\n\n";
      hProj_azure_elab->SetBinError(i, error);
    } else {
      hProj_azure_elab->SetBinError(i, 0); // Pas d'erreur si le contenu est nul
    }
  }

  // Set the output filename dynamically
  std::string outputFilename_proje =
      std::string(TString::Format("ACTAR_angle_%f_%f_%d.dat", angle_min_azure,
                                  angle_max_azure, pressure)
                      .Data());

  ExportToAZURE_1D_EProj(hProj_azure_elab, outputFilename_proje, angle_min_azure);

  ///////////////////////////////////////////
  ///////////Projection ThetaCM ////////////
  //////////////////////////////////////////
  auto hProjTheta =
      hEx2_range->ProjectionY("hProjX_norm", bin_min_theta, bin_max_theta);
  // Calcul des erreurs normalisées
  for (int i = 1; i <= hProjTheta->GetNbinsX(); ++i) {
    double bin_content_nonorm = hProjTheta_nonorm->GetBinContent(i); // Contenu
    // non normalisé
    double bin_content_norm = hProjTheta->GetBinContent(i) / 10.; //
    // Contenu normalisé

    // Vérification pour éviter division par zéro
    if (bin_content_nonorm > 0 && bin_content_norm > 0) {
      hProjTheta->SetBinContent(i, bin_content_norm);
      double norm_factor = bin_content_norm / bin_content_nonorm;
      // std::cout << "bin content no norm = " << bin_content_nonorm << "\n";
      // std::cout << "norm factor = " << norm_factor << "\n";
      double error = TMath::Sqrt(bin_content_nonorm) * norm_factor;
      // std::cout << "error no norm = " << TMath::Sqrt(bin_content_nonorm)  <<
      // "\n"; std::cout << "error norm = " << error  << "\n\n";
      hProjTheta->SetBinError(i, error);
    } else {
      hProjTheta->SetBinError(i, 0); // Pas d'erreur si le contenu est nul
    }
  }
  
  /////////////////////////////////////////
  /////Projection Thetalab for azure2 /////
  /////////////////////////////////////////
  auto hProj_azure_theta =
      hEx2_azure->ProjectionY("hProjX_azure_norm", bin_min_theta_AZURE, bin_max_theta_AZURE);
  // Calcul des erreurs normalisées
  for (int i = 1; i <= hProj_azure_theta->GetNbinsX(); ++i) {
    double bin_content_nonorm = hProj_azure_Theta_nonorm->GetBinContent(i); // Contenu
    // non normalisé
    double bin_content_norm = hProj_azure_theta->GetBinContent(i) / 10.; //
    // Contenu normalisé

    // Vérification pour éviter division par zéro
    if (bin_content_nonorm > 0 && bin_content_norm > 0) {
      hProj_azure_theta->SetBinContent(i, bin_content_norm);
      double norm_factor = bin_content_norm / bin_content_nonorm;
      // std::cout << "bin content no norm = " << bin_content_nonorm << "\n";
      // std::cout << "norm factor = " << norm_factor << "\n";
      double error = TMath::Sqrt(bin_content_nonorm) * norm_factor;
      // std::cout << "error no norm = " << TMath::Sqrt(bin_content_nonorm)  <<
      // "\n"; std::cout << "error norm = " << error  << "\n\n";
      hProj_azure_theta->SetBinError(i, error);
    } else {
      hProj_azure_theta->SetBinError(i, 0); // Pas d'erreur si le contenu est nul
    }
  }
  
  // Set the output filename dynamically
  std::string outputFilename_projtheta =
      std::string(TString::Format("ACTAR_angle_%f_%f_%d.dat", energy_min,
                                  energy_max, pressure)
                      .Data());

  ExportToAZURE_1D_ThetaProj(hProj_azure_elab, outputFilename_projtheta, energy_min);

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
      // std::cout<<"norm = "<<bin_content_2<<"\n";
      double cross_section_factor =
          1e24 * 1e3 /
          (2 * TMath::Pi() *
           TMath::Sin(((angle_max + angle_min) / 2) * TMath::DegToRad()) *
           0.02390941 * 2.58e9 * 1.29e19);
      if (bin_content_2 != 0 && bin_content_1 != 0) {
        hEx_proj->SetBinContent(ix, cross_section_factor *
                                        (bin_content_1 / bin_content_2));
        double norm_factor_2 = cross_section_factor / bin_content_2;
        double error_2 = TMath::Sqrt(bin_content_1) * norm_factor_2;
        hEx_proj->SetBinError(ix, error_2);
      } else {
        hEx_proj->SetBinContent(ix, 0);
        hEx_proj->SetBinError(ix, 0); // Pas d'erreur si le contenu est nul
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
  hEx2_azure->DrawClone("colz");
  c20->cd(3);
  // hProj->DrawClone("");
  // hEx_proj_nonorm->SetLineColor(kAzure-2); // Couleur des lignes et barres
  // d'erreurs
  hProjTheta->SetTitle(
      TString::Format("Projection of Ex(12Be) between %.2f and %.2f #circ",
                      angle_min, angle_max));
  // Créer le graphe avec erreurs
  int nBins_theta = hProjTheta->GetNbinsX();
  TGraphErrors *graph_theta = new TGraphErrors(nBins_theta);

  for (int i = 1; i <= nBins_theta; ++i) {
    double x = hProjTheta->GetBinCenter(i);  // Coordonnée X (centre de la bin)
    double y = hProjTheta->GetBinContent(i); // Coordonnée Y (contenu de la bin)
    double ex = 0;                           // Pas d'erreur sur X
    double ey = hProjTheta->GetBinError(i);  // Erreur sur Y
    graph_theta->SetPoint(i - 1, x, y);      // Ajouter le point
    graph_theta->SetPointError(i - 1, ex, ey); // Ajouter les barres d'erreurs
  }

  // Personnalisation des barres d'erreurs et des points
  graph_theta->SetMarkerStyle(20); // Style des points (cercles pleins)
  graph_theta->SetMarkerSize(0.5); // Taille des points
  graph_theta->SetMarkerColor(
      kBlack); // Couleur des points et des barres d'erreurs
  graph_theta->SetLineColor(kBlack); // Couleur des barres d'erreurs
  graph_theta->SetLineWidth(2);
  graph_theta->SetTitle(
      TString::Format("Proj of #theta_{CM} between %.2f and %.2f#circ",
                      energy_min, energy_max));
  graph_theta->GetXaxis()->SetTitle("#theta_{CM} (^{#circ})");
  graph_theta->GetYaxis()->SetTitle("#frac{d#sigma}{d#Omega} (mb/sr)");

  // Dessiner uniquement les barres d'erreurs et les points
  graph_theta->Draw("AP"); // "A" pour axes, "P" pour points et barres d'erreurs
  // hEx_proj_nonorm->DrawClone("E1 SAME");
  // Utiliser TGraphSmooth pour ajouter une courbe lissée
  TGraphSmooth *graphSmooth_theta = new TGraphSmooth("smooth");
  TGraph *smoothGraph_theta = graphSmooth_theta->SmoothKern(
      graph_theta, "normal", 0.1); // méthode Kernel

  // Personnalisation de la courbe lissée
  smoothGraph_theta->SetLineColor(kRed); // Couleur de la courbe (bleu)
  smoothGraph_theta->SetLineWidth(2);    // Épaisseur de la courbe

  // Dessiner la courbe lissée par-dessus
  smoothGraph_theta->Draw(
      "C SAME"); // "C" pour courbe lisse, "SAME" pour superposer
    
    //0+
    TF1 *line10 = new TF1("line10","(200*456518.17*(x-100)^(-4.0849579))",100, 160);
    line10->SetLineColor(kRed);
    line10->SetLineWidth(4);
    line10->Draw("same");

    //1-
    TF1 *line20 = new TF1("line20","(200*32523.259*(x-100)^(-3.2306928))",100, 160);
    line20->SetLineColor(kOrange);
    line20->SetLineWidth(4);
    line20->Draw("same");

    //2+
    TF1 *line30 = new TF1("line30","(200*13.182879*exp(-0.08076059*(x-100)))",100, 160);
    line30->SetLineColor(kBlue);
    line30->SetLineWidth(4);
    line30->Draw("same");
    
    //3-
    TF1 *line40 = new TF1("line40","150*(-0.63635095*(x-100)+29.04588)",100, 160);
    line40->SetLineColor(kGreen);
    line40->SetLineWidth(4);
    line40->Draw("same");

    //4+
    TF1 *line50 = new TF1("line50","200*130608.17*(x-100)^(-3.624432600)",100, 160);
    line50->SetLineColor(kViolet);
    line50->SetLineWidth(4);
    line50->Draw("same");
  
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
  hProj_azure_elab->SetTitle(
      TString::Format("Projection of Ex(12Be) between %.2f and %.2f #circ",
                      angle_min, angle_max));
  // Créer le graphe avec erreurs
  int nBins = hProj_azure_elab->GetNbinsX();
  TGraphErrors *graph = new TGraphErrors(nBins);

  for (int i = 1; i <= nBins; ++i) {
    double x = hProj_azure_elab->GetBinCenter(i);  // Coordonnée X (centre de la bin)
    double y = hProj_azure_elab->GetBinContent(i); // Coordonnée Y (contenu de la bin)
    double ex = 0;                       // Pas d'erreur sur X
    double ey = hProj_azure_elab->GetBinError(i);  // Erreur sur Y
    graph->SetPoint(i - 1, x, y);        // Ajouter le point
    graph->SetPointError(i - 1, ex, ey); // Ajouter les barres d'erreurs
  }

  // Personnalisation des barres d'erreurs et des points
  graph->SetMarkerStyle(20);     // Style des points (cercles pleins)
  graph->SetMarkerSize(0.5);     // Taille des points
  graph->SetMarkerColor(kBlack); // Couleur des points et des barres d'erreurs
  graph->SetLineColor(kBlack);   // Couleur des barres d'erreurs
  graph->SetLineWidth(2);
  graph->SetTitle(TString::Format("Proj of Ex between %.2f and %.2f#circ",
                                  angle_min, angle_max));
  graph->GetXaxis()->SetTitle("Ex_{^{12}Be} (MeV)");
  graph->GetYaxis()->SetTitle("#frac{d#sigma}{d#Omega} (mb/sr)");

  // Dessiner uniquement les barres d'erreurs et les points
  graph->Draw("AP"); // "A" pour axes, "P" pour points et barres d'erreurs

  // Utiliser TGraphSmooth pour ajouter une courbe lissée
  TGraphSmooth *graphSmooth = new TGraphSmooth("smooth");
  TGraph *smoothGraph =
      graphSmooth->SmoothKern(graph, "normal", 0.1); // méthode Kernel

  // Personnalisation de la courbe lissée
  smoothGraph->SetLineColor(kRed); // Couleur de la courbe (bleu)
  smoothGraph->SetLineWidth(2);    // Épaisseur de la courbe

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
