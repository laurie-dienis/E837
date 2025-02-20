#include "ActColors.h"
#include "TROOT.h"
#include "TString.h"

#include <iostream>
#include <string>

void Print(const std::string &beam, const std::string &target,
           const std::string &light, double e_beam_ini, int pressure,
           const std::string &what = "") {
  std::cout << BOLDGREEN << "···· Runner ····" << '\n';
  std::cout << "-> Beam   : " << beam << '\n';
  std::cout << "-> Target : " << target << '\n';
  std::cout << "-> Light  : " << light << '\n';
  std::cout << "-> EBeam  : " << e_beam_ini << " MeV / u" << '\n';
  std::cout << "-> P      : " << pressure << " mbar" << '\n';
  std::cout << "-> What   : " << what << '\n';
  std::cout << "······························" << RESET << '\n';
}

void Runner(TString what = "plot") {
  std::string beam{"8He"};
  std::string target{"4He"};
  std::string light{"4He"};
  double ebeam_i{};
  int pressure{700}; // mbar
  if (pressure == 700) {
    ebeam_i = 1.486; // MeV/u //1.486
  } else if (pressure == 900) {
    ebeam_i = 1.455; // MeV/u  from run 88 ebeam_i = 1.477 MeV before 1.455
  } else {
    throw std::runtime_error(
        "Runner : cannot set inital beam energy for given pressure ");
  }
  Print(beam, target, light, ebeam_i, pressure, what.Data());

  auto args{TString::Format("(\"%s\", \"%s\", \"%s\")", beam.c_str(),
                            target.c_str(), light.c_str())};
  auto extargs{TString::Format("(\"%s\", \"%s\", \"%s\", %f, %d)", beam.c_str(),
                               target.c_str(), light.c_str(), ebeam_i,
                               pressure)};
  TString path{"./Pipelines/"};
  TString func{};
  TString ext{".cxx"};
  if (what.Contains("0")) {
    func = "Pipe0_Beam";
    gROOT->LoadMacro(path + func + ext);
    gROOT->ProcessLine(func + "()");
  }

  if (what.Contains("1")) {
    func = "Pipe1_PID";
    gROOT->LoadMacro(path + func + ext);
    gROOT->ProcessLine(func + extargs);
  }

  if (what.Contains("2")) {
    func = "Pipe2_Ex";
    gROOT->LoadMacro(path + func + ext);
    gROOT->ProcessLine(func + extargs);
  }

  if (what.Contains("3")) {
    func = "Pipe3_Kine";
    gROOT->LoadMacro(path + func + ext);
    gROOT->ProcessLine(func + extargs);
  }

  if (what.Contains("4")) {
    func = "Pipe4_NN";
    gROOT->LoadMacro(path + func + ext);
    gROOT->ProcessLine(func + extargs);
  }

  if (what.Contains("plot")) {
    func = "Plotter";
    gROOT->LoadMacro("./" + func + ext);
    gROOT->ProcessLine(func + args);
  }
  if (what.Contains("write")) {
    func = "WriteEntries";
    gROOT->LoadMacro("./" + func + ext);
    gROOT->ProcessLine(func + args);
  }
}
