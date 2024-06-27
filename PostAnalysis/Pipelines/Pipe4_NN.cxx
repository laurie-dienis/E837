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

#include "ActMergerData.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "../HistConfig.h"
#include "../Utils.cxx"

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

void Pipe4_NN(const std::string &beam, const std::string &target,
              const std::string &light, double ebeam_i, int pressure) {
  // ROOT::EnableImplicitMT();
  //  Read data
  ROOT::RDataFrame d{"PID_Tree",
                     E837Utils::GetFileName(1, pressure, beam, target, light)};
  auto df = d.Filter("ESil>0");

  // Ereac, energy of the beam at the reaction point
  df = df.Define("Data",
                 [](ActRoot::MergerData &data) {
                   auto *p{data.fQprojX.GetArray()};
                   return std::vector<double>(p, p + data.fQprojX.GetNbinsX());
                 },
                 {"MergerData"});

  df.Snapshot("NN_Tree",
              TString::Format("/home/laurie/Analysis_e837/E837/PostAnalysis/"
                              "RootFiles/Pipe4/Tree_%s_%dmbar.root",
                              light.c_str(), pressure),
              {"Data"});

}
#endif
