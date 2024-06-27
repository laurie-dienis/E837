#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"
#include "Rtypes.h"

#include "TCanvas.h"
#include "THStack.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TString.h"
#include "TVirtualPad.h"

#include <iostream>

#include "./Profile.cxx"
#include "./Profile.h"

void TreatProfile(bool draw = false)
{
    // Read the data
    // 1-> Merge contains the ProfileX itself
    ActRoot::DataManager dataman {"../../configs/data.conf", ActRoot::ModeType::EMerge};
    auto chain {dataman.GetJoinedData()};
    // 2-> We need to friend the chain with Filter output to gate on cluster multiplicity (once processed)
    auto other {dataman.GetJoinedData(ActRoot::ModeType::EFilter)};
    chain->AddFriend(other.get());

    // Build the RDataFrame
    ROOT::RDataFrame d {*chain};

    // Gate on cluster multiplicity
    // We must have either 1 (beam + 11B fully contained in beam region) or 2 (part of 11B getting out of beam region)
    auto df {d.Filter([](const ActRoot::TPCData& tpc)
                      { return (1 <= tpc.fClusters.size()) && (tpc.fClusters.size() <= 2); },
                      {"TPCData"})};
    std::cout << "Number of events with OK multiplicity : " << df.Count().GetValue() << '\n';

    // Init canvas to visual inspect
    TCanvas* ci {};
    if(draw)
        ci = new TCanvas {"ci", "Inspection canvas"};
    // And process!
    // 1-> Buils lambda function to execute all the hard work
    auto process {[&](ActRoot::MergerData& merger)
                  {
                      // Get pointer to profile X stored in merger data
                      auto* hProf {&merger.fQprojX};
                      // Init return class
                      ChargeProfile::Data qdata {};

                      // 1. Compute derivative
                      auto* hDer {ChargeProfile::GetDerivative(hProf)};

                      // 2. Compute the vertex and qmax
                      // Set width for the Qmax search around the given
                      double width {10}; // pad or mm
                      auto [xmax, ymax] {ChargeProfile::GetVertexAndQmax(hDer, hProf, width)};
                      // Store results
                      qdata.fVertexP = {xmax, 0, 0};
                      qdata.fQVertex = ymax;

                      // 3. Compute range of recoils
                      auto xstop {ChargeProfile::GetStoppingPoint(hProf)};
                      // Store it
                      qdata.fStopP = {xstop, 0, 0};

                      // Visual inspection
                      if(draw)
                      {
                          ci->DivideSquare(4);
                          ci->cd(1);
                          // create a stack to both histograms in the same pad
                          auto* stack {new THStack};
                          stack->SetTitle(TString::Format("Run %d entry %d Qprof;X [pad or mm];Q and derivative [au]",
                                                          merger.fRun, merger.fEntry));
                          stack->Add(hProf, "hist");
                          hDer->SetLineColor(kMagenta);
                          stack->Add(hDer);
                          stack->Draw("nostack");

                          // Draw also a TLine
                          ci->Update();
                          auto* lv {new TLine {xmax, gPad->GetUymin(), xmax, gPad->GetUymax()}};
                          lv->SetLineWidth(2);
                          lv->SetLineColor(kRed);
                          lv->Draw();
                          auto* ls {new TLine {xstop, gPad->GetUymin(), xstop, gPad->GetUymax()}};
                          ls->SetLineWidth(2);
                          ls->SetLineColor(kOrange);
                          ls->Draw();

                          // And box with info
                          auto* box {new TPaveText {0.7, 0.7, 0.9, 0.9, "ndc"}};
                          box->SetBorderSize(0);
                          box->SetFillStyle(0);
                          box->AddText(TString::Format("X_{Vertex} = %.2f pad or mm", qdata.fVertexP.X()));
                          box->AddText(TString::Format("X_{Stop} = %.0f pad or mm", qdata.fStopP.X()));
                          box->AddText(TString::Format("Q_{Vertex} = %.0f au", qdata.fQVertex));
                          box->Draw();

                          ci->Update();
                          ci->WaitPrimitive("dummy", "");
                          ci->Clear();
                          ci->Update();
                          gROOT->SetSelectedPad(nullptr);
                          delete stack;
                          delete lv;
                          delete ls;
                          delete box;
                      }


                      // Delete news
                      delete hDer;
                      // And return
                      return qdata;
                  }};

    // Execute
    auto def {df.Define("ProfileData", process, {"MergerData"})};

    // Save
    // def.Snapshot("Test", "./test.root");
}
