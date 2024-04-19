#include "ActCutsManager.h"
#include "ActLine.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RResultPtr.hxx"
#include "Rtypes.h"

#include "TCanvas.h"
#include "TChain.h"
#include "TEllipse.h"
#include "TF1.h"
#include "TH2.h"
#include "TMath.h"
#include "TMathBase.h"
#include "TROOT.h"
#include "TString.h"

#include "Math/Point3Dfwd.h"

#include <fstream>
#include <iostream>
#include <vector>

const double padSide {2.};      // mm
const double timeFactor {0.08}; // tb to us
const int binFactor {4};
const ROOT::Math::XYZPointF sourcePoint {-22.5 * padSide, 39.5 * padSide, 0};

std::map<int, std::vector<double>> ReadFile(int run)
{
    std::map<int, std::vector<double>> ret;
    std::ifstream streamer {TString::Format("./Inputs/Drift/ell_run%d.dat", run)};
    if(!streamer)
        throw std::runtime_error("You must provide an ellipse conf file for this run");
    int alpha {};
    double rx1 {};
    double rx2 {};
    double ry1 {};
    double ry2 {};
    while(streamer >> alpha >> rx1 >> ry1 >> rx2 >> ry2)
        ret[alpha] = {rx1, ry1, rx2, ry2};
    return ret;
}

bool IsInsideEllipse(const std::vector<double>& pars, double deltaT, double Rxy)
{
    double ybig {TMath::Sqrt(TMath::Abs((pars[1] * pars[1]) * (1 - deltaT * deltaT / (pars[0] * pars[0]))))};
    double xbig {TMath::Sqrt(TMath::Abs((pars[0] * pars[0]) * (1 - Rxy * Rxy / (pars[1] * pars[1]))))};

    double ysmall {TMath::Sqrt(TMath::Abs((pars[3] * pars[3]) * (1 - deltaT * deltaT / (pars[2] * pars[2]))))};
    double xsmall {TMath::Sqrt(TMath::Abs((pars[2] * pars[2]) * (1 - Rxy * Rxy / (pars[3] * pars[3]))))};
    bool condR {ysmall <= Rxy && Rxy <= ybig};
    bool condT {xsmall <= TMath::Abs(deltaT) && TMath::Abs(deltaT) <= xbig};
    if(condR || condT)
        return true;
    else
        return false;
}

void Drift()
{
    // Enable MT
    ROOT::EnableImplicitMT();
    // Reading the data
    const int run {28};
    // Ellipses
    // Read data
    auto confs {ReadFile(run)};
    std::cout << "run number = " << run << '\n';
    auto* chain {new TChain("ACTAR_Merged")};
    chain->Add(TString::Format("./../RootFiles/Merger/Merged_Run_%04d.root", run));

    ROOT::RDataFrame df {*chain};
    // Only events with multiplicity 1
    auto gated {df.Filter("fBraggP.fCoordinates.fX != -1")};

    // Define varibales
    auto def {gated
                  .Define("Rxy",
                          [](const ROOT::Math::XYZPointF& braggP, const ROOT::Math::XYZPointF& wP)
                          {
                              // Just projections along 2D
                              ROOT::Math::XYZPointF bragg2D {braggP.X(), braggP.Y(), 0};
                              bragg2D *= padSide;
                              ROOT::Math::XYZPointF w2D {wP.X(), wP.Y(), 0};
                              w2D *= padSide;
                              double rangeChamber {(w2D - bragg2D).R()};
                              auto distToSource {(w2D - sourcePoint).R()};
                              return distToSource + rangeChamber;
                          },
                          {"fBraggP", "fWP"})
                  .Define("DeltaT",
                          [&](const ROOT::Math::XYZPointF& braggP, const ROOT::Math::XYZPointF& wP)
                          {
                              // Propagate to source point
                              ActPhysics::Line line {wP, braggP};
                              auto prop {line.MoveToX(sourcePoint.X() / padSide)};
                              // Just Z values!
                              return (double)(braggP.Z() - prop.Z()) * timeFactor * binFactor;
                          },
                          {"fBraggP", "fWP"})};
    def = def.Define("DeltaTSqr", "TMath::Power(DeltaT, 2)")
              .Define("RxySqr", "TMath::Power(Rxy, 2)")
              .Define("DeltaZ", "68.2477*DeltaT");

    // Plot!
    auto hDrift {
        def.Histo2D({"hDrift", "Drift;#Delta t [us];R_{xy} [mm]", 600, -30, 30, 200, 0, 250}, "DeltaT", "Rxy")};

    // Read and process cuts
    ActRoot::CutsManager<int> cuts;
    // cuts.ReadCut(1, TString::Format("./cut/alpha%d_run%d.root", 1, run));
    // cuts.ReadCut(2, TString::Format("./cut/alpha%d_run%d.root", 2, run));
    // cuts.ReadCut(3, TString::Format("./cut/alpha%d_run%d.root", 3, run));
    // Get the TGraphs
    std::vector<ROOT::RDF::RResultPtr<TGraph>> graphs;
    std::vector<double> vdrifts;
    for(int c = 1; c <= 3; c++)
    {
        auto node {def.Filter([&, c](double deltaT, double Rxy)
                              { return IsInsideEllipse(confs[c], deltaT, Rxy) && TMath::Abs(deltaT) <= 6; },
                              {"DeltaT", "Rxy"})};
        graphs.push_back(node.Graph("DeltaTSqr", "RxySqr"));
        graphs.back()->SetTitle(TString::Format("Gated linear drift %d", c));
        graphs.back()->Fit("pol1", "0Q+");
        auto* f {graphs.back()->GetFunction("pol1")};
        if(!f)
            continue;
        auto p1 {f->GetParameter(1) / 100.};
        auto up1 {f->GetParError(1) / 100.};
        auto vdrift {TMath::Sqrt(TMath::Abs(p1))};
        vdrifts.push_back(vdrift);
        std::cout << "-> Alpha " << c << " vdrift = " << vdrift << '\n';
        std::cout << "-> Alpha " << c << " vdrift error = " << up1 / TMath::Sqrt(4 * vdrift) << '\n';
    }
    std::cout << "-> Mean : " << TMath::Mean(vdrifts.begin(), vdrifts.end()) * 10 << " mm / us" << '\n';
    std::vector<TEllipse*> ells;
    for(auto& [alpha, vals] : confs)
    {
        // Lower ellipse
        auto* low {new TEllipse {0, 0, vals[0], vals[1], 60, 120}};
        // Upper ellipse
        auto* up {new TEllipse {0, 0, vals[2], vals[3], 60, 120}};
        // Push
        ells.push_back(low);
        ells.push_back(up);
    }

    auto* c0 {new TCanvas {"c0", "Drift canvas"}};
    c0->DivideSquare(graphs.size() + 1);
    c0->cd(1);
    hDrift->DrawClone("colz");
    // Draw ellipses
    for(auto& e : ells)
    {
        e->SetLineWidth(2);
        e->SetNoEdges(true);
        e->SetFillStyle(0);
        e->Draw("same");
    }
    cuts.DrawAll();
    for(int c = 0; c < graphs.size(); c++)
    {
        c0->cd(c + 2);
        graphs[c]->SetMarkerStyle(6);
        graphs[c]->DrawClone("ap");
        for(auto* o : *graphs[c]->GetListOfFunctions())
            if(o)
                o->DrawClone("same");
    }
    // hLinear->DrawClone("colz");
    // GLinear->DrawClone("ap");
    // hProf->Draw();
}
