#include "ActLine.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TChain.h"
#include "TROOT.h"

#include "Math/Point3D.h"
#include "Math/Point3Dfwd.h"
#include "Math/Vector3D.h"

#include <iostream>

const double kPadSide {2.};      // mm
const double kTimeFactor {0.08}; // tb to us
const int kBinFactor {4};
const ROOT::Math::XYZPointF sourcePoint {-22.5 * kPadSide, 39.5 * kPadSide, 0};

void ScalePoint(ROOT::Math::XYZPointF& p, double padSide, double driftFactor)
{
    p = {p.X() * padSide, p.Y() * padSide, p.Z() * driftFactor};
}

void Drift3D()
{
    // Reading the data
    const int run {28};
    // Setting the drift velocity
    const double vdrift {3.83}; // mm / us
    // Getting the drift factor
    const double driftFactor {kBinFactor * kTimeFactor * vdrift};
    std::cout << "---- 3D range computation ----" << '\n';
    std::cout << "-> Run         : " << run << '\n';
    std::cout << "-> TimeFactor  : " << kTimeFactor << '\n';
    std::cout << "-> Vdrift      : " << vdrift << '\n';
    std::cout << "-> DriftFactor : " << driftFactor << '\n';

    // Get the data
    auto* chain {new TChain("ACTAR_Merged")};
    chain->Add(TString::Format("./../RootFiles/Merger/Merged_Run_%04d.root", run));
    // Enable MT
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};
    // Only events with multiplicity 1
    auto gated {df.Filter("fBraggP.fCoordinates.fX != -1")};

    // Define range in 3D
    auto def {df.Define("R3D",
                        [&](ROOT::Math::XYZPointF& wp, ROOT::Math::XYZPointF& braggp)
                        {
                            // Convert all points to mm
                            // 1-> WP
                            ScalePoint(wp, kPadSide, driftFactor);
                            // 2-> BraggP
                            ScalePoint(braggp, kPadSide, driftFactor);
                            // 3-> Define line
                            ActPhysics::Line line {wp, braggp};
                            auto source {line.MoveToX(sourcePoint.X())};
                            return (source - braggp).R();
                        },
                        {"fWP", "fBraggP"})};
    // Histograms
    auto hR {def.Histo1D({"hR", "Range in 3D;R [mm]", 400, 0, 300}, "R3D")};
    // Plotting
    auto* c0 {new TCanvas {"c0", "3D range"}};
    hR->DrawClone();
}
