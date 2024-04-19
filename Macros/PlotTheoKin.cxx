#include "ActKinematics.h"
#include "ActParticle.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"

void PlotTheoKin()
{
    // Set initial energy at the pad plane
    const double Tini {1.486}; // MeV / u
    // Get beam particle
    ActPhysics::Particle pb {"8He"};

    // (4He, 8He)
    ActPhysics::Kinematics k4he {"8He", "4He", "4He", "8He", Tini * pb.GetAMU()};
    auto* g4he {k4he.GetKinematicLine3()};
    g4he->SetTitle("(#alpha, ^{8}He)");
    // (6He, 6He)
    ActPhysics::Kinematics k6he {"8He", "4He", "6He", "6He", Tini * pb.GetAMU()};
    auto* g6he {k6he.GetKinematicLine3()};
    g6he->SetTitle("(^{6}He, ^{6}He)");
    // Carbon reaction
    ActPhysics::Kinematics kc {"8He", "12C", "8He", "12C", Tini * pb.GetAMU()};
    auto* gc {kc.GetKinematicLine3()};
    gc->SetTitle("(^{8}He, ^{12}C)");

    // Multigraph
    auto* mg {new TMultiGraph};
    mg->SetTitle("Kinematics;#theta_{Lab} [#circ];E_{Lab} [MeV]");
    mg->Add(g4he);
    mg->Add(g6he);
    mg->Add(gc);

    // Plot
    auto* c0 {new TCanvas {"c0", "theoretical kinematics"}};
    mg->Draw("al plc pmc");
    c0->BuildLegend();
}
