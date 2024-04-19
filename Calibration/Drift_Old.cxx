#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActSRIM.h"
#include "ActLine.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TROOT.h"
#include "TChain.h"

#include "Math/Point3Dfwd.h"

#include <fstream>

const double padSide {2.};      // mm
const double timeFactor {0.08}; // tb to us
const int binFactor {4};
const ROOT::Math::XYZPointF sourcePoint {-22.5 * padSide, 39.5 * padSide, 0};

void Drift()
{
    // Enable MT
    ROOT::EnableImplicitMT();
    // Reading the data
    const int run {17};
    std::cout<<"run number = "<< run<<'\n';
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
    auto hDrift {def.Histo2D({"hDrift", "Drift;#Delta t [us];R_{xy} [mm]", 600, -30, 30, 200, 0, 250}, "DeltaT", "Rxy")};

    // Read and process cuts
    ActRoot::CutsManager<int> cuts;
    cuts.ReadCut(1, TString::Format("./cut/alpha%d_run%d.root", 1, run));
    cuts.ReadCut(2, TString::Format("./cut/alpha%d_run%d.root", 2, run));
    cuts.ReadCut(3, TString::Format("./cut/alpha%d_run%d.root", 3, run));
    // Get the TGraphs
    std::vector<ROOT::RDF::RResultPtr<TGraph>> graphs;
    if(cuts.GetCut(1) && cuts.GetCut(2) && cuts.GetCut(3)){
      for(int c = 1; c <= 3; c++)
	{
	  auto node {def.Filter([&, c](double deltaT, double Rxy){return cuts.IsInside(c, deltaT, Rxy);}, {"DeltaT","Rxy"})};
	  graphs.push_back(node.Graph("DeltaTSqr","RxySqr"));
	  graphs.back()->SetTitle(TString::Format("Gated linear drift %d",c));
	  graphs.back()->Fit("pol1", "0Q+");
	  auto* f {graphs.back()->GetFunction("pol1")};
	  auto p1 {f->GetParameter(1)/100.};
	  auto up1 {f->GetParError(1)/100.};
	  auto vdrift {TMath::Sqrt(TMath::Abs(p1))};
	  std::cout<<"-> Alpha "<<c<<" vdrift = "<<vdrift<<'\n';
	  std::cout<<"-> Alpha "<<c<<" vdrift error = "<<up1/TMath::Sqrt(4*vdrift)<<'\n';
	}
    }
    // auto cut {def.Filter("DeltaT>=-0.8 && DeltaT <= 0.8")};
    // auto hLinear {
    //   cut.Histo2D({"hLinear", "Linear drift;R_{xy}^{2} [mm^{2}];#Delta t^{2} [us^{2}]", 625, 0, 8, 600, 0, 50000},
    // 		  "DeltaTSqr", "RxySqr")};
    // auto GLinear {cut.Graph("RxySqr","DeltaTSqr")};
    // auto hProf {hLinear->ProfileY()};
    //
    // Debug
    // ActRoot::CutsManager<int> cuts;
    // cuts.ReadCut(0, "./debug.root");
    // {
    //     std::ofstream streamer {"debug.dat"};
    //     def.Foreach(
    //         [&](const ActRoot::MergerData& d, double rxy, double deltaT)
    //         {
    // 	      if(cuts.IsInside(0, deltaT, rxy))
    //                 streamer << d.fRun << " " << d.fEntry << '\n';
    //         },
    //         {"MergerData", "Rxy", "DeltaT"});
    //     streamer.close();
    // }
    auto* c0 {new TCanvas {"c0", "Drift canvas"}};
    c0->DivideSquare(graphs.size() + 1);
    c0->cd(1);
    hDrift->DrawClone("colz");
    cuts.DrawAll();
    for(int c = 0; c < graphs.size(); c++){
      c0->cd(c + 2);
      graphs[c]->DrawClone("ap");
      for(auto* o : *graphs[c]->GetListOfFunctions())
	if(o)
	  o->DrawClone("same");
    }
    //hLinear->DrawClone("colz");
    //GLinear->DrawClone("ap");
    //hProf->Draw();
}
