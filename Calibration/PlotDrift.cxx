#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TMultiGraph.h"

#include <string>
#include <unordered_map>
#include <vector>

TGraphErrors* GetMean(std::vector<TGraphErrors*> gs)
{
    auto* ret {new TGraphErrors};
    ret->SetNameTitle("Mean", "Mean");
    ret->SetLineWidth(2);
    ret->SetMarkerStyle(24);
    // Run for each point
    for(int p = 0; p < gs.front()->GetN(); p++)
    {
        // Fill vectors
        std::vector<double> vdrifts;
        for(auto& g : gs)
            vdrifts.push_back(g->GetPointY(p));
        // Mean
        auto mean {TMath::Mean(vdrifts.begin(), vdrifts.end())};
        auto std {TMath::StdDev(vdrifts.begin(), vdrifts.end())};
        // Set points
        ret->SetPoint(p, gs.front()->GetPointX(p), mean);
        auto umean {std / TMath::Sqrt(gs.size())};
        ret->SetPointError(p, 0, umean);
    }
    return ret;
}
void PlotDrift()
{
    std::unordered_map<std::string, std::string> files {{"#alpha_{1}", "./Inputs/drift_alpha1.dat"},
                                                        {"#alpha_{2}", "./Inputs/drift_alpha2.dat"},
                                                        {"#alpha_{3}", "./Inputs/drift_alpha3.dat"}};

    // Create multigraphs
    // 1-> Alpha 1 2 3 to mean comparison
    auto* mg {new TMultiGraph};
    mg->SetTitle(";E [V / cm / mbar];v_{drift} [cm / #mus]");
    // 2-> Mean to simu comp
    auto* ms {new TMultiGraph};
    ms->SetTitle(";E [V / cm / mbar];v_{drift} [cm / #mus]");
    std::vector<TGraphErrors*> gs;
    for(const auto& [name, file] : files)
    {
        auto* g {new TGraphErrors {file.c_str(), "%lg %lg %lg"}};
        g->SetNameTitle(name.c_str(), name.c_str());
        g->SetLineWidth(2);
        // g->SetMarkerStyle(24);
        mg->Add(g, "lp");
        gs.push_back(g);
    }
    // Get mean
    auto* gmean {GetMean(gs)};
    mg->Add(gmean, "lp");
    ms->Add(gmean, "lp");
    // Simulation
    auto* gsimu {new TGraphErrors {"./Inputs/drift_simu.dat", "%lg %lg"}};
    gsimu->SetNameTitle("Simu", "Simu");
    gsimu->SetLineWidth(2);
    gsimu->SetMarkerStyle(25);
    ms->Add(gsimu, "lp");

    // Plot
    auto* c0 {new TCanvas {"c0", "Drift velocity canvas"}};
    c0->DivideSquare(2);
    c0->cd(1);
    mg->Draw("apl plc pmc");
    c0->cd(1)->BuildLegend();
    c0->cd(2);
    ms->Draw("apl plc pmc");
    c0->cd(2)->BuildLegend();
}
