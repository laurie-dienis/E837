#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TSpectrum.h"
#include "TString.h"
#include "TVirtualPad.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

typedef std::vector<std::vector<double>> Peaks;
typedef std::vector<TGraph*> Graphs;

std::vector<double> SortPeaks(TSpectrum* spe)
{
    std::vector<double> ret;
    auto n {spe->GetNPeaks()};
    for(int p = 0; p < n; p++)
        ret.push_back(spe->GetPositionX()[p]);
    // Sort
    std::sort(ret.begin(), ret.end());
    // Delete baseline
    if(ret.size() > 0)
        ret.erase(ret.begin());
    return ret;
}

std::vector<double> FitPeaks(TH1D* p, const std::vector<double>& peaks)
{
    std::vector<double> ret;
    double width {50};
    for(const auto& peak : peaks)
    {
        auto xmin {peak - width};
        auto xmax {peak + width};
        // auto* f {new TF1 {TString::Format("f%s", p->GetName()), "gaus", 0, 5000}};
        p->Fit("gaus", "0QR", "", xmin, xmax);
        ret.push_back(p->GetFunction("gaus")->GetParameter("Mean"));
    }
    return ret;
}

void FillFitGraph(TGraph* g, const std::vector<double>& x, const std::vector<double>& y, std::ofstream& ofs)
{
    // Assert same size
    if(!(x.size() == y.size()))
    {
        std::cout << "Different sizes for " << g->GetTitle() << '\n';
        // Dummy matching in case of wrong channel
        ofs << 0 << " " << 0 << " " << 0 << '\n';
        return;
    }
    for(int p = 0; p < x.size(); p++)
        g->SetPoint(g->GetN(), x[p], y[p]);
    // Fit
    TString func {"pol2"};
    g->Fit(func, "0Q+");
    auto* f {g->GetFunction(func)};
    // Write parameters
    ofs << f->GetParameter(0) << " " << f->GetParameter(1) << " " << f->GetParameter(2) << '\n';
}

void DoGainMatching()
{
    // Read histogram
    auto* f {new TFile {"./Inputs/gain.root"}};
    auto* h {f->Get<TH2D>("h")};
    if(!h)
        throw std::runtime_error("Gain histogram could be read");

    const int nchannels {17408};
    auto* gpeaks {new TGraph};
    gpeaks->SetTitle(";Channel;Peaks");
    // 1-> Find peaks
    std::vector<TH1D*> projs;
    Peaks peaks;
    for(int c = 0; c < nchannels; c++)
    {
        // Get projection
        // bin starts at 1!
        auto p {h->ProjectionY(TString::Format("p%d", c), c + 1, c + 1)};
        // Find peaks
        TSpectrum spe {11};
        spe.Search(p, 2, "nodraw", 0.05);
        // Sort and clean
        auto sorted {SortPeaks(&spe)};
        // Fit with gaussian
        peaks.push_back(FitPeaks(p, sorted));
        // Fill graph
        for(const auto& peak : peaks.back())
            gpeaks->SetPoint(gpeaks->GetN(), c, peak);
        // save projections
        projs.push_back(p);
    }

    // 2-> Set reference pad and fit
    const int ref {4196};
    Graphs gfits;
    std::ofstream streamer {"./Outputs/gain_matching_v0.dat"};
    for(int c = 0; c < nchannels; c++)
    {
        // Init graph
        gfits.push_back(new TGraph);
        auto& g {gfits.back()};
        g->SetTitle(TString::Format("Fit for %d", c));
        // Fill graph
        FillFitGraph(g, peaks[c], peaks[ref], streamer);
    }
    streamer.close();


    // Plot
    auto* c0 {new TCanvas {"c0", "Gain matching"}};
    h->Draw("colz");

    int toPlot {15};
    auto* c1 {new TCanvas {"c1", "Inspect projs"}};
    c1->DivideSquare(toPlot);
    for(int p = 0; p < toPlot; p++)
    {
        c1->cd(p + 1);
        projs[p]->Draw();
        for(auto* o : *(projs[p]->GetListOfFunctions()))
            if(o)
                o->Draw("same");
    }

    auto* c2 {new TCanvas {"c2", "Peaks by channel"}};
    gpeaks->SetMarkerStyle(6);
    gpeaks->Draw("ap");

    auto* c3 {new TCanvas {"c3", "Inspect fits"}};
    c3->DivideSquare(toPlot);
    for(int p = 0; p < toPlot; p++)
    {
        c3->cd(p + 1);
        gfits[p]->SetMarkerStyle(24);
        gfits[p]->Draw("ap");
        for(auto* o : *(gfits[p]->GetListOfFunctions()))
            if(o)
                o->Draw("same");
    }
}
