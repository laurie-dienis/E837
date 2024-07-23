#ifndef Profile_cxx
#define Profile_cxx
#include "Profile.h"

#include "TF1.h"
#include "TGraph.h"
#include "TH1.h"
#include "TSpline.h"

#include <utility>

// Implement functions
TH1D* ChargeProfile::GetDerivative(TH1* h)
{
    // Build the derivative
    // Using EXACTLY same binning as input profiel
    auto* hd {new TH1D {"hd", "Derivative of profile;X [pad];Derivative", h->GetNbinsX(), h->GetXaxis()->GetXmin(),
                        h->GetXaxis()->GetXmax()}};
    // Bin 0 is underflow and bin max + 1 is overflow, hence data is in [1, max]
    for(auto bin = 1; bin <= h->GetNbinsX(); bin++)
    {
        double diff {};
        if(bin > 1)
        {
            auto y1 {h->GetBinContent(bin)};
            auto y0 {h->GetBinContent(bin - 1)};
            diff = (y1 - y0);
        }
        hd->SetBinContent(bin, diff);
    }
    return hd;
}

std::pair<float, float> ChargeProfile::GetVertexAndQmax(TH1* der, TH1* prof, double width)
{
    // der is the Derivative
    // prof is the Profile of charge
    // 1-> Find the maximum
    auto bmax {der->GetMaximumBin()};
    auto xmax {der->GetBinCenter(bmax)}; // xmax should be the vertex!
    // 2-> Find the Qmax within a given width
    prof->GetXaxis()->SetRangeUser(xmax - width, xmax + width);
    auto bmax2 {prof->GetMaximumBin()};
    auto ymax {
        prof->GetBinContent(bmax2)}; // ymax should be the Qmax (ofc computed in the profile hist, not derivative)
    // 3-> Reset the range
    prof->GetXaxis()->SetRangeUser(prof->GetXaxis()->GetXmin(), prof->GetXaxis()->GetXmax());
    return {(float)xmax, (float)ymax};
}

float ChargeProfile::GetStoppingPoint(TH1* prof)
{
    // 1-> Find the GLOBAL maximum of the profile histogram
    auto bmax {prof->GetMaximumBin()};
    auto xmax {prof->GetBinCenter(bmax)};
    auto ymax {prof->GetBinContent(bmax)};
    // 2-> Set charge to seek for
    auto range {ymax / 5};
    // 3-> Create the Spline to interpolate
    auto spe {std::make_unique<TSpline3>(prof)};
    // And now function
    auto func {std::make_unique<TF1>(
        "func", [&](double* x, double* p) { return spe->Eval(x[0]); }, 0, prof->GetXaxis()->GetXmax(), 1)};
    // Find maximum in the range [xMax, xRangeMax of histogram]
    auto ret {func->GetX(range, xmax, prof->GetXaxis()->GetXmax())};
    return ret;
}

#endif
