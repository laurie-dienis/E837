#ifndef CalibrationRunner_h
#define CalibrationRunner_h

#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1D.h"
#include "TList.h"
#include "TSpectrum.h"

#include "CalibrationSource.h"

#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Calibration
{
class Runner
{
    typedef std::unordered_map<std::string, std::vector<std::shared_ptr<TF1>>> SatelliteCont;
    typedef std::unordered_map<std::string, std::shared_ptr<TF1>> GaussCont;

private:
    Source* fSource {}; // Pointer to previously defined Source obj
    TH1D* fData {};     // Pointer to hist containing channel data (could be rebinned)
    TH1D* fRawData {};  // Same but always unrebinned!
    bool fDebug {};     // Print info of calibration procedure

    // Settings for data
    std::pair<double, double> fRange {}; // Channel range to search for peaks
    // Settings for TSpectrum
    double fSpeSigma {2};
    double fSpeThresh {0.1};
    // Precalibration
    double fPreGaussWidth {15};
    GaussCont fGaussPre {};
    std::shared_ptr<TGraphErrors> fGraphPre {};
    std::shared_ptr<TF1> fCalibPre {};
    std::pair<double, double> fHistOpts {10., 0.015};
    std::shared_ptr<TH1D> fHistPre {};
    // Fine calibration
    GaussCont fGaussFine {};
    std::shared_ptr<TGraphErrors> fGraphFine {};
    std::shared_ptr<TF1> fCalibFine {};
    std::unordered_map<std::string, std::vector<std::string>> fSatelliteStr {};
    SatelliteCont fFineSat {};
    // Final: combination of both calibrations
    std::shared_ptr<TF1> fCalibFinal {};
    std::shared_ptr<TH1D> fHistFinal {};
    GaussCont fGaussFinal {};
    SatelliteCont fFinalSat {};
    // Fit settings
    // Graphs
    // Range is not important as we use the whole graph range
    // (had to delete R bc it caused wrong plotting in some case, altough fit was done perfectly)
    std::string fFitOptsGraph {"0QM+"}; // range determined by graph itself
    std::string fFitOptsGraphDebug {"0M+"};
    // Histograms
    // Using range (R) of each function is important (to fit subranges of histogram)
    std::string fFitOpts {"0QRM+"};
    std::string fFitOptsDebug {"0RM+"};

public:
    Runner(Source* source, TH1D* data, TH1D* originalData = nullptr, bool debug = false)
        : fSource(source),
          fData(data),
          fRawData(originalData),
          fDebug(debug)
    {
        TH1::AddDirectory(false);
    }

    // Setters
    void SetRange(double min, double max)
    {
        fRange = {min, max};
        ApplyRange();
    }
    void SetGaussPreWidth(double w) { fPreGaussWidth = w; }

    // Main execution function
    bool DoIt();

    // Draw
    void Draw(TCanvas* c = nullptr) const { Debug(c); }

    // Method to retrieve resolution
    void PrintRes() const;

private:
    void ApplyRange();
    void DoPreCalibration();
    std::vector<std::pair<double, double>> FilterPeaks(const TSpectrum& spe);
    void FillHistPre();
    void DoFineCalibration();
    std::pair<double, double> GetAmpMeanInRange(TH1D* h, double min, double max);
    void AddSatellite(SatelliteCont& map, const GaussCont& funcs);
    void InitFinalCalib();
    void FillHistFinal();
    void DoFinalPlots();
    void Debug(TCanvas* c = nullptr) const;
    inline void DrawAll(TList* l, bool clone = false) const
    {
        for(auto* o : *l)
            if(o)
            {
                if(clone)
                    o->DrawClone("same");
                else
                    o->Draw("same");
            }
    }
    inline void DrawSat(const SatelliteCont& sats) const
    {
        for(const auto& [_, funcs] : sats)
            for(const auto& f : funcs)
            {
                f->DrawClone("same");
            }
    }
};
} // namespace Calibration

#endif // !CalibrationRunner_h
