#ifndef CalibrationSource_h
#define CalibrationSource_h

#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>
namespace Calibration
{
class Source
{
public:
    typedef std::vector<std::vector<double>> Data;

private:
    std::string fName {};                                                  // Name of source
    Data fEnergies {};                                                     // Energies
    Data fSigmas {};                                                       // Sigma of each peak
    Data fBR {};                                                           // Branching ratio
    std::vector<std::string> fLabels {};                                   // Names of isotope
    std::unordered_map<std::string, std::pair<double, double>> fLimits {}; // Limits of fits per source

public:
    Source(const std::string& name = "3AlphaGanil");

    // Getters
    std::tuple<Data, Data, Data> GetComponents() const { return {fEnergies, fSigmas, fBR}; }
    std::vector<std::string> GetLabels() const { return fLabels; }
    std::unordered_map<std::string, double> GetMajorPeaks() const; // Returns peak of isotope with the largest BR
    std::pair<double, double> GetLimits(const std::string& peak) const { return fLimits.at(peak); }

    void Print() const;

private:
    void Add3AlphaGanil();
    void Add3AlphaUSC();
};
} // namespace Calibration

#endif // !CalibrationSource_h
