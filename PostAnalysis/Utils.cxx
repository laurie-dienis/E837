#ifndef Utils_cxx
#define Utils_cxx

#include "ActColors.h"
#include "ActSilMatrix.h"

#include "TString.h"

#include <iostream>
#include <stdexcept>
#include <string>

namespace E837Utils
{
TString GetFileName(int pipe, int pressure,const std::string& beam, const std::string& target, const std::string& light,
                    const std::string& type = "tree")
{
    auto path {TString::Format("/home/dienis/Analysis_e837/PostAnalysis/RootFiles/Pipe%d/", pipe)};
    TString name {};
    name = TString::Format("%s_beam_%s_target_%s_light_%s_%d_mbar.root", type.c_str(), beam.c_str(), target.c_str(),
                               light.c_str(),pressure);
    return path + name;
}

std::string GetParticleName(const std::string& filename, const std::string& what)
{
    auto base {filename.find(what)};
    auto init {filename.find_first_of('_', base) + 1};
    auto end {filename.find_first_of('_', init)};
    return filename.substr(init, (end - init));
}

struct Signature
{
    std::string beam {};
    std::string target {};
    std::string light {};
    bool isSide {};
};

Signature ExtractSignature(const std::string& file)
{
    auto beam {GetParticleName(file, "beam")};
    auto target {GetParticleName(file, "target")};
    auto light {GetParticleName(file, "light")};
    auto side {TString(file).Contains("side")};
    return {beam, target, light, side};
}

} // namespace E796Utils
#endif // !Utils_cxx
