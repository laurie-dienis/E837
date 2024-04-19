#ifndef HistConfig_h
#define HistConfig_h

#include "ROOT/RDF/HistoModels.hxx"
namespace HistConfig
{

inline ROOT::RDF::TH2DModel ESilR {"hESilR", "Range vs E_{sil};BSP [mm];E_{Sil} [MeV]", 128, 0, 256, 200, 0, 15};

inline ROOT::RDF::TH2DModel Kin {"hKin", "Kinematics;#theta_{Lab} [#circ];E_{Vertex} [MeV]", 200, 0, 60, 200, 0, 30};

inline ROOT::RDF::TH2DModel ThetaCMECM {
    "hCM", "#theta_{CM} vs E_{CM};#theta_{CM} [#circ];E_{CM} [MeV]", 200, 0, 60, 200, 19500, 19600};

inline ROOT::RDF::TH2DModel ThetaLabEReac {
    "hLab", "#theta_{Lab} vs E*_{^{21}Na,elastic};#theta_{Lab} [#circ];E*_{^{21}Na,elastic} [MeV]", 200, 0, 60, 200, 2, 6};

inline ROOT::RDF::TH1DModel TBeam {"hBeam", "T_{beam};T_{beam} [MeV]", 200, 0, 120};

inline ROOT::RDF::TH2DModel TBeamRPx {"hTBeamX", "T_{Beam} vs RP.X;RP.X() [mm];T_{Beam} [MeV]", 200, 0, 260, 200, 0,
                                      120};

inline ROOT::RDF::TH1DModel Ex {"hEx", "Excitation energy;E_{x} [MeV]", 200, 0, 30};

inline ROOT::RDF::TH2DModel ESilEx {"hESilEx", "E_{x} vs E_{Sil};E_{Sil} [MeV];E_{x} [MeV]", 200, 0, 15, 200, 0, 20};
} // namespace HistConfig

#endif // !HistConfig_h
