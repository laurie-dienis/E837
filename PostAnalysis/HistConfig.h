#ifndef HistConfig_h
#define HistConfig_h

#include "ROOT/RDF/HistoModels.hxx"

#include "TString.h"

namespace HistConfig
{
using namespace ROOT::RDF;

const TH2DModel PID {"hPID", "PID;E_{Sil} [MeV];Q_{ave} [mm^{-1}]", 300, 0, 15, 800, 0, 2000};

const TH2DModel PID_selected {
    "hPID_selected", "PID_selected;E_{Sil} [MeV];Q_{ave} [mm^{-1}]", 300, 0, 15, 800, 0, 2000};

const TH2DModel Kine {"hKine", "Kine;E_{Sil} [MeV];#theta_{Lab} [#circ]", 300, 0, 100, 180, 0, 40};

const TH2DModel Angle {
    "hAngle", "Angle;#theta_{3 = hits sil} [#circ];#theta_{4 = does not hit sil} [#circ]]", 400, 0, 180, 400, 0, 180};

const TH1DModel ProjQ {"ProjQ", "ProjQ;QProjX", 300, 0, 400};

const TH2DModel Ex21Nalab {
    "hExlab", "#theta_{Lab} vs E*_{^{12}Be,elastic};#theta_{Lab} [#circ];E*_{^{12}Be,elastic} [MeV]",
    200,      0,
    60,       200,
    5,        25};

const TH2DModel Ex21Na {
    "hEx", "#theta_{CM} vs E*_{^{12}Be,elastic};#theta_{CM} [#circ];E*_{^{12}Be,elastic} [MeV]", 200, 50, 180, 200, 10,
    15};

const TH2DModel Ex21Na2 {
    "hEx2", "#theta_{Lab} vs E*_{^{12}Be,elastic};E*_{^{12}Be,elastic} [MeV];#theta_{CM} [#circ]", 200, 10, 15, 200, 90,
    180};

const TH2DModel EexSiVsRange {"hEex_Si_vs_range", "Ex(8He);E_{x}(8He) computed with range (MeV);E_{x}(8He) computed with Si (MeV);", 500, -2, 4, 500, -2, 4};

const TH2DModel EbeamSiVsRange {"hEbeam_Si_vs_range", "Ebeam;E_{beam} computed with range (MeV);E_{beam} computed with Si (MeV);", 500, 5, 13, 500, 5, 13};

const TH2DModel PIDTwo {"hPIDUncal", "PIDTWO;E_{Sil} [Channel];Q_{ave} [mm^{-1}];", 500, 0, 15000, 1000, 0, 3000};

const TH2DModel SP {"hSP", "SP;Y [mm];Z [mm]", 200, -10, 300, 200, -10, 500};

const TH2DModel RP {"hRP", "RP;X [mm];Y [mm]", 200, -10, 300, 200, -10, 300};

const TH2DModel TPC {"hTPC", "TPC;X [mm];Y [mm]", 200, -10, 300, 200, -10, 300};

const TH1DModel TL {"hTL", "Track length; TL [mm]", 300, 0, 600};

const TH1DModel Ereac {"hEreac", "Ereac; Ereac [MeV]", 100, 5, 20};

const TH1DModel Ediff {"hEdiff", "Ediff; Ereac_diff [MeV]", 50, 0, 10};

const TH1DModel Elight {"hElight", "Elight; Elight [MeV]", 50, 0, 100};

const TH2DModel Kin {"hKin", "Kinematics;#theta_{Lab} [#circ];E_{Vertex} [MeV]", 250, 0, 60, 250, 0, 20};

const TH2DModel KinEl {"hKinEl", "Lab kinematics;#theta_{Lab} [#circ];E_{Vertex} [MeV]", 600, 0, 180, 300, 0, 20};

const TH2DModel KinSimu {"hKin", "Simulation kinematics;#theta_{Lab} [#circ];E_{Vertex} [MeV]", 400, 0, 90, 300, 0, 40};

const TH2DModel KinCM {"hKinCM", "CM kinematics;#theta_{CM} [#circ];E_{Vertex} [MeV]", 400, 0, 60, 400, 0, 20};

const TH1DModel Ex {"hEx", "Excitation energy of 8He;E_{x} [MeV];Counts", 300, -5, 10};

const TH1DModel Ex2 {"hEx2", "Excitation energy of 8He;E_{x} [MeV];Counts", 300, -0.01, 0.01};

const TH1DModel ThetaCM {"hThetaCM", "ThetaCM;#theta_{CM} [#circ]", 600, 0, 180};

const TH2DModel ZThetaZ {"hZThetaZ", "Emittance along Z;Z [mm];#theta_{Z} [#circ]", 600, 0, 270, 600, -10, 10};

const TH2DModel YPhiY {"hYPhiY", "Emittance along Y;Y [mm];#phi_{Y} [#circ]", 600, 0, 270, 600, -10, 10};

const TH2DModel ThetaBeam {
    "hThetaBeam", "#theta_{Beam} against RP.X;RP.X() [mm];#theta_{Beam} [#circ]", 200, -5, 270, 200, -1, 10};

const TH2DModel ExZ {"hExZ", "E_{x} dependence on SP.Z();SP.Z() [mm];E_{x} [MeV]", 200, -10, 300, 200, -10, 20};

const TH2DModel ExThetaCM {
    "hExThetaCM", "E_{x} vs #theta_{CM};#theta_{CM} [#circ];E_{x} [MeV]", 400, 0, 60, 200, -10, 20};

const TH2DModel ExThetaLab {
    "hExThetaLab", "E_{x} vs #theta_{Lab};#theta_{Lab} [#circ];E_{x} [MeV]", 400, 0, 60, 200, -10, 20};

const TH2DModel ExRPZ {"hExRPZ", "E_{x} vs RP.Z;RP.Z() [mm];E_{x} [MeV]", 200, -10, 300, 200, -10, 20};

const TH2DModel ThetaHeavyLight {
    "hThetaHL", "#theta heavy vs light;#theta_{Light} [#circ];#theta_{Heavy} [#circ]", 400, 0, 60, 400, 0, 60};

const TH2DModel ThetaCMLab {
    "hThetaCMLab", "CM vs Lab correlations;#theta_{Lab} [#circ];#theta_{CM} [#circ]", 400, 0, 60, 400, 0, 60};

const TH2DModel RPxThetaCM {
    "hRPxThetaCM", "RP.X vs #theta_{CM} correlations;RP.X [mm];#theta_{CM} [#circ]", 200, 0, 300, 100, 0, 60};

template <typename T>
T ChangeTitle(T model, const TString& title, const TString& label = "");
} // namespace HistConfig

template <typename T>
T HistConfig::ChangeTitle(T model, const TString& title, const TString& label)
{
    auto ret {model};
    if(label.Length() > 0)
        ret.fName = model.fName + label;
    TString old {model.fTitle};
    auto end {old.First(';')};
    TString nt {title + old(end, old.Length() - end)};
    ret.fTitle = nt;
    return ret;
}

#endif
