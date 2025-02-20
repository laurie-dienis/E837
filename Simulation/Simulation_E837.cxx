#include "ActColors.h"
#include "ActGeometry.h"
#include "ActKinematicGenerator.h"
#include "ActKinematics.h"
#include "ActParticle.h"
#include "ActRunner.h"
#include "ActSRIM.h"

#include "Rtypes.h"

#include "TCanvas.h"
#include "TEfficiency.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TProfile2D.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TTree.h"

#include "Math/Point3Dfwd.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <utility>

#include "HistConfig.h"

void set_plot_style() {
  const Int_t NRGBs = 5;
  const Int_t NCont = 455;

  Double_t stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
  Double_t red[NRGBs] = {0.00, 0.00, 0.87, 1.00, 0.51};
  Double_t green[NRGBs] = {0.00, 0.81, 1.00, 0.20, 0.00};
  Double_t blue[NRGBs] = {0.51, 1.00, 0.12, 0.00, 0.00};
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}

void Simulation_E837(const std::string &beam, const std::string &target,
                     const std::string &light, const std::string &heavy,
                     int neutronPS, int protonPS, double T1, double Ex,
                     int pressure, bool standalone) {
  // set batch mode if not an independent function
  if (!standalone)
    gROOT->SetBatch(true);

  double angle_min{130.};
  double angle_max{140.};

  TRandom random;

  // SIGMAS
  const double sigmaSil{0};
  const double sigmaPercentBeam{0};
  const double sigmaAngleLight{0.95 / 2.355};
  // Parameters of beam in mm
  // Beam has to be manually placed in the simulation
  // Centered in Z and Y with a width of 4 mm
  // Center in Z
  const double zVertexMean{128. + 18.}; // beam not centered in chamber, upwards
  const double zVertexSigma{4};
  // Center in Y
  const double yVertexMean{128.};
  const double yVertexSigma{4};
  // const double zVertexMean {83.59};
  // const double zVertexSigma {3.79};

  // THRESHOLDS FOR SILICONS
  const double thresholdSi0{1.};
  const double thresholdSi1{1.};

  // number of iterations
  const int iterations{static_cast<int>(16000000)};

  // ACTIVATE STRAGGLING OR NOT
  bool stragglingInGas{false};
  bool stragglingInSil{false};
  bool silResolution{false};
  bool thetaResolution{true};
  bool slowDownBeam{true};

  // CUTS ON SILICON ENERGY, depending on particle
  std::pair<double, double> eLoss0Cut{0, 1000}; //{6.5, 27.};

  //---- SIMULATION STARTS HERE
  // ROOT::EnableImplicitMT();

  // timer
  TStopwatch timer{};
  timer.Start();

  // Init particles
  ActPhysics::Particle p1{beam};
  ActPhysics::Particle p2{target};
  ActPhysics::Particle p3{light};
  ActPhysics::Particle p4{heavy};
  // Init kinematics generator
  ActSim::KinematicGenerator kingen{p1, p2, p3, p4, protonPS, neutronPS};
  kingen.Print();
  // Get threshold
  auto T1Thresh{ActPhysics::Kinematics(p1, p2, p3, p4, -1, Ex).GetT1Thresh()};

  // Histograms
  // To compute a fine-grain efficiency, we require at least a binning width of
  // 0.25 degrees!
  auto hThetaCM{HistConfig::ThetaCM.GetHistogram()};
  auto hThetaCMAll{
      HistConfig::ChangeTitle(HistConfig::ThetaCM, "ThetaCM all", "All")
          .GetHistogram()};

  auto hDistL0{
      HistConfig::ChangeTitle(HistConfig::TL, "Distance to L0").GetHistogram()};

  auto hKinVertex{
      HistConfig::ChangeTitle(HistConfig::KinSimu, "Kinematics at vertex")
          .GetHistogram()};

  auto hSP{HistConfig::SP.GetHistogram()};
  auto hSPCut{(TH2D *)hSP->Clone("hSPCut")};

  auto hEexBefore{
      HistConfig::ChangeTitle(HistConfig::Ex, "Ex before resolutions", "Bef")
          .GetHistogram()};
  auto hEexAfter{
      HistConfig::ChangeTitle(HistConfig::Ex, "Ex after resolutions", "After")
          .GetHistogram()};

  auto hSPTheta{std::make_unique<TProfile2D>(
      "hSPTheta", "SP vs #theta_{CM};Y [mm];Z [mm];#theta_{CM} [#circ]", 75, 0,
      300, 75, 0, 300)};

  auto hRP{HistConfig::RP.GetHistogram()};

  auto *hRPxSimu{HistConfig::RP.GetHistogram()->ProjectionX("hRPxSimu")};

  auto *hRPxRange{new TH2D{
      "hRPxRange",
      "Light range vs dist to sil wall;Dist to wall [mm];Light range [mm]", 300,
      0, 500, 300, 0, 500}};

  auto *hEBeam{new TH2D{"hEBeam", "Range angle", 300, 0, 180, 300, 0, 500}};

  auto *hPID{
      new TH2D{"hPID", "PID;E[MeV];Eloss [MeV]", 100, 0, 15, 100, 0, 0.1}};

  auto *hEx2{
      new TH2D{"hEx2",
               "#theta_{Lab} vs E*_{^{12}Be,elastic};E*_{^{12}Be,elastic} "
               "[MeV];#theta_{CM} [#circ]",
               200, 10, 15, 200, 90, 180}};

  auto *hnorm{new TH2D{
      "hnorm",
      "geometric efficiency;E*_{^{12}Be,elastic} [MeV];#theta_{CM} [#circ]", 90,
      10, 15, 90, 90, 180}};

  auto *hnorm_phi{new TH2D{
      "hnorm_phi",
      "geometric efficiency;E*_{^{12}Be,elastic} [MeV];#theta_{CM} [#circ]", 90,
      10, 15, 90, 0, 200}};

  auto *hECM{
      new TH1D{"hECM", "ECM + threshold;E_{CM} + Thresh [MeV]", 100, 8, 18}};

  auto *hprojection{new TH1D{"hprojection", "Ex12Be_norm", 800, 12.2, 12.4}};

  auto *hEx_8He{
      new TH1D{"hEx8He", "Ex_He; E_{x} [MeV]; Counts", 300, -0.1, 0.1}};

  auto *hKin_6He{new TH2D{"hkin6He",
                          "kinematics;#theta_{lab} [#circ];E_{lab} [MeV]", 90,
                          0, 50, 90, 0, 9}};

  auto *htrack{new TH2D{
      "htrack", "trackslength;#theta_{heavy lab} [#circ];track length [mm]", 90,
      0, 100, 90, 0, 250}};

  // Load SRIM tables
  // The name of the file sets particle + medium
  auto *srim{new ActPhysics::SRIM()};
  srim->ReadTable("light",
                  TString::Format("../Inputs/SRIM/%s_He97_butane_%dmbar.txt",
                                  light.c_str(), pressure)
                      .Data());
  srim->ReadTable("beam",
                  TString::Format("../Inputs/SRIM/%s_He97_butane_%dmbar.txt",
                                  beam.c_str(), pressure)
                      .Data());
  srim->ReadTable(
      "lightInSil",
      TString::Format("../Inputs/SRIM/%s_silicon.txt", light.c_str()).Data());

  // Load geometry
  auto *geometry{new ActSim::Geometry()};
  if (pressure == 900) {
    geometry->ReadGeometry("Geometry/", "e837");
  } else {
    geometry->ReadGeometry("Geometry/", "e837");
  }

  // Random generator
  auto *rand{new TRandom3()};
  rand->SetSeed(); // random path in each execution of macro

  // Runner: contains utility funcstions to execute multiple actions
  ActSim::Runner runner(srim, geometry, rand, sigmaSil);

  // Output from simulation!
  // We only store a few things in the TTree
  // 1-> Excitation energy
  // 2-> Theta in CM frame
  // 3-> Weight of the generator: for three-body reactions (phase spaces) the
  // other two variables need to be weighted by this value. For binary
  // reactions, weight = 1 4-> Energy at vertex 5-> Theta in Lab frame
  auto *outFile{new TFile(
      TString::Format("../PostAnalysis/Input/Norm_%s_%dmbar_%.0f-%.0fdeg.root",
                      light.c_str(), pressure, angle_min, angle_max),

      "recreate")};
  auto *outTree{
      new TTree("SimulationTTree",
                "A TTree containing only our Eex obtained by simulation")};
  double theta3CM_tree{};
  outTree->Branch("theta3CM", &theta3CM_tree);
  double Eex_tree{};
  outTree->Branch("Eex", &Eex_tree);
  double weight_tree{};
  outTree->Branch("weight", &weight_tree);
  double EVertex_tree{};
  outTree->Branch("EVertex", &EVertex_tree);
  double theta3Lab_tree{};
  outTree->Branch("theta3Lab", &theta3Lab_tree);
  double rpx_tree{};
  outTree->Branch("RPx", &rpx_tree);

  // RUN!
  // print fancy info
  std::cout << BOLDMAGENTA << "Running for Ex = " << Ex << " MeV" << RESET
            << '\n';
  std::cout << BOLDGREEN;
  const int percentPrint{5};
  int step{iterations / (100 / percentPrint)};
  int nextPrint{step};
  int percent{};
  double maxBinContent = -1;
  for (long int reaction = 0; reaction < iterations; reaction++) {
    // Print progress
    if (reaction >= nextPrint) {
      percent = 100 * (reaction + 1) / iterations;
      int nchar{percent / percentPrint};
      std::cout << "\r" << std::string((int)(percent / percentPrint), '|')
                << percent << "%";
      std::cout.flush();
      nextPrint += step;
    }
    // 1-> Sample vertex and apply same cut as in analysis
    auto vertex{runner.SampleVertex(yVertexMean, yVertexSigma, zVertexMean,
                                    zVertexSigma, nullptr)};
    if (vertex.second.X() < 127.8 || vertex.second.X() > 128.3) {
      continue;
    }
    // 2-> Beam energy according to its sigma
    auto TBeam_i = T1 * p1.GetAMU(); // T1 in Mev / u * mass of beam
    double TBeam;
    double TBeam_u = -11;
    // in u = total kinetic energy
    // Slow down it according to vertex position
    std::cout << "vertex = " << vertex.second.X() << "\n";
    double vertex_gauss = random.Gaus(vertex.second.X(), 0.3);

    std::cout << "vertex_gauss = " << vertex_gauss << "\n";

    if (slowDownBeam) {
      TBeam = srim->Slow("beam", TBeam_i, vertex.second.X(), 0.);
      if (vertex_gauss < 0) {
        TBeam_u = srim->Slow("beam", TBeam_i, 0., 0.);
      } else {
        TBeam_u = srim->Slow("beam", TBeam_i, vertex_gauss, 0.);
      }
    }
    std::cout << "tbeam= " << TBeam * 1e3 << "\n";
    TBeam = TBeam_u;
    std::cout << "tbeam u  = " << TBeam_u * 1e3 << "\n";
    // runner energy functions return std::nan when the particle is stopped in
    // the gas! if nan (aka stopped in gas, continue) if not stopped but beam
    // energy below kinematic threshold, continue
    if (std::isnan(TBeam) || TBeam < T1Thresh)
      continue;

    // 3-> Run kinematics!
    kingen.SetBeamAndExEnergies(TBeam, Ex);
    double weight{kingen.Generate()};
    // focus on recoil 3 (light)
    // obtain thetas and energies
    auto *PLight{kingen.GetLorentzVector(0)};
    // auto *PHeavy{kingen.GetLorentzVector(1)};
    auto theta3Lab{PLight->Theta()};
    std::cout << "theta3 = " << theta3Lab << "\n";
    std::cout << "ebeam = " << TBeam << "\n";
    // auto theta4Lab{PHeavy->Theta()};
    // auto phi3Lab{PLight->Phi()};
    auto T3Lab{PLight->Energy() - p3.GetMass()};
    std::cout << "talpha= " << T3Lab << "\n";
    std::cout << "\n";
    // auto T4Lab{PHeavy->Energy() - p4.GetMass()};
    auto EexBefore{kingen.GetBinaryKinematics().ReconstructExcitationEnergy(
        T3Lab, theta3Lab)};
    // to compute geometric efficiency by CM interval and with our set reference
    // direction
    double theta3CMBefore{
        kingen.GetBinaryKinematics().ReconstructTheta3CMFromLab(T3Lab,
                                                                theta3Lab)};
    hThetaCMAll->Fill(theta3CMBefore * TMath::RadToDeg());
    // 4-> Include thetaLab resolution to compute thetaCM and Ex
    if (thetaResolution) // resolution in
      theta3Lab = runner.GetRand()->Gaus(theta3Lab,
                                         sigmaAngleLight * TMath::DegToRad());
    auto phi3Lab{runner.GetRand()->Uniform(0., 2 * TMath::Pi())};
    // Eval range of light particle
    auto lightRange{srim->EvalDirect("light", T3Lab)};
    // auto heavyRange{srim->EvalDirect("heavy", T4Lab)};
    // Ecm from formula using masses
    auto Ecm{(p2.GetAMU() / (p2.GetAMU() + p1.GetAMU())) * TBeam};
    // double theta_gauss = random.Gaus(theta3Lab, 0.1);

    auto Ex_calc =
        T3Lab *
        ((2 * TMath::Cos(theta3Lab) *
          TMath::Sqrt((TBeam_u * p3.GetAMU()) / (T3Lab * p4.GetAMU()))) -
         ((p3.GetAMU() / p4.GetAMU()) + 1));

    // random Ecm for geometric efficiency
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distri(10.0, 15.0);
    double Ecm_norm = 0;
    // std::cout << "ecm_norm = "<< Ecm_norm<<"\n";

    // 5-> Propagate track from vertex to silicon wall using Geometry class
    ROOT::Math::XYZVector direction{
        TMath::Cos(theta3Lab), TMath::Sin(theta3Lab) * TMath::Sin(phi3Lab),
        TMath::Sin(theta3Lab) * TMath::Cos(phi3Lab)};
    auto vertexInGeoFrame{runner.DisplacePointToTGeometryFrame(vertex.second)};
    ROOT::Math::XYZPoint silPoint0{};
    int silType0{};
    int silIndex0{};
    double distance0{};
    bool side0{};
    // Assembly 0
    int assemblyIndex{0};
    runner.GetGeo()->PropagateTrackToSiliconArray(
        vertexInGeoFrame, direction, assemblyIndex, side0, distance0, silType0,
        silIndex0, silPoint0);
    // convert to mm (Geometry::PropagateTracksToSiliconArray works in cm but we
    // need mm to use in SRIM)
    distance0 *= 10.;

    // skip tracks that doesn't reach silicons or are in silicon index cut
    if (silIndex0 == -1)
      continue;
    auto silPoint0InMM{runner.DisplacePointToStandardFrame(silPoint0)};

    // Range-distance evaluation for light particle
    hRPxRange->Fill(distance0, lightRange);
    hEBeam->Fill(theta3Lab * TMath::RadToDeg(), lightRange);
    auto eloss{srim->EvalInverse("light", distance0)};
    // if particle doest reach the silicon
    if (lightRange < distance0)
      hSPCut->Fill(silPoint0InMM.Y(), silPoint0InMM.Z());
    if (lightRange > distance0)
      Ecm_norm = distri(gen);
    auto T3EnteringSil{
        runner.EnergyAfterGas(T3Lab, distance0, "light", stragglingInGas)};
    // nan if stopped in gas
    if (!std::isfinite(T3EnteringSil))
      continue;

    // SILICON0
    auto [eLoss0, T3AfterSil0]{runner.EnergyAfterSilicons(
        T3EnteringSil, geometry->GetAssemblyUnitWidth(0) * 10., thresholdSi0,
        "lightInSil", silResolution, stragglingInSil)};
    // nan if bellow threshold
    if (!std::isfinite(eLoss0))
      continue;
    // 7->
    // we are ready to reconstruct Eex with all resolutions implemented
    //(d,light) is investigated gating on Esil1 = 0!
    bool cutEAfterSil0{T3AfterSil0 == 0.};
    bool cutELoss0{eLoss0Cut.first <= eLoss0 && eLoss0 <= eLoss0Cut.second};
    if (cutEAfterSil0 && cutELoss0) // fill histograms
    {
      auto T3Recon{runner.EnergyBeforeGas(eLoss0, distance0, "light")};
      auto EexAfter{kingen.GetBinaryKinematics().ReconstructExcitationEnergy(
          T3Recon, theta3Lab)};
      auto theta3CM{kingen.GetBinaryKinematics().ReconstructTheta3CMFromLab(
          T3Recon, theta3Lab)};

      // fill histograms
      // hThetaCM->Fill(theta3CM * TMath::RadToDeg());
      // TF1 *fit = new TF1("fit","[0]*(x^[1])",10, -0.5);
      // hPID->Fill(T3Lab, eloss/distance0);
      // hThetaCM->Fill(theta3CMBefore * TMath::RadToDeg());
      hEexBefore->Fill(
          EexBefore,
          weight); // with the weight from each TGenPhaseSpace::Generate()
      hDistL0->Fill(distance0);
      hKinVertex->Fill(theta3Lab * TMath::RadToDeg(), Ecm + 8.96); // T3Recon);
      hEx2->Fill(Ecm + 8.96, theta3CM * TMath::RadToDeg());
      hnorm->Fill(Ecm + 8.96, theta3CM * TMath::RadToDeg());
      hnorm_phi->Fill(Ecm + 8.96, phi3Lab * TMath::RadToDeg());
      hEx_8He->Fill(Ex_calc);
      if (Ecm + 8.96 > 0 && Ecm + 8.96 < 130) {
        // htrack->Fill(theta4Lab* TMath::RadToDeg(), heavyRange);
        hKin_6He->Fill(theta3Lab * TMath::RadToDeg(), T3EnteringSil);
      }
      if (theta3CM * TMath::RadToDeg() >= angle_min &&
          theta3CM * TMath::RadToDeg() <= angle_max) {
        hprojection->Fill(Ecm + 8.96);
        // compute maxBinContent
        if (hprojection->GetBinContent(hprojection->GetMaximumBin()) >
            maxBinContent) {
          maxBinContent =
              hprojection->GetBinContent(hprojection->GetMaximumBin());
        }
      }
      // std::cout<<"maxbincontent"<<maxBinContent<<std::endl;
      hEexAfter->Fill(EexAfter, weight);
      hSP->Fill(silPoint0InMM.Y(), silPoint0InMM.Z());
      // Fill histogram of SP with thetaCM as weight
      hSPTheta->Fill(silPoint0InMM.Y(), silPoint0InMM.Z(),
                     theta3CM * TMath::RadToDeg());
      // RP histogram
      hRP->Fill(vertex.second.X(), vertex.second.Y());
      // Beam energy at RP
      // hEBeam->Fill(theta3Lab*TMath::RadToDeg(),lightRange);
      // CM energy
      hECM->Fill(Ecm + 8.96);

      // write to TTree
      Eex_tree = EexAfter;
      weight_tree = weight;
      theta3CM_tree = theta3CM * TMath::RadToDeg();
      EVertex_tree = T3Recon;
      theta3Lab_tree = theta3Lab * TMath::RadToDeg();
      rpx_tree = vertex.second.X();
      hRPxSimu->Fill(vertex.second.X());
      outTree->Fill();
    }
    // hProj = hptoection->ProjectionX("", 50, 60);
  }

  // Obtenir le nombre de bins
  int nBinsX = hnorm->GetNbinsX();
  int nBinsY = hnorm->GetNbinsY();

  for (int i = 1; i <= nBinsX; i++) {
    for (int j = 1; j <= nBinsY; j++) {

      double n_observe = hnorm->GetBinContent(i, j);
      // std::cout << "n_data" << n_observe << "\n";

      // double n_theorique = hnorm_uniform->GetBinContent(i, j);
      double n_theorique = (iterations / 8100.) * 3.;
      // std::cout << "n_norm" << n_theorique << "\n";
      // std::cout << "\n";

      // Normalise le bin
      double normalise = 0;
      if (n_theorique > 0) {
        normalise = n_observe / n_theorique;
      }

      // Mets à jour le contenu du bin normalisé
      hnorm->SetBinContent(i, j, normalise);
    }
  }

  // // Obtenir le nombre de bins
  // int nBinsX_azure = hnorm_azure->GetNbinsX();
  // int nBinsY_azure = hnorm_azure->GetNbinsY();

  // for (int i = 1; i <= nBinsX_azure; i++) {
  //   for (int j = 1; j <= nBinsY_azure; j++) {

  //     double n_observe_azure = hnorm_azure->GetBinContent(i, j);
  //     std::cout<<"n_data"<<n_observe_azure<<"\n";

  //     //double n_theorique = hnorm_uniform->GetBinContent(i, j);
  //     double n_theorique_azure = (iterations / 8100.)*3.;
  //     std::cout<<"n_norm"<<n_theorique_azure<<"\n";
  //     std::cout<<"\n";

  //     // Normalise le bin
  //     double normalise_azure = 0;
  //     if (n_theorique_azure > 0) {
  //       normalise_azure = n_observe_azure / n_theorique_azure;
  //     }

  //     // Mets à jour le contenu du bin normalisé
  //     hnorm_azure->SetBinContent(i, j, normalise_azure);
  //   }
  // }

  // Efficiencies as quotient of histograms in TEfficiency class
  auto *eff{new TEfficiency(*hThetaCM, *hThetaCMAll)};
  eff->SetNameTitle("eff",
                    TString::Format("#theta_{CM} eff E_{x} = %.2f MeV", Ex));

  // auto hProj = GetProjectionX(hnorm.GetPtr(), 167, 176);
  // auto hprojection = hnorm->ProjectionX("", 0, 1);
  if (maxBinContent > 0.0) {
    hprojection->Scale(1.0 / maxBinContent);
    hprojection->Sumw2(kFALSE);
  }

  // SAVING
  outFile->cd();
  hnorm->Write();
  hnorm_phi->Write();
  hprojection->Write();
  outTree->Write();
  eff->Write();
  outFile->Close();
  delete outFile;
  outFile = nullptr;

  // plotting
  if (standalone) {
    // auto* c0 {new TCanvas("c0", "Canvas for inspection 0")};
    // c0->DivideSquare(4);
    // c0->cd(1);
    // hThetaCM->DrawClone();
    // c0->cd(2);
    // hThetaCMAll->DrawClone();
    // c0->cd(3);
    // hEexBefore->DrawClone("hist");
    // c0->cd(4);
    // // hThetaCMAll->DrawClone();
    // hRP->DrawClone("colz");

    // draw theoretical kinematics
    ActPhysics::Kinematics theokin{p1, p2, p3, p4, T1 * p1.GetAMU(), Ex};
    auto *gtheo{theokin.GetKinematicLine3()};

    set_plot_style();
    auto *c1{new TCanvas("cAfter", "Canvas for inspection 1")};
    c1->DivideSquare(6);
    c1->cd(1);
    hKin_6He->DrawClone("colz");
    // hKinVertex->DrawClone("colz");
    // gtheo->Draw("same");
    c1->cd(2);
    htrack->DrawClone("col");
    // c1->cd(3);
    // // hEexAfter->DrawClone("hist");
    // hSPCut->SetTitle("SP for particles not reaching sils");
    // hSPCut->DrawClone("colz");
    c1->cd(3);
    // eff->Draw("apl");
    hEx_8He->DrawClone();
    c1->cd(4);
    hSPTheta->DrawClone("colz");
    hPID->Fit("fit");
    hPID->DrawClone("colz");
    c1->cd(5);
    auto *fdummy{new TF1{"fdummy", "x", 0, 500}};
    hRPxRange->DrawClone("colz");
    fdummy->SetLineWidth(2);
    fdummy->SetLineColor(kMagenta);
    fdummy->Draw("same");
    c1->cd(6);
    // hEx2->DrawClone("col");
    hprojection->DrawClone("hist");
    std::cout << "maxbincontent" << maxBinContent << std::endl;
    // c1->cd(8);
    // hnorm_phi->DrawClone("col");

    // hRPx->DrawNormalized();
    // hRPxSimu->SetLineColor(kRed);
    // hRPxSimu->DrawNormalized("same");
  }

  // deleting news
  delete geometry;
  delete srim;
  delete rand;

  timer.Stop();
  timer.Print();
}
