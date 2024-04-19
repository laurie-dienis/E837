#include "TString.h"

#include <stdexcept>
#include <string>
#include <vector>

// #include "./Plotter.cpp"
#include "./Simulation_E837.cxx"
// what
// if simu = runs simulation
// if plot = plots results
// standalone (only applies to simu setting):
// if true, runs only first item in Ex vector and plots in-simulation results
// if false, runs all Ex simulations but doesn't plot
void Runner(TString what = "plot", bool standalone = true)
{
    // Settings
    // Names of particles
    std::string beam {"8He"};
    std::string target {"4He"};
    std::string light {"6He"};
    std::string heavy {"6He"};
    // Phase space reactions: when the heavy decays by proton or neutron emission
    // So we have something like: 4He + n + 17N (needs to be simulated to be
    // included as background in fits)
    int neutronPS {0}; // number of neutrons in final state
    int protonPS {0};  // number of protons in final state
    double T1 {1.51};  // Beam energy at entrance of pad plane

    std::vector<double> Eexs;
    if(neutronPS == 0 && protonPS == 0)
        Eexs = {1.79};
    else if(neutronPS > 0 && protonPS == 0)
        Eexs = {0}; // only gs for n phase space
    else if(neutronPS == 0 && protonPS > 0)
        Eexs = {0};
    else
        throw std::runtime_error("No confs with neutronPS and protonPS enabled at the same time");

    if(what.Contains("simu"))
    {
        for(const auto& Eex : Eexs)
        {
            Simulation_E837(beam, target, light, heavy, neutronPS, protonPS, T1, Eex, standalone);
            if(standalone)
                break;
        }
    }
    if(what.Contains("plot"))
    {
        throw std::runtime_error("Runner for E837 hasn't implemented yet a Plotter function");
        // Plotter(Eexs, beam, target, light, heavy, T1, neutronPS, protonPS);
    }
}
