#include "ActGeometry.h"

#include <map>
#include <string>
#include <utility>

void BuildGeometry_nocentralSi(bool draw = true)
{
    // Define parameters
    // Remember that we work with HALF LENGTHS
    // Drift cage for ACTAR
    double driftX {25.6 / 2}; // cm
    double driftY {25.6 / 2};
    double driftZ {25.6 / 2};
    ActSim::DriftChamber actar(driftX, driftY, driftZ);
    // unit silicon size
    double silicon1X {5.0E-2 / 2}; // cm
    double silicon1Y {8. / 2};
    double silicon1Z {5.0 / 2};
    ActSim::SilUnit silUnit(0, silicon1X, silicon1Y, silicon1Z);
    // set placements for front L0
    std::map<int, std::pair<double, double>> l0Placements {
        {0, {+2 * silicon1Y, -3 * silicon1Z}}, {1, {0, -3 * silicon1Z}},  {2, {-2 * silicon1Y, -3 * silicon1Z}},
        {3, {+2 * silicon1Y, -1 * silicon1Z}},  {5, {-2 * silicon1Y, -1 * silicon1Z}},
        {9, {+2 * silicon1Y, +3 * silicon1Z}}, {11, {-2 * silicon1Y, +3 * silicon1Z}}};
    ActSim::SilAssembly l0Assembly(0, silUnit, true, false);
    // offset from flange of ACTAR
    double l0offset {7.9}; // cm
    l0Assembly.SetOffsets(l0offset);
    l0Assembly.SetAssemblyPlacements(l0Placements);

    // BUILD GEOMETRY
    ActSim::Geometry geo {};
    geo.SetDrift(actar);
    geo.AddAssemblyData(l0Assembly);
    geo.Construct();
    geo.Print();

    // SAVE GEO
    std::string path {"./"};
    geo.WriteGeometry(path, "e837_noHe8");

    // and draw it if necessary
    if(draw)
        geo.Draw();
}
