#ifndef Profile_h
#define Profile_h

// Header containing utility functions to deal with profiles

#include "TGraph.h"
#include "TH1.h"

#include "Math/Point3D.h"
#include "Math/Vector3D.h"

using XYZPoint = ROOT::Math::XYZPointF;
using XYZVector = ROOT::Math::XYZVectorF;

namespace ChargeProfile
{

class Data
{
public:
    XYZPoint fVertexP {-1, -1, -1};
    XYZPoint fStopP {-1, -1, -1};
    float fQVertex {-1};
};

// Functions
TH1D* GetDerivative(TH1* h);
std::pair<float, float> GetVertexAndQmax(TH1* der, TH1* prof, double width = 10);
float GetStoppingPoint(TH1* prof);
} // namespace ChargeProfile

#endif
