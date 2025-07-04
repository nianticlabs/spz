#include "load-spz.h"
#include <cassert>
#include <iostream>

using namespace spz;

Quat4f getQuaternion(const GaussianCloud& gaussians, size_t i)
{
  return normalized(
    Quat4f{gaussians.rotations[4*i], gaussians.rotations[4*i +1], gaussians.rotations[4*i + 2], gaussians.rotations[4*i + 3]}
  );
}

void computeAngularError(const GaussianCloud& a, const GaussianCloud& b, float& maxError, float& avgError)
{
  assert(a.numPoints == b.numPoints);
  maxError = avgError = 0;
  for (size_t i = 0; i < a.numPoints; i++)
  {    
    const float error = quaternionGeodesicDistanceDegree(getQuaternion(a, i), getQuaternion(b, i));
    if (error > maxError)
    {
      maxError = error;
    }
    avgError += error;
  }
  avgError /= a.numPoints;
}

int main()
{
    GaussianCloud before = loadSplatFromPly("D:\\tests\\spz\\before.ply", UnpackOptions());

    saveSpz(before, PackOptions(), "D:\\tests\\spz\\compressed.spz");

    GaussianCloud after = loadSpz("D:\\tests\\spz\\compressed.spz", UnpackOptions());
    saveSplatToPly(after, PackOptions(), "D:\\tests\\spz\\after.ply");

    float maxError, avgError;
    computeAngularError(before, after, maxError, avgError);
    std::cout << "Max error in deg: " << maxError << std::endl;
    std::cout << "Avg error in deg: " << avgError << std::endl;

    return 0;
}
