#ifndef Utils_H
#define Utils_H

#include <Eigen/Eigen>
#include <vector>
#include "Fractures.hpp"

using namespace std;
namespace Geometry{
bool readFractures(const string& fileName, vector<Fracture>& vec, double tol);
vector<Trace> findTraces(vector<Fracture> fractures, double tol);
}

namespace Algebra{
Vector3d findPlaneEquation(const vector<Vector3d>& points, double& constantTerm);
}

#endif
