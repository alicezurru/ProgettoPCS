#ifndef Utils_H
#define Utils_H

#include <Eigen/Eigen>
#include <vector>
#include "Fractures.hpp"


using namespace std;
using namespace Geometry;

namespace Geometry{
bool readFractures(const string& fileName, vector<Fracture>& vec, double tol);
vector<Trace> findTraces(vector<Fracture> fractures, double tol);
}

namespace Algebra{
inline Vector3d findPlaneEquation(const vector<Vector3d>& points, double& constantTerm);
inline Vector3d intersectionPlaneLine(const Vector3d& coeff, const double d, const Vector3d& p1, const Vector3d& p2);
inline bool findIntersectionPoints(Fracture& f1, Fracture& f2, array<Vector3d,4>& intPoints, double tol);
}

#endif
