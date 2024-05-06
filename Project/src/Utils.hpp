#ifndef Utils_H
#define Utils_H

#include <Eigen/Eigen>
#include <vector>
#include "Fractures.hpp"


using namespace std;
using namespace Geometry;

namespace Geometry{
bool readFractures(const string& fileName, vector<Fracture>& vec, double tol);
vector<Trace> findTraces(vector<Fracture>& fractures, double tol);
void printGlobalResults (const string& fileName, vector<Trace>& traces);
void printLocalResults (const string& fileName,const vector<Fracture>&fractures, const vector<Trace>& traces);
}

namespace Algebra{
inline Vector3d findPlaneEquation(vector<Vector3d>& points, double& constantTerm);
inline Vector3d intersectionPlaneLine(const Vector3d& coeff, const double d, const Vector3d& p1, const Vector3d& p2);
inline bool findIntersectionPoints(Fracture& f1, Fracture& f2, array<Vector3d,4>& intPoints, double tol);
}

namespace detail{
void merge(vector<unsigned int>& vecIdTraces, const vector<Trace>& traces, size_t left, size_t center, size_t right);
void mergesort(vector<unsigned int>& vecIdTraces, const vector<Trace>& traces, size_t left, size_t right);
}
void mergesort(vector<unsigned int>& data, const vector<Trace>& traces);

#endif
