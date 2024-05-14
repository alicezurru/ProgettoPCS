#ifndef Utils_H
#define Utils_H

#include <Eigen/Eigen>
#include <vector>
#include "Fractures.hpp"
#include "PolygonalMesh.hpp"
#include <cmath>


using namespace std;
using namespace Geometry;

namespace Geometry{
bool readFractures(const string& fileName, vector<Fracture>& vec, double tol);//T
vector<Trace> findTraces(vector<Fracture>& fractures, double tol);
void printGlobalResults (const string& fileName, vector<Trace>& traces);
void printLocalResults (const string& fileName,const vector<Fracture>&fractures, const vector<Trace>& traces);
}

namespace Algebra{
Vector3d findPlaneEquation(vector<Vector3d>& points, double& constantTerm); //T
Vector3d intersectionPlaneLine(const Vector3d& coeff, const double d, const Vector3d& p1, const Vector3d& p2); //T
bool findIntersectionPoints(Fracture& f1, Fracture& f2, array<Vector3d,4>& intPoints, double tol, array<bool,2>& onThePlane); //T
bool passBoundingBox(Fracture& f1, Fracture& f2); //T
bool findInternalPoints(array<Vector3d,4>& intPoints, double tol, array<Vector3d,2>& extremities, array<bool,2>& tips);//T
inline bool areVectorsEqual (const Vector3d& v1, const Vector3d&v2, double tol2){ //passo già la tol2 in modo da non doverla ricavare ogni volta
    if((v1-v2).squaredNorm()<= tol2)
        return true;
    return false;
}
}

namespace detail{
void merge(vector<unsigned int>& vecIdTraces, const vector<Trace>& traces, size_t left, size_t center, size_t right);//T
void mergesort(vector<unsigned int>& vecIdTraces, const vector<Trace>& traces, size_t left, size_t right);//T
void mergesort(vector<unsigned int>& data, const vector<Trace>& traces); //T
}

namespace PolygonalMeshLibrary{
vector<PolygonalMesh> cutFractures(const vector<Fracture>& fractures, const vector <Trace>& traces);
void makeCuts (const Fracture& fr, vector <Trace>& traces);
}

#endif
