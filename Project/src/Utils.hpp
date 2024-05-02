#ifndef Utils_H
#define Utils_H

#include <Eigen/Eigen>
#include <vector>
#include "Fractures.hpp"

using namespace std;
namespace Geometry{
bool readFractures(const string& fileName, vector<Fracture>& vec);
}

#endif
