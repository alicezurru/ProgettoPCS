#ifndef Fractures_H
#define Fractures_H

#include <Eigen/Eigen>
#include <vector>
#include <array>

using namespace std;
using namespace Eigen;

namespace Geometry{
struct Fracture{
    unsigned int idFrac;
    unsigned int numVertices;
    vector <Vector3d> vertices; //vettore con le coordinate dei vertici della frattura
    vector <unsigned int> passingTraces; //vettore con gli id delle tracce passanti per la frattura corrente
    vector <unsigned int> notPassingTraces;

};

struct Trace{
    unsigned int idTr;
    array <Vector3d, 2> extremitiesCoord; //array delle coordinate dei punti estremi della traccia
    array <unsigned int, 2> fracturesIds; //array degli id delle due fratture che formano la traccia
    array<bool,2> Tips; //memorizza se la traccia è passante (F) o no (V) in ognuna delle due fratture coinvolte
    double length;
};

}

#endif
