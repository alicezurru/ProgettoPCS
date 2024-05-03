#include "Utils.hpp"
#include "Fractures.hpp"
#include <iostream>
#include <Eigen/Eigen>
#include <vector>
#include <algorithm>

using namespace std;
using namespace Geometry;

int main(int argc, char **argv) //passo la tolleranza
{
    double tolInput=stod(argv[1]);
    double tol=max(10*numeric_limits<double>::epsilon(), tolInput);
    vector<Fracture> vec;
    string path="./DFN";
    bool flag=readFractures("./DFN/FR3_data.txt",vec, tol);
    if (!flag){ //ci son stati problemi nella lettura file
        return 1;
    }

    return 0;
}

    /*for (Fracture fr : vec){
        cout << "numero fratture " << vec.size() << endl;
        cout << "frattura1 " << " id " << fr.idFrac << " numero vertici " << fr.numVertices;
        for (unsigned int i=0; i<fr.numVertices; i++){
            cout << " vertice " << i << ": " << fr.vertices[i][0] << " " << fr.vertices[i][1] << " " << fr.vertices[i][2] << endl;
        }
    }*/

