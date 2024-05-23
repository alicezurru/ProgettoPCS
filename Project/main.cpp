#include "Utils.hpp"
#include "Fractures.hpp"
#include <iostream>
#include <Eigen/Eigen>
#include <vector>
#include <algorithm>
#include <chrono>

using namespace std;
using namespace Geometry;
using namespace PolygonalMeshLibrary;

int main(int argc, char **argv) //passo la tolleranza
{
    chrono::steady_clock::time_point t_begin = std::chrono::steady_clock::now();
    double tolInput=stod(argv[1]);
    double tol=max(10*numeric_limits<double>::epsilon(), tolInput);
    vector<Fracture> vec;
    string path="./DFN";
    bool flag=readFractures(path+"/FR3_data.txt",vec, tol);
    if (!flag){ //ci son stati problemi nella lettura file
        return 1;
    }
    vector<Trace> vecTraces=findTraces(vec,tol);
    printGlobalResults("results", vecTraces);
    printLocalResults("lresults",vec,vecTraces);

    //parte 2
    vector<PolygonalMesh> vecMesh = cutFractures(vec, vecTraces,tol);
    for (unsigned int i=0; i<vecMesh.size(); i++){
        printPolygonalMesh(vecMesh[i], "printMesh"+to_string(i));
    }

    chrono::steady_clock::time_point t_end = std::chrono::steady_clock::now();
    double duration=std::chrono::duration_cast<std::chrono::milliseconds>(t_end-t_begin).count();
    cout<<"Tempo impiegato: "<<duration<<endl;



    return 0;
}



