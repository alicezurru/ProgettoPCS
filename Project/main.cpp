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

int main()
{

    cout<<"Inserire la tolleranza: "<<endl;
    double tolInput;
    cin>>tolInput; //passo la tolleranza
    double tol=max(10*numeric_limits<double>::epsilon(), tolInput);
    chrono::steady_clock::time_point t_begin = std::chrono::steady_clock::now();
    vector<Fracture> vec;
    string path="./DFN";
    bool flag=readFractures(path+"/FR10_data.txt",vec, tol);
    if (!flag){ //ci son stati problemi nella lettura file
        return 1;
    }
    vector<Trace> vecTraces=findTraces(vec,tol);
    printGlobalResults("results", vecTraces);
    printLocalResults("lresults",vec,vecTraces);

    //parte 2
    vector<PolygonalMesh> vecMesh = cutFractures(vec, vecTraces,tol);
    printPolygonalMesh(vecMesh, "printMesh");

    //paraview
    Export::exportMesh(vec[0],vecMesh);


    chrono::steady_clock::time_point t_end = std::chrono::steady_clock::now();
    double duration=std::chrono::duration_cast<std::chrono::milliseconds>(t_end-t_begin).count();
    cout<<"Tempo impiegato: "<<duration<<endl;


    return 0;
}



