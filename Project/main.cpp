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
    vector<Fracture> vec;
    string path="./DFN";
    bool flag=readFractures(path+"/FR362_data.txt",vec, tol);
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


    return 0;
}



