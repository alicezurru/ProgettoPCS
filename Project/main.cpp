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

    vector<Vector3d> points = {Vector3d(1.0, 2.5, -0.5), Vector3d(2.0, -1.0, 0.0), Vector3d(0.0, 2.0, 1.0)};
    //non importa quanti punti gli do, userà sempre solo i primi 3, quindi posso dargli direttamente 3 punti
    double constantTerm=0;
    Vector3d n = Algebra::findPlaneEquation(points, constantTerm);

    double tolInput=stod(argv[1]);
    double tol=max(10*numeric_limits<double>::epsilon(), tolInput);
    vector<Fracture> vec;
    string path="./DFN";
    bool flag=readFractures(path+"/FR3_data.txt",vec, tol);
    if (!flag){ //ci son stati problemi nella lettura file
        return 1;
    }
    vector<Trace> vecTraces=findTraces(vec,tol);
    //printGlobalResults("results", vecTraces);
    //printLocalResults("lresults",vec,vecTraces);


    //prova findTraces
    /*vector<Trace> vecTraces=findTraces(vec,tol);
    cout << "Numero tracce: " << vecTraces.size() << endl;
    for (Trace& tr : vecTraces){
        cout << "Traccia: " << " id " << tr.idTr << ", Frattura1 "<<tr.fracturesIds[0]<<", Frattura2 "<<tr.fracturesIds[1]<<endl;
        cout<<tr.Tips[0]<<tr.Tips[1]<<endl;

    }*/


    //prova mergesort
    /*vector<Trace> vecProvaMS;
    Trace t1;
    t1.idTr=0;
    t1.length=3.0;
    Trace t2;
    t2.idTr=1;
    t2.length=5.0;
    Trace t3;
    t3.idTr=2;
    t3.length=2.0;
    Trace t4;
    t3.idTr=3;
    t3.length=1.0;
    Trace t5;
    t3.idTr=4;
    t3.length=1.5;
    vector<unsigned int> idProvaMS ={0,1,2,3,4};
    vecProvaMS.push_back(t1);
    vecProvaMS.push_back(t2);
    vecProvaMS.push_back(t3);
    vecProvaMS.push_back(t4);
    vecProvaMS.push_back(t5);*/
    /*detail::mergesort(idProvaMS, vecProvaMS, 0, 4);
    for(unsigned int i=0; i<5; i++)
        cout << idProvaMS[i] << endl;*/

    /*vector<Vector3d> points = {Vector3d(1.0, 2.5, -0.5), Vector3d(2.0, -1.0, 0.0), Vector3d(0.0, 2.0, 1.0)};
    double constantTerm=0;
    Vector3d n = Algebra::findPlaneEquation(points, constantTerm);*/


    return 0;
}

    /*for (Fracture fr : vec){
        cout << "numero fratture " << vec.size() << endl;
        cout << "frattura1 " << " id " << fr.idFrac << " numero vertici " << fr.numVertices;
        for (unsigned int i=0; i<fr.numVertices; i++){
            cout << " vertice " << i << ": " << fr.vertices[i][0] << " " << fr.vertices[i][1] << " " << fr.vertices[i][2] << endl;
        }
    }*/

