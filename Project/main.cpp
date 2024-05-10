#include "Utils.hpp"
#include "Fractures.hpp"
#include <iostream>
#include <Eigen/Eigen>
#include <vector>
#include <algorithm>
#include <chrono>

using namespace std;
using namespace Geometry;

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


    //prova findTraces
    /*vector<Trace> vecTraces=findTraces(vec,tol);
    cout << "Numero tracce: " << vecTraces.size() << endl;
    for (Trace& tr : vecTraces){
        cout << "Traccia: " << " id " << tr.idTr << ", Frattura1 "<<tr.fracturesIds[0]<<", Frattura2 "<<tr.fracturesIds[1]<<endl;
        cout<<tr.Tips[0]<<tr.Tips[1]<<endl;

    }*/


    chrono::steady_clock::time_point t_end = std::chrono::steady_clock::now();
    double duration=std::chrono::duration_cast<std::chrono::milliseconds>(t_end-t_begin).count();
    cout<<"Tempo impiegato: "<<duration<<endl;

    vector<Trace> vecProvaMS(5);
    Trace t1;
    t1.idTr=0;
    t1.length=5.0;
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
    vector<unsigned int> idProvaMS ={0,1,3,4};
    //ordine:1403
    vecProvaMS.push_back(t1);
    vecProvaMS.push_back(t2);
    vecProvaMS.push_back(t3);
    vecProvaMS.push_back(t4);
    vecProvaMS.push_back(t5);
    detail::mergesort(idProvaMS, vecProvaMS,0,3);
    cout<<idProvaMS[0];
    cout<<idProvaMS[1];
    cout<<idProvaMS[2];
    cout<<idProvaMS[3];

    return 0;
}

    /*for (Fracture fr : vec){
        cout << "numero fratture " << vec.size() << endl;
        cout << "frattura1 " << " id " << fr.idFrac << " numero vertici " << fr.numVertices;
        for (unsigned int i=0; i<fr.numVertices; i++){
            cout << " vertice " << i << ": " << fr.vertices[i][0] << " " << fr.vertices[i][1] << " " << fr.vertices[i][2] << endl;
        }
    }*/

