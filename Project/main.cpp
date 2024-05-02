#include "Utils.hpp"
#include "Fractures.hpp"
#include <iostream>
#include <Eigen/Eigen>
#include <vector>

using namespace std;
using namespace Geometry;

int main()
{
    vector<Fracture> vec;
    string path="./DFN";
    bool flag=readFractures("C:/Users/sophi/Desktop/ProgettoPCS/Project/DFN/FR3_data.txt",vec);
    //"C:\Users\sophi\Desktop\ProgettoPCS\Project\DFN\FR10_data.txt"

    for (Fracture fr : vec){
        cout << "numero fratture " << vec.size() << endl;
        cout << "frattura1 " << " id " << fr.idFrac << " numero vertici " << fr.numVertices;
        for (unsigned int i=0; i<fr.numVertices; i++){
            cout << " vertice " << i << ": " << fr.vertices[i][0] << " " << fr.vertices[i][1] << " " << fr.vertices[i][2] << endl;
        }
    }

    return 0;
}
