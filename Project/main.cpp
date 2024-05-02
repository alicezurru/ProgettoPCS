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
    bool flag=readFractures("C:/Users/sophi/Desktop/ProgettoPCS/Project/DFN/FR10_data.txt",vec);
    //"C:\Users\sophi\Desktop\ProgettoPCS\Project\DFN\FR10_data.txt"

    return 0;
}
