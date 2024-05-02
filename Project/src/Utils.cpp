#include "Utils.hpp"
#include "Fractures.hpp"

#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>


using namespace std;
using namespace Eigen;

namespace Geometry{
bool readFractures(const string& fileName, vector<Fracture>& vec){
    ifstream ifstr(fileName);
    if(ifstr.fail()){
        cerr << "errore nell'apertura del file" << endl;
        return false;
    }
    cout << "prova" << endl;
    string header; //da ignorare
    getline(ifstr, header);

    string line;
    char c; //lo uso dopo per togliere i ';'
    getline(ifstr, line);
    istringstream convert(line);
    unsigned int numFractures;
    convert >> numFractures;
    vec.resize(numFractures);

    while(getline(ifstr,line)){ //toglie già la prima riga con #
        Fracture frac;
        getline(ifstr, line);
        istringstream convert(line);
        convert >> frac.idFrac >> c >> frac.numVertices;

        getline(ifstr, line); //da ignorare
        bool firstTime; //per gestire i ';'
        for (unsigned int i=0; i<3; i++){ //3 dimensioni (i=0 componente x)
            getline(ifstr, line);
            istringstream convert(line);
            firstTime=true;
            for (unsigned int j=0; j<frac.numVertices; j++){ //numero di vertici
                if (! firstTime){
                    convert >> c; //tolgo il ';' prima di ogni componente tranne la prima
                    firstTime=false;
                }
                convert >> (frac.vertices[j])[i];
            }
        }

    }


    return true;
}
}

