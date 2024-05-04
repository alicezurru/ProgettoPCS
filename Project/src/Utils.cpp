#include "Utils.hpp"
#include "Fractures.hpp"

#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>
#include <Eigen/Eigen>
#include <functional> //per le lambda function
#include <list>


using namespace std;
using namespace Eigen;
using namespace Algebra;

namespace Geometry{
bool readFractures(const string& fileName, vector<Fracture>& vec, double tol){
    ifstream ifstr(fileName);
    if(ifstr.fail()){
        cerr << "errore nell'apertura del file" << endl;
        return false;
    }
    string header; //da ignorare
    getline(ifstr, header);

    string line;
    char c; //lo uso dopo per togliere i ';'
    getline(ifstr, line);
    istringstream convert(line);
    unsigned int numFractures;
    convert >> numFractures;
    vec.reserve(numFractures);

    while(getline(ifstr,line)){ //toglie già la prima riga con #
        Fracture frac;
        getline(ifstr, line);
        istringstream convert(line);
        convert >> frac.idFrac >> c >> frac.numVertices;
        frac.vertices.resize(frac.numVertices);

        getline(ifstr, line); //da ignorare
        bool firstTime; //per gestire i ';'
        for (unsigned int i=0; i<3; i++){ //3 dimensioni (i=0 componente x)
            getline(ifstr, line);
            istringstream convert2(line);
            firstTime=true;
            for (unsigned int j=0; j<frac.numVertices; j++){ //numero di vertici
                if (! firstTime){
                    convert2 >> c; //tolgo il ';' prima di ogni componente tranne la prima
                }
                convert2 >> ((frac.vertices)[j])[i];
                firstTime=false;
            }
        }

        for (unsigned int k=0; k<frac.numVertices-1; k++){ //controllo che non ci siano lati di lunghezza nulla.
            //In questo caso tolgo la frattura e metto nel vettore una frattura con id NULL
            Vector3d edge = (frac.vertices)[k]-(frac.vertices)[k+1]; //differenza tra due vertici consecutivi
            if (edge[0]*edge[0]+edge[1]*edge[1]+edge[2]*edge[2]<tol*tol){//distanza al quadrato in modo da non dover valutare la radice quadrata
                cerr << "la frattura " << frac.idFrac << " ha lati di lunghezza nulla" << endl;
                frac.idFrac=-1; //come se fosse null
                break;
            }
        }
        Vector3d edgeF = (frac.vertices)[0]-(frac.vertices)[frac.numVertices-1]; //faccio lo stesso per il primo e l'ultimo vertice
        if (edgeF[0]*edgeF[0]+edgeF[1]*edgeF[1]+edgeF[2]*edgeF[2]<tol*tol){
            cerr << "la frattura " << frac.idFrac << " ha lati di lunghezza nulla" << endl;
            frac.idFrac=-1;
        }
        vec.push_back(frac);
    }
    return true;
}

vector<Trace> findTraces(vector<Fracture> fractures, double tol){ //date tutte le fratture, trova tutte le tracce e le restituisce in un vettore
    list<Trace> listTraces; //metto prima in una lista per efficienza nell'aggiungere le nuove tracce, poi in un vettore per efficienza nell'accesso casuale
    vector<Trace> vectorTraces;
    for (unsigned int i=0; i<fractures.size(); i++){//devo controllare ogni coppia di fratture possibile
        for(unsigned int j=i+1; j<fractures.size(); j++){
                //controllo se i due poligoni sono molto lontani (bounding box) e in quel caso passo alla coppia successiva
                //prima cerco centro (facendo la media dei vertici) e raggio (massima distanza tra il centro e i vertici) delle bounding box
                Vector3d centreBB1 = Vector3d::Zero();
                for(unsigned int k1=0; k1<(fractures[i]).numVertices; k1++){
                    centreBB1 = centreBB1 + (fractures[i]).vertices[k1];
                }
                centreBB1= centreBB1/(fractures[i]).numVertices;
                double radiusBB1 = 0.0;
                for(unsigned int k1=0; k1<(fractures[i]).numVertices; k1++){
                    radiusBB1=max(radiusBB1,((fractures[i]).vertices[k1]-centreBB1).squaredNorm());
                }

                Vector3d centreBB2 = Vector3d::Zero();
                for(unsigned int k2=0; k2<(fractures[j]).numVertices; k2++){
                    centreBB2 = centreBB2 + (fractures[j]).vertices[k2];
                }
                double radiusBB2 = 0.0;
                centreBB2= centreBB2/(fractures[j]).numVertices;

                for(unsigned int k2=0; k2<(fractures[j]).numVertices; k2++){
                    radiusBB2=max(radiusBB2,((fractures[j]).vertices[k2]-centreBB2).squaredNorm());
                }
                if((radiusBB1+radiusBB2)*(radiusBB1+radiusBB2) > (centreBB1-centreBB2).squaredNorm()){ //vado avanti solo se il controllo della
//                    bounding box non è passato (se le due probabilmente si intersecano)
                    //ora vedo se c'è effettivamente intersezione
                    array<Vector3d,4> intPoints;//qui metterò i potenziali punti di intersezione
                    bool intersection = findIntersectionPoints(fractures[i],fractures[j],intPoints,tol);
                    if(intersection){
                        //ora stabiliamo tra i 4 potenziali chi sono i punti di intersezione
                        //QUAAAA

                    }



                }
            }
        }
    }
    //convertire la lista di tracce in vettore
    //return
}


namespace Algebra{
//per trovare l'equazione del piano che contiene i vertici di un poligono

inline Vector3d findPlaneEquation(const vector<Vector3d>& points, double& constantTerm){ //restituisce la normale e modifica il dato in input che corrisponde al termine noto
    //assumiamo che non ci possano essere 3 punti allineati e che le fratture siano planari
    //calcolo della normale:
    Vector3d v1=points[1]-points[0];
    Vector3d v2= points[2]-points[0];
    Vector3d n= v1.cross(v2); //le componenti della normale definiscono i primi 3 coefficienti del piano
    constantTerm=-(n[0]*(points[0])[0]+n[1]*(points[0])[1]+n[2]*(points[0])[2]);

    return n;
}


inline Vector3d intersectionPlaneLine(const Vector3d& coeff, const double d, const Vector3d& p1, const Vector3d& p2 ){
    // eq retta: s=p1+t(p2-p1)
    // eq piano: ax+by+cz+d=0
    double t=-(p1.dot(coeff)+d)/(coeff.dot(p2-p1));
    Vector3d inter = p1+t*(p2-p1);
    return inter;
}

inline bool findIntersectionPoints(Fracture& f1, Fracture& f2, array<Vector3d,4>& intPoints, double tol){
    bool intersection=false;
    //controllo se i piani che contengono le due fratture sono parallelli (non possono intersecarsi)
    double d1; //termine noto piano 1
    double d2;
    Vector3d coeff1 = findPlaneEquation(f1.vertices, d1); //coefficienti del piano che contengono la frattura 1
    Vector3d coeff2 = findPlaneEquation(f2.vertices, d2);
    if((coeff1.cross(coeff2)).squaredNorm()>tol*tol){ //i piani sono paralleli se hanno normali parallele
        //posizione dei punti del poligono 2 rispetto al piano 1
        bool positive=false; //per segnalare se il vertice analizzato in questo momento sta "sopra" (true) o "sotto"(false) il piano
        bool previous=false;
        if ((f2.vertices[0]).dot(coeff1)+d1>0){
            previous=true;//vedo se si comincia "sopra" o "sotto"
        }
        unsigned int count1=0; //conto quanti vertici stanno dall'altra parte
        unsigned int firstVertexOtherSide2=0;
        for (unsigned int i=1; i<f2.numVertices; i++){
            if ((f2.vertices[i]).dot(coeff1)+d1>0){ //vedo da che parte stanno i vertici
                positive=true;
            }
            else {
                positive=false;
            }
            if((previous!=positive) && (count1!=0)){
                break; //non è necessario andare avanti: saranno tutti di nuovo dalla prima parte
            }
            if ((previous!=positive) && (count1==0)){ //si incontra un vertice dall'altra parte per la prima volta
                firstVertexOtherSide2=i;
                count1++;
            }
            if((previous=positive) && (count1!=0)){ //conto quanti stanno dall'altra parte
                count1++;
            }

            previous=positive;
        }
        if(count1!=0){ //c'è almeno un vertice dall'altra parte
            intersection=true;
        }
        //ora posizione dei punti del poligono 1 rispetto al piano 2
        if(intersection){ //solo se già risulta accaduto per il piano 1
            previous=false;
            if ((f1.vertices[0]).dot(coeff2)+d1>0){
                previous=true;//vedo se si comincia "sopra" o "sotto"
            }
            unsigned int count2=0;
            unsigned int firstVertexOtherSide1=0;
            for (unsigned int i=1; i<f2.numVertices; i++){
                if ((f2.vertices[i]).dot(coeff1)+d1>0){ //vedo da che parte stanno i vertici
                    positive=true;
                }
                else {
                    positive=false;
                }
                if((previous!=positive) && (count2!=0)){
                    break; //non è necessario andare avanti: saranno tutti di nuovo dalla prima parte
                }
                if ((previous!=positive) && (count2==0)){ //si incontra un vertice dall'altra parte per la prima volta
                    firstVertexOtherSide1=i;
                    count2++;
                }
                if((previous=positive) && (count2!=0)){ //conto quanti stanno dall'altra parte
                    count2++;
                }

                previous=positive;
            }
            if(count2!=0){ //c'è intersezione
                //ora individuo i punti di intersezione
                intPoints[0]=intersectionPlaneLine(coeff2, d2,f1.vertices[firstVertexOtherSide1-1],f1.vertices[firstVertexOtherSide1]);
                intPoints[1]=intersectionPlaneLine(coeff2, d2,f1.vertices[firstVertexOtherSide1+count2-1],f1.vertices[firstVertexOtherSide1+count2]);
                //i primi due punti sono del poligono 1, i successivi 2 del poligono 2
                intPoints[2]=intersectionPlaneLine(coeff1, d1,f2.vertices[firstVertexOtherSide2-1],f2.vertices[firstVertexOtherSide2]);
                intPoints[3]=intersectionPlaneLine(coeff1, d1,f2.vertices[firstVertexOtherSide2+count1-1],f2.vertices[firstVertexOtherSide2+count1]);
            }
            else{
                intersection=false;
            }
        }

        }
        return intersection;
}

}



