
#include "Fractures.hpp"
#include "Utils.hpp"
#include "PolygonalMesh.hpp"

#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>
#include <Eigen/Eigen>
#include <functional> //per le lambda function
#include <list>
#include <cassert> //per mergesort


using namespace std;
using namespace Eigen;
using namespace Algebra;
using namespace detail;


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
            //In questo caso tolgo la frattura e metto nel vettore una frattura con id -1
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
    ifstr.close();
    return true;
}

vector<Trace> findTraces(vector<Fracture>& fractures, double tol){ //date tutte le fratture, trova tutte le tracce e le restituisce in un vettore
    list<Trace> listTraces; //metto prima in una lista per efficienza nell'aggiungere le nuove tracce, poi in un vettore per efficienza nell'accesso casuale
    vector<Trace> vectorTraces;
    for (unsigned int i=0; i<fractures.size(); i++){//devo controllare ogni coppia di fratture possibile
        for(unsigned int j=i+1; j<fractures.size(); j++){
                //controllo se i due poligoni sono molto lontani (bounding box) e in quel caso passo alla coppia successiva
                //prima cerco centro (facendo la media dei vertici) e raggio (massima distanza tra il centro e i vertici) delle bounding box
                bool pass = passBoundingBox(fractures[i],fractures[j]);
                if(pass){ //vado avanti solo se il controllo della
                //bounding box non è passato (se le due probabilmente si intersecano)
                    array<Vector3d,4> intPoints;//qui metterò i potenziali punti di intersezione.
                    //Li tengo comunque memorizzati perché mi serviranno nella parte 2
                    array<bool,2> onThePlane;
                    bool passIntPlane = findIntersectionPoints(fractures[i],fractures[j],intPoints,tol, onThePlane);
                    if(passIntPlane){//ora vedo se c'è effettivamente intersezione
                        //e stabiliamo tra i 4 potenziali chi sono i punti di intersezione
                        array<Vector3d,2> extremities; //qui salverò i due punti estremi della traccia
                        array<bool,2> tips = {true,true}; //se resta così è non passante per entrambi
                        bool intersection = findInternalPoints(intPoints,tol,extremities,tips);
                        //modifico tips in caso di onThePlane
                        if(onThePlane[0]){
                            tips[0]=false;
                        }
                        if(onThePlane[1]){
                            tips[1]=false;
                        }

                        //calcolo la lunghezza ed escludo il caso in cui sia un unico punto a toccare il poligono
                        double len =(extremities[0]-extremities[1]).norm();
                        if (len>tol&&intersection){
                            //creo finalmente la traccia
                            Trace tr;
                            tr.idTr=listTraces.size();
                            tr.extremitiesCoord = extremities;
                            tr.fracturesIds = {fractures[i].idFrac,fractures[j].idFrac};
                            tr.length=len;
                            tr.Tips=tips;
                            tr.onThePlane=onThePlane;
                            tr.intPoints=intPoints; //per la parte 2
                            listTraces.push_back(tr);//inserisco la traccia ora creata nella lista
                            if(tips[0]){ //e nel vettore corrispondente nella frattura
                                fractures[i].notPassingTraces.push_back(tr.idTr);
                            }
                            else{
                                fractures[i].passingTraces.push_back(tr.idTr);

                            }
                            if(tips[1]){
                                fractures[j].notPassingTraces.push_back(tr.idTr);

                            }
                            else{
                                fractures[j].passingTraces.push_back(tr.idTr);

                            }//vedi poi se c'è modo migliore per farlo
                        }
                    }
                }
            }
        }
    //ora ho una lista di tracce: popolo il vettore
        vectorTraces.reserve(listTraces.size());
        for(Trace& tr:listTraces){
            vectorTraces.push_back(tr);
        }

        return vectorTraces;
    }

void printGlobalResults (const string& fileName, vector<Trace>& traces){ //primo file di ouput, con le informazioni sulle tracce
    ofstream ofstr(fileName); //se il file non esiste, lo crea
    ofstr << "# Number of Traces" << endl;
    ofstr << traces.size() << endl;
    ofstr << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2" << endl;
    for (Trace& tr:traces){
        ofstr << tr.idTr << "; " << tr.fracturesIds[0] << "; " << tr.fracturesIds[1] << "; " << (tr.extremitiesCoord[0])[0] << "; "
              << (tr.extremitiesCoord[0])[1] << "; " << (tr.extremitiesCoord[0])[2] << "; " << (tr.extremitiesCoord[1])[0] << "; "
              << (tr.extremitiesCoord[1])[1] << "; " << (tr.extremitiesCoord[1])[2] << endl;
    }
    ofstr.close();
}

void printLocalResults (const string& fileName, vector<Fracture>& fractures, const vector<Trace>& traces){
    //secondo file di ouput, con le informazioni sulle fratture e sulle tracce corrispondenti
    ofstream ofstr(fileName);
    bool firstTime;
    for (Fracture& fr:fractures){
        if (fr.idFrac!=-1){ //non gestisco quelle eventualmente problematiche
            firstTime=true;
            ofstr << "# FractureId; NumTraces" << endl;
            ofstr << fr.idFrac << "; " << (fr.passingTraces.size())+(fr.notPassingTraces.size()) << endl;
            mergesort(fr.passingTraces, traces); //ordino le tracce passanti per lunghezza decrescente (chiamo quella fuori dal namespace details così da ordinarlo tutto)
            mergesort(fr.notPassingTraces, traces); //ordino le tracce non passanti per lunghezza decrescente
            for (unsigned int trId:fr.passingTraces){
                if (firstTime){
                     ofstr << "# TraceId; Tips; Length" << endl;
                    firstTime=false;
                }
                ofstr << trId << "; " << "false; " << traces[trId].length << endl; //sto stampando prima tutte quelle passanti, quindi avranno tutte tips=false

            }
            for (unsigned int trId:fr.notPassingTraces){
                if (firstTime){
                    ofstr << "# TraceId; Tips; Length" << endl;
                    firstTime=false;
                }
                ofstr << trId << "; " << "true; " << traces[trId].length << endl; //sto stampando prima tutte quelle non passanti, quindi avranno tutte tips=true
            }
        }
    }
    ofstr.close();
}
}


namespace Algebra{

//per trovare l'equazione del piano che contiene i vertici di un poligono
Vector3d findPlaneEquation(vector<Vector3d>& points, double& constantTerm){ //restituisce la normale e modifica il dato in input che corrisponde al termine noto
    //assumiamo che non ci possano essere 3 punti allineati e che le fratture siano planari
    //calcolo della normale:
    Vector3d v1=points[1]-points[0];
    Vector3d v2= points[2]-points[0];
    Vector3d n= v1.cross(v2); //le componenti della normale definiscono i primi 3 coefficienti del piano
    constantTerm=-(n[0]*(points[0])[0]+n[1]*(points[0])[1]+n[2]*(points[0])[2]);

    return n;
}


Vector3d intersectionPlaneLine(const Vector3d& coeff, const double d, const Vector3d& p1, const Vector3d& p2 ){
    // eq retta: s=p1+t(p2-p1)
    // eq piano: ax+by+cz+d=0
    double t=-(p1.dot(coeff)+d)/(coeff.dot(p2-p1));
    Vector3d inter = p1+t*(p2-p1);
    return inter;
}

bool passBoundingBox(Fracture& f1, Fracture& f2){
    bool pass = false;
    //controllo se i due poligoni sono molto lontani (bounding box)
    //prima cerco centro (facendo la media dei vertici) e raggio (massima distanza tra il centro e i vertici) delle bounding box

    Vector3d centreBB1 = Vector3d::Zero();
    for(unsigned int k1=0; k1<(f1).numVertices; k1++){
        centreBB1 = centreBB1 + (f1).vertices[k1];
    }
    centreBB1= centreBB1/(f1).numVertices;
    double radiusBB1 = 0.0;
    for(unsigned int k1=0; k1<(f1).numVertices; k1++){
        radiusBB1=max(radiusBB1,((f1).vertices[k1]-centreBB1).norm());
    }

    Vector3d centreBB2 = Vector3d::Zero();
    for(unsigned int k2=0; k2<(f2).numVertices; k2++){
        centreBB2 = centreBB2 + (f2).vertices[k2];
    }
    double radiusBB2 = 0.0;
    centreBB2= centreBB2/(f2).numVertices;

    for(unsigned int k2=0; k2<(f2).numVertices; k2++){
        radiusBB2=max(radiusBB2,((f2).vertices[k2]-centreBB2).norm());
    }


    if((radiusBB1+radiusBB2)*(radiusBB1+radiusBB2) > (centreBB1-centreBB2).squaredNorm()){
        pass=true;
    }
    return pass;
}

bool findIntersectionPoints(Fracture& f1, Fracture& f2, array<Vector3d,4>& intPoints, double tol, array<bool,2>& onThePlane){
    bool passIntPoint=false; //CAMBIA: è tipo probabile intersezione, caso Andrea
    //controllo se i piani che contengono le due fratture sono parallelli (non possono intersecarsi)
    double d1; //termine noto piano 1
    double d2;

    Vector3d coeff1 = findPlaneEquation(f1.vertices, d1); //coefficienti del piano che contengono la frattura 1
    Vector3d coeff2 = findPlaneEquation(f2.vertices, d2);
    double tol2=max(10*numeric_limits<double>::epsilon(), tol*tol);
    if((coeff1.cross(coeff2)).squaredNorm()>tol2){ //i piani sono paralleli se hanno normali parallele
        //posizione dei punti del poligono 2 rispetto al piano 1
        bool positive=false; //per segnalare se il vertice analizzato in questo momento sta "sopra" (true) o "sotto"(false) il piano
        bool previous=false;
        if ((f2.vertices[0]).dot(coeff1)+d1>tol){
            previous=true;//vedo se si comincia "sopra" o "sotto"
        }
        unsigned int count1=0; //conto quanti vertici stanno dall'altra parte
        unsigned int firstVertexOtherSide2=0;
        onThePlane[1] = false; //metto anche il caso in cui si toccano senza che uno passi attraverso l'altro
        for (unsigned int i=1; i<f2.numVertices; i++){
            if(abs((f2.vertices[i]).dot(coeff1)+d1)<tol){
                //se ne ho due di fila sul piano, ho traccia particolare in cui uno non passa attraverso l'altro
                if(abs((f2.vertices[(i+1)%f2.numVertices]).dot(coeff1)+d1)<tol){
                    onThePlane[1]=true;
                    intPoints[2]=f2.vertices[i];
                    intPoints[3]=f2.vertices[(i+1)%f2.numVertices];
                    break;
                }
                else if(abs((f2.vertices[(i-1)%f2.numVertices]).dot(coeff1)+d1)<tol){
                    onThePlane[1]=true;
                    intPoints[2]=f2.vertices[(i-1)%f2.numVertices];
                    intPoints[3]=f2.vertices[i];
                    break;
                }

                }
            if ((f2.vertices[i]).dot(coeff1)+d1>tol){ //vedo da che parte stanno i vertici
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
            if((previous==positive) && (count1!=0)){ //conto quanti stanno dall'altra parte
                count1++;
            }

            previous=positive;

        }
        if(count1!=0||onThePlane[1]){ //c'è almeno un vertice dall'altra parte
            passIntPoint=true;
        }
        //ora posizione dei punti del poligono 1 rispetto al piano 2
        if(passIntPoint){ //solo se già risulta accaduto per il piano 1
            previous=false;
            onThePlane[0] = false; //metto anche il caso in cui si toccano senza che uno passi attraverso l'altro
            if ((f1.vertices[0]).dot(coeff2)+d2>tol){
                previous=true;//vedo se si comincia "sopra" o "sotto"
            }
            unsigned int count2=0;
            unsigned int firstVertexOtherSide1=0;
            for (unsigned int i=1; i<f1.numVertices; i++){
                if(abs((f1.vertices[i]).dot(coeff2)+d2)<tol){
                    //se ne ho due di fila sul piano, ho traccia particolare in cui uno non passa attraverso l'altro
                    if(abs((f1.vertices[(i+1)%f1.numVertices]).dot(coeff2)+d2)<tol){
                        onThePlane[0]=true;
                        intPoints[0]=f1.vertices[i];
                        intPoints[1]=f1.vertices[(i+1)%f1.numVertices];
                        break;
                    }
                    else if(abs((f1.vertices[(i-1)%f1.numVertices]).dot(coeff2)+d2)<tol){
                        onThePlane[0]=true;
                        intPoints[0]=f1.vertices[(i-1)%f1.numVertices];
                        intPoints[1]=f1.vertices[i];
                        break;
                    }

                }
                if ((f1.vertices[i]).dot(coeff2)+d2>tol){ //vedo da che parte stanno i vertici
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
                if((previous==positive) && (count2!=0)){ //conto quanti stanno dall'altra parte
                    count2++;
                }

                previous=positive;
            }

            if(count2!=0||onThePlane[0]){ //c'è intersezione
                //ora individuo i punti di intersezione
                if(!onThePlane[0]){//se non so gia che i due punti sono sul piano
                intPoints[0]=intersectionPlaneLine(coeff2, d2,f1.vertices[(firstVertexOtherSide1-1)],f1.vertices[firstVertexOtherSide1]);
                intPoints[1]=intersectionPlaneLine(coeff2, d2,f1.vertices[(firstVertexOtherSide1+count2-1)%f1.numVertices],f1.vertices[(firstVertexOtherSide1+count2)%f1.numVertices]);
                }
                //i primi due punti sono del poligono 1, i successivi 2 del poligono 2
                if(!onThePlane[1]){
                intPoints[2]=intersectionPlaneLine(coeff1, d1,f2.vertices[(firstVertexOtherSide2-1)],f2.vertices[firstVertexOtherSide2]);
                intPoints[3]=intersectionPlaneLine(coeff1, d1,f2.vertices[(firstVertexOtherSide2+count1-1)%f2.numVertices],f2.vertices[(firstVertexOtherSide2+count1)%f2.numVertices]);
                }
            }
            else{
                passIntPoint=false;
            }
        }

        }
        return passIntPoint;
}

bool findInternalPoints(array<Vector3d,4>& intPoints, double tol, array<Vector3d,2>& extremities, array<bool,2>& tips)
{//trovo se effettivamente si interseca e i punti che delimitano la traccia
    //vedo la posizione reciproca dei punti per stabilire i due più interni: saranno gli estremi della traccia
    //inoltre così stabilisco anche che tipo di traccia è:
    //se i due interni sono dello stesso poligono, la traccia è passante per quel poligono
    //se sono uno di un poligono e uno di un altro è non passante per entrambi
    //se inoltre i punti interni e esterni coincidono a due a due è passante per entrambi
    //ricordando che in intPoints ci sono prima due punti del primo poligono (0,1), poi due punti del secondo poligono (2,3)
    //escludo casi di punti coincidenti:
    bool done=false;
    bool intersection=true; //vede se c'è effettivamente intersezione
    double tol2=max(10*numeric_limits<double>::epsilon(), tol*tol);
    if(((intPoints[0]-intPoints[2]).squaredNorm()<tol2 && (intPoints[1]-intPoints[3]).squaredNorm()<tol2)||
        ((intPoints[1]-intPoints[2]).squaredNorm()<tol2 && (intPoints[0]-intPoints[3]).squaredNorm()<tol2)){
        tips={false,false};
        done=true;
        extremities = {intPoints[0],intPoints[1]};//di uno dei due poligoni
    }//passante per entrambi se i punti coincidono
    //vedo i casi in cui solo uno coincide
    else if ((intPoints[0]-intPoints[2]).squaredNorm()<tol2){
        done=true;
        if((intPoints[0]-intPoints[1]).dot(intPoints[0]-intPoints[3])>tol){
            if((intPoints[0]-intPoints[1]).dot(intPoints[3]-intPoints[1])>tol){
                //031
                tips={true,false};
                extremities={intPoints[0],intPoints[3]};}
            else{
                //013
                tips={false,true};
                extremities={intPoints[0],intPoints[1]};
            }
        }
        else{
            //103
            extremities={intPoints[0],intPoints[2]};//si vedrà dopo che ha lunghezza nulla
        }
    }
    else if ((intPoints[1]-intPoints[2]).squaredNorm()<tol2){
        done=true;
        if((intPoints[1]-intPoints[0]).dot(intPoints[1]-intPoints[3])>tol){
            if((intPoints[1]-intPoints[0]).dot(intPoints[3]-intPoints[0])>tol){
                //130
                tips={true,false};
                extremities={intPoints[1],intPoints[3]};}
            else{
                //103
                tips={false,true};
                extremities={intPoints[1],intPoints[0]};
            }
        }
        else{
            //013
            extremities={intPoints[1],intPoints[2]};//si vedrà dopo che ha lunghezza nulla
        }
    }
    else if ((intPoints[0]-intPoints[3]).squaredNorm()<tol2){
        done=true;
        if((intPoints[0]-intPoints[1]).dot(intPoints[0]-intPoints[2])>tol){
            if((intPoints[0]-intPoints[1]).dot(intPoints[2]-intPoints[1])>tol){
                //021
                tips={true,false};
                extremities={intPoints[0],intPoints[2]};}
            else{
                //012
                tips={false,true};
                extremities={intPoints[0],intPoints[1]};
            }
        }
        else{
            //102
            extremities={intPoints[0],intPoints[3]};//si vedrà dopo che ha lunghezza nulla
        }
    }
    else if ((intPoints[1]-intPoints[3]).squaredNorm()<tol2){
        done=true;
        if((intPoints[1]-intPoints[0]).dot(intPoints[1]-intPoints[2])>tol){
            if((intPoints[1]-intPoints[0]).dot(intPoints[2]-intPoints[0])>tol){
                //120
                tips={true,false};
                extremities={intPoints[1],intPoints[2]};}
            else{
                //102
                tips={false,true};
                extremities={intPoints[1],intPoints[0]};
            }
        }
        else{
            //012
            extremities={intPoints[1],intPoints[3]};//si vedrà dopo che ha lunghezza nulla
        }
    }

    if(!done){
        //se non è passante per entrambi, vedo la posizione reciproca con i prodotti scalari:
        if((intPoints[1]-intPoints[0]).dot(intPoints[2]-intPoints[0])>tol){ //confronto posizione di 2 rispetto a 0
            if((intPoints[0]-intPoints[1]).dot(intPoints[2]-intPoints[1])>tol){//2 rispetto a 1
                if((intPoints[0]-intPoints[2]).dot(intPoints[3]-intPoints[2])>tol){//3 rispetto a 2
                    if((intPoints[1]-intPoints[0]).dot(intPoints[3]-intPoints[0])>tol){//3 rispetto a 0
                        //0321
                        extremities = {intPoints[3],intPoints[2]};
                        tips[1]=false;
                        tips[0]=true; //S
                    }
                    else{
                        //3021
                        extremities = {intPoints[0],intPoints[2]};
                        tips[0]=true; //S
                        tips[1]=true; //S
                    }
                }
                else{
                    if((intPoints[0]-intPoints[1]).dot(intPoints[3]-intPoints[1])>tol){//3 rispetto a 1
                        //0231
                        extremities = {intPoints[3],intPoints[2]};
                        tips[1]=false;
                        tips[0]=true; //S
                    }
                    else{
                        //0213
                        extremities = {intPoints[1],intPoints[2]};
                        tips[0]=true; //S
                        tips[1]=true; //S
                    }
                }
            }
            else{
                if((intPoints[0]-intPoints[1]).dot(intPoints[3]-intPoints[1])>tol){//3 rispetto a 1
                    if((intPoints[1]-intPoints[0]).dot(intPoints[3]-intPoints[0])>tol){//3 rispetto a 0
                        //0312
                        extremities = {intPoints[3],intPoints[1]};
                        tips[0]=true; //S
                        tips[1]=true; //S
                    }
                    else{
                        //3012
                        extremities = {intPoints[0],intPoints[1]};
                        tips[0]=false;
                        tips[1]=true; //S
                    }
                }
                else{
                    if((intPoints[0]-intPoints[2]).dot(intPoints[3]-intPoints[2])>tol){//3 rispetto a 2
                        //0132
                        extremities = {intPoints[3],intPoints[1]};
                        intersection=false;
                    }
                    else{
                        //0123
                        extremities = {intPoints[1],intPoints[2]};
                        intersection=false;
                    }
                }
            }
        }
        else{
            if((intPoints[1]-intPoints[0]).dot(intPoints[3]-intPoints[0])>tol){ //3 rispetto a 0
                if((intPoints[0]-intPoints[1]).dot(intPoints[3]-intPoints[1])>tol){ //3 rispetto a 1
                    //2031
                    extremities = {intPoints[3],intPoints[0]};
                    tips[0]=true; //S
                    tips[1]=true; //S
                }
                else{
                    //2013
                    extremities = {intPoints[0],intPoints[1]};
                    tips[1]=true; //S (ho cambiato, era scritto false)
                    tips[0]=false;
                }
            }
            else{
                if((intPoints[0]-intPoints[2]).dot(intPoints[3]-intPoints[2])>tol){ //3 rispetto a 2
                    //2301
                    extremities = {intPoints[3],intPoints[0]};
                    intersection=false;
                }
                else{
                    //3201
                    extremities = {intPoints[0],intPoints[2]};
                    intersection=false;
                }
            }
        }

    }
    return intersection;

}

Vector3d intersectionLines(array<Vector3d,2>& line1, array<Vector3d,2>& line2){
    //y=mx+q
    double m1=(line1[1][1]-line1[0][1])/(line1[1][0]-line1[0][0]);
    double q1=line1[0][1]-m1*line1[0][0];
    double m2=(line2[1][1]-line2[0][1])/(line2[1][0]-line2[0][0]);
    double q2=line2[0][1]-m2*line2[0][0];


}

}

namespace detail { //per ordinare
//modifico la funzione già esistente di mergesort in modo da far ordinare le tracce per lunghezza

//template<typename T>
void merge(vector<unsigned int>& vecIdTraces, const vector<Trace>& traces, size_t left, size_t center, size_t right)
//devo dare in input non solo il vettore degli id delle tracce (passanti o non passanti) che è quello che verrà ordinato,
//ma anche il vettore contenente tutte le tracce in modo da poter risalire dall'id alla lunghezza della traccia
{
    assert(right >= left);
    size_t i = left;
    size_t j = center+1;
    size_t k = 0;

    vector<unsigned int> tmp(right - left + 1);

    while (i <= center && j <= right) {
        if (traces[vecIdTraces[i]].length >= traces[vecIdTraces[j]].length) { //cambio <= con >= per ordinare in modo decrescente e invece di prendere solo l'elemento ne guardo la lunghezza
            assert(k < tmp.size());
            assert(i < vecIdTraces.size());
            tmp[k++] = vecIdTraces[i++];
        }
        else {
            assert(k < tmp.size());
            assert(j < vecIdTraces.size());
            tmp[k++] = vecIdTraces[j++];
        }
    }

    while (i <= center) {
        assert(k < tmp.size());
        assert(i < vecIdTraces.size());
        tmp[k++] = vecIdTraces[i++];
    }

    while (j <= right) {
        assert(k < tmp.size());
        assert(j < vecIdTraces.size());
        tmp[k++] = vecIdTraces[j++];
    }

    assert(k == (right - left + 1));

    for (size_t h = left; h <= right; h++){
        vecIdTraces[h] = tmp[h-left];
    }
}

//template<typename T>
void mergesort(vector<unsigned int>& vecIdTraces, const vector<Trace>& traces, size_t left, size_t right)
{
    assert(left <= vecIdTraces.size());
    assert(right <= vecIdTraces.size());

    if (left < right) {
        size_t center = (left + right)/2;
        mergesort(vecIdTraces, traces, left, center);
        mergesort(vecIdTraces, traces, center+1, right);
        /* Ipotesi induttiva: [left, center] e
         * [center+1, right] sono ordinati.
         * Assumo merge() corretta */
        merge(vecIdTraces, traces, left, center, right);
        /* Qui [left,right] Ã¨ ordinato */
    }
    /* Caso base: il vettore di un solo elemento Ã¨ ordinato */
}

void mergesort(vector<unsigned int>& data, const vector<Trace>& traces)
{
    if (data.size()>0){ //va aggiunto il controllo perché se data.size()=0 avrei come right -1 (ma size_t non può essere negativo)
        mergesort(data, traces, 0, data.size()-1);
    }
}

} // namespace detail

namespace PolygonalMeshLibrary{
//vettore di PolygonalMesh dove a ogni frattura ne corrisponde una (nella stessa posizione)
vector<PolygonalMesh> cutFractures(vector<Fracture>& fractures, const vector <Trace>& traces, double tol){
    vector <PolygonalMesh> vec (fractures.size());
    for (Fracture& fr:fractures){
        unsigned int countIdV=0;
        PolygonalMesh mesh;
        list<Vector3d> verticesMesh; //a priori non so quanti saranno. Faccio una lista che poi copierò in un vettore
        list<unsigned int> idVerticesMesh;
        if (fr.idFrac!=-1){//escludo quelle problematiche
            list<Trace> allTraces;
            for (auto it=fr.passingTraces.begin(); it!=fr.passingTraces.end(); it++){
                allTraces.push_back(traces[*it]);
            }
            for (auto it=fr.notPassingTraces.begin(); it!=fr.passingTraces.end(); it++){
                allTraces.push_back(traces[*it]);
            }
        }
        //inserisco nella mesh i vertici della frattura (poi aggiungerò man mano quelli che trovo con i tagli)
        for (unsigned int i=0; i<fr.numVertices; i++){
            verticesMesh.push_back(fr.vertices[i]);
            idVerticesMesh.push_back(countIdV);
            countIdV++;
        }

    }
    return vec;
}
                                                           //gli passo anche la mesh?
void makeCuts (list<Vector3d>& vertices, list<Trace>& traces, double tol, int idFrac, PolygonalMesh& mesh, unsigned int& countIdV,
              list<Vector3d>& verticesMesh, list<unsigned int>& idVerticesMesh){
    //prende in input la lista con i vertici correnti del sottopoligono che deve ancora essere tagliato dalle tracce memorizzate in traces
    //scelgo di usare una lista perché, man mano che taglio il sottopoligono, si creano altri vertici che possono essere in mezzo a quelli
    //precedenti; inoltre non so a priori quanti saranno. Anche per le tracce uso una lista perché dovrò poi toglierle dalla testa man mano
    //che taglio il sottopoligono.
    //Vediamo da che parte della retta passante per la traccia sta il primo vertice
    bool firstVertex;
    bool previous;
    bool current;
    bool firstFracture;
    Vector3d previousVertex; //per capire quale dei due intPoints relativi alla frattura prendere (vedo quale sta sul lato tra i due vertici successivi)
    list<Vector3d> subvertices1; //ci memorizzo i vertici del primo sottopoligono dato dal taglio della prima traccia
    list<Vector3d> subvertices2;
    list<Trace> ;
    //PRIMA DI QUESTO FAI IL CONTROLLO CHE IL PROD VETTORIALE NON SIA 0 (CASO PARTICOLARE)
    Vector3d vecOnTheTrace=traces.front().extremitiesCoord[0]-traces.front().extremitiesCoord[1];
    previousVertex=vertices.front();
    firstVertex = findSideOfTheLine(vecOnTheTrace, previousVertex-traces.front().extremitiesCoord[0], tol);
    subvertices1.push_back(previousVertex);
    previous=firstVertex;
    //guardo per quella traccia se la frattura corrente è la prima o la seconda (per sapere quali intPoints prendere)
    if(idFrac==traces.front().fracturesIds[0]){
        firstFracture=true;
    }
    for (Vector3d& v:vertices){
        current = findSideOfTheLine(vecOnTheTrace,v-traces.front().extremitiesCoord[0],tol);
        if (current != previous){//aggiungo il vertice preso da intPoints a entrambi i sottopoligoni
            //cerco il vertice giusto da aggiungere
            if(firstFracture){ //so che è uno dei primi due punti di intPoints (quello sul lato corrente)
                if(abs((traces.front().intPoints[0]-previousVertex).dot(v-previousVertex))<tol){
                    verticesMesh.push_back(traces.front().intPoints[0]);
                    idVerticesMesh.push_back(countIdV);
                    countIdV++;
                    subvertices1.push_back(traces.front().intPoints[0]);
                    subvertices2.push_back(traces.front().intPoints[0]);
                }
                if(abs((traces.front().intPoints[1]-previousVertex).dot(v-previousVertex))<tol){
                    verticesMesh.push_back(traces.front().intPoints[1]);
                    idVerticesMesh.push_back(countIdV);
                    countIdV++;
                    subvertices1.push_back(traces.front().intPoints[1]);
                    subvertices2.push_back(traces.front().intPoints[1]);
                }

            }
            else{
                if(abs((traces.front().intPoints[2]-previousVertex).dot(v-previousVertex))<tol){
                    verticesMesh.push_back(traces.front().intPoints[2]);
                    idVerticesMesh.push_back(countIdV);
                    countIdV++;
                    subvertices1.push_back(traces.front().intPoints[2]);
                    subvertices2.push_back(traces.front().intPoints[2]);
                }
                if(abs((traces.front().intPoints[3]-previousVertex).dot(v-previousVertex))<tol){
                    verticesMesh.push_back(traces.front().intPoints[3]);
                    idVerticesMesh.push_back(countIdV);
                    countIdV++;
                    subvertices1.push_back(traces.front().intPoints[3]);
                    subvertices2.push_back(traces.front().intPoints[3]);
                }

            }
        }
        //controllo da che parte sta il vertice corrente (se dalla stessa del primo o dall'altra) e lo salvo nella lista di subvertices corrispondente
        //questo va fatto in entrambi i casi (sia current == previous sia !=)
        if(current==firstVertex){
            subvertices1.push_back(v);
        }
        else{
            subvertices2.push_back(v);
        }

        previousVertex=v;
    }

    /*if((coeff1.cross(coeff2)).squaredNorm()>tol2){ //i piani sono paralleli se hanno normali parallele
        //posizione dei punti del poligono 2 rispetto al piano 1
        bool positive=false; //per segnalare se il vertice analizzato in questo momento sta "sopra" (true) o "sotto"(false) il piano
        bool previous=false;
        if ((f2.vertices[0]).dot(coeff1)+d1>tol){
            previous=true;//vedo se si comincia "sopra" o "sotto"
        }
        unsigned int count1=0; //conto quanti vertici stanno dall'altra parte
        unsigned int firstVertexOtherSide2=0;
        onThePlane[1] = false; //metto anche il caso in cui si toccano senza che uno passi attraverso l'altro
        for (unsigned int i=1; i<f2.numVertices; i++){
            if(abs((f2.vertices[i]).dot(coeff1)+d1)<tol){
                //se ne ho due di fila sul piano, ho traccia particolare in cui uno non passa attraverso l'altro
                if(abs((f2.vertices[(i+1)%f2.numVertices]).dot(coeff1)+d1)<tol){
                    onThePlane[1]=true;
                    intPoints[2]=f2.vertices[i];
                    intPoints[3]=f2.vertices[(i+1)%f2.numVertices];
                    break;
                }
                else if(abs((f2.vertices[(i-1)%f2.numVertices]).dot(coeff1)+d1)<tol){
                    onThePlane[1]=true;
                    intPoints[2]=f2.vertices[(i-1)%f2.numVertices];
                    intPoints[3]=f2.vertices[i];
                    break;
                }

            }
            if ((f2.vertices[i]).dot(coeff1)+d1>tol){ //vedo da che parte stanno i vertici
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
            if((previous==positive) && (count1!=0)){ //conto quanti stanno dall'altra parte
                count1++;
            }

            previous=positive;
        }
    }*/

}

}






