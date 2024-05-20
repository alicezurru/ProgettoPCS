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
#include <queue>

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
   //r: p=p1+t(p2-p1)
    //[V1|V2](t1;t2)=(P2-P1)
    MatrixXd M(3, 2);
    M.col(0)=line1[1]-line1[0];
    M.col(1)=-(line2[1]-line2[0]);
    Vector2d t=M.householderQr().solve(line2[0]-line1[0]);//trovo t
    Vector3d intersectionPoint=line1[0]+t[0]*(line1[1]-line1[0]);
    return intersectionPoint;

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
        unsigned int countIdE=0;
        PolygonalMesh mesh;
        list<Vector3d> verticesMesh; //a priori non so quanti saranno. Faccio una lista che poi copierò in un vettore
        list<unsigned int> idVerticesMesh;
        list<array<unsigned int,2>> edgesMesh;
        list<unsigned int> idEdgesMesh;

        queue<Vector3d> vertices;//coda perchè mi serve solo inserire e riprendere nello stesso ordine (senso antiorario)
        queue<unsigned int> verticesId;
        queue<array<unsigned int,2>> edges;
        queue<unsigned int> edgesId;
        if (fr.idFrac!=-1){//escludo quelle problematiche
            queue<Trace> allTraces;
            for (auto it=fr.passingTraces.begin(); it!=fr.passingTraces.end(); it++){
                allTraces.push(traces[*it]);
            }
            for (auto it=fr.notPassingTraces.begin(); it!=fr.passingTraces.end(); it++){
                allTraces.push(traces[*it]);
            }
            //inserisco nella mesh i vertici e i lati della frattura (poi aggiungerò man mano quelli che trovo con i tagli)
            for (unsigned int i=0; i<fr.numVertices; i++){
                //vertici
                vertices.push(fr.vertices[i]);
                verticesId.push(countIdV);
                verticesMesh.push_back(fr.vertices[i]);
                idVerticesMesh.push_back(countIdV);
                //lati
                array<unsigned int,2> edge={countIdV,(countIdV+1)%fr.numVertices};
                edges.push(edge);
                edgesMesh.push_back(edge);//li aggiungo già da qua anche se poi magari li taglio? SECONDO ME NO -> RISOLVERE
                edgesId.push(countIdE);
                idEdgesMesh.push_back(countIdE);
                countIdV++;
                countIdE++;
            }

            //do il via ai tagli ricorsivi
            makeCuts(vertices,verticesId, allTraces,tol,mesh,countIdV,countIdE,edges,edgesId,verticesMesh,idVerticesMesh,edgesMesh,idEdgesMesh,fr.idFrac);


            //trasformo le liste di verticesMesh e idVerticesMesh in vettori e le aggiungo alla mesh
            vector<Vector3d> verticesMeshV;
            vector<unsigned int> idVerticesMeshV;
            idVerticesMeshV.reserve(countIdV);
            verticesMeshV.reserve(countIdV);
            for(unsigned int i=0;i<countIdV;i++){
                verticesMeshV.push_back(verticesMesh.front());
                verticesMesh.pop_front();
                idVerticesMeshV.push_back(idVerticesMesh.front());
                idVerticesMesh.pop_front();
            }
            mesh.numVertices=countIdV;
            mesh.idVertices=idVerticesMeshV;
            mesh.coordVertices=verticesMeshV;

            //anche per i lati
            vector<array<unsigned int,2>> edgesMeshV;
            vector<unsigned int> idEdgesMeshV;
            idEdgesMeshV.reserve(countIdE);
            edgesMeshV.reserve(countIdE);
            for(unsigned int i=0;i<countIdE;i++){
                edgesMeshV.push_back(edgesMesh.front());
                edgesMesh.pop_front();
                idEdgesMeshV.push_back(idEdgesMesh.front());
                idEdgesMesh.pop_front();
            }
            mesh.numEdges=countIdE;
            mesh.idEdges=idEdgesMeshV;
            mesh.extremitiesEdges=edgesMeshV;
        }

    }
    return vec;
}
void makeCuts (queue<Vector3d>& vertices, queue<unsigned int>& verticesId, queue<Trace>& traces, double tol, PolygonalMesh& mesh, unsigned int& countIdV, unsigned int& countIdE,
              queue<array<unsigned int,2>>& edges, queue<unsigned int>& edgesId, list<Vector3d>& verticesMesh, list<unsigned int>& idVerticesMesh,
              list<array<unsigned int,2>> edgesMesh,list<unsigned int> idEdgesMesh, int idFrac){
    //prende in input la coda con i vertici correnti del sottopoligono che deve ancora essere tagliato dalle tracce memorizzate in traces
    //scelgo di usare una coda perché, man mano che taglio il sottopoligono, si creano altri vertici che andranno sempre aggiunti in coda a quelli
    //del nuovo sottopoligono; inoltre non so a priori quanti saranno. Anche per le tracce uso una coda perché dovrò poi toglierle dalla testa man mano
    //che taglio il sottopoligono.
    if(traces.size()!=0){
        //Vediamo da che parte della retta passante per la traccia sta il primo vertice (1 se sta da una parte, 0 dall'altra e -1 se è sulla retta)
        int first;
        int previous;
        int current;
        bool onThePlaneCase; //per distinguere il caso particolare in cui la traccia cade su un lato
        bool previousSide1; //per sapere da che parte ho messo il vertice precedente (serve nel caso current sia -1)
        //CON QUESTO RIESCO ANCHE A TRATTARE IL CASO first==-1?
        queue<Vector3d> subvertices1; //ci memorizzo i vertici del primo sottopoligono dato dal taglio della prima traccia
        queue<Vector3d> subvertices2; //poi per metterle nella mesh TRASFORMA IN VECTOR
        queue<unsigned int> subverticesId1;//poi nella mesh memorizzo gli id per ogni poligono
        queue<unsigned int> subverticesId2;
        queue<array<unsigned int,2>> subedges1;
        queue<array<unsigned int,2>> subedges2;
        queue<unsigned int> subedgesId1;
        queue<unsigned int> subedgesId2;
        queue<Trace> subtraces1;
        queue<Trace> subtraces2;
        //PRIMA DI QUESTO FAI IL CONTROLLO CHE IL PROD VETTORIALE NON SIA 0 (CASI PARTICOLARI)

        Vector3d vecOnTheTrace=traces.front().extremitiesCoord[0]-traces.front().extremitiesCoord[1]; //cambia qualcosa per findSideOfTheLine se faccio 1-0??
        Vector3d firstVertex=vertices.front();
        vertices.pop();
        first = findSideOfTheLine(vecOnTheTrace, firstVertex-traces.front().extremitiesCoord[0], tol);
        subvertices1.push(firstVertex);
        subverticesId1.push(verticesId.front());
        if (first==-1){//se il primo sta sulla retta della traccia, va messo in entrambi i sottopoligoni
            //e scelgo arbitrariamente da che parte iniziare a salvare i vertici successivi
            subvertices2.push(firstVertex);
            subverticesId2.push(verticesId.front());
            previousSide1=false; //SERVE FALSE PER COME VIENE SALVATO L'ULTIMO LATO (?)
        }
        verticesId.pop();
        previous=first;
        Vector3d previousVertex=firstVertex;
        Vector3d lineTrace = traces.front().extremitiesCoord[0]-traces.front().extremitiesCoord[1];//o 1-0??? Forse si (controllare come funziona side of the line)
        //linetrace=vecOnTheTrace?

        //caso particolare onThePlane=true: non si generano due nuovi sottopoligoni ma solo uno che ha dei vertici in più sui lati già esistenti
        //salvo tutto nel sottopoligono1
        if((traces.front().onThePlane[0] && traces.front().fracturesIds[0]==idFrac) || (traces.front().onThePlane[1] && traces.front().fracturesIds[1]==idFrac)){
            onThePlaneCase=true;
            previousSide1=true;
        }
        for (unsigned int i=0;i<vertices.size();i++){
            Vector3d v=vertices.front();
            current = findSideOfTheLine(vecOnTheTrace,v-traces.front().extremitiesCoord[0],tol);
            if(current==-1 && previous==-1){//ovvero onThePlaneCase e sono sul lato interessato dalla traccia
                //non va bene fare come i casi precedenti perché non posso fare l'intersezione tra rette (sono parallele)
                //devo vedere se uno o entrambi gli estremi della traccia cadono dentro al lato

            }
            if ((current==0 && previous==1) || (current==1 && previous==0)){//aggiungo il nuovo vertice ad entrambi i sottopoligoni
                //trovo il nuovo vertice: intersezione tra prolungamento della traccia e lato
                array<Vector3d,2> lineEdge = {v,previousVertex};
                Vector3d intersection = intersectionLines(traces.front().extremitiesCoord,lineEdge);
                verticesMesh.push_back(intersection);
                idVerticesMesh.push_back(countIdV);
                subverticesId1.push(countIdV);
                subverticesId2.push(countIdV);
                countIdV++;
                subvertices1.push(intersection);
                subvertices2.push(intersection);
                //creo due nuovi lati a partire dal precedente
                array<unsigned int,2> newEdge1={edges.front()[0],countIdE}; //PERCHé countIdE? (quello dovrebbe essere l'id di tutto il lato) -> countIdV (prima di incrementarlo)??
                array<unsigned int,2> newEdge2={countIdE,edges.front()[1]};
                if(previous==first){ //vedo dove mettere i nuovi lati
                    //CAMBIARE
                    //aggiungo ogni volta i nuovi lati quindi conta come lato uno intero e poi anche le sue suddivisioni //L
                    edgesMesh.push_back(newEdge1); //L
                    idEdgesMesh.push_back(countIdE);//L
                    subedges1.push(newEdge1);
                    subedgesId1.push(countIdE);
                    countIdE++;
                    edgesMesh.push_back(newEdge2);//L
                    idEdgesMesh.push_back(countIdE);//L
                    subedges2.push(newEdge2);
                    subedgesId2.push(countIdE);
                }
                else{
                    edgesMesh.push_back(newEdge2);//L
                    idEdgesMesh.push_back(countIdE);//L
                    subedges1.push(newEdge2);
                    subedgesId2.push(countIdE);
                    countIdE++;
                    edgesMesh.push_back(newEdge1);//L
                    idEdgesMesh.push_back(countIdE);//L
                    subedges2.push(newEdge1);
                    subedgesId1.push(countIdE);
                }
                edgesId.pop();
                edges.pop();
                countIdE++;
                if(first==-1){//aggiungo qui i nuovi vertici nel caso in cui first==-1 perché mi serve sapere se ho "cambiato
                    //lato della traccia" o no (gli altri casi li tratto dopo)
                    //sono nel caso in cui ho cambiato lato
                    if(previousSide1){
                        subvertices2.push(v);
                        subverticesId2.push(verticesId.front());
                        previousSide1=false;
                    }
                    else{
                        subvertices1.push(v);
                        subverticesId1.push(verticesId.front());
                        previousSide1=true;
                    }
                }
            }
                                                                                                                      //caso current ==-1 lo tratto dopo (?)
            else if ((current==0 && previous ==0) || (current==1 && previous ==1) || (current!=-1 && previous==-1)){
                //se il lato non è tagliato dalla traccia (o se la traccia lo taglia solo su un estremo) devo aggiungerlo
                if(current==first || previousSide1){ //vedo dove mettere il lato
                    //va nell'1 se il vertice corrente è dalla stessa parte del primo. Nel caso in cui first==-1 non tratto qui il caso current==first
                    //ma guardo da che parte avevo salvato il vertice precedente e lo metto lì
                    subedges1.push(edges.front());
                    subedgesId1.push(edgesId.front());
                }
                else{
                    subedges2.push(edges.front());
                    subedgesId2.push(edgesId.front());
                }
                edges.pop();
                edgesId.pop();
                //vertici nel caso first==-1
                if(first==-1){//aggiungo qui i nuovi vertici nel caso in cui first==-1 perché mi serve sapere se ho "cambiato
                    //lato della traccia" o no (gli altri casi li tratto dopo)
                    //sono nel caso in cui non ho cambiato lato
                    if(previousSide1){
                        subvertices1.push(v);
                        subverticesId1.push(verticesId.front());
                        previousSide1=true;
                    }
                    else{
                        subvertices2.push(v);
                        subverticesId2.push(verticesId.front());
                        previousSide1=false;
                    }
                }
            }
            else{ //current==-1
                //vertice da entrambe le parti
                subvertices1.push(v);
                subverticesId1.push(verticesId.front());
                subvertices2.push(v);
                subverticesId2.push(verticesId.front());
                //lato dalla parte di quello precedente e non è tagliato
                if (previousSide1){
                    subedges1.push(edges.front());
                    subedgesId1.push(edgesId.front());
                }
                else{
                    subedges2.push(edges.front());
                    subedgesId2.push(edgesId.front());
                }
                if(!onThePlaneCase){//se sono nel caso onThePlane voglio continuare a salvare tutto nel sottopoligono 1
                    //(ed è l'unico caso in cui posso avere previous==-1 e current==-1)
                    previousSide1= !previousSide1; //lo inverto perché moralmente sto "passando dall'altra parte" della traccia
                }
            }
            //controllo da che parte sta il vertice corrente (se dalla stessa del primo o dall'altra) e lo salvo nella lista di subvertices corrispondente
            //questo va fatto in entrambi i casi (sia current == previous sia !=)
            //il caso current==-1 l'ho già trattato e anche quello first==-1
            if(current==first && current !=-1){
                //vertici
                subvertices1.push(v);
                subverticesId1.push(verticesId.front());
                previousSide1=true;
            }
            else if (current != first && current !=-1 && first!=-1){ //&& first!=-1??
                //vertici
                subvertices2.push(v);
                subverticesId2.push(verticesId.front());
                previousSide1=false;
                //verticesId.pop(); va bene metterlo una volta sola dopo?
            }

            previousVertex=v;
            verticesId.pop();
            vertices.pop();
            previous=current;
        }
        //manca ancora da controllare se c'è intersezione nell'ultimo lato
        if (first != previous && previous!=-1 && first!=-1){//aggiungo il nuovo vertice ad entrambi i sottopoligoni
            //trovo il nuovo vertice: intersezione tra prolungamento della traccia e lato
            array<Vector3d,2> lineEdge = {firstVertex,previousVertex};
            Vector3d intersection = intersectionLines(traces.front().extremitiesCoord,lineEdge);
            verticesMesh.push_back(intersection);
            idVerticesMesh.push_back(countIdV);
            subverticesId1.push(countIdV);
            subverticesId2.push(countIdV);
            countIdV++;
            subvertices1.push(intersection);
            subvertices2.push(intersection);
            //creo due nuovi lati a partire dal precedente
            array<unsigned int,2> newEdge1={edges.front()[0],countIdE};
            array<unsigned int,2> newEdge2={countIdE,edges.front()[1]};
            subedges1.push(newEdge2);
            subedgesId2.push(countIdE);
            idEdgesMesh.push_back(countIdE);//L
            edgesMesh.push_back(newEdge2);//L
            countIdE++;
            idEdgesMesh.push_back(countIdE);//L
            edgesMesh.push_back(newEdge1);//L
            subedges2.push(newEdge1);
            subedgesId1.push(countIdE);

            edgesId.pop();
            edges.pop();
            countIdE++;
        }
        else{//se l'ultimo lato non è tagliato dalla traccia devo aggiungerlo
            subedges1.push(edges.front());
            subedgesId1.push(edgesId.front());
            edges.pop();
            edgesId.pop();
        }

        //ora controllo da che parte stanno le tracce rimanenti
        traces.pop(); //tolgo la prima traccia: ho già tagliato
        for(unsigned int i=0;i<traces.size();i++){
            Trace t=traces.front();
            if((findSideOfTheLine(lineTrace,traces.front().extremitiesCoord[0],tol)==first)&&(findSideOfTheLine(lineTrace,traces.front().extremitiesCoord[1],tol)==first)){
                subtraces1.push(t);
            }
            else if((findSideOfTheLine(lineTrace,t.extremitiesCoord[0],tol)!=first)&&(findSideOfTheLine(lineTrace,t.extremitiesCoord[1],tol)!=first)){
                subtraces2.push(t);
            }
            else{//la traccia ha tagliato quest'altra traccia
                subtraces1.push(t);
                subtraces2.push(t);
            }
            traces.pop();
        }
        //ricorsione:
        makeCuts(subvertices1,subverticesId1,subtraces1,tol,mesh,countIdV,countIdE,subedges1,subedgesId1,verticesMesh,idVerticesMesh,edgesMesh,idEdgesMesh,idFrac);
        if (!onThePlaneCase){ //se onThePlane ho un solo sottopoligono
            makeCuts(subvertices2,subverticesId2,subtraces2,tol,mesh,countIdV,countIdE,subedges1,subedgesId1,verticesMesh,idVerticesMesh,edgesMesh,idEdgesMesh,idFrac);
        }

    }
    else{//vuol dire che per questo sottopoligono ho finito di tagliare
        //trasformo la coda in vettore e la aggiungo alla mesh
        vector<unsigned int> v1;
        v1.reserve(verticesId.size());
        for(unsigned int i=0;i<verticesId.size();i++){
            v1.push_back(verticesId.front());
            verticesId.pop();
        }
        mesh.verticesPolygons.push_back(v1);
        //stessa cosa per i lati
        vector<unsigned int> v2;
        v2.reserve(edgesId.size());
        for(unsigned int i=0;i<edgesId.size();i++){
            v2.push_back(edgesId.front());
            edgesId.pop();
        }
        mesh.edgesPolygons.push_back(v2);

    }
}

}






