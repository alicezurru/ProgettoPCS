#include "Utils.hpp"
#include "Fractures.hpp"

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
    ifstr.close();

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
                        array<Vector3d,2> extremities; //qui salverò i due punti estremi della traccia
                        array<bool,2> tips = {true,true}; //se resta così è non passante per entrambi
                        //vedo la posizione reciproca dei punti per stabilire i due più interni: saranno gli estremi della traccia
                        //inoltre così stabilisco anche che tipo di traccia è:
                        //se i due interni sono dello stesso poligono, la traccia è passante per quel poligono
                        //se sono uno di un poligono e uno di un altro è non passante per entrambi
                        //se inoltre i punti interni e esterni coincidono a due a due è passante per entrambi
                        //ricordando che in intPoints ci sono prima due punti del primo poligono (0,1), poi due punti del secondo poligono (2,3)
                        if(((intPoints[0]-intPoints[2]).squaredNorm()<tol*tol && (intPoints[1]-intPoints[3]).squaredNorm()<tol*tol)||
                            ((intPoints[1]-intPoints[2]).squaredNorm()<tol*tol && (intPoints[0]-intPoints[3]).squaredNorm()<tol*tol)){
                            tips={false,false};
                            extremities = {intPoints[0],intPoints[1]};//di uno dei due poligoni
                        }//passante per entrambi se i punti coincidono
                        if(tips[0]){
                            //se non è passante per entrambi, vedo la posizione reciproca con i prodotti scalari:
                            if((intPoints[1]-intPoints[0]).dot(intPoints[2]-intPoints[0])>0){ //confronto posizione di 2 rispetto a 0
                                if((intPoints[0]-intPoints[1]).dot(intPoints[2]-intPoints[1])>0){//2 rispetto a 1
                                    if((intPoints[0]-intPoints[2]).dot(intPoints[3]-intPoints[2])>0){//3 rispetto a 2
                                        if((intPoints[1]-intPoints[0]).dot(intPoints[3]-intPoints[0])>0){//3 rispetto a 0
                                            //0321
                                            extremities = {intPoints[3],intPoints[2]};
                                            tips[1]=false;
                                        }
                                        else{
                                            //3021
                                            extremities = {intPoints[0],intPoints[2]};
                                        }
                                    }
                                    else{
                                        if((intPoints[0]-intPoints[1]).dot(intPoints[3]-intPoints[1])>0){//3 rispetto a 1
                                            //0231
                                            extremities = {intPoints[3],intPoints[2]};
                                            tips[1]=false;
                                        }
                                        else{
                                            //0213
                                            extremities = {intPoints[1],intPoints[2]};
                                        }
                                    }
                                }
                                else{
                                    if((intPoints[0]-intPoints[1]).dot(intPoints[3]-intPoints[1])>0){//3 rispetto a 1
                                        if((intPoints[1]-intPoints[0]).dot(intPoints[3]-intPoints[0])>0){//3 rispetto a 0
                                            //0312
                                            extremities = {intPoints[3],intPoints[1]};
                                        }
                                        else{
                                            //3012
                                            extremities = {intPoints[0],intPoints[1]};
                                            tips[0]=false;
                                        }
                                    }
                                    else{
                                        if((intPoints[0]-intPoints[2]).dot(intPoints[3]-intPoints[2])>0){//3 rispetto a 2
                                            //0132
                                            extremities = {intPoints[3],intPoints[1]};
                                        }
                                        else{
                                            //0123
                                            extremities = {intPoints[1],intPoints[2]};
                                        }
                                    }
                                }
                            }
                            else{
                                if((intPoints[1]-intPoints[0]).dot(intPoints[3]-intPoints[0])>0){ //3 rispetto a 0
                                    if((intPoints[0]-intPoints[1]).dot(intPoints[3]-intPoints[1])>0){ //3 rispetto a 1
                                        //2031
                                        extremities = {intPoints[3],intPoints[0]};
                                    }
                                    else{
                                        //2013
                                        extremities = {intPoints[0],intPoints[1]};
                                        tips[1]=false;
                                    }
                                }
                                else{
                                    if((intPoints[0]-intPoints[2]).dot(intPoints[3]-intPoints[2])>0){ //3 rispetto a 2
                                        //2301
                                        extremities = {intPoints[3],intPoints[0]};
                                    }
                                    else{
                                        //3201
                                        extremities = {intPoints[0],intPoints[2]};
                                    }
                                }
                            }
                        }
                        //calcolo la lunghezza ed escludo il caso in cui sia un unico punto a toccare il poligono
                        double len =(extremities[0]-extremities[1]).norm();
                        if (len>tol){
                            //creo finalmente la traccia
                            Trace tr;
                            tr.idTr=listTraces.size();
                            tr.extremitiesCoord = extremities;
                            tr.fracturesIds = {fractures[i].idFrac,fractures[j].idFrac};
                            tr.length=len;
                            tr.Tips=tips;
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
    for (Trace& tr:traces){
        ofstr << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2" << endl;
        ofstr << tr.idTr << "; " << tr.fracturesIds[0] << "; " << tr.fracturesIds[1] << "; " << (tr.extremitiesCoord[0])[0] << "; "
              << (tr.extremitiesCoord[0])[1] << "; " << (tr.extremitiesCoord[0])[2] << "; " << (tr.extremitiesCoord[1])[0] << "; "
              << (tr.extremitiesCoord[1])[1] << "; " << (tr.extremitiesCoord[1])[2] << endl;
    }
    ofstr.close();
}

void printLocalResults (const string& fileName,const vector<Fracture>& fractures, const vector<Trace>& traces){
    //secondo file di ouput, con le informazioni sulle fratture e sulle tracce corrispondenti
    ofstream ofstr(fileName);
    for (Fracture fr:fractures){
        if (fr.idFrac!=-1){ //non gestisco quelle eventualmente problematiche
            ofstr << "# FractureId; NumTraces" << endl;
            ofstr << fr.idFrac << "; " << (fr.passingTraces.size())+(fr.notPassingTraces.size()) << endl;
            detail::mergesort(fr.passingTraces, traces, 0, fr.passingTraces.size()-1); //ordino le tracce passanti per lunghezza decrescente
            //??size-1 è giusto??  si:)
            detail::mergesort(fr.notPassingTraces, traces, 0, fr.notPassingTraces.size()-1); //ordino le tracce non passanti per lunghezza decrescente
            for (unsigned int trId:fr.passingTraces){
                if(traces[trId].length >0){ //non gestisco quelle eventualmente problematiche
                    ofstr << "# TraceId; Tips; Length" << endl;
                    ofstr << trId << "; " << "false; " << traces[trId].length << endl; //sto stampando prima tutte quelle passanti, quindi avranno tutte tips=false
                }
            }
            for (unsigned int trId:fr.notPassingTraces){
                if(traces[trId].length >0){ //non gestisco quelle eventualmente problematiche
                    ofstr << "# TraceId; Tips; Length" << endl;
                    ofstr << trId << "; " << "true; " << traces[trId].length << endl; //sto stampando prima tutte quelle non passanti, quindi avranno tutte tips=true
                }
            }
        }

    }
    ofstr.close();
}
}


namespace Algebra{

//per trovare l'equazione del piano che contiene i vertici di un poligono
inline Vector3d findPlaneEquation(vector<Vector3d>& points, double& constantTerm){ //restituisce la normale e modifica il dato in input che corrisponde al termine noto
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
            if((previous==positive) && (count1!=0)){ //conto quanti stanno dall'altra parte
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
            if ((f1.vertices[0]).dot(coeff2)+d2>0){
                previous=true;//vedo se si comincia "sopra" o "sotto"
            }
            unsigned int count2=0;
            unsigned int firstVertexOtherSide1=0;
            for (unsigned int i=1; i<f2.numVertices; i++){
                if ((f1.vertices[i]).dot(coeff2)+d2>0){ //vedo da che parte stanno i vertici
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

            if(count2!=0){ //c'è intersezione
                //ora individuo i punti di intersezione
                intPoints[0]=intersectionPlaneLine(coeff2, d2,f1.vertices[(firstVertexOtherSide1-1)],f1.vertices[firstVertexOtherSide1]);
                intPoints[1]=intersectionPlaneLine(coeff2, d2,f1.vertices[(firstVertexOtherSide1+count2-1)%f1.numVertices],f1.vertices[(firstVertexOtherSide1+count2)%f1.numVertices]);
                //i primi due punti sono del poligono 1, i successivi 2 del poligono 2
                intPoints[2]=intersectionPlaneLine(coeff1, d1,f2.vertices[(firstVertexOtherSide2-1)],f2.vertices[firstVertexOtherSide2]);
                intPoints[3]=intersectionPlaneLine(coeff1, d1,f2.vertices[(firstVertexOtherSide2+count1-1)%f2.numVertices],f2.vertices[(firstVertexOtherSide2+count1)%f2.numVertices]);
            }
            else{
                intersection=false;
            }
        }

        }
        return intersection;
}

}

/*template<typename T>
void print(const std::vector<T>& vs)
{
    for (auto & v : vs)
        std::cout << v << " ";
    std::cout << std::endl;
}*/


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

    for (size_t h = left; h <= right; h++)
        vecIdTraces[h] = tmp[h-left];
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

} // namespace detail

/*template<typename T>
void
mergesort(std::vector<T>& data)
{
    detail::mergesort(data, 0, data.size()-1);
}*/ //-> cos'é?





