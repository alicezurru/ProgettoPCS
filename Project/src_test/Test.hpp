#ifndef __TEST_HPP
#define __TEST_HPP

#include <gtest/gtest.h>
#include "Fractures.hpp"
#include "Utils.hpp"
#include "Eigen/Eigen"
#include <cmath>
#include <iostream>
#include "PolygonalMesh.hpp"

using namespace Eigen;
using namespace std;
using namespace Algebra;
using namespace PolygonalMeshLibrary;

namespace Geometry {

TEST(FRACTURESTEST, TestReadingFractures){
    //test per la lettura da file
    //bool readFractures(const string& fileName, vector<Fracture>& vec, double tol);
    double tol=10*numeric_limits<double>::epsilon();
    vector<Fracture> vec;
    string path="./DFN";
    bool flag=readFractures(path+"/FR3_data.txt",vec, tol);
    EXPECT_TRUE(flag);
    EXPECT_EQ(vec[0].idFrac,0);
    EXPECT_EQ(vec[1].idFrac,1);
    EXPECT_EQ(vec[2].idFrac,2);
    //controllo i vertici del primo
    EXPECT_TRUE(areVectorsEqual(vec[0].vertices[0],Vector3d(0,0,0),tol));
    EXPECT_TRUE(areVectorsEqual(vec[0].vertices[1],Vector3d(1,0,0),tol));
    EXPECT_TRUE(areVectorsEqual(vec[0].vertices[2],Vector3d(1,1,0),tol));
    EXPECT_TRUE(areVectorsEqual(vec[0].vertices[3],Vector3d(0,1,0),tol));
}

TEST(FRACTURETEST, TestFindTraces){
    //uso un file creato da noi (FR5_data_prova) e testo i risultati con quelli che conosciamo a priori
    double tol=10*numeric_limits<double>::epsilon();
    vector<Fracture> vec;
    readFractures("./DFN/FR5_data_prova.txt",vec, tol);
    vector<Trace> vecTraces=findTraces(vec,tol);
    array<bool, 2> traces_ok;
    EXPECT_EQ(vecTraces.size(), 4);
    traces_ok={false, false};
    EXPECT_EQ(vecTraces[0].Tips, traces_ok);
    traces_ok={true, false};
    EXPECT_EQ(vecTraces[1].Tips, traces_ok);
    EXPECT_EQ(vecTraces[2].Tips, traces_ok);
    EXPECT_EQ(vecTraces[3].Tips, traces_ok);
    //così controlliamo i tips. Le estremità le avevamo già controllate in TestFindIntersectionPoints e TestFindInternalPoints
    //abbiamo comunque controllato a mano che fossero giuste.

}

}

namespace Algebra{

TEST(ALGEBRATEST, TestFindPlaneEquation){
    vector<Vector3d> points = {Vector3d(1.0, 2.5, -0.5), Vector3d(2.0, -1.0, 0.0), Vector3d(0.0, 2.0, 1.0)};
    //non importa quanti punti gli do, userà sempre solo i primi 3, quindi posso dargli direttamente 3 punti
    double constantTerm=0;
    Vector3d n = findPlaneEquation(points, constantTerm);
    ASSERT_EQ(n, Vector3d(-5, -2, -4));
    ASSERT_EQ(constantTerm, 8);

    points = {Vector3d(1.0, 1.0, 0.0), Vector3d(0.5, 1.0, 1.0), Vector3d(2.0, 0.0, -4.0)};
    n=findPlaneEquation(points,constantTerm);
    ASSERT_EQ(n, Vector3d(1, -1, 0.5));
    ASSERT_EQ(constantTerm, 0);

}

//Vector3d intersectionPlaneLine(const Vector3d& coeff, const double d, const Vector3d& p1, const Vector3d& p2);
TEST(ALGEBRATEST, TestIntersectionPlaneLine){//non testo casi particolari come retta contenuta nel piano o assenza di intersezione perchè nel
    //programma si controlla e si fa in modo non vengano passati alla funzione questi casi
    Vector3d coeff = {1.0, -1.0, 1.0};
    double d= 1.0;
    Vector3d p1 =Vector3d(1.0, 0.0, 2.0);
    Vector3d p2 =Vector3d(-0.5, 1.5, -1.0);
    Vector3d intersection = intersectionPlaneLine(coeff, d, p1, p2);
    EXPECT_EQ(intersection, Vector3d(0,1,0));

}
//bool findIntersectionPoints(Fracture& f1, Fracture& f2, array<Vector3d,4>& intPoints, double tol, array<bool,2> onThePlane);
TEST(ALGEBRATEST, TestFindIntersectionPoints){
    double tol=10*numeric_limits<double>::epsilon();
    //testo un caso in cui so che c'è intersezione
    Fracture f1;
    f1.numVertices=4;
    f1.vertices={Vector3d(0,0,0),Vector3d(2,0,0),Vector3d(2,0,1),Vector3d(0,0,1)};
    Fracture f2;
    f2.numVertices=4;
    f2.vertices={Vector3d(1,-0.5,0.5),Vector3d(3,-0.5,0.5),Vector3d(3,0.5,0.5),Vector3d(1,0.5,0.5)};
    array<Vector3d,4> intPoints={Vector3d(0,0,0),Vector3d(0,0,0),Vector3d(0,0,0),Vector3d(0,0,0)};
    array<bool,2> onThePlane;
    bool pass = findIntersectionPoints(f1,f2,intPoints,tol,onThePlane);
    EXPECT_TRUE(pass);
    EXPECT_FALSE(onThePlane[0]);
    EXPECT_FALSE(onThePlane[1]);
    array<Vector3d,4> ok={Vector3d(2,0,0.5),Vector3d(0,0,0.5),Vector3d(3,0,0.5),Vector3d(1,0,0.5)};
    EXPECT_EQ(intPoints, ok);

    //testo un caso in cui non c'è intersezione perchè i piani sono paralleli
    f1.numVertices=5;
    f1.vertices={Vector3d(0,0,1),Vector3d(1,-1,1),Vector3d(0,-2,-1),Vector3d(-1,0,0),Vector3d(0.5,0.5,2)};
    f2.numVertices=3;
    f2.vertices={Vector3d(0,1,1),Vector3d(1,0,1),Vector3d(0.5,0.5,1)};
    pass=true;
    pass = findIntersectionPoints(f1,f2,intPoints,tol,onThePlane);
    EXPECT_FALSE(pass);
    EXPECT_FALSE(onThePlane[0]);
    EXPECT_FALSE(onThePlane[1]);

    //testo un caso limite: la traccia coincide con un lato di una frattura
    f1.numVertices=4;
    f1.vertices={Vector3d(1,-0.5,0.5),Vector3d(3,-0.5,0.5),Vector3d(3,0.5,0.5),Vector3d(1,0.5,0.5)};
    f2.numVertices=3;
    f2.vertices={Vector3d(2,0,0.5),Vector3d(2,1,1),Vector3d(1.5,0,0.5)};
    pass=false;
    pass = findIntersectionPoints(f1,f2,intPoints,tol,onThePlane);
    EXPECT_TRUE(pass);
    EXPECT_TRUE(onThePlane[1]);
    EXPECT_FALSE(onThePlane[0]);
    ok={Vector3d(3,0,0.5),Vector3d(1,0,0.5),Vector3d(1.5,0,0.5),Vector3d(2,0,0.5)};
    EXPECT_EQ(intPoints, ok);

    //testo un altro caso limite: un vertice è sul piano dell'altro ma la traccia non coincide con il lato di una frattura
    //f1 resta uguale
    f2.numVertices=4;
    f2.vertices={Vector3d(1.6,0,1),Vector3d(1.5,0,0.5),Vector3d(1.6,0,0),Vector3d(2,0,0.5)};
    pass=false;
    pass = findIntersectionPoints(f1,f2,intPoints,tol,onThePlane);
    EXPECT_TRUE(pass);
    EXPECT_FALSE(onThePlane[1]);
    EXPECT_FALSE(onThePlane[0]);
    ok={Vector3d(3,0,0.5),Vector3d(1,0,0.5),Vector3d(1.5,0,0.5),Vector3d(2,0,0.5)};
    EXPECT_EQ(intPoints, ok);


}
TEST(ALGEBRATEST, TestPassBoundingBox){
    //testo due poligoni vicini ma che non si intersecano
    Fracture f1;
    f1.numVertices=4;
    f1.vertices={Vector3d(0,0,0),Vector3d(2,0,0),Vector3d(2,0,1),Vector3d(0,0,1)};
    Fracture f2;
    f2.numVertices=5;
    f2.vertices={Vector3d(0,1,0),Vector3d(1.5,1,0),Vector3d(2,1,0.5),Vector3d(2,1,1),Vector3d(0,1,1)};
    bool pass=passBoundingBox(f1,f2);
    EXPECT_TRUE(pass);

    //testo due poligoni lontani
    f2.vertices={Vector3d(0,5,0),Vector3d(1.5,5,0),Vector3d(2,5,0.5),Vector3d(2,5,1),Vector3d(0,5,1)};
    pass=passBoundingBox(f1,f2);
    EXPECT_FALSE(pass);

}

TEST(ALGEBRATEST, TestFindInternalPoints){
    //testo un caso in cui non c'è intersezione ma i piani non sono paralleli
    double tol=10*numeric_limits<double>::epsilon();
    double tol2=max(10*numeric_limits<double>::epsilon(), tol*tol);
    Fracture f1;
    f1.numVertices=4;
    f1.vertices={Vector3d(0,0,0),Vector3d(2,0,0),Vector3d(2,0,1),Vector3d(0,0,1)};
    Fracture f2;
    f2.numVertices=4;
    f2.vertices={Vector3d(3,-0.5,0.5),Vector3d(5,-0.5,0.5),Vector3d(5,0.5,0.5),Vector3d(3,0.5,0.5)};
    array<bool,2> onThePlane = {false, false}; //nessuno dei vertici sta sul piano dell'altra frattura
    //DOMANDA: devo usare findIntersectionPoints per trovarla??
    array<Vector3d,4> intPoints; //= {Vector3d(0,0,0.5), Vector3d(2,0,0.5), Vector3d(3,0,0.5), Vector3d(5,0,0.5)};
    findIntersectionPoints(f1, f2, intPoints, tol, onThePlane);
    array<Vector3d, 2> extremities;
    array<bool, 2> tips;
    bool intersection = findInternalPoints(intPoints, tol, extremities, tips);
    EXPECT_FALSE(intersection);

    //testo un caso in cui c'è intersezione e la traccia è non passante per entrambe le fratture e f1 è poggiato sopra il piano di f2
    f1.numVertices=3;
    f1.vertices = {Vector3d(4,0,0.5), Vector3d(6,0,0.5), Vector3d(5,0.5,1)};
    //come frattura f2 tengo la stessa di prima
    findIntersectionPoints(f1, f2, intPoints, tol, onThePlane); //aggiorno i valori di intPoints da passare a findInternalpoints
    intersection = findInternalPoints(intPoints, tol, extremities, tips);
    array<Vector3d,2> extremities_ok = {Vector3d(4,0,0.5), Vector3d(5,0,0.5)};
    EXPECT_TRUE(intersection);
    EXPECT_EQ(extremities, extremities_ok);

    //testo un caso in cui c'è intersezione e la traccia è non passante per entrambe le fratture
    f1.numVertices=3;
    f1.vertices = {Vector3d(4,0,0.4), Vector3d(6,0,0.4), Vector3d(5,0.5,1)};
    //come frattura f2 tengo la stessa di prima
    findIntersectionPoints(f1, f2, intPoints, tol, onThePlane); //aggiorno i valori di intPoints da passare a findInternalpoints
    intersection = findInternalPoints(intPoints, tol, extremities, tips);
    extremities_ok = {Vector3d(4.166666666666667,0.083333333333333,0.5), Vector3d(5,0.083333333333333,0.5)};
    array <bool, 2> tips_ok = {true, true}; //non passante per entrambe
    EXPECT_TRUE(intersection);
    EXPECT_TRUE(areVectorsEqual(extremities[0], extremities_ok[0], tol));
    EXPECT_TRUE(areVectorsEqual(extremities[1], extremities_ok[1], tol));
    EXPECT_EQ(tips, tips_ok);

    //testo un caso in cui c'è intersezione e la traccia è passante per entrambe le fratture
    f1.numVertices=4;
    f1.vertices = {Vector3d(4,-0.5,0), Vector3d(4,0.5,0), Vector3d(4,0.5,1), Vector3d(4,-0.5,1)};
    //come frattura f2 tengo la stessa di prima
    findIntersectionPoints(f1, f2, intPoints, tol, onThePlane); //aggiorno i valori di intPoints da passare a findInternalpoints
    intersection = findInternalPoints(intPoints, tol, extremities, tips);
    extremities_ok = {Vector3d(4,0.5,0.5),Vector3d(4,-0.5,0.5)};
    tips_ok = {false, false}; //passante per entrambe
    EXPECT_TRUE(intersection);
    EXPECT_TRUE(areVectorsEqual(extremities[0], extremities_ok[0], tol));
    EXPECT_TRUE(areVectorsEqual(extremities[1], extremities_ok[1], tol));
    EXPECT_EQ(tips, tips_ok);

    //testo un caso in cui c'è intersezione e la traccia è passante per una frattura (triangolo f1) e non per l'altra (rettangolo f2)
    f1.numVertices=3;
    f1.vertices = {Vector3d(4,-0.5,0), Vector3d(4,0.5,0), Vector3d(4,0,1)};
    //come frattura f2 tengo la stessa di prima
    findIntersectionPoints(f1, f2, intPoints, tol, onThePlane); //aggiorno i valori di intPoints da passare a findInternalpoints
    intersection = findInternalPoints(intPoints, tol, extremities, tips);
    extremities_ok = {Vector3d(4,0.25,0.5),Vector3d(4,-0.25,0.5)};
    tips_ok = {false, true}; //passante per entrambe
    EXPECT_TRUE(intersection);
    EXPECT_TRUE(areVectorsEqual(extremities[0], extremities_ok[0], tol));
    EXPECT_TRUE(areVectorsEqual(extremities[1], extremities_ok[1], tol));
    EXPECT_EQ(tips, tips_ok);

    //testo un caso particolare in cui c'è un solo punto di intersezione (traccia di lunghezza nulla)
    f1.vertices = {Vector3d(4,0,0.5), Vector3d(4,1,1), Vector3d(3.5,-1,1)};
    //come frattura f2 tengo la stessa di prima
    findIntersectionPoints(f1, f2, intPoints, tol, onThePlane); //aggiorno i valori di intPoints da passare a findInternalpoints
    intersection = findInternalPoints(intPoints, tol, extremities, tips);
    EXPECT_FALSE(intersection);
}

//int findSideOfTheLine (const Vector3d &vecLine, const Vector3d &vecToTest, const Vector3d &n, double tol);
TEST(ALGEBRATEST, TestFindSideOfTheLine){
    //caso base
    double tol=10*numeric_limits<double>::epsilon();
    Vector3d n = Vector3d(0,1,0);
    Vector3d vecLine= Vector3d(-4,0,0);
    Vector3d vecToTest= Vector3d(0,0,2);
    int i =findSideOfTheLine(vecLine,vecToTest,n,tol);
    EXPECT_EQ(i,0);

    vecToTest= Vector3d(0,0,-2);
    i =findSideOfTheLine(vecLine,vecToTest,n,tol);
    EXPECT_EQ(i,1);

    //considero un punto molto vicino alla retta
    vecToTest= Vector3d(0,0,2*0.00000000000001);
    i =findSideOfTheLine(vecLine,vecToTest,n,tol);
    EXPECT_EQ(i,0);

    //considero un punto molto vicino alla retta: così vicino che è sul piano
    vecToTest= Vector3d(0,0,2*0.0000000000000000001);
    i =findSideOfTheLine(vecLine,vecToTest,n,tol);
    EXPECT_EQ(i,-1);

}







}

namespace detail{

TEST(SORTINGTEST, TestMergeSort){
    vector<Trace> vecProvaMS(5);
    Trace t0;
    t0.idTr=0;
    t0.length=1.0;
    Trace t1;
    t1.idTr=1;
    t1.length=5.0;
    Trace t2;
    t2.idTr=2;
    t2.length=2.0;
    Trace t3;
    t3.idTr=3;
    t3.length=1.0;
    Trace t4;
    t4.idTr=4;
    t4.length=1.5;
    vector<unsigned int> idProvaMS ={0,1,2,3,4};
    //ordine:124(03)
    vecProvaMS[0]=t0;
    vecProvaMS[1]=t1;
    vecProvaMS[2]=t2;
    vecProvaMS[3]=t3;
    vecProvaMS[4]=t4;
    mergesort(idProvaMS, vecProvaMS);
    EXPECT_EQ(idProvaMS[0],1);
    EXPECT_EQ(idProvaMS[1],2);
    EXPECT_EQ(idProvaMS[2],4);
    EXPECT_EQ(idProvaMS[3],0);
    EXPECT_EQ(idProvaMS[4],3);


}

}

namespace PolygonalMeshLibrary{

/*void makeCuts (queue<Vector3d>& vertices, queue<unsigned int>& verticesId, queue<Trace>& traces, double tol, PolygonalMesh& mesh, unsigned int& countIdV, unsigned int& countIdE,
              list<Vector3d>& verticesMesh, list<unsigned int>& idVerticesMesh,
              list<array<unsigned int,2>>& edgesMesh,list<unsigned int>& idEdgesMesh, int idFrac, Vector3d& n,map<array<unsigned int,2>,unsigned int>& mapEdges);*/

TEST(POLYGONTEST,TestMakeCuts){
    //caso base: 1 traccia passante
    double tol=10*numeric_limits<double>::epsilon();
    Vector3d n = Vector3d(0,0,1);
    unsigned int countIdV=4;
    unsigned int countIdE=0;
    PolygonalMesh mesh;
    Vector3d a = Vector3d(0,0,0);
    Vector3d b = Vector3d(2,0,0);
    Vector3d c = Vector3d(2,2,0);
    Vector3d d = Vector3d(0,2,0);
    list<Vector3d> verticesMesh={a,b,c,d};
    list<unsigned int> idVerticesMesh={0,1,2,3};
    list<unsigned int> idEdgesMesh={0,1,2,3};
    list<array<unsigned int,2>> edgesMesh;
    edgesMesh.push_back({0,1});
    edgesMesh.push_back({1,2});
    edgesMesh.push_back({2,3});
    edgesMesh.push_back({3,0});
    map<array<unsigned int,2>,unsigned int> mapEdges;
    queue<unsigned int> verticesId;
    verticesId.push(0);
    verticesId.push(1);
    verticesId.push(2);
    verticesId.push(3);
    queue<Vector3d> vertices;
    vertices.push(a);
    vertices.push(b);
    vertices.push(c);
    vertices.push(d);
    queue<Trace> allTraces;
    Trace tr;

    tr.idTr=0;
    tr.extremitiesCoord={Vector3d(1,0,0),Vector3d(1,2,0)};
    tr.length=2;
    tr.fracturesIds={0,1};
    tr.Tips={false,false};
    tr.onThePlane={false,false};
    allTraces.push(tr);
    makeCuts(vertices,verticesId, allTraces,tol,mesh,countIdV,countIdE,verticesMesh,
             idVerticesMesh,edgesMesh,idEdgesMesh,0,n,mapEdges);
    EXPECT_EQ(countIdV,6);
    EXPECT_EQ(countIdE,7);
    EXPECT_EQ(mesh.verticesPolygons[0].size(),4);
    EXPECT_EQ(mesh.verticesPolygons[1].size(),4);
    EXPECT_EQ(mesh.edgesPolygons[0].size(),4);
    EXPECT_EQ(mesh.edgesPolygons[1].size(),4);

    //poligono1
    EXPECT_EQ(mesh.verticesPolygons[0][0],0);
    EXPECT_EQ(mesh.verticesPolygons[0][1],4);
    EXPECT_EQ(mesh.verticesPolygons[0][2],5);
    EXPECT_EQ(mesh.verticesPolygons[0][3],3);

    EXPECT_EQ(mesh.edgesPolygons[0][0],0);
    EXPECT_EQ(mesh.edgesPolygons[0][1],1);
    EXPECT_EQ(mesh.edgesPolygons[0][2],2);
    EXPECT_EQ(mesh.edgesPolygons[0][3],3);

    //poligono2
    EXPECT_EQ(mesh.verticesPolygons[1][0],4);
    EXPECT_EQ(mesh.verticesPolygons[1][1],1);
    EXPECT_EQ(mesh.verticesPolygons[1][2],2);
    EXPECT_EQ(mesh.verticesPolygons[1][3],5);

    EXPECT_EQ(mesh.edgesPolygons[1][0],4);
    EXPECT_EQ(mesh.edgesPolygons[1][1],5);
    EXPECT_EQ(mesh.edgesPolygons[1][2],6);
    EXPECT_EQ(mesh.edgesPolygons[1][3],1);

    //vertici
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),a,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),b,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),c,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),d,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),Vector3d(1,0,0),tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),Vector3d(1,2,0),tol));
    verticesMesh.pop_front();

    //lati
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),{0,4},tol));
    edgesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),{4,5},tol));
    edgesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),{5,3},tol));
    edgesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),{3,0},tol));
    edgesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),{4,1},tol));
    edgesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),{1,2},tol));
    edgesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),{2,5},tol));
    edgesMesh.pop_front();

    //caso base: 1 traccia non passante
    //stessi dati di prima ma traccia non passante
    countIdV=4;
    countIdE=0;
    tr.extremitiesCoord={Vector3d(1,0.5,0),Vector3d(1,1.5,0)};
    tr.length=1;
    allTraces.push(tr);
    verticesMesh={a,b,c,d};
    edgesMesh.push_back({0,1});
    edgesMesh.push_back({1,2});
    edgesMesh.push_back({2,3});
    edgesMesh.push_back({3,0});
    mapEdges={};
    verticesId={};
    verticesId.push(0);
    verticesId.push(1);
    verticesId.push(2);
    verticesId.push(3);
    vertices={};
    vertices.push(a);
    vertices.push(b);
    vertices.push(c);
    vertices.push(d);

    makeCuts(vertices,verticesId, allTraces,tol,mesh,countIdV,countIdE,verticesMesh,
             idVerticesMesh,edgesMesh,idEdgesMesh,0,n,mapEdges);
    EXPECT_EQ(countIdV,6);
    EXPECT_EQ(countIdE,7);
    EXPECT_EQ(mesh.verticesPolygons[0].size(),4);
    EXPECT_EQ(mesh.verticesPolygons[1].size(),4);
    EXPECT_EQ(mesh.edgesPolygons[0].size(),4);
    EXPECT_EQ(mesh.edgesPolygons[1].size(),4);

    //poligono1
    EXPECT_EQ(mesh.verticesPolygons[0][0],0);
    EXPECT_EQ(mesh.verticesPolygons[0][1],4);
    EXPECT_EQ(mesh.verticesPolygons[0][2],5);
    EXPECT_EQ(mesh.verticesPolygons[0][3],3);

    EXPECT_EQ(mesh.edgesPolygons[0][0],0);
    EXPECT_EQ(mesh.edgesPolygons[0][1],1);
    EXPECT_EQ(mesh.edgesPolygons[0][2],2);
    EXPECT_EQ(mesh.edgesPolygons[0][3],3);

    //poligono2
    EXPECT_EQ(mesh.verticesPolygons[1][0],4);
    EXPECT_EQ(mesh.verticesPolygons[1][1],1);
    EXPECT_EQ(mesh.verticesPolygons[1][2],2);
    EXPECT_EQ(mesh.verticesPolygons[1][3],5);

    EXPECT_EQ(mesh.edgesPolygons[1][0],4);
    EXPECT_EQ(mesh.edgesPolygons[1][1],5);
    EXPECT_EQ(mesh.edgesPolygons[1][2],6);
    EXPECT_EQ(mesh.edgesPolygons[1][3],1);

    //vertici
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),a,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),b,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),c,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),d,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),Vector3d(1,0,0),tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),Vector3d(1,2,0),tol));
    verticesMesh.pop_front();

    //lati
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),{0,4},tol));
    edgesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),{4,5},tol));
    edgesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),{5,3},tol));
    edgesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),{3,0},tol));
    edgesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),{4,1},tol));
    edgesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),{1,2},tol));
    edgesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),{2,5},tol));
    edgesMesh.pop_front();



}







}


#endif
