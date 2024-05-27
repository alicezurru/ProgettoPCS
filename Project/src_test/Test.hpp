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

TEST(FRACTURESTEST, TestReadFractures){
    //test per la lettura da file
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
    //così controlliamo i tips. Le estremità le controlliamo in TestFindIntersectionPoints e TestFindInternalPoints
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

TEST(ALGEBRATEST, TestIntersectionPlaneLine){//non testo casi particolari come retta contenuta nel piano o assenza di intersezione perchè nel
    //programma si controlla e si fa in modo non vengano passati alla funzione questi casi
    Vector3d coeff = {1.0, -1.0, 1.0};
    double d= 1.0;
    Vector3d p1 =Vector3d(1.0, 0.0, 2.0);
    Vector3d p2 =Vector3d(-0.5, 1.5, -1.0);
    Vector3d intersection = intersectionPlaneLine(coeff, d, p1, p2);
    EXPECT_EQ(intersection, Vector3d(0,1,0));

}

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
    Fracture f1;
    f1.numVertices=4;
    f1.vertices={Vector3d(0,0,0),Vector3d(2,0,0),Vector3d(2,0,1),Vector3d(0,0,1)};
    Fracture f2;
    f2.numVertices=4;
    f2.vertices={Vector3d(3,-0.5,0.5),Vector3d(5,-0.5,0.5),Vector3d(5,0.5,0.5),Vector3d(3,0.5,0.5)};
    array<bool,2> onThePlane = {false, false}; //nessuno dei vertici sta sul piano dell'altra frattura
    array<Vector3d,4> intPoints;
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

TEST(ALGEBRATEST, TestIntersectionLines){
    double tol=10*numeric_limits<double>::epsilon();
    //test 1
    //punti retta 1
    Vector3d A=Vector3d(0,0,0);
    Vector3d B=Vector3d(1,0,0);
    array<Vector3d,2> line1={A,B};
    //punti retta 2
    Vector3d C=Vector3d(0,1,0);
    Vector3d D=Vector3d(0,-1,0);
    array<Vector3d,2> line2={C,D};
    Vector3d intersection =intersectionLines(line1,line2);
    EXPECT_TRUE(areVectorsEqual(A,intersection,tol));

    //test 2
    //punti retta 1
    A=Vector3d(0,0,0);
    B=Vector3d(1,1,0);
    line1={A,B};
    //punti retta 2
    C=Vector3d(0,2,0);
    D=Vector3d(2,0,0);
    line2={C,D};
    intersection =intersectionLines(line1,line2);
    EXPECT_TRUE(areVectorsEqual(Vector3d(1,1,0),intersection,tol));
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
    list<Trace> allTraces;
    Trace tr;
    tr.idTr=0;
    tr.extremitiesCoord={Vector3d(1,0,0),Vector3d(1,2,0)};
    tr.length=2;
    tr.fracturesIds={0,1};
    tr.Tips={false,false};
    tr.onThePlane={false,false};
    allTraces.push_back(tr);
    list<reference_wrapper<Trace>> traceRefs(allTraces.begin(), allTraces.end());
    makeCuts(vertices,verticesId, traceRefs,tol,mesh,countIdV,countIdE,verticesMesh,
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
    array<unsigned int,2> ok={0,4};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={4,5};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={5,3};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={3,0};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={4,1};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={1,2};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={2,5};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();

    //caso base: 1 traccia non passante
    //stessi dati di prima ma traccia non passante
    countIdV=4;
    countIdE=0;
    tr.extremitiesCoord={Vector3d(1,0.5,0),Vector3d(1,1.5,0)};
    tr.length=1;
    allTraces={};
    allTraces.push_back(tr);
    verticesMesh={a,b,c,d};
    edgesMesh={};
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
    mesh.verticesPolygons={};
    mesh.edgesPolygons={};
    traceRefs={};
    traceRefs.assign(allTraces.begin(), allTraces.end());
    makeCuts(vertices,verticesId, traceRefs,tol,mesh,countIdV,countIdE,verticesMesh,
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
    ok={0,4};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={4,5};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={5,3};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={3,0};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={4,1};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={1,2};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={2,5};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();


    //caso base: 1 traccia non passante
    //stessi dati di prima ma traccia che incontra l'ultimo lato (che abbiamo gestito diversamente)
    countIdV=4;
    countIdE=0;
    tr.extremitiesCoord={Vector3d(0.5,0.5,0),Vector3d(1.5,0.5,0)};
    tr.length=1;
    allTraces={};
    allTraces.push_back(tr);
    verticesMesh={a,b,c,d};
    edgesMesh={};
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
    mesh.verticesPolygons={};
    mesh.edgesPolygons={};
    traceRefs={};
    traceRefs.assign(allTraces.begin(), allTraces.end());
    makeCuts(vertices,verticesId, traceRefs,tol,mesh,countIdV,countIdE,verticesMesh,
             idVerticesMesh,edgesMesh,idEdgesMesh,0,n,mapEdges);
    EXPECT_EQ(countIdV,6);
    EXPECT_EQ(countIdE,7);
    EXPECT_EQ(mesh.verticesPolygons[0].size(),4);
    EXPECT_EQ(mesh.verticesPolygons[1].size(),4);
    EXPECT_EQ(mesh.edgesPolygons[0].size(),4);
    EXPECT_EQ(mesh.edgesPolygons[1].size(),4);

    //poligono1
    EXPECT_EQ(mesh.verticesPolygons[0][0],0);
    EXPECT_EQ(mesh.verticesPolygons[0][1],1);
    EXPECT_EQ(mesh.verticesPolygons[0][2],4);
    EXPECT_EQ(mesh.verticesPolygons[0][3],5);

    EXPECT_EQ(mesh.edgesPolygons[0][0],0);
    EXPECT_EQ(mesh.edgesPolygons[0][1],1);
    EXPECT_EQ(mesh.edgesPolygons[0][2],2);
    EXPECT_EQ(mesh.edgesPolygons[0][3],3);

    //poligono2
    EXPECT_EQ(mesh.verticesPolygons[1][0],4);
    EXPECT_EQ(mesh.verticesPolygons[1][1],2);
    EXPECT_EQ(mesh.verticesPolygons[1][2],3);
    EXPECT_EQ(mesh.verticesPolygons[1][3],5);

    EXPECT_EQ(mesh.edgesPolygons[1][0],4);
    EXPECT_EQ(mesh.edgesPolygons[1][1],5);
    EXPECT_EQ(mesh.edgesPolygons[1][2],6);
    EXPECT_EQ(mesh.edgesPolygons[1][3],2);

    //vertici
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),a,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),b,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),c,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),d,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),Vector3d(2,0.5,0),tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),Vector3d(0,0.5,0),tol));
    verticesMesh.pop_front();

    //lati
    ok={0,1};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={1,4};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={4,5};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={5,0};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={4,2};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={2,3};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={3,5};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();

    //traccia sul primo vertice (first=-1)
    countIdV=4;
    countIdE=0;
    tr.extremitiesCoord={Vector3d(0,0,0),Vector3d(1,2,0)};
    tr.length=sqrt(5);
    allTraces={};
    allTraces.push_back(tr);
    verticesMesh={a,b,c,d};
    edgesMesh={};
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
    mesh.verticesPolygons={};
    mesh.edgesPolygons={};
    traceRefs={};
    traceRefs.assign(allTraces.begin(), allTraces.end());
    makeCuts(vertices,verticesId, traceRefs,tol,mesh,countIdV,countIdE,verticesMesh,
             idVerticesMesh,edgesMesh,idEdgesMesh,0,n,mapEdges);
    EXPECT_EQ(countIdV,5);
    EXPECT_EQ(countIdE,6);
    EXPECT_EQ(mesh.verticesPolygons[0].size(),4);
    EXPECT_EQ(mesh.verticesPolygons[1].size(),3);
    EXPECT_EQ(mesh.edgesPolygons[0].size(),4);
    EXPECT_EQ(mesh.edgesPolygons[1].size(),3);

    //poligono1
    EXPECT_EQ(mesh.verticesPolygons[0][0],0);
    EXPECT_EQ(mesh.verticesPolygons[0][1],1);
    EXPECT_EQ(mesh.verticesPolygons[0][2],2);
    EXPECT_EQ(mesh.verticesPolygons[0][3],4);

    EXPECT_EQ(mesh.edgesPolygons[0][0],0);
    EXPECT_EQ(mesh.edgesPolygons[0][1],1);
    EXPECT_EQ(mesh.edgesPolygons[0][2],2);
    EXPECT_EQ(mesh.edgesPolygons[0][3],3);

    //poligono2
    EXPECT_EQ(mesh.verticesPolygons[1][0],0);
    EXPECT_EQ(mesh.verticesPolygons[1][1],4);
    EXPECT_EQ(mesh.verticesPolygons[1][2],3);


    EXPECT_EQ(mesh.edgesPolygons[1][0],3);
    EXPECT_EQ(mesh.edgesPolygons[1][1],4);
    EXPECT_EQ(mesh.edgesPolygons[1][2],5);


    //vertici
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),a,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),b,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),c,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),d,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),Vector3d(1,2,0),tol));
    verticesMesh.pop_front();


    //lati
    ok={0,1};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={1,2};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={2,4};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={4,0};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={4,3};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={3,0};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();

    //traccia sul primo vertice (first=-1) e sul terzo
    countIdV=4;
    countIdE=0;
    tr.extremitiesCoord={Vector3d(0,0,0),Vector3d(2,2,0)};
    tr.length=sqrt(8);
    allTraces={};
    allTraces.push_back(tr);
    verticesMesh={a,b,c,d};
    edgesMesh={};
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
    mesh.verticesPolygons={};
    mesh.edgesPolygons={};
    traceRefs={};
    traceRefs.assign(allTraces.begin(), allTraces.end());
    makeCuts(vertices,verticesId, traceRefs,tol,mesh,countIdV,countIdE,verticesMesh,
             idVerticesMesh,edgesMesh,idEdgesMesh,0,n,mapEdges);
    EXPECT_EQ(countIdV,4);
    EXPECT_EQ(countIdE,5);
    EXPECT_EQ(mesh.verticesPolygons[0].size(),3);
    EXPECT_EQ(mesh.verticesPolygons[1].size(),3);
    EXPECT_EQ(mesh.edgesPolygons[0].size(),3);
    EXPECT_EQ(mesh.edgesPolygons[1].size(),3);

    //poligono1
    EXPECT_EQ(mesh.verticesPolygons[0][0],0);
    EXPECT_EQ(mesh.verticesPolygons[0][1],1);
    EXPECT_EQ(mesh.verticesPolygons[0][2],2);

    EXPECT_EQ(mesh.edgesPolygons[0][0],0);
    EXPECT_EQ(mesh.edgesPolygons[0][1],1);
    EXPECT_EQ(mesh.edgesPolygons[0][2],2);

    //poligono2
    EXPECT_EQ(mesh.verticesPolygons[1][0],0);
    EXPECT_EQ(mesh.verticesPolygons[1][1],2);
    EXPECT_EQ(mesh.verticesPolygons[1][2],3);


    EXPECT_EQ(mesh.edgesPolygons[1][0],2);
    EXPECT_EQ(mesh.edgesPolygons[1][1],3);
    EXPECT_EQ(mesh.edgesPolygons[1][2],4);


    //vertici
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),a,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),b,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),c,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),d,tol));
    verticesMesh.pop_front();


    //lati
    ok={0,1};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={1,2};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={2,0};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={2,3};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={3,0};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();


    //traccia su un vertice che non è il primo
    countIdV=4;
    countIdE=0;
    tr.extremitiesCoord={Vector3d(0,2,0),Vector3d(2,1,0)};
    tr.length=sqrt(5);
    allTraces={};
    allTraces.push_back(tr);
    verticesMesh={a,b,c,d};
    edgesMesh={};
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
    mesh.verticesPolygons={};
    mesh.edgesPolygons={};
    traceRefs={};
    traceRefs.assign(allTraces.begin(), allTraces.end());
    makeCuts(vertices,verticesId, traceRefs,tol,mesh,countIdV,countIdE,verticesMesh,
             idVerticesMesh,edgesMesh,idEdgesMesh,0,n,mapEdges);
    EXPECT_EQ(countIdV,5);
    EXPECT_EQ(countIdE,6);
    EXPECT_EQ(mesh.verticesPolygons[0].size(),4);
    EXPECT_EQ(mesh.verticesPolygons[1].size(),3);
    EXPECT_EQ(mesh.edgesPolygons[0].size(),4);
    EXPECT_EQ(mesh.edgesPolygons[1].size(),3);

    //poligono1
    EXPECT_EQ(mesh.verticesPolygons[0][0],0);
    EXPECT_EQ(mesh.verticesPolygons[0][1],1);
    EXPECT_EQ(mesh.verticesPolygons[0][2],4);
    EXPECT_EQ(mesh.verticesPolygons[0][3],3);

    EXPECT_EQ(mesh.edgesPolygons[0][0],0);
    EXPECT_EQ(mesh.edgesPolygons[0][1],1);
    EXPECT_EQ(mesh.edgesPolygons[0][2],2);
    EXPECT_EQ(mesh.edgesPolygons[0][3],3);

    //poligono2
    EXPECT_EQ(mesh.verticesPolygons[1][0],4);
    EXPECT_EQ(mesh.verticesPolygons[1][1],2);
    EXPECT_EQ(mesh.verticesPolygons[1][2],3);


    EXPECT_EQ(mesh.edgesPolygons[1][0],4);
    EXPECT_EQ(mesh.edgesPolygons[1][1],5);
    EXPECT_EQ(mesh.edgesPolygons[1][2],2);


    //vertici
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),a,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),b,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),c,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),d,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),Vector3d(2,1,0),tol));
    verticesMesh.pop_front();


    //lati
    ok={0,1};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={1,4};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={4,3};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={3,0};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={4,2};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={2,3};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();


    //traccia su 2 vertici (non il primo)
    countIdV=4;
    countIdE=0;
    tr.extremitiesCoord={Vector3d(0,2,0),Vector3d(2,0,0)};
    tr.length=sqrt(8);
    allTraces={};
    allTraces.push_back(tr);
    verticesMesh={a,b,c,d};
    edgesMesh={};
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
    mesh.verticesPolygons={};
    mesh.edgesPolygons={};
    traceRefs={};
    traceRefs.assign(allTraces.begin(), allTraces.end());
    makeCuts(vertices,verticesId, traceRefs,tol,mesh,countIdV,countIdE,verticesMesh,
             idVerticesMesh,edgesMesh,idEdgesMesh,0,n,mapEdges);
    EXPECT_EQ(countIdV,4);
    EXPECT_EQ(countIdE,5);
    EXPECT_EQ(mesh.verticesPolygons[0].size(),3);
    EXPECT_EQ(mesh.verticesPolygons[1].size(),3);
    EXPECT_EQ(mesh.edgesPolygons[0].size(),3);
    EXPECT_EQ(mesh.edgesPolygons[1].size(),3);

    //poligono1
    EXPECT_EQ(mesh.verticesPolygons[0][0],0);
    EXPECT_EQ(mesh.verticesPolygons[0][1],1);
    EXPECT_EQ(mesh.verticesPolygons[0][2],3);

    EXPECT_EQ(mesh.edgesPolygons[0][0],0);
    EXPECT_EQ(mesh.edgesPolygons[0][1],1);
    EXPECT_EQ(mesh.edgesPolygons[0][2],2);

    //poligono2
    EXPECT_EQ(mesh.verticesPolygons[1][0],1);
    EXPECT_EQ(mesh.verticesPolygons[1][1],2);
    EXPECT_EQ(mesh.verticesPolygons[1][2],3);


    EXPECT_EQ(mesh.edgesPolygons[1][0],3);
    EXPECT_EQ(mesh.edgesPolygons[1][1],4);
    EXPECT_EQ(mesh.edgesPolygons[1][2],1);


    //vertici
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),a,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),b,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),c,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),d,tol));
    verticesMesh.pop_front();


    //lati
    ok={0,1};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={1,3};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={3,0};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={1,2};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={2,3};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();


    //caso OnThePlane: primo lato
    countIdV=4;
    countIdE=0;
    tr.extremitiesCoord={Vector3d(0.5,0,0),Vector3d(1.5,0,0)};
    tr.length=1;
    tr.onThePlane={true,false};
    allTraces={};
    allTraces.push_back(tr);
    verticesMesh={a,b,c,d};
    edgesMesh={};
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
    mesh.verticesPolygons={};
    mesh.edgesPolygons={};
    traceRefs={};
    traceRefs.assign(allTraces.begin(), allTraces.end());
    makeCuts(vertices,verticesId, traceRefs,tol,mesh,countIdV,countIdE,verticesMesh,
             idVerticesMesh,edgesMesh,idEdgesMesh,0,n,mapEdges);
    EXPECT_EQ(countIdV,6);
    EXPECT_EQ(countIdE,6);
    EXPECT_EQ(mesh.verticesPolygons[0].size(),6);
    EXPECT_EQ(mesh.edgesPolygons[0].size(),6);
    EXPECT_EQ(mesh.verticesPolygons.size(),1);
    EXPECT_EQ(mesh.edgesPolygons.size(),1);

    //poligono
    EXPECT_EQ(mesh.verticesPolygons[0][0],0);
    EXPECT_EQ(mesh.verticesPolygons[0][1],4);
    EXPECT_EQ(mesh.verticesPolygons[0][2],5);
    EXPECT_EQ(mesh.verticesPolygons[0][3],1);
    EXPECT_EQ(mesh.verticesPolygons[0][4],2);
    EXPECT_EQ(mesh.verticesPolygons[0][5],3);

    EXPECT_EQ(mesh.edgesPolygons[0][0],0);
    EXPECT_EQ(mesh.edgesPolygons[0][1],1);
    EXPECT_EQ(mesh.edgesPolygons[0][2],2);
    EXPECT_EQ(mesh.edgesPolygons[0][3],3);
    EXPECT_EQ(mesh.edgesPolygons[0][4],4);
    EXPECT_EQ(mesh.edgesPolygons[0][5],5);

    //vertici
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),a,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),b,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),c,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),d,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),Vector3d(0.5,0,0),tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),Vector3d(1.5,0,0),tol));
    verticesMesh.pop_front();

    //lati
    ok={0,4};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={4,5};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={5,1};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={1,2};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={2,3};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={3,0};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();


    //caso OnThePlane: la traccia è in un vertice
    countIdV=4;
    countIdE=0;
    tr.extremitiesCoord={Vector3d(1,0,0),Vector3d(2,0,0)};
    tr.length=1;
    tr.onThePlane={true,false};
    allTraces={};
    allTraces.push_back(tr);
    verticesMesh={a,b,c,d};
    edgesMesh={};
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
    mesh.verticesPolygons={};
    mesh.edgesPolygons={};
    traceRefs={};
    traceRefs.assign(allTraces.begin(), allTraces.end());
    makeCuts(vertices,verticesId, traceRefs,tol,mesh,countIdV,countIdE,verticesMesh,
             idVerticesMesh,edgesMesh,idEdgesMesh,0,n,mapEdges);
    EXPECT_EQ(countIdV,5);
    EXPECT_EQ(countIdE,5);
    EXPECT_EQ(mesh.verticesPolygons[0].size(),5);
    EXPECT_EQ(mesh.edgesPolygons[0].size(),5);
    EXPECT_EQ(mesh.verticesPolygons.size(),1);
    EXPECT_EQ(mesh.edgesPolygons.size(),1);

    //poligono
    EXPECT_EQ(mesh.verticesPolygons[0][0],0);
    EXPECT_EQ(mesh.verticesPolygons[0][1],4);
    EXPECT_EQ(mesh.verticesPolygons[0][2],1);
    EXPECT_EQ(mesh.verticesPolygons[0][3],2);
    EXPECT_EQ(mesh.verticesPolygons[0][4],3);

    EXPECT_EQ(mesh.edgesPolygons[0][0],0);
    EXPECT_EQ(mesh.edgesPolygons[0][1],1);
    EXPECT_EQ(mesh.edgesPolygons[0][2],2);
    EXPECT_EQ(mesh.edgesPolygons[0][3],3);
    EXPECT_EQ(mesh.edgesPolygons[0][4],4);

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

    //lati
    ok={0,4};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={4,1};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={1,2};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={2,3};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={3,0};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();

    //caso OnThePlane: la traccia è sui primi due vertici
    countIdV=4;
    countIdE=0;
    tr.extremitiesCoord={Vector3d(0,0,0),Vector3d(2,0,0)};
    tr.length=2;
    tr.onThePlane={true,false};
    allTraces={};
    allTraces.push_back(tr);
    verticesMesh={a,b,c,d};
    edgesMesh={};
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
    mesh.verticesPolygons={};
    mesh.edgesPolygons={};
    traceRefs={};
    traceRefs.assign(allTraces.begin(), allTraces.end());
    makeCuts(vertices,verticesId, traceRefs,tol,mesh,countIdV,countIdE,verticesMesh,
             idVerticesMesh,edgesMesh,idEdgesMesh,0,n,mapEdges);
    EXPECT_EQ(countIdV,4);
    EXPECT_EQ(countIdE,4);
    EXPECT_EQ(mesh.verticesPolygons[0].size(),4);
    EXPECT_EQ(mesh.edgesPolygons[0].size(),4);
    EXPECT_EQ(mesh.verticesPolygons.size(),1);
    EXPECT_EQ(mesh.edgesPolygons.size(),1);

    //poligono
    EXPECT_EQ(mesh.verticesPolygons[0][0],0);
    EXPECT_EQ(mesh.verticesPolygons[0][1],1);
    EXPECT_EQ(mesh.verticesPolygons[0][2],2);
    EXPECT_EQ(mesh.verticesPolygons[0][3],3);

    EXPECT_EQ(mesh.edgesPolygons[0][0],0);
    EXPECT_EQ(mesh.edgesPolygons[0][1],1);
    EXPECT_EQ(mesh.edgesPolygons[0][2],2);
    EXPECT_EQ(mesh.edgesPolygons[0][3],3);

    //vertici
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),a,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),b,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),c,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),d,tol));
    verticesMesh.pop_front();

    //lati
    ok={0,1};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={1,2};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={2,3};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={3,0};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();

    //caso OnThePlane: la traccia è su due vertici (non sui primi due)
    countIdV=4;
    countIdE=0;
    tr.extremitiesCoord={Vector3d(2,0,0),Vector3d(2,2,0)};
    tr.length=2;
    tr.onThePlane={true,false};
    allTraces={};
    allTraces.push_back(tr);
    verticesMesh={a,b,c,d};
    edgesMesh={};
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
    mesh.verticesPolygons={};
    mesh.edgesPolygons={};
    traceRefs={};
    traceRefs.assign(allTraces.begin(), allTraces.end());
    makeCuts(vertices,verticesId, traceRefs,tol,mesh,countIdV,countIdE,verticesMesh,
             idVerticesMesh,edgesMesh,idEdgesMesh,0,n,mapEdges);
    EXPECT_EQ(countIdV,4);
    EXPECT_EQ(countIdE,4);
    EXPECT_EQ(mesh.verticesPolygons[0].size(),4);
    EXPECT_EQ(mesh.edgesPolygons[0].size(),4);
    EXPECT_EQ(mesh.verticesPolygons.size(),1);
    EXPECT_EQ(mesh.edgesPolygons.size(),1);

    //poligono
    EXPECT_EQ(mesh.verticesPolygons[0][0],0);
    EXPECT_EQ(mesh.verticesPolygons[0][1],1);
    EXPECT_EQ(mesh.verticesPolygons[0][2],2);
    EXPECT_EQ(mesh.verticesPolygons[0][3],3);

    EXPECT_EQ(mesh.edgesPolygons[0][0],0);
    EXPECT_EQ(mesh.edgesPolygons[0][1],1);
    EXPECT_EQ(mesh.edgesPolygons[0][2],2);
    EXPECT_EQ(mesh.edgesPolygons[0][3],3);

    //vertici
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),a,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),b,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),c,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),d,tol));
    verticesMesh.pop_front();

    //lati
    ok={0,1};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={1,2};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={2,3};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={3,0};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();


    //caso OnThePlane: la traccia è su due vertici ed è più lunga del lato
    countIdV=4;
    countIdE=0;
    tr.extremitiesCoord={Vector3d(2,-1,0),Vector3d(2,2,0)};
    tr.length=3;
    tr.onThePlane={true,false};
    allTraces={};
    allTraces.push_back(tr);
    verticesMesh={a,b,c,d};
    edgesMesh={};
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
    mesh.verticesPolygons={};
    mesh.edgesPolygons={};
    traceRefs={};
    traceRefs.assign(allTraces.begin(), allTraces.end());
    makeCuts(vertices,verticesId, traceRefs,tol,mesh,countIdV,countIdE,verticesMesh,
             idVerticesMesh,edgesMesh,idEdgesMesh,0,n,mapEdges);
    EXPECT_EQ(countIdV,4);
    EXPECT_EQ(countIdE,4);
    EXPECT_EQ(mesh.verticesPolygons[0].size(),4);
    EXPECT_EQ(mesh.edgesPolygons[0].size(),4);
    EXPECT_EQ(mesh.verticesPolygons.size(),1);
    EXPECT_EQ(mesh.edgesPolygons.size(),1);

    //poligono
    EXPECT_EQ(mesh.verticesPolygons[0][0],0);
    EXPECT_EQ(mesh.verticesPolygons[0][1],1);
    EXPECT_EQ(mesh.verticesPolygons[0][2],2);
    EXPECT_EQ(mesh.verticesPolygons[0][3],3);

    EXPECT_EQ(mesh.edgesPolygons[0][0],0);
    EXPECT_EQ(mesh.edgesPolygons[0][1],1);
    EXPECT_EQ(mesh.edgesPolygons[0][2],2);
    EXPECT_EQ(mesh.edgesPolygons[0][3],3);

    //vertici
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),a,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),b,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),c,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),d,tol));
    verticesMesh.pop_front();

    //lati
    ok={0,1};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={1,2};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={2,3};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={3,0};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();

    //caso OnThePlane: la traccia è su due vertici ed è più lunga del lato da entrambe le parti
    countIdV=4;
    countIdE=0;
    tr.extremitiesCoord={Vector3d(-1,0,0),Vector3d(3,0,0)};
    tr.length=4;
    tr.onThePlane={true,false};
    allTraces={};
    allTraces.push_back(tr);
    verticesMesh={a,b,c,d};
    edgesMesh={};
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
    mesh.verticesPolygons={};
    mesh.edgesPolygons={};
    traceRefs={};
    traceRefs.assign(allTraces.begin(), allTraces.end());
    makeCuts(vertices,verticesId, traceRefs,tol,mesh,countIdV,countIdE,verticesMesh,
             idVerticesMesh,edgesMesh,idEdgesMesh,0,n,mapEdges);
    EXPECT_EQ(countIdV,4);
    EXPECT_EQ(countIdE,4);
    EXPECT_EQ(mesh.verticesPolygons[0].size(),4);
    EXPECT_EQ(mesh.edgesPolygons[0].size(),4);
    EXPECT_EQ(mesh.verticesPolygons.size(),1);
    EXPECT_EQ(mesh.edgesPolygons.size(),1);

    //poligono
    EXPECT_EQ(mesh.verticesPolygons[0][0],0);
    EXPECT_EQ(mesh.verticesPolygons[0][1],1);
    EXPECT_EQ(mesh.verticesPolygons[0][2],2);
    EXPECT_EQ(mesh.verticesPolygons[0][3],3);

    EXPECT_EQ(mesh.edgesPolygons[0][0],0);
    EXPECT_EQ(mesh.edgesPolygons[0][1],1);
    EXPECT_EQ(mesh.edgesPolygons[0][2],2);
    EXPECT_EQ(mesh.edgesPolygons[0][3],3);

    //vertici
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),a,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),b,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),c,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),d,tol));
    verticesMesh.pop_front();

    //lati
    ok={0,1};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={1,2};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={2,3};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={3,0};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();



    //ORA TESTIAMO LA RICORSIONE: FRATTURA CON PIU' TRACCE
    //caso1
    countIdV=4;
    countIdE=0;
    tr.extremitiesCoord={Vector3d(1,0,0),Vector3d(1,2,0)};
    tr.length=2;
    tr.onThePlane={false,false};
    allTraces={};
    allTraces.push_back(tr);
    Trace tr2;
    tr2.idTr=1;
    tr2.extremitiesCoord={Vector3d(0.5,0,0),Vector3d(1.5,0,0)};
    tr2.length=1;
    tr2.fracturesIds={0,2};
    tr2.Tips={false,false};
    tr2.onThePlane={true,false};
    allTraces.push_back(tr2);
    verticesMesh={a,b,c,d};
    edgesMesh={};
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
    mesh.verticesPolygons={};
    mesh.edgesPolygons={};
    traceRefs={};
    traceRefs.assign(allTraces.begin(), allTraces.end());
    makeCuts(vertices,verticesId, traceRefs,tol,mesh,countIdV,countIdE,verticesMesh,
             idVerticesMesh,edgesMesh,idEdgesMesh,0,n,mapEdges);
    EXPECT_EQ(countIdV,8);
    EXPECT_EQ(countIdE,9);
    EXPECT_EQ(mesh.verticesPolygons[0].size(),5);
    EXPECT_EQ(mesh.edgesPolygons[0].size(),5);
    EXPECT_EQ(mesh.verticesPolygons[1].size(),5);
    EXPECT_EQ(mesh.edgesPolygons[1].size(),5);

    //poligono1
    EXPECT_EQ(mesh.verticesPolygons[0][0],0);
    EXPECT_EQ(mesh.verticesPolygons[0][1],6);
    EXPECT_EQ(mesh.verticesPolygons[0][2],4);
    EXPECT_EQ(mesh.verticesPolygons[0][3],5);
    EXPECT_EQ(mesh.verticesPolygons[0][4],3);


    EXPECT_EQ(mesh.edgesPolygons[0][0],0);
    EXPECT_EQ(mesh.edgesPolygons[0][1],1);
    EXPECT_EQ(mesh.edgesPolygons[0][2],2);
    EXPECT_EQ(mesh.edgesPolygons[0][3],3);
    EXPECT_EQ(mesh.edgesPolygons[0][4],4);

    //poligono2
    EXPECT_EQ(mesh.verticesPolygons[1][0],4);
    EXPECT_EQ(mesh.verticesPolygons[1][1],7);
    EXPECT_EQ(mesh.verticesPolygons[1][2],1);
    EXPECT_EQ(mesh.verticesPolygons[1][3],2);
    EXPECT_EQ(mesh.verticesPolygons[1][4],5);

    EXPECT_EQ(mesh.edgesPolygons[1][0],5);
    EXPECT_EQ(mesh.edgesPolygons[1][1],6);
    EXPECT_EQ(mesh.edgesPolygons[1][2],7);
    EXPECT_EQ(mesh.edgesPolygons[1][3],8);
    EXPECT_EQ(mesh.edgesPolygons[1][4],2);

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
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),Vector3d(0.5,0,0),tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),Vector3d(1.5,0,0),tol));
    verticesMesh.pop_front();

    //lati
    ok={0,6};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={6,4};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={4,5};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={5,3};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={3,0};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={4,7};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={7,1};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={1,2};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={2,5};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();

    //caso2
    countIdV=4;
    countIdE=0;
    tr.extremitiesCoord={Vector3d(0,2,0),Vector3d(2,0,0)};
    tr.length=sqrt(8);
    tr.onThePlane={false,false};
    allTraces={};
    allTraces.push_back(tr);
    tr2.idTr=1;
    tr2.extremitiesCoord={Vector3d(0,0,0),Vector3d(1,1,0)};
    tr2.length=sqrt(8)/2;
    tr2.fracturesIds={0,2};
    tr2.Tips={false,false};
    tr2.onThePlane={false,false};
    allTraces.push_back(tr2);
    verticesMesh={a,b,c,d};
    edgesMesh={};
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
    mesh.verticesPolygons={};
    mesh.edgesPolygons={};
    traceRefs={};
    traceRefs.assign(allTraces.begin(), allTraces.end());
    makeCuts(vertices,verticesId, traceRefs,tol,mesh,countIdV,countIdE,verticesMesh,
             idVerticesMesh,edgesMesh,idEdgesMesh,0,n,mapEdges);
    EXPECT_EQ(countIdV,5);
    EXPECT_EQ(countIdE,8);
    EXPECT_EQ(mesh.verticesPolygons[0].size(),3);
    EXPECT_EQ(mesh.edgesPolygons[0].size(),3);
    EXPECT_EQ(mesh.verticesPolygons[1].size(),3);
    EXPECT_EQ(mesh.edgesPolygons[1].size(),3);
    EXPECT_EQ(mesh.verticesPolygons[2].size(),3);
    EXPECT_EQ(mesh.edgesPolygons[2].size(),3);


    //poligono1
    EXPECT_EQ(mesh.verticesPolygons[0][0],0);
    EXPECT_EQ(mesh.verticesPolygons[0][1],1);
    EXPECT_EQ(mesh.verticesPolygons[0][2],4);

    EXPECT_EQ(mesh.edgesPolygons[0][0],0);
    EXPECT_EQ(mesh.edgesPolygons[0][1],1);
    EXPECT_EQ(mesh.edgesPolygons[0][2],2);

    //poligono2
    EXPECT_EQ(mesh.verticesPolygons[1][0],0);
    EXPECT_EQ(mesh.verticesPolygons[1][1],4);
    EXPECT_EQ(mesh.verticesPolygons[1][2],3);

    EXPECT_EQ(mesh.edgesPolygons[1][0],2);
    EXPECT_EQ(mesh.edgesPolygons[1][1],3);
    EXPECT_EQ(mesh.edgesPolygons[1][2],4);

    //poligono3
    EXPECT_EQ(mesh.verticesPolygons[2][0],1);
    EXPECT_EQ(mesh.verticesPolygons[2][1],2);
    EXPECT_EQ(mesh.verticesPolygons[2][2],3);

    EXPECT_EQ(mesh.edgesPolygons[2][0],5);
    EXPECT_EQ(mesh.edgesPolygons[2][1],6);
    EXPECT_EQ(mesh.edgesPolygons[2][2],7);

    //vertici
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),a,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),b,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),c,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),d,tol));
    verticesMesh.pop_front();
    EXPECT_TRUE(areVectorsEqual(verticesMesh.front(),Vector3d(1,1,0),tol));
    verticesMesh.pop_front();

    //lati
    ok={0,1};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={1,4};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={4,0};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={4,3};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={3,0};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={1,2};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={2,3};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={3,1};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();



    //caso4
    countIdV=4;
    countIdE=0;
    tr.extremitiesCoord={Vector3d(1,0,0),Vector3d(1,2,0)};
    tr.length=2;
    tr.onThePlane={false,false};
    allTraces={};
    allTraces.push_back(tr);
    tr2.idTr=1;
    tr2.extremitiesCoord={Vector3d(0,1,0),Vector3d(2,1,0)};
    tr2.length=2;
    tr2.fracturesIds={0,2};
    tr2.Tips={false,false};
    tr2.onThePlane={false,false};
    allTraces.push_back(tr2);
    verticesMesh={a,b,c,d};
    edgesMesh={};
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
    mesh.verticesPolygons={};
    mesh.edgesPolygons={};
    traceRefs={};
    traceRefs.assign(allTraces.begin(), allTraces.end());
    makeCuts(vertices,verticesId, traceRefs,tol,mesh,countIdV,countIdE,verticesMesh,
             idVerticesMesh,edgesMesh,idEdgesMesh,0,n,mapEdges);
    EXPECT_EQ(countIdV,9);
    EXPECT_EQ(countIdE,12);
    EXPECT_EQ(mesh.verticesPolygons.size(),4);
    EXPECT_EQ(mesh.edgesPolygons.size(),4);
    EXPECT_EQ(mesh.verticesPolygons[0].size(),4);
    EXPECT_EQ(mesh.edgesPolygons[0].size(),4);
    EXPECT_EQ(mesh.verticesPolygons[1].size(),4);
    EXPECT_EQ(mesh.edgesPolygons[1].size(),4);
    EXPECT_EQ(mesh.verticesPolygons[2].size(),4);
    EXPECT_EQ(mesh.edgesPolygons[2].size(),4);



    //poligono1
    EXPECT_EQ(mesh.verticesPolygons[0][0],0);
    EXPECT_EQ(mesh.verticesPolygons[0][1],4);
    EXPECT_EQ(mesh.verticesPolygons[0][2],6);
    EXPECT_EQ(mesh.verticesPolygons[0][3],7);

    EXPECT_EQ(mesh.edgesPolygons[0][0],0);
    EXPECT_EQ(mesh.edgesPolygons[0][1],1);
    EXPECT_EQ(mesh.edgesPolygons[0][2],2);
    EXPECT_EQ(mesh.edgesPolygons[0][3],3);

    //poligono2
    EXPECT_EQ(mesh.verticesPolygons[1][0],6);
    EXPECT_EQ(mesh.verticesPolygons[1][1],5);
    EXPECT_EQ(mesh.verticesPolygons[1][2],3);
    EXPECT_EQ(mesh.verticesPolygons[1][3],7);

    EXPECT_EQ(mesh.edgesPolygons[1][0],4);
    EXPECT_EQ(mesh.edgesPolygons[1][1],5);
    EXPECT_EQ(mesh.edgesPolygons[1][2],6);
    EXPECT_EQ(mesh.edgesPolygons[1][3],2);

    //poligono3
    EXPECT_EQ(mesh.verticesPolygons[2][0],4);
    EXPECT_EQ(mesh.verticesPolygons[2][1],1);
    EXPECT_EQ(mesh.verticesPolygons[2][2],8);
    EXPECT_EQ(mesh.verticesPolygons[2][3],6);

    EXPECT_EQ(mesh.edgesPolygons[2][0],7);
    EXPECT_EQ(mesh.edgesPolygons[2][1],8);
    EXPECT_EQ(mesh.edgesPolygons[2][2],9);
    EXPECT_EQ(mesh.edgesPolygons[2][3],1);

    //poligono4
    EXPECT_EQ(mesh.verticesPolygons[3][0],8);
    EXPECT_EQ(mesh.verticesPolygons[3][1],2);
    EXPECT_EQ(mesh.verticesPolygons[3][2],5);
    EXPECT_EQ(mesh.verticesPolygons[3][3],6);

    EXPECT_EQ(mesh.edgesPolygons[3][0],10);
    EXPECT_EQ(mesh.edgesPolygons[3][1],11);
    EXPECT_EQ(mesh.edgesPolygons[3][2],4);
    EXPECT_EQ(mesh.edgesPolygons[3][3],9);



    //lati
    ok={0,4};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={4,6};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={6,7};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={7,0};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={6,5};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={5,3};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={3,7};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={4,1};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={1,8};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={8,6};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={8,2};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={2,5};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();


    //caso4
    countIdV=4;
    countIdE=0;
    tr.extremitiesCoord={Vector3d(1,0,0),Vector3d(1,2,0)};
    tr.length=2;
    tr.onThePlane={false,false};
    allTraces={};
    allTraces.push_back(tr);
    tr2.idTr=1;
    tr2.extremitiesCoord={Vector3d(0,1,0),Vector3d(2,1,0)};
    tr2.length=sqrt(8)/2;
    tr2.fracturesIds={0,2};
    tr2.Tips={false,false};
    tr2.onThePlane={false,false};
    allTraces.push_back(tr2);
    Trace tr3;
    tr3.idTr=2;
    tr3.extremitiesCoord={Vector3d(0.5,2,0),Vector3d(2,0,0)};
    tr3.length=1;
    tr3.fracturesIds={0,3};
    tr3.Tips={false,false};
    tr3.onThePlane={false,false};
    allTraces.push_back(tr3);
    verticesMesh={a,b,c,d};
    edgesMesh={};
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
    mesh.verticesPolygons={};
    mesh.edgesPolygons={};
    traceRefs={};
    traceRefs.assign(allTraces.begin(), allTraces.end());
    makeCuts(vertices,verticesId, traceRefs,tol,mesh,countIdV,countIdE,verticesMesh,
             idVerticesMesh,edgesMesh,idEdgesMesh,0,n,mapEdges);
    EXPECT_EQ(countIdV,12);
    EXPECT_EQ(countIdE,18);
    EXPECT_EQ(mesh.verticesPolygons.size(),7);
    EXPECT_EQ(mesh.edgesPolygons.size(),7);
    EXPECT_EQ(mesh.verticesPolygons[0].size(),4);
    EXPECT_EQ(mesh.edgesPolygons[0].size(),4);
    EXPECT_EQ(mesh.verticesPolygons[1].size(),5);
    EXPECT_EQ(mesh.edgesPolygons[1].size(),5);
    EXPECT_EQ(mesh.verticesPolygons[2].size(),3);
    EXPECT_EQ(mesh.edgesPolygons[2].size(),3);
    EXPECT_EQ(mesh.verticesPolygons[3].size(),4);
    EXPECT_EQ(mesh.edgesPolygons[3].size(),4);
    EXPECT_EQ(mesh.verticesPolygons[4].size(),3);
    EXPECT_EQ(mesh.edgesPolygons[4].size(),3);
    EXPECT_EQ(mesh.verticesPolygons[5].size(),5);
    EXPECT_EQ(mesh.edgesPolygons[5].size(),5);
    EXPECT_EQ(mesh.verticesPolygons[6].size(),3);
    EXPECT_EQ(mesh.edgesPolygons[6].size(),3);



    //poligono1
    EXPECT_EQ(mesh.verticesPolygons[0][0],0);
    EXPECT_EQ(mesh.verticesPolygons[0][1],4);
    EXPECT_EQ(mesh.verticesPolygons[0][2],6);
    EXPECT_EQ(mesh.verticesPolygons[0][3],7);

    EXPECT_EQ(mesh.edgesPolygons[0][0],0);
    EXPECT_EQ(mesh.edgesPolygons[0][1],1);
    EXPECT_EQ(mesh.edgesPolygons[0][2],2);
    EXPECT_EQ(mesh.edgesPolygons[0][3],3);

    //poligono2
    EXPECT_EQ(mesh.verticesPolygons[1][0],6);
    EXPECT_EQ(mesh.verticesPolygons[1][1],8);
    EXPECT_EQ(mesh.verticesPolygons[1][2],9);
    EXPECT_EQ(mesh.verticesPolygons[1][3],3);
    EXPECT_EQ(mesh.verticesPolygons[1][4],7);

    EXPECT_EQ(mesh.edgesPolygons[1][0],4);
    EXPECT_EQ(mesh.edgesPolygons[1][1],5);
    EXPECT_EQ(mesh.edgesPolygons[1][2],6);
    EXPECT_EQ(mesh.edgesPolygons[1][3],7);
    EXPECT_EQ(mesh.edgesPolygons[1][4],2);

    //poligono3
    EXPECT_EQ(mesh.verticesPolygons[2][0],8);
    EXPECT_EQ(mesh.verticesPolygons[2][1],5);
    EXPECT_EQ(mesh.verticesPolygons[2][2],9);

    EXPECT_EQ(mesh.edgesPolygons[2][0],8);
    EXPECT_EQ(mesh.edgesPolygons[2][1],9);
    EXPECT_EQ(mesh.edgesPolygons[2][2],5);

    //poligono4
    EXPECT_EQ(mesh.verticesPolygons[3][0],4);
    EXPECT_EQ(mesh.verticesPolygons[3][1],1);
    EXPECT_EQ(mesh.verticesPolygons[3][2],11);
    EXPECT_EQ(mesh.verticesPolygons[3][3],6);

    EXPECT_EQ(mesh.edgesPolygons[3][0],10);
    EXPECT_EQ(mesh.edgesPolygons[3][1],11);
    EXPECT_EQ(mesh.edgesPolygons[3][2],12);
    EXPECT_EQ(mesh.edgesPolygons[3][3],1);

    //poligono5
    EXPECT_EQ(mesh.verticesPolygons[4][0],1);
    EXPECT_EQ(mesh.verticesPolygons[4][1],10);
    EXPECT_EQ(mesh.verticesPolygons[4][2],11);

    EXPECT_EQ(mesh.edgesPolygons[4][0],13);
    EXPECT_EQ(mesh.edgesPolygons[4][1],14);
    EXPECT_EQ(mesh.edgesPolygons[4][2],11);

    //poligono6
    EXPECT_EQ(mesh.verticesPolygons[5][0],10);
    EXPECT_EQ(mesh.verticesPolygons[5][1],2);
    EXPECT_EQ(mesh.verticesPolygons[5][2],5);
    EXPECT_EQ(mesh.verticesPolygons[5][3],8);
    EXPECT_EQ(mesh.verticesPolygons[5][4],11);

    EXPECT_EQ(mesh.edgesPolygons[5][0],15);
    EXPECT_EQ(mesh.edgesPolygons[5][1],16);
    EXPECT_EQ(mesh.edgesPolygons[5][2],8);
    EXPECT_EQ(mesh.edgesPolygons[5][3],17);
    EXPECT_EQ(mesh.edgesPolygons[5][4],14);

    //poligono7
    EXPECT_EQ(mesh.verticesPolygons[6][0],8);
    EXPECT_EQ(mesh.verticesPolygons[6][1],6);
    EXPECT_EQ(mesh.verticesPolygons[6][2],11);

    EXPECT_EQ(mesh.edgesPolygons[6][0],4);
    EXPECT_EQ(mesh.edgesPolygons[6][1],12);
    EXPECT_EQ(mesh.edgesPolygons[6][2],17);


    //lati
    ok={0,4};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={4,6};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={6,7};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={7,0};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={6,8};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={8,9};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={9,3};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={3,7};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={8,5};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={5,9};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={4,1};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={1,11};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={11,6};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={1,10};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={10,11};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={10,2};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={2,5};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();
    ok={8,11};
    EXPECT_EQ(edgesMesh.front(),ok);
    edgesMesh.pop_front();




}



TEST(POLYGONTEST,TestCutFractures){
    //test con il nostro file con casi particolari e onThePlane
    //controlliamo le fratture più particolari, le altre le abbiamo controllate a mano
    double tol=10*numeric_limits<double>::epsilon();
    vector<Fracture> vecFrac;
    string path="./DFN";
    readFractures(path+"/FR5_data_prova.txt",vecFrac, tol);
    vector<Trace> vecTraces={};
    vecTraces=findTraces(vecFrac,tol);
    printLocalResults("lresults",vecFrac,vecTraces); //la chiamo perché così ordina le tracce
    vector<PolygonalMesh> vecMesh = cutFractures(vecFrac, vecTraces,tol);

    EXPECT_EQ(vecMesh[0].numVertices,8);
    EXPECT_EQ(vecMesh[1].numVertices,3);
    EXPECT_EQ(vecMesh[2].numVertices,6);
    EXPECT_EQ(vecMesh[3].numVertices,6);
    EXPECT_EQ(vecMesh[4].numVertices,6);

    //frattura0, sottopoligono0
    EXPECT_EQ(vecMesh[0].verticesPolygons[0][0],0);
    EXPECT_EQ(vecMesh[0].verticesPolygons[0][1],4);
    EXPECT_EQ(vecMesh[0].verticesPolygons[0][2],6);
    EXPECT_EQ(vecMesh[0].verticesPolygons[0][3],7);

    EXPECT_EQ(vecMesh[0].edgesPolygons[0][0],0);
    EXPECT_EQ(vecMesh[0].edgesPolygons[0][1],1);
    EXPECT_EQ(vecMesh[0].edgesPolygons[0][2],2);
    EXPECT_EQ(vecMesh[0].edgesPolygons[0][3],3);

    //frattura0, sottopoligono1
    EXPECT_EQ(vecMesh[0].verticesPolygons[1][0],6);
    EXPECT_EQ(vecMesh[0].verticesPolygons[1][1],5);
    EXPECT_EQ(vecMesh[0].verticesPolygons[1][2],3);
    EXPECT_EQ(vecMesh[0].verticesPolygons[1][3],7);

    EXPECT_EQ(vecMesh[0].edgesPolygons[1][0],4);
    EXPECT_EQ(vecMesh[0].edgesPolygons[1][1],5);
    EXPECT_EQ(vecMesh[0].edgesPolygons[1][2],6);
    EXPECT_EQ(vecMesh[0].edgesPolygons[1][3],2);

    //frattura0, sottopoligono2
    EXPECT_EQ(vecMesh[0].verticesPolygons[2][0],4);
    EXPECT_EQ(vecMesh[0].verticesPolygons[2][1],1);
    EXPECT_EQ(vecMesh[0].verticesPolygons[2][2],2);
    EXPECT_EQ(vecMesh[0].verticesPolygons[2][3],5);

    EXPECT_EQ(vecMesh[0].edgesPolygons[2][0],7);
    EXPECT_EQ(vecMesh[0].edgesPolygons[2][1],8);
    EXPECT_EQ(vecMesh[0].edgesPolygons[2][2],9);
    EXPECT_EQ(vecMesh[0].edgesPolygons[2][3],10);

    //frattura1, sottopoligono0
    EXPECT_EQ(vecMesh[1].verticesPolygons[0][0],0);
    EXPECT_EQ(vecMesh[1].verticesPolygons[0][1],1);
    EXPECT_EQ(vecMesh[1].verticesPolygons[0][2],2);

    EXPECT_EQ(vecMesh[1].edgesPolygons[0][0],0);
    EXPECT_EQ(vecMesh[1].edgesPolygons[0][1],1);
    EXPECT_EQ(vecMesh[1].edgesPolygons[0][2],2);

    //frattura3, sottopoligono0
    EXPECT_EQ(vecMesh[3].verticesPolygons[0][0],0);
    EXPECT_EQ(vecMesh[3].verticesPolygons[0][1],1);
    EXPECT_EQ(vecMesh[3].verticesPolygons[0][2],4);
    EXPECT_EQ(vecMesh[3].verticesPolygons[0][3],5);

    EXPECT_EQ(vecMesh[3].edgesPolygons[0][0],0);
    EXPECT_EQ(vecMesh[3].edgesPolygons[0][1],1);
    EXPECT_EQ(vecMesh[3].edgesPolygons[0][2],2);
    EXPECT_EQ(vecMesh[3].edgesPolygons[0][3],3);

    //frattura3, sottopoligono1
    EXPECT_EQ(vecMesh[3].verticesPolygons[1][0],4);
    EXPECT_EQ(vecMesh[3].verticesPolygons[1][1],2);
    EXPECT_EQ(vecMesh[3].verticesPolygons[1][2],3);
    EXPECT_EQ(vecMesh[3].verticesPolygons[1][3],5);

    EXPECT_EQ(vecMesh[3].edgesPolygons[1][0],4);
    EXPECT_EQ(vecMesh[3].edgesPolygons[1][1],5);
    EXPECT_EQ(vecMesh[3].edgesPolygons[1][2],6);
    EXPECT_EQ(vecMesh[3].edgesPolygons[1][3],2);
}

TEST(POLYGONTEST,TestAddVerticesOnThePlane){

    double tol=10*numeric_limits<double>::epsilon();
    queue<Vector3d> subvertices1={};
    queue<unsigned int> subverticesId1={};
    list<Vector3d> verticesMesh={};
    list<unsigned int> idVerticesMesh={};
    unsigned int countIdV=1;
    Vector3d p=Vector3d(0,0,0);
    Vector3d c=Vector3d(5,0,0);
    Vector3d e0=Vector3d(4,0,0);
    Vector3d e1=Vector3d(7,0,0);

    addVerticesOnThePlane(subvertices1,subverticesId1,verticesMesh,idVerticesMesh,countIdV,c,p,e0,e1,tol);
    EXPECT_TRUE(areVectorsEqual(e0,subvertices1.front(),tol));
    EXPECT_EQ(countIdV,2);

    subvertices1={};
    countIdV=1;
    e1=Vector3d(2,0,0);
    addVerticesOnThePlane(subvertices1,subverticesId1,verticesMesh,idVerticesMesh,countIdV,c,p,e0,e1,tol);
    EXPECT_TRUE(areVectorsEqual(e1,subvertices1.front(),tol));
    subvertices1.pop();
    EXPECT_TRUE(areVectorsEqual(e0,subvertices1.front(),tol));
    EXPECT_EQ(countIdV,3);

    subvertices1={};
    countIdV=1;
    e0=Vector3d(0,0,0);
    e1=Vector3d(7,0,0);
    addVerticesOnThePlane(subvertices1,subverticesId1,verticesMesh,idVerticesMesh,countIdV,c,p,e0,e1,tol);
    EXPECT_EQ(subvertices1.size(),0);
    EXPECT_EQ(countIdV,1);

    subvertices1={};
    countIdV=1;
    e0=Vector3d(-1,0,0);
    e1=Vector3d(5,0,0);
    addVerticesOnThePlane(subvertices1,subverticesId1,verticesMesh,idVerticesMesh,countIdV,c,p,e0,e1,tol);
    EXPECT_EQ(subvertices1.size(),0);
    EXPECT_EQ(countIdV,1);

    subvertices1={};
    countIdV=1;
    e0=Vector3d(0,0,0);
    e1=Vector3d(5,0,0);
    addVerticesOnThePlane(subvertices1,subverticesId1,verticesMesh,idVerticesMesh,countIdV,c,p,e0,e1,tol);
    EXPECT_EQ(subvertices1.size(),0);
    EXPECT_EQ(countIdV,1);
}
}

#endif
