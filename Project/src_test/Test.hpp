#ifndef __TEST_HPP
#define __TEST_HPP

#include <gtest/gtest.h>
#include "Fractures.hpp"
#include "Utils.hpp"
#include "Eigen/Eigen"
#include <cmath>
#include <iostream>

using namespace Eigen;
using namespace std;
using namespace Algebra;

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



}

namespace detail{

TEST(SORTINGTEST, TestMergeSort){
    vector<Trace> vecProvaMS(5);
    Trace t1;
    t1.idTr=0;
    t1.length=1.0;
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
    vector<unsigned int> idProvaMS ={0,1,2,3,4};
    //ordine:1403
    vecProvaMS.push_back(t1);
    vecProvaMS.push_back(t2);
    vecProvaMS.push_back(t3);
    vecProvaMS.push_back(t4);
    vecProvaMS.push_back(t5);
    detail::mergesort(idProvaMS, vecProvaMS,0,4);
    EXPECT_EQ(idProvaMS[0],1);
    EXPECT_EQ(idProvaMS[1],2);
    EXPECT_EQ(idProvaMS[2],0);
    EXPECT_EQ(idProvaMS[3],4);

    //provo con lunghezza nulla
    idProvaMS ={};



}

}




#endif
