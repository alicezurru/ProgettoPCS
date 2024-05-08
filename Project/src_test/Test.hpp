#ifndef __TEST_HPP
#define __TEST_HPP

#include <gtest/gtest.h>
#include "Fractures.hpp"
#include "Utils.hpp"
#include "Eigen/Eigen"
#include <iostream>

using namespace Eigen;
using namespace std;
//using namespace Algebra;

namespace Geometry {
/*TEST(FRACTURESTEST, TestReadingFractures){
    //facciamo un test per la lettura da file?
    //EXPECT_EQ();
}*/
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
    ASSERT_EQ(intersection, Vector3d(0,1,0));

}
//bool findIntersectionPoints(Fracture& f1, Fracture& f2, array<Vector3d,4>& intPoints, double tol);
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
    bool intersection = findIntersectionPoints(f1,f2,intPoints,tol,onThePlane);
    EXPECT_TRUE(intersection);
    array<Vector3d,4> ok={Vector3d(2,0,0.5),Vector3d(0,0,0.5),Vector3d(3,0,0.5),Vector3d(1,0,0.5)};
    EXPECT_EQ(intPoints, ok);

    //testo un caso in cui non c'è intersezione ma i piani non sono paralleli
    //tengo f1 come prima
    f2.vertices={Vector3d(10,-0.5,0.5),Vector3d(12,-0.5,0.5),Vector3d(12,0.5,0.5),Vector3d(10,0.5,0.5)};
    intersection=true;
    intersection = findIntersectionPoints(f1,f2,intPoints,tol,onThePlane);
    EXPECT_FALSE(intersection);//DA CAMBIARE


    //testo un caso in cui non c'è intersezione perchè i piani sono paralleli
    f1.numVertices=5;
    f1.vertices={Vector3d(0,0,1),Vector3d(1,-1,1),Vector3d(0,-2,-1),Vector3d(-1,0,0),Vector3d(0.5,0.5,2)};
    f2.numVertices=3;
    f2.vertices={Vector3d(0,1,1),Vector3d(1,0,1),Vector3d(0.5,0.5,1)};
    intersection=true;
    intersection = findIntersectionPoints(f1,f2,intPoints,tol,onThePlane);
    EXPECT_FALSE(intersection);

    //




}

}

namespace detail{
/*TEST(SORTINGTEST, TestMerge){

}

TEST(SORTINGTEST, TestMergeSort){

}*/

}




#endif
