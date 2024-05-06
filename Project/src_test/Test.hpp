#ifndef __TESTPOLYGONS_H
#define __TESTPOLYGONS_H

#include <gtest/gtest.h>
#include "Fractures.hpp"
#include "Utils.hpp"
#include "Eigen/Eigen"
#include <iostream>
//#include "UCDUtilities.hpp"

using namespace Eigen;
using namespace std;
//using namespace Algebra;

namespace Geometry {
/*TEST(FRACTURESTEST, TestReadingFractures){
    //facciamo un test per la lettura da file?
    //EXPECT_EQ();
}*/
};

namespace Algebra{
TEST(ALGEBRATEST, TestFindPlaneEquation){
    vector<Vector3d> points = {Vector3d(1.0, 2.5, -0.5), Vector3d(2.0, -1.0, 0.0), Vector3d(0.0, 2.0, 1.0)};
    //non importa quanti punti gli do, userà sempre solo i primi 3, quindi posso dargli direttamente 3 punti
    double constantTerm=0;
    Vector3d n = findPlaneEquation(points, constantTerm);
    EXPECT_EQ(n, Vector3d(-5, -2, -4));
    EXPECT_EQ(constantTerm, 8);
}

/*TEST(ALGEBRATEST, TestIntersectionPlaneLine){

}

TEST(ALGEBRATEST, TestFindIntersectionPoints){

}*/

};

namespace detail{
/*TEST(SORTINGTEST, TestMerge){

}

TEST(SORTINGTEST, TestMergeSort){

}*/

}




#endif
