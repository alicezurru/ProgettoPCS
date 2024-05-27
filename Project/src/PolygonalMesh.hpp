#ifndef PolygonalMesh_H
#define PolygonalMesh_H

#include <iostream>
#include <vector>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

namespace PolygonalMeshLibrary{
struct PolygonalMesh{

    unsigned int numVertices;
    vector<unsigned int> idVertices;
    vector<Vector3d> coordVertices; //ci metto dentro un vector per farci operazioni matematiche

    unsigned int numEdges;
    vector<unsigned int> idEdges;
    vector<array<unsigned int, 2>> extremitiesEdges; //sono id, non serve fare operazioni matematiche

    unsigned int numPolygons;
    vector<vector<unsigned int>> verticesPolygons; //non so a priori quanti sono e non mi serve accedere per indice (non hanno id)
    vector<vector<unsigned int>> edgesPolygons;
/*
    //Triangolazione per usare paraview: da esercitazione
    vector<vector<vector<unsigned int>>> Triangulation()
    {
        const unsigned int NumPolygons = verticesPolygons.size();
        vector<vector<vector<unsigned int>>> triangleList(NumPolygons); //vettore di poligoni che ha un vettore di triangoli che ha un vettore di id di vertici

        for(unsigned int p = 0; p < NumPolygons; p++)
        {
            const unsigned int numPolygonVertices = verticesPolygons[p].size();

            for (unsigned int v = 0; v < numPolygonVertices; v++)
            {
                const unsigned int nextVertex = verticesPolygons[p][(v + 1) % numPolygonVertices];
                const unsigned int nextNextVertex = verticesPolygons[p][(v + 2) % numPolygonVertices];

                if ((v + 2) % numPolygonVertices == 0)
                    break;

                vector<unsigned int> triangle_vertices = {verticesPolygons[p][0], nextVertex, nextNextVertex};

                triangleList[p].push_back(triangle_vertices);
            }
        }
        return triangleList;
    }

    void GedimInterface(vector<vector<unsigned int>>& triangles,
                        VectorXi& materials);*/
};

}


#endif
