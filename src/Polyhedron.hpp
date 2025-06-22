#pragma once 

#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;


namespace PolyhedronLibrary {

struct Polyhedron
{   
    // VERTICI
    unsigned int NumCell0Ds = 0; // NumCell0Ds contiene il numero di vertici
    vector<unsigned int> Cell0DsId = {}; // Cell0DsId contiene gli id dei vertici
    MatrixXd Cell0DsCoordinates = {}; // Cell0DCoordinates è la matrice contenente le coordinate dei vertici (x, y, z per ogni vertice)

    // LATI
    unsigned int NumCell1Ds = 0; // NumCell1Ds contiene il numero di lati
    vector<unsigned int> Cell1DsId = {}; // Cell1DsId contiene gli id dei lati
    MatrixXi Cell1DsExtrema = {}; // Cell1DExtrema è la matrice contenente gli id dei vertici iniziale e finale che caratterizzano il lato (id_v1, id_v2)
    
    // FACCE
    unsigned int NumCell2Ds = 0; // NumCell2Ds contiene il numero di facce
    vector<unsigned int> Cell2DsId = {}; // Cell2DsId contiene gli id delle facce
    vector<vector<unsigned int>> Cell2DsVertices = {}; // Cell2DsVertices è la matrice contenente gli id dei vertici che caratterizzano la faccia
    vector<vector<unsigned int>> Cell2DsEdges = {}; // Cell2DsEdges è la matrice contenente gli id dei lati che caratterizzano la faccia

    // POLIEDRO
    unsigned int NumCell3Ds = 1; // NumCell3Ds è impostato a 0 perchè considriamo un solo poliedro dato in input
    unsigned int Cell3DsId = 0; // Cell3DsId contiene l'id del poliedro
    vector<unsigned int> Cell3DsVertices = {}; // Cell3DsVertices contiene gli id dei vertici del poliedro
    vector<unsigned int> Cell3DsEdges = {}; // Cell3DsEdges contiene gli id dei lati del poliedro
    vector<unsigned int> Cell3DsFaces = {}; // Cell3DsFaces contiene gli id delle facce del poliedro
};

}

