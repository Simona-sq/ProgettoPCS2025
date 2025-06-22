#pragma once

#include <iostream>
#include "Eigen/Eigen"

// creo la struttura Polyhedron
using namespace std;
using namespace Eigen;


namespace PolyhedronLibrary {

struct Polyhedron
{
    unsigned int NumCell0Ds = 0; // Numero di vertici
    vector<unsigned int> Cell0DsId = {}; // Cell0DsId contiene gli id dei vertici
    MatrixXd Cell0DsCoordinates = {}; // Cell0DCoordinates è la matrice delle coordinate dei vertici

    unsigned int NumCell1Ds = 0; // Numero di lati
    vector<unsigned int> Cell1DsId = {}; // Cell1DsId contiene gli id dei lati
    MatrixXi Cell1DsExtrema = {}; // Cell1DExtrema è la matrice contenente gli id dei vertici iniziale e finale che caratterizzano il lato
    
    unsigned int NumCell2Ds = 0; // Numero di facce
    vector<unsigned int> Cell2DsId = {}; // Cell2DsId contiene gli id delle facce
    vector<vector<unsigned int>> Cell2DsVertices = {}; // Cell2DsVertices è la matrice contenente gli id dei vertici che caratterizzano la faccia
    vector<vector<unsigned int>> Cell2DsEdges = {}; // Cell2DsEdges è la matrice contenente gli id dei lati che caratterizzano la faccia

    unsigned int NumCell3Ds = 1; // Un solo poliedro
    unsigned int Cell3DsId = 0; // E' l'ID del poliedtro 
    vector<vector<unsigned int>> Cell3DsVertices = {}; // Cell3DsVertices contiene gli id dei vertici del poliedro
    vector<vector<unsigned int>> Cell3DsEdges = {}; // Cell3DsEdges contiene gli id dei lati del poliedro
    vector<vector<unsigned int>> Cell3DsFaces = {}; // Cell3DsFaces contiene gli id delle facce del poliedro

};

}
