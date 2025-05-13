#pragma once

#include <iostream>
#include "Eigen/Eigen"

// creo la struttura Polyhedron

using namespace std;
using namespace Eigen;


namespace PolyhedronLibrary {

struct Polyhedron
{
    unsigned int NumCell0Ds = 0; // numero di vertici
    std::vector<unsigned int> Cell0DsId = {}; // Cell0DsId è il vettore di ID dei vertici (dimensione 1 x NumCell0Ds)
    Eigen::MatrixXd Cell0DsCoordinates = {}; // Cell0DCoordinates è la matrice delle coordinate (dimensione 3 x NumCell0Ds)
    std::map<unsigned int, vector<double>> IdCell0Ds = {}; // IdCell0Ds è un dizionario di id_punto - coordinate_punto


    unsigned int NumCell1Ds = 0; // numero di lati
    std::vector<unsigned int> Cell1DsId = {}; // Cell1DsId è il vettore di ID dei lati (dimensione 1 x NumCell1Ds)
    Eigen::MatrixXi Cell1DsExtrema = {}; // Cell1DExtrema è la matrice contenente gli ID dei vertici iiniziale e finale che caratterizzano il lato (dimensione 2 x NumCell1Ds)
    
    unsigned int NumCell2Ds = 0; // numero di facce
    std::vector<unsigned int> Cell2DsId = {}; // Cell2DsId è il vettore di ID delle facce (dimensione 1 x NumCell2Ds)
    //visto che il numero di lati è variabile, uso un vettore dinamico
    std::vector<std::vector<unsigned int>> Cell2DsVertices = {}; // Cell2DsVertices è una matrice con dimensione 1 x NumberCell2DVertices
    std::vector<std::vector<unsigned int>> Cell2DsEdges = {}; // Cell2DsEdges è una matrice con dimensione 1 x NumberCell2DEdges

    unsigned int NumCell3Ds = 1; // 1 run ==> un solo poliedro
    std::vector<unsigned int> Cell3DsId = {}; // Cell3DsId è il vettore di ID delle facce (dimensione 1 x NumCell2Ds)
    //visto che il numero di lati è variabile, uso un vettore dinamico
    std::vector<std::vector<unsigned int>> Cell3DsVertices = {}; // Cell2DsVertices è una matrice con dimensione 1 x NumberCell2DVertices
    std::vector<std::vector<unsigned int>> Cell3DsEdges = {}; // Cell2DsEdges è una matrice con dimensione 1 x NumberCell2DEdges
    std::vector<std::vector<unsigned int>> Cell3DsFaces = {}; // Cell2DsFaces è una matrice con dimensione 1 x NumberCell2DFaces

};

}
