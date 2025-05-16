#pragma once

#include <iostream>
#include "Polyhedron.hpp"

using namespace std;

namespace PolyhedronLibrary
{
/// Import the triangular mesh and test if the mesh is correct
/// mesh: a TriangularMesh struct
/// return the result of the reading, true if is success, false otherwise
bool ImportMesh(Polyhedron& mesh);

/// Import the Cell0D properties from Cell0Ds.csv file
/// mesh: a TriangularMesh struct
/// return the result of the reading, true if is success, false otherwise
bool ImportCell0Ds(Polyhedron& mesh);

/// Import the Cell1D properties from Cell1Ds.csv file
/// mesh: a TriangularMesh struct
/// return the result of the reading, true if is success, false otherwise
bool ImportCell1Ds(Polyhedron& mesh);

/// Import the Cell2D properties from Cell2Ds.csv file
/// mesh: a TriangularMesh struct
/// return the result of the reading, true if is success, false otherwise
bool ImportCell2Ds(Polyhedron& mesh);

Vector3d normalize(Vector3d v);

pair<vector<Vector3d>, vector<vector<unsigned int>>> getSolidData(unsigned int& q);

Polyhedron buildPlatonicSolid(unsigned int& p, unsigned int& q, unsigned int& b, unsigned int& c);

}
