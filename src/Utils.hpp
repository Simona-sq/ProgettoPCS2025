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

pair<vector<Vector3d>, vector<vector<unsigned int>>> getSolidData(unsigned int& q);

Polyhedron buildPlatonicSolid(unsigned int& q);

Polyhedron triangulateClass1(PolyhedronLibrary::Polyhedron& P, unsigned int& t_value);

//Polyhedron triangulateClass2(PolyhedronLibrary::Polyhedron& P, unsigned int& t_value);

vector<unsigned int> cycled_face_for_dual(vector<unsigned int>& face_new, const Eigen::MatrixXd& coord);

Polyhedron Dualize(const Polyhedron& P_original);

}
