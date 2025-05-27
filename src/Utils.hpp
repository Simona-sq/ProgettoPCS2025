#pragma once

#include "UCDUtilities.hpp"
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

Polyhedron triangulateClass1(const Polyhedron& P, const unsigned int& t_value);

Polyhedron Dualize(const Polyhedron& P_normale);

Polyhedron projectPolyhedronOnSphere(const Polyhedron& P);

void ExportPolyhedron(const Polyhedron& P,
                    const vector<Gedim::UCDProperty<double>>& points_properties = {},
                    const vector<Gedim::UCDProperty<double>>& segments_properties = {});

void Esporta_file(const Polyhedron& P);

vector<unsigned int> Cammini_minimi(const Polyhedron& P, const int& v1, const int& v2);

}
