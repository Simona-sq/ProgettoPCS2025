#pragma once

#include "UCDUtilities.hpp"
#include <iostream>
#include "Polyhedron.hpp"

using namespace std;

namespace PolyhedronLibrary
{
pair<vector<Vector3d>, vector<vector<unsigned int>>> getSolidData(unsigned int& q);

Polyhedron buildPlatonicSolid(unsigned int& q);

Polyhedron triangulateClass1(const Polyhedron& P, const unsigned int& t_value);

Polyhedron triangulateClass2(const Polyhedron& P, const unsigned int& b);

Polyhedron Dualize(const Polyhedron& P_normale);

Polyhedron projectPolyhedronOnSphere(const Polyhedron& P);

void ExportPolyhedron(const Polyhedron& P,
                        const vector<Gedim::UCDProperty<double>>& points_properties = {},
                        const vector<Gedim::UCDProperty<double>>& segments_properties = {});

void Esporta_file(const Polyhedron& P);

vector<unsigned int> Cammini_minimi(const Polyhedron& P, const unsigned int& v1, const unsigned int& v2);
}
