#pragma once

#include "UCDUtilities.hpp"
#include <iostream>
#include "Polyhedron.hpp"

using namespace std;

namespace PolyhedronLibrary
{
pair<vector<Vector3d>, vector<vector<unsigned int>>> SalvataggioDati(unsigned int& q);

Polyhedron CreazioneSolidoPlatonico(unsigned int& q);

Polyhedron TriangolazioneClasse1(const Polyhedron& P, const unsigned int& t_value);

Polyhedron TriangolazioneClasse2(const Polyhedron& P, const unsigned int& b);

Polyhedron Dualizzazione(const Polyhedron& P);

Polyhedron ProiezioneSullaSfera(const Polyhedron& P);

void EsportazionePoliedro(const Polyhedron& P,
                        const vector<Gedim::UCDProperty<double>>& points_properties = {},
                        const vector<Gedim::UCDProperty<double>>& segments_properties = {});

void EsportazioneFile(const Polyhedron& P);

vector<unsigned int> CalcoloCamminoMinimo(const Polyhedron& P, const unsigned int& v1, const unsigned int& v2);
}
